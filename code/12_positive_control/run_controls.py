#!/usr/bin/env python3
"""Sanity controls: positive (Tabula Muris x height), negative (Smillie x SCZ),
and MHC sensitivity (Smillie x de Lange MHC-included).

Spec: see ./README.md.

Invokes scDRS via subprocess using the same conventions as
scripts/slurm/03_scdrs_compute.slurm, then runs the group analysis.
Warn-don't-fail semantics for the assertions (per README sanity-track
caveats).
"""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
import sys
import time
from pathlib import Path

LOGGER = logging.getLogger(Path(__file__).stem)

_REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO / "code"))

from _shared.constants import SCDRS_N_CTRL, SEED  # noqa: E402


def require_path(p: Path, descr: str) -> None:
    if not p.exists():
        LOGGER.error("Missing input: %s (%s)", p, descr)
        raise SystemExit(2)


def run_scdrs(h5ad_path: Path, gene_set: Path, out_dir: Path,
              cov_file: Path | None, seed: int) -> Path:
    """Run scdrs compute-score; return the *.full_score.gz path."""
    if shutil.which("scdrs") is None:
        LOGGER.error("scdrs CLI not on PATH; install per pyproject [atlas].")
        raise SystemExit(2)
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "scdrs", "compute-score",
        "--h5ad-file", str(h5ad_path),
        "--gs-file", str(gene_set),
        "--flag-filter-data", "True",
        "--flag-raw-count", "False",
        "--n-ctrl", str(SCDRS_N_CTRL),
        "--random-seed", str(seed),
        "--out-folder", str(out_dir) + "/",
    ]
    if cov_file is not None and cov_file.exists():
        cmd.extend(["--cov-file", str(cov_file)])
    LOGGER.info("scdrs compute-score: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)
    scores = sorted(out_dir.glob("*.full_score.gz"))
    if not scores:
        raise RuntimeError(f"scdrs compute-score produced no scores in {out_dir}")
    return scores[0]


def run_group_analysis(h5ad_path: Path, score_file: Path,
                       tier: str, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "scdrs", "perform-downstream",
        "--h5ad-file", str(h5ad_path),
        "--score-file", str(score_file),
        "--group-analysis", f"cell_type_{tier}",
        "--out-folder", str(out_dir) + "/",
    ]
    LOGGER.info("scdrs perform-downstream: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)
    groups = sorted(out_dir.glob("*.scdrs_group"))
    if not groups:
        raise RuntimeError(f"perform-downstream produced no group file in {out_dir}")
    return groups[0]


def assert_at_least_one_significant(
    group_file: Path, fdr: float, predicate, predicate_name: str,
) -> bool:
    """Log a warning (do not fail) when predicate-matching rows lack signal."""
    import pandas as pd
    df = pd.read_csv(group_file, sep="\t")
    fdr_col = "assoc_mcq" if "assoc_mcq" in df.columns else "fdr"
    name_col = "group" if "group" in df.columns else "cell_type"
    match = df[df[name_col].astype(str).map(predicate)]
    sig = match[match[fdr_col] < fdr]
    LOGGER.info("Positive control match=%d sig=%d (predicate=%s)",
                len(match), len(sig), predicate_name)
    if sig.empty:
        LOGGER.warning(
            "Positive control: no %s cell types pass FDR < %.2f. "
            "Informational; review pipeline.",
            predicate_name, fdr,
        )
        return False
    return True


def assert_no_colon_significant(
    group_file: Path, fdr: float, colon_predicate,
) -> bool:
    """Negative control: warn if any colon cell type has FDR < threshold."""
    import pandas as pd
    df = pd.read_csv(group_file, sep="\t")
    fdr_col = "assoc_mcq" if "assoc_mcq" in df.columns else "fdr"
    name_col = "group" if "group" in df.columns else "cell_type"
    match = df[df[name_col].astype(str).map(colon_predicate)]
    sig = match[match[fdr_col] < fdr]
    if not sig.empty:
        LOGGER.warning(
            "Negative control: %d colon cell types pass FDR < %.2f for SCZ. "
            "Possible MHC / complement leakage — see README caveat.",
            len(sig), fdr,
        )
        return False
    return True


def mhc_sensitivity(group_file_mhc_included: Path,
                    group_file_baseline: Path) -> None:
    """Compare cell-type rankings; log Spearman rho and top-5 Jaccard."""
    import pandas as pd
    from scipy.stats import spearmanr

    def _load(p):
        df = pd.read_csv(p, sep="\t")
        z_col = "assoc_mcz" if "assoc_mcz" in df.columns else "z_mean"
        name_col = "group" if "group" in df.columns else "cell_type"
        return df[[name_col, z_col]].rename(
            columns={name_col: "cell_type", z_col: "z"}
        )

    a = _load(group_file_mhc_included)
    b = _load(group_file_baseline)
    merged = a.merge(b, on="cell_type", suffixes=("_mhc", "_base"))
    if len(merged) < 3:
        LOGGER.warning("MHC sensitivity: only %d shared cell types", len(merged))
        return
    rho, _ = spearmanr(merged["z_mhc"], merged["z_base"])
    top5_mhc = set(merged.nlargest(5, "z_mhc")["cell_type"])
    top5_base = set(merged.nlargest(5, "z_base")["cell_type"])
    union = top5_mhc | top5_base
    jaccard = len(top5_mhc & top5_base) / len(union) if union else float("nan")
    LOGGER.info("MHC sensitivity: spearman_rho=%.3f top5_jaccard=%.3f",
                rho, jaccard)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--control", required=True,
                        choices=["positive", "negative", "mhc-sensitivity"])
    parser.add_argument("--atlas", required=True)
    parser.add_argument("--h5ad-path", type=Path, required=True)
    parser.add_argument("--gwas-gene-set", type=Path, required=True,
                        help="Path to scDRS-format gene set (.gs)")
    parser.add_argument("--baseline-group-file", type=Path,
                        help="For --control mhc-sensitivity: the MHC-excluded "
                             "scDRS group file to compare against.")
    parser.add_argument("--cov-file", type=Path,
                        help="Optional covariate TSV.")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--tier", default="broad", choices=["broad", "fine"])
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s: %(message)s",
    )

    require_path(args.h5ad_path, f"atlas h5ad ({args.atlas})")
    require_path(args.gwas_gene_set, "GWAS gene set (.gs)")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    score_file = run_scdrs(
        args.h5ad_path, args.gwas_gene_set, args.out_dir,
        args.cov_file, args.seed,
    )
    group_file = run_group_analysis(
        args.h5ad_path, score_file, args.tier, args.out_dir,
    )

    if args.control == "positive":
        def is_skel(name: str) -> bool:
            n = name.lower()
            return any(t in n for t in (
                "chondrocyt", "osteoblast", "skeletal", "musculoskeletal",
                "osteo", "muscle",
            ))
        # TODO: positive control requires mouse->human ortholog mapping
        # for the Tabula Muris atlas. The actual mapping happens in MAGMA
        # gene-set generation upstream (or in a thin wrapper here). Defer
        # to upstream gene-set construction; this driver assumes the .gs
        # file is already human-ortholog-mapped.
        assert_at_least_one_significant(
            group_file, fdr=0.05, predicate=is_skel,
            predicate_name="skeletal/musculoskeletal",
        )
    elif args.control == "negative":
        def is_colon(name: str) -> bool:
            n = name.lower()
            return any(t in n for t in (
                "colon", "epithel", "enterocyte", "goblet", "paneth", "tuft",
                "stromal", "fibroblast", "endothel", "myeloid", "t cell",
                "plasma",
            ))
        assert_no_colon_significant(group_file, fdr=0.05, colon_predicate=is_colon)
    elif args.control == "mhc-sensitivity":
        if args.baseline_group_file is None:
            LOGGER.error("--baseline-group-file required for mhc-sensitivity")
            raise SystemExit(2)
        require_path(args.baseline_group_file, "MHC-excluded baseline group")
        mhc_sensitivity(group_file, args.baseline_group_file)

    LOGGER.info("Done in %.1fs", time.time() - t0)
    return 0


# TODO (stretch #3): LDL x Tabula Sapiens hepatocytes positive control.
# Stub only; do not implement until M6 stretches are green-lit.
def stretch3_ldl_hepatocytes_stub() -> None:
    raise NotImplementedError(
        "Stretch #3 (LDL x Tabula Sapiens hepatocytes) deferred to M6."
    )


if __name__ == "__main__":
    sys.exit(main())
