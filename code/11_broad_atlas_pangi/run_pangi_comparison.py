#!/usr/bin/env python3
"""Pan-GI comparator: full + no-Elmentaite + no-Smillie sensitivities + donor audit.

Spec: see ./README.md and DECISIONS.md correction 2026-05-20 (3/7).

Loads pre-existing per-atlas method results (scDRS, seismicGWAS) for Pan-GI
(under each of three result roots: full, no-Elmentaite, no-Smillie) and
each UC trio atlas, then computes concordance with the min-50-cells filter.

Donor audit: scans h5ad obs columns for cross-atlas donor presence; runs
independently from the concordance pipelines and skips gracefully when an
atlas h5ad is missing.
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
import time
from pathlib import Path

LOGGER = logging.getLogger(Path(__file__).stem)

_REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO / "code"))
sys.path.insert(0, str(_REPO / "code" / "06_concordance"))

from _shared.constants import (  # noqa: E402
    BOOTSTRAP_N,
    METHODS,
    MIN_CELLS_PER_TYPE,
    SEED,
    TIERS,
    UC_ATLASES,
    UC_GWAS,
)
from _shared.git import git_sha  # noqa: E402
from _shared.result_loading import (  # noqa: E402
    load_scdrs_group,
    load_seismic_results,
    shared_cell_types,
    to_lookup,
)
from metrics import concordance  # noqa: E402

SMILLIE_DONOR_RE = re.compile(r"^(N|UC)\d+$")


def _load_one(method, results_base, atlas, gwas, tier):
    if method == "scdrs":
        df = load_scdrs_group(results_base, atlas, gwas, tier)
        df = df.rename(columns={"z_mean": "score"})
        return df
    df = load_seismic_results(results_base, atlas, gwas, tier)
    from statsmodels.stats.multitest import multipletests
    _, fdr, _, _ = multipletests(df["pvalue"], method="fdr_bh")
    df = df.rename(columns={"coefficient": "score"})
    df["fdr"] = fdr
    return df


def _compare(pangi_results_base, uc_atlas_results_dirs, method, tier, gwas,
             uc_atlas, args):
    try:
        pangi_df = _load_one(
            method, pangi_results_base, "pangi", gwas, tier,
        )
    except FileNotFoundError as exc:
        LOGGER.warning("skip Pan-GI load %s/%s/%s: %s",
                       method, gwas, tier, exc)
        return None
    scdrs_dir, seismic_dir = uc_atlas_results_dirs
    base = scdrs_dir if method == "scdrs" else seismic_dir
    try:
        uc_df = _load_one(method, base, uc_atlas, gwas, tier)
    except FileNotFoundError as exc:
        LOGGER.warning("skip UC load %s/%s/%s/%s: %s",
                       uc_atlas, method, gwas, tier, exc)
        return None

    shared = shared_cell_types(pangi_df, uc_df, min_cells=args.min_cells)
    if len(shared) < 3:
        LOGGER.warning(
            "skip %s/%s/%s/%s: %d shared cell types",
            uc_atlas, method, gwas, tier, len(shared),
        )
        return None
    pangi_df = pangi_df[pangi_df["cell_type"].isin(shared)]
    uc_df = uc_df[uc_df["cell_type"].isin(shared)]

    cr = concordance(
        scores_a=to_lookup(pangi_df, "cell_type", "score"),
        scores_b=to_lookup(uc_df, "cell_type", "score"),
        qvals_a=to_lookup(pangi_df, "cell_type", "fdr"),
        qvals_b=to_lookup(uc_df, "cell_type", "fdr"),
        cell_counts_a=dict(zip(
            pangi_df["cell_type"].astype(str),
            pangi_df["n_cells"].astype(int),
        )),
        cell_counts_b=dict(zip(
            uc_df["cell_type"].astype(str),
            uc_df["n_cells"].astype(int),
        )),
        min_cells=args.min_cells,
        n_bootstrap=args.bootstrap_n,
        seed=args.seed,
        is_fine_tier=(tier == "fine"),
    )
    return {
        "uc_atlas": uc_atlas,
        "method": method,
        "tier": tier,
        "gwas": gwas,
        "spearman_rho": cr.spearman_rho,
        "ci_lo": cr.spearman_ci_lo,
        "ci_hi": cr.spearman_ci_hi,
        "jaccard_top5": cr.jaccard_top5,
        "jaccard_top10": cr.jaccard_top10,
        "jaccard_top20": cr.jaccard_top20,
        "kappa": cr.kappa,
        "kappa_threshold": cr.kappa_threshold,
        "kappa_saturated": cr.kappa_threshold_used_due_to_saturation,
        "n_sig_pangi": cr.n_sig_a,
        "n_sig_uc": cr.n_sig_b,
        "n_common": cr.n_common,
    }


def run_donor_audit(pangi_h5ad, uc_atlas_h5ads):
    """Build a cross-atlas donor table.

    Output schema: donor_id, atlas, n_cells, present_in_pangi (bool),
    matches_smillie_pattern (bool)
    """
    import anndata as ad
    import pandas as pd

    pangi_donors: set[str] = set()
    if pangi_h5ad is not None and pangi_h5ad.exists():
        pangi = ad.read_h5ad(pangi_h5ad, backed="r")
        col = "donorID_unified" if "donorID_unified" in pangi.obs.columns else "donor"
        pangi_donors = set(pangi.obs[col].astype(str).unique())
        LOGGER.info("Pan-GI: %d unique donors", len(pangi_donors))
    else:
        LOGGER.warning("Pan-GI h5ad missing; donor audit will mark "
                       "present_in_pangi=False for all rows.")

    rows = []
    for atlas_name, atlas_path in uc_atlas_h5ads.items():
        if atlas_path is None or not Path(atlas_path).exists():
            LOGGER.warning("Skip %s donor audit: h5ad not on disk", atlas_name)
            continue
        a = ad.read_h5ad(atlas_path, backed="r")
        col = "donor" if "donor" in a.obs.columns else "donor_id"
        if col not in a.obs.columns:
            LOGGER.warning("%s: no donor column; skipping", atlas_name)
            continue
        donors = a.obs[col].astype(str)
        for d in donors.unique():
            n_cells = int((donors == d).sum())
            rows.append({
                "donor_id": d,
                "atlas": atlas_name,
                "n_cells": n_cells,
                "present_in_pangi": d in pangi_donors,
                "matches_smillie_pattern": bool(SMILLIE_DONOR_RE.match(d)),
            })
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--pangi-results-base", type=Path, required=True)
    parser.add_argument("--pangi-no-elmentaite-results-base", type=Path)
    parser.add_argument("--pangi-no-smillie-results-base", type=Path)
    parser.add_argument("--uc-atlases", nargs="+", default=list(UC_ATLASES))
    parser.add_argument("--gwas", nargs="+", default=list(UC_GWAS))
    parser.add_argument("--methods", nargs="+", default=list(METHODS),
                        choices=list(METHODS))
    parser.add_argument("--tiers", nargs="+", default=list(TIERS),
                        choices=["broad", "fine"])
    parser.add_argument("--scdrs-dir", type=Path, required=True)
    parser.add_argument("--seismic-dir", type=Path, required=True)
    parser.add_argument("--out-full", type=Path, required=True)
    parser.add_argument("--out-no-elmentaite", type=Path)
    parser.add_argument("--out-no-smillie", type=Path)
    parser.add_argument("--out-donor-audit", type=Path)
    parser.add_argument("--pangi-h5ad", type=Path,
                        help="Pan-GI h5ad path for donor audit (optional)")
    parser.add_argument("--uc-h5ad-dir", type=Path,
                        default=Path("data/atlases"),
                        help="Dir containing {atlas}.h5ad for the UC atlases")
    parser.add_argument("--min-cells", type=int, default=MIN_CELLS_PER_TYPE)
    parser.add_argument("--bootstrap-n", type=int, default=BOOTSTRAP_N)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s: %(message)s",
    )

    import pandas as pd

    t0 = time.time()
    sha = git_sha()
    uc_dirs = (args.scdrs_dir, args.seismic_dir)

    def _sweep(base):
        rows = []
        for uc_atlas in args.uc_atlases:
            for method in args.methods:
                for tier in args.tiers:
                    for gwas in args.gwas:
                        r = _compare(
                            base, uc_dirs, method, tier, gwas, uc_atlas, args,
                        )
                        if r is not None:
                            rows.append(r)
        return rows

    # Full Pan-GI
    rows_full = _sweep(args.pangi_results_base)
    if not rows_full:
        LOGGER.error("No full-Pan-GI rows produced.")
        raise SystemExit(2)
    df_full = pd.DataFrame(rows_full)
    df_full["sensitivity"] = "full"
    df_full["tool_version"] = "uc-cross-atlas-v1"
    df_full["git_sha"] = sha
    args.out_full.parent.mkdir(parents=True, exist_ok=True)
    df_full.to_csv(args.out_full, sep="\t", index=False)
    LOGGER.info("Wrote %s (%d rows)", args.out_full, len(df_full))

    if args.out_no_elmentaite and args.pangi_no_elmentaite_results_base:
        rows = _sweep(args.pangi_no_elmentaite_results_base)
        if rows:
            df = pd.DataFrame(rows)
            df["sensitivity"] = "no_elmentaite"
            df["tool_version"] = "uc-cross-atlas-v1"
            df["git_sha"] = sha
            args.out_no_elmentaite.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(args.out_no_elmentaite, sep="\t", index=False)
            LOGGER.info("Wrote %s (%d rows)", args.out_no_elmentaite, len(df))

    if args.out_no_smillie and args.pangi_no_smillie_results_base:
        rows = _sweep(args.pangi_no_smillie_results_base)
        if rows:
            df = pd.DataFrame(rows)
            df["sensitivity"] = "no_smillie"
            df["tool_version"] = "uc-cross-atlas-v1"
            df["git_sha"] = sha
            args.out_no_smillie.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(args.out_no_smillie, sep="\t", index=False)
            LOGGER.info("Wrote %s (%d rows)", args.out_no_smillie, len(df))

    if args.out_donor_audit:
        uc_h5ads = {
            atlas: args.uc_h5ad_dir / f"{atlas}.h5ad"
            for atlas in args.uc_atlases
        }
        audit = run_donor_audit(args.pangi_h5ad, uc_h5ads)
        audit["tool_version"] = "uc-cross-atlas-v1"
        audit["git_sha"] = sha
        args.out_donor_audit.parent.mkdir(parents=True, exist_ok=True)
        audit.to_csv(args.out_donor_audit, sep="\t", index=False)
        LOGGER.info("Wrote %s (%d donors)",
                    args.out_donor_audit, len(audit))

    LOGGER.info("Done in %.1fs", time.time() - t0)
    return 0


if __name__ == "__main__":
    sys.exit(main())
