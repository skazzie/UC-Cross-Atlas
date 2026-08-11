#!/usr/bin/env python3
"""Cross-method concordance: scDRS vs seismic on real (atlas, gwas, tier) pairs.

For each pair present on disk:
  - Read scDRS group file: uses `group` as cell type, `assoc_mcz` as the
    signed cell-type-level Z-score, `assoc_mcp` as the p-value.
  - Read seismic tsv (1.0.0 API): columns cell_type, pvalue, FDR, n_cells.
  - Inner-join on cell_type; compute Spearman(assoc_mcz, -log10(pvalue)).

Interpretation caveat
    scDRS `assoc_mcz` is SIGNED (direction of enrichment). Seismic
    -log10(pvalue) is UNSIGNED (1.0.0 strips coefficient/direction).
    The correlation therefore measures enrichment-rank concordance
    against an unsigned seismic axis: a cell type that scDRS scores as
    depleted (large-negative z) but that seismic flags as significant
    (large -log10 p) drags the correlation down even though both
    methods flagged it.

MHC pairing
    Seismic's gene-Z is MHC-EXCLUDED (see code/04_seismic manifest);
    only MHC-excluded scDRS runs are paired against seismic. MHC-included
    scDRS runs are intentionally omitted — the two methods would be
    operating on different gene sets.

Pair construction
    Pairs are built from UC_ATLASES x UC_GWAS x TIERS (see
    _shared.constants), so adding an atlas to UC_ATLASES automatically
    picks it up here. Per-(atlas, gwas) scDRS group-dir naming
    conventions differ, so SCDRS_GROUP_DIRS below is an explicit
    atlas -> gwas -> dirname map. Seismic path is the locked convention
    results/seismic/{atlas}_{gwas}_{tier}.tsv.

    Every missing file logs an explicit skip so a silently-absent atlas
    is impossible — the failure mode that previously hid TAURUS.

Output columns
    atlas, gwas, tier, n_cell_types_matched, spearman_rho, spearman_pval

Usage
    python code/08_cross_method/run_scdrs_seismic_concordance.py \\
        --out results/concordance/scdrs_vs_seismic.csv
"""

from __future__ import annotations

import argparse
import logging
import sys
from itertools import product
from pathlib import Path
from typing import NamedTuple

import numpy as np
import pandas as pd

LOGGER = logging.getLogger(Path(__file__).stem)

_REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO / "code"))
sys.path.insert(0, str(_REPO / "code" / "06_concordance"))
from _shared.constants import TIERS, UC_ATLASES, UC_GWAS  # noqa: E402
from metrics import spearman  # noqa: E402


# atlas -> gwas -> scDRS group dir under results/scdrs/.
# All entries reference MHC-EXCLUDED scDRS runs (matches seismic's MHC-excluded
# gene-Z). Conventions genuinely differ per (atlas, gwas) — see results/MANIFEST.md
# for provenance:
#   garrido_trigo x delange: smoketest run (VM 2026-07-07, depth+donor cov)
#   garrido_trigo x liu:     non-smoketest run
#   smillie:                 fullcov MHC-excluded (locked-primary per DECISIONS.md)
#   taurus:                  MHC-excluded; delange only exists today, liu will
#                            skip cleanly until the run lands.
# Entries whose files don't exist yet skip via the file-existence check in
# _compute — they must still appear here or the (atlas, gwas) is silently dropped.
SCDRS_GROUP_DIRS: dict[str, dict[str, str]] = {
    "garrido_trigo": {
        "delange": "garrido_delange_smoketest_group",
        "liu": "garrido_liu_group",
    },
    "smillie": {
        "delange": "smillie_delange_excl_fullcov_group",
        "liu": "smillie_liu_excl_fullcov_group",
    },
    "taurus": {
        "delange": "taurus_delange_excl_group",
        "liu": "taurus_liu_excl_group",
    },
}

SCDRS_GROUP_FILE_TEMPLATE = "UC.scdrs_group.cell_type_{tier}"
SEISMIC_PATH_TEMPLATE = "results/seismic/{atlas}_{gwas}_{tier}.tsv"


class Pair(NamedTuple):
    atlas: str
    gwas: str
    tier: str
    scdrs: str  # repo-relative path
    seismic: str  # repo-relative path


def build_manifest() -> list[Pair]:
    """Cartesian product of UC_ATLASES x UC_GWAS x TIERS, resolved to paths.

    An (atlas, gwas) with no SCDRS_GROUP_DIRS entry is skipped with a warning
    so an atlas added to UC_ATLASES but not wired here is loud, not silent.
    """
    pairs: list[Pair] = []
    for atlas, gwas, tier in product(UC_ATLASES, UC_GWAS, TIERS):
        group_dir = SCDRS_GROUP_DIRS.get(atlas, {}).get(gwas)
        if group_dir is None:
            LOGGER.warning(
                "skip %s/%s/%s: no SCDRS_GROUP_DIRS entry for (%s, %s) — "
                "add to the map in this file",
                atlas, gwas, tier, atlas, gwas,
            )
            continue
        scdrs = f"results/scdrs/{group_dir}/{SCDRS_GROUP_FILE_TEMPLATE.format(tier=tier)}"
        seismic = SEISMIC_PATH_TEMPLATE.format(atlas=atlas, gwas=gwas, tier=tier)
        pairs.append(Pair(atlas, gwas, tier, scdrs, seismic))
    return pairs


def _load_scdrs(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    rename = {"group": "cell_type", "assoc_mcz": "z", "assoc_mcp": "pvalue"}
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
    for col in ("cell_type", "z"):
        if col not in df.columns:
            raise KeyError(f"{path}: missing required column '{col}'. "
                           f"Got: {list(df.columns)}")
    return df[["cell_type", "z"]].copy()


def _load_seismic(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    for col in ("cell_type", "pvalue"):
        if col not in df.columns:
            raise KeyError(f"{path}: missing required column '{col}'. "
                           f"Got: {list(df.columns)}")
    return df[["cell_type", "pvalue"]].copy()


def _compute(pair: Pair, repo: Path) -> dict | None:
    sd_path = repo / pair.scdrs
    ss_path = repo / pair.seismic
    missing_sd = not sd_path.exists()
    missing_ss = not ss_path.exists()
    if missing_sd or missing_ss:
        if missing_sd:
            LOGGER.warning("skip %s/%s/%s: scDRS missing %s",
                           pair.atlas, pair.gwas, pair.tier, sd_path)
        if missing_ss:
            LOGGER.warning("skip %s/%s/%s: seismic missing %s",
                           pair.atlas, pair.gwas, pair.tier, ss_path)
        return None

    sd = _load_scdrs(sd_path)
    ss = _load_seismic(ss_path)

    # cell_type is broad/fine label in both files (colonocyte, goblet, ...).
    sd["cell_type"] = sd["cell_type"].astype(str)
    ss["cell_type"] = ss["cell_type"].astype(str)
    merged = sd.merge(ss, on="cell_type", how="inner").dropna(subset=["z", "pvalue"])
    n_matched = len(merged)
    LOGGER.info("%s/%s/%s: %d matched cell types (scDRS n=%d, seismic n=%d)",
                pair.atlas, pair.gwas, pair.tier, n_matched, len(sd), len(ss))

    if n_matched < 3:
        return {"atlas": pair.atlas, "gwas": pair.gwas, "tier": pair.tier,
                "n_cell_types_matched": n_matched,
                "spearman_rho": float("nan"), "spearman_pval": float("nan")}

    # Clip to avoid -log10(0) = inf on p==0 (shouldn't happen but be safe).
    neg_log10_p = -np.log10(np.clip(merged["pvalue"].to_numpy(dtype=float),
                                    1e-300, 1.0))
    scores_a = dict(zip(merged["cell_type"], merged["z"].astype(float)))
    scores_b = dict(zip(merged["cell_type"], neg_log10_p))

    # No min-cell filter here (both counts=None -> filter disabled); the
    # inner-join is already the shared cell-type set.
    rho, pval, n_kept, _ = spearman(scores_a, scores_b)
    return {"atlas": pair.atlas, "gwas": pair.gwas, "tier": pair.tier,
            "n_cell_types_matched": n_kept,
            "spearman_rho": rho, "spearman_pval": pval}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--repo-dir", type=Path, default=_REPO,
                        help="Repo root (default: auto-detect).")
    parser.add_argument("--out", type=Path,
                        default=_REPO / "results/concordance/scdrs_vs_seismic.csv",
                        help="Output CSV.")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s: %(message)s",
    )

    manifest = build_manifest()
    LOGGER.info("Attempting %d (atlas, gwas, tier) pairs from UC_ATLASES x UC_GWAS x TIERS",
                len(manifest))

    rows = []
    for pair in manifest:
        row = _compute(pair, args.repo_dir)
        if row is not None:
            rows.append(row)

    if not rows:
        LOGGER.error("No (atlas, gwas, tier) pairs found on disk.")
        return 2

    out = pd.DataFrame(rows, columns=[
        "atlas", "gwas", "tier",
        "n_cell_types_matched", "spearman_rho", "spearman_pval",
    ])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    LOGGER.info("Wrote %s (%d rows of %d attempted)",
                args.out, len(out), len(manifest))
    print(out.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
