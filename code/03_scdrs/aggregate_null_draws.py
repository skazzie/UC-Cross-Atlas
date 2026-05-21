"""Aggregate scDRS per-cell null-draw z-scores into per-cell-type tensors.

Brown's method (used downstream in ``code/06_concordance/`` and
``code/07_regime2_meta/``) needs the empirical correlation matrix of
cell-type-level statistics across the null distribution. The scDRS
output ships per-cell scores for ``n_ctrl=1000`` control gene sets;
this script aggregates them to per-cell-type means so that the
correlation matrix can be estimated.

Output shape: ``(n_cell_types, n_ctrl)`` per tier, written as a feather
file (and optionally also as a ``.npz``) under
``$UCC_RESULTS/null_draws/<atlas>_nulls.npz`` plus the SLURM-local
``$OUT/null_aggregations_<tier>.feather`` consumed by
``scripts/slurm/03_scdrs_compute.slurm``.

CLI matches the SLURM wrapper invocation:

    python aggregate_null_draws.py \
        --scdrs-out-folder $OUT \
        --h5ad-file $H5AD \
        --tier broad \
        --out $OUT/null_aggregations_broad.feather
"""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def _read_scores_file(path: Path) -> pd.DataFrame:
    """Read a scDRS ``*.full_score.gz`` file.

    Format: tab-separated, first column the cell index, columns include
    ``norm_score``, ``mc_pval``, ``pval``, ``zscore``, and one column per
    null draw named ``ctrl_norm_score_<i>``.
    """
    df = pd.read_csv(path, sep="\t", index_col=0, compression="gzip")
    return df


def _find_score_file(scdrs_out_folder: Path) -> Path:
    candidates = sorted(scdrs_out_folder.glob("*.full_score.gz"))
    if not candidates:
        raise FileNotFoundError(
            f"No *.full_score.gz file in {scdrs_out_folder}. "
            "Run scdrs compute-score first."
        )
    if len(candidates) > 1:
        logger.warning(
            "Multiple full_score.gz files in %s; using %s",
            scdrs_out_folder,
            candidates[0].name,
        )
    return candidates[0]


def _read_cell_type_labels(h5ad_path: Path, tier: str) -> pd.Series:
    """Read ``obs['cell_type_<tier>']`` lazily.

    Uses ``backed='r'`` so we don't load X into memory.
    """
    import anndata as ad

    obs_col = f"cell_type_{tier}"
    adata = ad.read_h5ad(h5ad_path, backed="r")
    if obs_col not in adata.obs.columns:
        raise KeyError(
            f"{h5ad_path} has no obs column '{obs_col}'. "
            f"Available: {list(adata.obs.columns)}"
        )
    labels = adata.obs[obs_col].astype("category").copy()
    return labels


def aggregate(
    scdrs_out_folder: Path,
    h5ad_path: Path,
    tier: str,
) -> pd.DataFrame:
    """Return a DataFrame of shape (n_cell_types, n_ctrl).

    Rows indexed by cell type, columns are ``ctrl_<i>``.
    """
    scores = _read_scores_file(_find_score_file(scdrs_out_folder))
    ctrl_cols = [c for c in scores.columns if c.startswith("ctrl_norm_score_")]
    if not ctrl_cols:
        raise ValueError(
            f"No ctrl_norm_score_* columns in score file at {scdrs_out_folder}. "
            "Did scdrs compute-score run with --n-ctrl > 0?"
        )
    logger.info("Found %d control draws", len(ctrl_cols))

    labels = _read_cell_type_labels(h5ad_path, tier)
    # Align: scDRS scores are indexed by cell name; subset to the cells
    # present in both indices (defensive — should match exactly).
    common = scores.index.intersection(labels.index)
    if len(common) < len(scores):
        logger.warning(
            "Cell intersection: %d in scores vs %d in h5ad; using %d common",
            len(scores),
            len(labels),
            len(common),
        )
    scores = scores.loc[common, ctrl_cols]
    labels = labels.loc[common]

    grouped = scores.groupby(labels.values, observed=True).mean()
    grouped.index.name = "cell_type"
    grouped.columns = [f"ctrl_{i}" for i in range(len(ctrl_cols))]
    return grouped


def _save_per_atlas_npz(df: pd.DataFrame, atlas: str, tier: str) -> Path | None:
    results_base = os.environ.get("UCC_RESULTS")
    if not results_base:
        return None
    out_dir = Path(results_base) / "null_draws"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{atlas}_{tier}_nulls.npz"
    np.savez_compressed(
        out_path,
        cell_types=np.asarray(df.index, dtype=object),
        nulls=df.to_numpy(dtype=np.float32),
    )
    logger.info("Wrote %s (shape %s)", out_path, df.shape)
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scdrs-out-folder", required=True, type=Path)
    parser.add_argument("--h5ad-file", required=True, type=Path)
    parser.add_argument("--tier", required=True, choices=["broad", "fine"])
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument(
        "--atlas",
        default=None,
        help="Atlas slug for the optional per-atlas npz mirror under "
        "$UCC_RESULTS/null_draws/. Defaults to the h5ad filename stem.",
    )
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")

    df = aggregate(args.scdrs_out_folder, args.h5ad_file, args.tier)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.reset_index().to_feather(args.out)
    logger.info("Wrote %s (shape %s)", args.out, df.shape)

    atlas = args.atlas or args.h5ad_file.stem
    _save_per_atlas_npz(df, atlas, args.tier)


if __name__ == "__main__":
    main()
