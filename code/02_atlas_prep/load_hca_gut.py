"""Loader for HCA Gut (Elmentaite 2021) broad atlas.

Source: CELLxGENE deposit f34d2b82-9265-4a73-bda4-852933bf2a8d.h5ad
Reference: ``code/02_atlas_prep/atlas_schemas.md``;
DECISIONS.md correction 2026-05-20 (3/7), (5/7), (6/7).
"""

from __future__ import annotations

import logging

import anndata as ad
from anndata import AnnData

from hgnc_remap import ensembl_to_hgnc

logger = logging.getLogger(__name__)

EXPECTED_OBS_COLS = (
    "Age_group",
    "tissue",
    "category",
    "author_cell_type",
    "disease",
    "donor_id",
    "sex",
    "assay",
    "batch",
    "Fraction",
)

AGE_KEEP = ("Adult", "Adult_MLN")

TISSUE_KEEP = (
    "ascending colon",
    "caecum",
    "colon",
    "descending colon",
    "large intestine",
    "rectum",
    "sigmoid colon",
    "transverse colon",
)


def load(
    h5ad_path: str,
    apply_v1_filter: bool = True,
    raw_count_mode: bool = False,
) -> AnnData:
    """Load HCA Gut, apply v1 filter (Adult + colon tissues), standardize obs."""
    if raw_count_mode:
        raise ValueError(
            "raw_count_mode=True is not supported for v1 "
            "(see DECISIONS.md correction 2026-05-20 (5/7))."
        )

    logger.info("Reading %s", h5ad_path)
    adata = ad.read_h5ad(h5ad_path)

    missing = [c for c in EXPECTED_OBS_COLS if c not in adata.obs.columns]
    if missing:
        raise KeyError(
            f"HCA Gut loader: missing expected obs columns {missing}. "
            f"Got: {list(adata.obs.columns)}"
        )

    if apply_v1_filter:
        m = (
            adata.obs["Age_group"].isin(AGE_KEEP)
            & adata.obs["tissue"].isin(TISSUE_KEEP)
        )
        n_drop = int((~m).sum())
        adata = adata[m].copy()
        logger.info(
            "HCA Gut v1 filter: dropped %d cells, kept %d", n_drop, adata.n_obs
        )

    obs = adata.obs
    obs["cell_type_broad"] = obs["category"].astype("category")
    obs["cell_type_fine"] = obs["author_cell_type"].astype("category")
    obs["donor"] = obs["donor_id"].astype("category")
    adata.obs = obs

    logger.info("Post-filter cell count: %d", adata.n_obs)
    logger.info("Per-Fraction cell counts:")
    for frac, count in adata.obs["Fraction"].value_counts().items():
        logger.info("  %s: %d cells", frac, count)

    adata = ensembl_to_hgnc(adata)
    return adata


def load_hca_gut_no_crohn(
    h5ad_path: str,
    apply_v1_filter: bool = True,
    raw_count_mode: bool = False,
) -> AnnData:
    """Sensitivity: HCA Gut with Crohn-disease cells removed.

    Tests whether residual Crohn signal in this broad reference atlas
    affects UC cell-type prioritization. See DECISIONS.md correction (6/7).
    """
    adata = load(
        h5ad_path,
        apply_v1_filter=apply_v1_filter,
        raw_count_mode=raw_count_mode,
    )
    keep = adata.obs["disease"] == "normal"
    n_drop = int((~keep).sum())
    adata = adata[keep].copy()
    logger.info(
        "load_hca_gut_no_crohn: dropped %d non-normal cells, kept %d",
        n_drop,
        adata.n_obs,
    )
    return adata
