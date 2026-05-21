"""Loader for Garrido-Trigo 2023 (Atlas 2, UC subset).

Source: CELLxGENE deposit b1a62801-f509-45f8-b55f-533fbb7e7800.h5ad
Reference: ``code/02_atlas_prep/atlas_schemas.md``;
DECISIONS.md correction 2026-05-20 (2/7), (4/7), (5/7).
"""

from __future__ import annotations

import logging

import anndata as ad
from anndata import AnnData
import pandas as pd

from hgnc_remap import ensembl_to_hgnc

logger = logging.getLogger(__name__)

EXPECTED_OBS_COLS = (
    "cell_type",
    "disease",
    "donor_id",
    "sex",
    "biospsy_or_surgical_resection_area",  # NOTE: typo "biospsy" is canonical
    "assay",
)

UC_DISEASES = ("normal", "ulcerative colitis")


def load(
    h5ad_path: str,
    apply_v1_filter: bool = True,
    raw_count_mode: bool = False,
) -> AnnData:
    """Load Garrido-Trigo, filter to UC subset, standardize obs columns.

    Parameters
    ----------
    h5ad_path
        Path to the downloaded CELLxGENE .h5ad.
    apply_v1_filter
        If True (default), keep only ``disease in {normal, ulcerative colitis}``.
    raw_count_mode
        Must remain False for v1 (DECISIONS.md correction 5/7).
    """
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
            f"Garrido-Trigo loader: missing expected obs columns {missing}. "
            f"Got: {list(adata.obs.columns)}"
        )

    if apply_v1_filter:
        keep = adata.obs["disease"].isin(UC_DISEASES)
        n_drop = int((~keep).sum())
        adata = adata[keep].copy()
        logger.info(
            "UC subset filter: dropped %d cells, kept %d", n_drop, adata.n_obs
        )

    obs = adata.obs
    obs["cell_type_broad"] = obs["cell_type"].astype("category")
    obs["cell_type_fine"] = pd.Series(pd.NA, index=obs.index, dtype="object")
    obs["donor"] = obs["donor_id"].astype("category")
    if "tissue" not in obs.columns:
        obs["tissue"] = "colonic mucosa"
    obs["batch"] = pd.Series(pd.NA, index=obs.index, dtype="object")
    adata.obs = obs

    logger.info("Post-filter cell count: %d", adata.n_obs)
    per_donor = adata.obs["donor"].value_counts()
    logger.info("Donors (n=%d):", per_donor.size)
    for donor, count in per_donor.items():
        logger.info("  %s: %d cells", donor, count)

    adata = ensembl_to_hgnc(adata)
    return adata
