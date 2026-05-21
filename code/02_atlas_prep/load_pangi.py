"""Loader for Pan-GI Extended+ broad atlas (Oliver 2024).

Source: CELLxGENE deposit 1dcf15ee-c103-4aaa-8b8c-0fc697fcccc8.h5ad
Reference: ``code/02_atlas_prep/atlas_schemas.md``;
DECISIONS.md correction 2026-05-20 (3/7), (5/7).
"""

from __future__ import annotations

import logging

import anndata as ad
from anndata import AnnData
import pandas as pd

from hgnc_remap import ensembl_to_hgnc

logger = logging.getLogger(__name__)

EXPECTED_OBS_COLS = (
    "disease",
    "organ_unified",
    "sample_type",
    "level_2_annot",
    "level_3_annot",
    "study",
    "donorID_unified",
    "sex",
    "assay",
    "cell_type",
)

DISEASE_KEEP = (
    "normal",
    "ulcerative colitis",
    "inflammatory bowel disease",
)

ORGAN_KEEP = (
    "ascending colon",
    "caecum",
    "colon",
    "descending colon",
    "rectum",
    "sigmoid colon",
    "transverse colon",
)

EXCLUDE_SAMPLE_TYPE = "Organ_donor_resection"

# Conservative pattern for Smillie 2019 donor IDs as documented in
# DECISIONS.md correction (3/7) point 4. Smillie's naming convention is
# ``N1..N12`` (healthy) and ``UC1..UC68`` (UC). The loader scans for any
# donorID_unified that exactly matches these patterns.
SMILLIE_DONOR_REGEX = r"^(N|UC)\d+$"


def _smillie_overlap(obs: pd.DataFrame) -> pd.Series:
    return (
        obs["donorID_unified"].astype(str).str.match(SMILLIE_DONOR_REGEX, na=False)
    )


def load(
    h5ad_path: str,
    apply_v1_filter: bool = True,
    raw_count_mode: bool = False,
) -> AnnData:
    """Load Pan-GI Extended+, apply v1 filter, standardize obs."""
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
            f"Pan-GI loader: missing expected obs columns {missing}. "
            f"Got: {list(adata.obs.columns)}"
        )

    if apply_v1_filter:
        m = (
            adata.obs["disease"].isin(DISEASE_KEEP)
            & adata.obs["organ_unified"].isin(ORGAN_KEEP)
            & (adata.obs["sample_type"] != EXCLUDE_SAMPLE_TYPE)
        )
        n_drop = int((~m).sum())
        adata = adata[m].copy()
        logger.info(
            "Pan-GI v1 filter: dropped %d cells, kept %d", n_drop, adata.n_obs
        )

    obs = adata.obs
    obs["cell_type_broad"] = obs["level_2_annot"].astype("category")
    obs["cell_type_fine"] = obs["level_3_annot"].astype("category")
    obs["donor"] = obs["donorID_unified"].astype("category")
    obs["tissue"] = obs["organ_unified"].astype("category")
    obs["batch"] = pd.Series(pd.NA, index=obs.index, dtype="object")
    adata.obs = obs

    logger.info("Post-filter cell count: %d", adata.n_obs)

    elm_count = int((adata.obs["study"] == "Elmentaite2021").sum())
    logger.info("Elmentaite2021 cells (HCA Gut overlap): %d", elm_count)

    smillie_mask = _smillie_overlap(adata.obs)
    smillie_count = int(smillie_mask.sum())
    if smillie_count == 0:
        logger.info("Smillie donor-ID scan: 0 matching cells (expected).")
    else:
        unique = (
            adata.obs.loc[smillie_mask, "donorID_unified"]
            .astype(str)
            .unique()
            .tolist()
        )
        logger.warning(
            "Smillie donor-ID scan: %d matching cells across donors %s. "
            "Produce removal sensitivity.",
            smillie_count,
            unique,
        )

    adata = ensembl_to_hgnc(adata)
    return adata


def load_pangi_no_elmentaite(
    h5ad_path: str,
    apply_v1_filter: bool = True,
    raw_count_mode: bool = False,
) -> AnnData:
    """Sensitivity: Pan-GI with Elmentaite2021 (HCA Gut) cells removed.

    Tests whether HCA Gut donor overlap drives Pan-GI cell-type
    prioritization. See DECISIONS.md correction (3/7).
    """
    adata = load(
        h5ad_path,
        apply_v1_filter=apply_v1_filter,
        raw_count_mode=raw_count_mode,
    )
    keep = adata.obs["study"] != "Elmentaite2021"
    n_drop = int((~keep).sum())
    adata = adata[keep].copy()
    logger.info(
        "load_pangi_no_elmentaite: dropped %d Elmentaite2021 cells, kept %d",
        n_drop,
        adata.n_obs,
    )
    return adata


def load_pangi_no_smillie(
    h5ad_path: str,
    apply_v1_filter: bool = True,
    raw_count_mode: bool = False,
) -> AnnData:
    """Sensitivity: Pan-GI with any Smillie-pattern donor IDs removed.

    Expected to be a no-op based on initial inspection (no overlap found);
    documents the empirical overlap rather than asserting absence.
    """
    adata = load(
        h5ad_path,
        apply_v1_filter=apply_v1_filter,
        raw_count_mode=raw_count_mode,
    )
    mask = _smillie_overlap(adata.obs)
    n_drop = int(mask.sum())
    if n_drop == 0:
        logger.info(
            "load_pangi_no_smillie: 0 cells matched Smillie donor pattern; no-op."
        )
        return adata
    adata = adata[~mask].copy()
    logger.info(
        "load_pangi_no_smillie: dropped %d Smillie-pattern cells, kept %d",
        n_drop,
        adata.n_obs,
    )
    return adata
