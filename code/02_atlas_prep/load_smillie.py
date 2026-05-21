"""Loader for Smillie 2019 UC core atlas (SKELETON — deferred to next session).

The CELLxGENE deposit (e373cf41-...) is unusable: 34,772 cells, all
healthy-epithelial. Canonical source is Single Cell Portal SCP259
(366,650 cells, 30 donors, 18 UC + 12 HC, all compartments).

Access requires SCP account creation, email verification, and a
browser-mediated consent click. See DECISIONS.md correction 2026-05-20 (7/7).

This loader will be filled in once SCP259 is on disk.
"""

from __future__ import annotations

from anndata import AnnData


def load(
    h5ad_path: str,
    apply_v1_filter: bool = True,
    raw_count_mode: bool = False,
) -> AnnData:
    """Not yet implemented; SCP259 download required first.

    TODO once SCP259 is on disk:
      - Download URL (after SCP authentication)
      - Discover obs schema (donor, disease, tissue, compartment columns)
      - Filter chain (likely: keep all 30 donors; the whole cohort is v1
        relevant)
      - Identify broad and fine tier columns
      - Add log1p(CP10k) normalization if the deposit ships raw counts
        (per Correction 5/7 cross-atlas-comparable input)
    """
    raise NotImplementedError(
        "Smillie loader is deferred. SCP259 must be downloaded from "
        "https://singlecell.broadinstitute.org/single_cell/study/SCP259 "
        "(requires account + browser consent) before this loader can be "
        "completed. See DECISIONS.md correction 2026-05-20 (7/7)."
    )
