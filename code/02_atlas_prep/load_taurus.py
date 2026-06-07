"""Loader for TAURUS-IBD (Atlas 3, UC core; replaces Mennillo per DECISIONS 16).

Source: Zenodo deposit by Thomas, Dendrou, Agarwal (Oxford), 2024.
  "A longitudinal single-cell atlas of anti-tumour necrosis factor
   treatment in inflammatory bowel disease"
  - DOI cited in DECISIONS swap directive: 10.5281/zenodo.13768607
  - DOI resolved on Zenodo 2026-06-06:        10.5281/zenodo.14007626
  - These may be related versioned records; the version to use must be
    pinned before loader implementation. Flagged in DECISIONS 16.

Files (Zenodo, current resolved record):
- ``TAURUS_raw_counts_annotated_final.h5ad`` (12.7 GB) — pooled, all
  lineages, ANNOTATED — the production input.
- Per-lineage h5ads (epithelial, CD4/CD8 T, B, plasma, myeloid, ILC,
  fibroblasts, vascular) for targeted analyses.
- Paired raw BAM/CellRanger outputs at GEO **GSE282122**.

Scope alignment with v1
-----------------------
TAURUS is an **IBD** atlas (UC + CD donors) sampled longitudinally
during anti-TNF treatment. For v1 cross-atlas UC concordance, the
loader MUST subset to:

- UC donors only (drop CD).
- A single biopsy time-point per donor — likely **pre-treatment
  baseline**, mirroring the Mennillo plan. Confirm the metadata field
  that distinguishes baseline vs on-treatment, and verify the per-donor
  cell count after subsetting (≥50 cells per fine cluster per donor).

Open questions before implementation (DECISIONS 16):
- Which Zenodo DOI version is canonical (13768607 vs 14007626)?
- What is the obs schema? (Thomas 2024 supplement lists annotation
  columns; capture at download time, same discipline as Smillie / Garrido.)
- Counts: ``raw_counts`` in the filename suggests raw integer counts ->
  apply ``log1p(CP10k)`` on load (DECISIONS 5/7), preserve raw in
  ``layers['counts']``.
- Cell-type vocab: TAURUS annotations differ from Salas (Garrido) and
  Broad (Smillie); a TAURUS-specific FINE_TO_BROAD will need writing,
  and the canonical broad set in `canonical_broad_DRAFT.md` is the
  target.

References: DECISIONS.md correction (16) for the Mennillo→TAURUS swap;
(5/7) for the normalization requirement; (10), (12) for the
loader-discipline pattern that this loader will follow.
"""

from __future__ import annotations

from anndata import AnnData


def load(
    h5ad_path: str,
    apply_v1_filter: bool = True,
    raw_count_mode: bool = False,
) -> AnnData:
    """Not yet implemented; Zenodo download required first.

    TODO once TAURUS h5ad is on disk:
      - Confirm Zenodo DOI version (13768607 vs 14007626 — see module
        docstring and DECISIONS 16).
      - Capture obs schema (donor, disease subtype, treatment timepoint,
        tissue, batch, fine cell-type column).
      - Subset to UC donors only; drop CD.
      - Subset to a single time-point per donor (pre-treatment baseline
        unless metadata argues otherwise).
      - Map TAURUS fine labels into the canonical broad vocab; verify
        the two-gate canonical-vocab assertion pattern from
        ``load_garrido_trigo.py`` and ``load_smillie.py``.
      - Apply ``log1p(CP10k)`` (raw counts -> log-normalized) to match
        the other atlases (DECISIONS 5/7); preserve raw in
        ``layers['counts']``.
      - Donor-structure invariant: confirm ≥8 UC donors survive the
        subset filter.
    """
    raise NotImplementedError(
        "TAURUS loader is deferred. Download from Zenodo "
        "(10.5281/zenodo.14007626 — confirm version vs 13768607 first) "
        "and capture the obs schema before this loader can be completed. "
        "See DECISIONS.md correction (16)."
    )
