"""Loader for TAURUS-IBD (Atlas 3, UC core; replaces Mennillo per DECISIONS 16).

**Source pin** (DECISIONS 16):

- Zenodo v3: `10.5281/zenodo.14007626`
- Pooled file: ``TAURUS_raw_counts_annotated_final.h5ad`` (12.7 GB)
- md5: ``c1bd13b92cacb164a401c6c4a4e7912c``
- Paper: Thomas et al., *Nat Immunol* 25:2152-2165 (2024),
  doi:10.1038/s41590-024-01994-8, PMC11519010.

**Cohort**: 41 subjects total (16 CD + 22 UC + 3 HC), ~1 million cells,
109 cell states. v1 subset is **UC × colonic region × pretreatment
baseline**, expected to yield 22 UC donors and ~50 inflamed baseline
samples per the paper's Fig. 2b. CD and HC arms are dropped at the
filter step.

**Annotation hierarchy** (verified against actual h5ad obs, 2026-08-09):
TAURUS ships only TWO cell-type tiers — ``major`` (broad; 8 classes
[Plasma, Endothelium, CD4_T, CD8_T, Non_ileal_epithelium, Pericyte, B,
Mono_macro]) and ``minor`` (fine; ~109 classes). There is NO
``compartment`` / ``low`` / ``intermediate`` sub-hierarchy like the
Garrido or Smillie atlases. ``minor`` is stored as
``obs['cell_type_fine']``; ``major`` is mapped to the canonical
broad vocab via ``MAJOR_TO_BROAD`` and stored as
``obs['cell_type_broad']``. The compartment / low / intermediate obs
columns are intentionally omitted for TAURUS.

**Filter chain (v1)** — three stages, NOT four:

1. Disease == UC (drop CD + HC).
2. Region: ``Ileum_vs_Colon`` in {Colon, Rectum} (drop Ileum). The
   finer ``Site`` column (Ascending_Colon / Descending_Colon / Sigmoid /
   Rectum / Terminal_Ileum) is preserved for cohort validation and
   biopsy-level covariate structure.
3. Timepoint == baseline / W0 / pretreatment (drop post-treatment).

The Zenodo deposit description suggests ``inflammation_score > 6.5``
as the cutoff for the paper's inflamed-baseline remission analysis;
this loader **does not apply that filter** because doing so would
(a) drop ~half of UC baseline samples (Fig. 2b: 50 inflamed vs 53
non-inflamed); (b) pre-empt OPEN_FLAGS F1, which must set one
inflamed/non-inflamed/pooled policy uniformly across all three UC
atlases (Smillie and Garrido are not inflamed-only); (c) skew the
cell-type composition toward inflammation-expanded subsets and
distort GWAS prioritization. ``inflammation_score`` is preserved as
obs metadata so F1 can apply its policy downstream of all three
loaders. See DECISIONS correction (18)(a).

**Validation gates** (mirror correction 9 / 12 / 16 patterns):

- Donor-structure hard invariant: 22 UC donors (``Patient`` column)
  after the full filter. Per-donor ``Site`` sets checked against Supp
  Table 1B (EXPECTED_UC_COHORT).
- Cell-count tripwire (soft): expect ~30-50k cells but TBD on first run;
  log only.
- Canonical-vocab assertions (gates 1 + 2): same two-gate pattern as
  the Garrido + Smillie loaders. ``MAJOR_TO_BROAD`` ships **empty** in
  this v0 — gate 2 will fail loud on first run listing the 8 major
  labels [Plasma, Endothelium, CD4_T, CD8_T, Non_ileal_epithelium,
  Pericyte, B, Mono_macro] so their canonical-broad targets can be
  populated (one commit per Muskaan biology pass).
- Counts pipeline: ``log1p(CP10k)`` on load per DECISIONS (5/7); raw
  integer counts preserved in ``layers['counts']``;
  ``raw_count_mode=True`` unsupported.

**Open before first run** (DECISIONS 16 + load gate will catch):

- The ``major``-tier label set → ``MAJOR_TO_BROAD`` map (8 labels
  from schema inspection 2026-08-09; canonical-broad targets TBD).

References: DECISIONS.md (16) [Mennillo→TAURUS swap]; (5/7)
[normalization]; (11) [HGNC pin]; (9) / (12) [loader discipline].
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import anndata as ad
from anndata import AnnData
import numpy as np
import pandas as pd
import scanpy as sc

from hgnc_remap import ensembl_to_hgnc

logger = logging.getLogger(__name__)

# ---- Source pin (DECISIONS 16) ------------------------------------------
ZENODO_DOI = "10.5281/zenodo.14007626"  # v3
POOLED_FILENAME = "TAURUS_raw_counts_annotated_final.h5ad"
POOLED_MD5 = "c1bd13b92cacb164a401c6c4a4e7912c"

# ---- Filter constants ---------------------------------------------------
# The Zenodo deposit description suggests ``inflammation_score > 6.5`` as
# the cutoff for the paper's inflamed-baseline analysis. This loader
# deliberately does NOT apply that filter — see DECISIONS (18)(a) and the
# module docstring. Constant retained as documentation for any future
# inflamed-only sensitivity that goes through the F1 cross-atlas
# inflammation policy.
PAPER_BASELINE_INFLAMMATION_MIN: float = 6.5  # NOT applied; see (18)(a)

# Disease values to keep / drop after canonicalization (UC + HC + CD all
# get canonicalized to short strings by _canonicalize_disease).
KEEP_DISEASE = ("UC",)  # v1 strict — drop CD and HC per user spec

# Region filter operates on the coarse `Ileum_vs_Colon` obs column,
# whose value set is exactly {Rectum, Colon, Ileum}. Keep Colon+Rectum,
# drop Ileum. The finer `Site` column (Ascending_Colon / Descending_Colon
# / Sigmoid / Rectum / Terminal_Ileum) is preserved separately for
# cohort validation and covariate structure — see load(): scol.
COLONIC_REGION_VALUES: frozenset[str] = frozenset({"colon", "rectum"})

# Timepoint values that count as baseline. Matched by EXACT equality
# after lower-case + whitespace-normalize. **Includes the bare token
# `"pre"`** — TAURUS Supp Table 1B uses "Pre"/"Post" verbatim, and a
# substring-style match in the wrong direction would have silently
# dropped TAURUS's "Pre" rows on the first run (bug caught 2026-06-06
# via Supp Table 1 dry-run; see DECISIONS 20).
BASELINE_TIMEPOINT_KEYS: frozenset[str] = frozenset({
    "pre", "baseline", "pretreatment", "pre-treatment", "pre_treatment",
    "w0", "week 0", "week_0", "wk0", "v1", "visit 1", "visit_1",
})

# ---- Donor + sample expectations (Supp Table 1B, dry-run 2026-06-06) ----
# Derived empirically by filtering Supplementary Table 1B
# (Disease == "UC" & Site in colonic & Treatment == "Pre"). All 22 UC
# patients contributed at least one Pre colonic sample (no patient is
# excluded by the v1 filter chain). Fig. 4c's smaller cohort (4 + 13 =
# 17) results from the paper's *additional* inflamed-baseline filter
# which (18)(a) deliberately does NOT apply.
EXPECTED_N_UC_DONORS_POST_FILTER: int = 22
# 52 = 39 Inflamed + 13 Non_Inflamed UC × colonic × Pre samples.
EXPECTED_N_UC_SAMPLES_POST_FILTER: int = 52

# Per-donor expected (region-set, sample-count) tuples from Supp Table
# 1B (v3 has metadata corrections vs prior versions, so a mismatch here
# is meaningful — see DECISIONS 21). Used in the donor-invariant
# failure path to dump a structured diff rather than "expected 22, got
# 21". Note: TAURUS Site values use underscores ("Ascending_Colon");
# the loader passes obs[Site] through unchanged into the comparison.
EXPECTED_UC_COHORT: dict[str, frozenset[str]] = {
    "UC1":  frozenset({"Ascending_Colon", "Rectum"}),
    "UC2":  frozenset({"Descending_Colon", "Rectum", "Sigmoid"}),
    "UC3":  frozenset({"Ascending_Colon", "Descending_Colon", "Rectum"}),
    "UC4":  frozenset({"Ascending_Colon", "Descending_Colon"}),
    "UC5":  frozenset({"Descending_Colon", "Rectum", "Sigmoid"}),
    "UC6":  frozenset({"Rectum", "Sigmoid"}),
    "UC7":  frozenset({"Descending_Colon", "Rectum", "Sigmoid"}),
    "UC8":  frozenset({"Rectum", "Sigmoid"}),
    "UC9":  frozenset({"Descending_Colon", "Rectum", "Sigmoid"}),
    "UC10": frozenset({"Descending_Colon", "Rectum", "Sigmoid"}),
    "UC11": frozenset({"Ascending_Colon", "Rectum"}),
    "UC12": frozenset({"Descending_Colon", "Rectum", "Sigmoid"}),
    "UC13": frozenset({"Descending_Colon", "Rectum", "Sigmoid"}),
    "UC14": frozenset({"Rectum", "Sigmoid"}),
    "UC15": frozenset({"Descending_Colon", "Rectum"}),
    "UC16": frozenset({"Ascending_Colon", "Descending_Colon"}),
    "UC17": frozenset({"Descending_Colon", "Rectum"}),
    "UC18": frozenset({"Rectum"}),
    "UC19": frozenset({"Descending_Colon", "Rectum"}),
    "UC20": frozenset({"Descending_Colon", "Rectum"}),
    "UC21": frozenset({"Descending_Colon", "Rectum"}),
    "UC22": frozenset({"Ascending_Colon", "Descending_Colon", "Rectum"}),
}

# Inflammation breakdown of the v1 cohort (informational; F1 governs
# whether the subset gets further split before concordance):
#   n(inflammation_score > 6.5)  ≈ 31
#   n(inflammation_score ≤ 6.5)  ≈ 21
#   inflammation_score range     ≈ [0.00, 8.98]

# Cell count post-filter is not pinned anywhere in the paper Methods —
# left as a soft tripwire on first run.
EXPECTED_N_CELLS_HINT: int | None = None  # set after first-run capture

# ---- Schema auto-detect candidates --------------------------------------
# Candidate name lists are TAURUS-Supp-Table-1B-confirmed first, then
# fallbacks. Auto-detect returns the first that exists in obs.columns.
_DONOR_COL_CANDIDATES = (
    "Patient",  # TAURUS Supp Table 1B
    "patient_id", "Patient_ID", "donor_id", "Donor_ID",
    "patient", "subject", "Subject", "subject_id",
)
_DISEASE_COL_CANDIDATES = (
    "Disease",  # TAURUS Supp Table 1B
    "diagnosis", "Diagnosis", "disease",
    "condition", "Condition", "disease_status", "group",
)
_REGION_FILTER_COL_CANDIDATES = (
    "Ileum_vs_Colon",  # TAURUS: coarse [Rectum, Colon, Ileum] — used for filter
    "ileum_vs_colon", "region_broad",
)
_SITE_COL_CANDIDATES = (
    "Site",  # TAURUS Supp Table 1B — fine [Ascending_Colon, Descending_Colon,
             # Rectum, Sigmoid, Terminal_Ileum]. Used for cohort + covariate.
    "site", "biopsy_site", "region", "Region", "tissue", "Tissue",
    "anatomical_region", "location", "Location", "tissue_region",
)
_TIMEPOINT_COL_CANDIDATES = (
    "Treatment",  # TAURUS Supp Table 1B ("Pre" / "Post")
    "timepoint", "Timepoint", "time_point", "visit", "Visit",
    "week", "Week", "treatment_status", "treatment_timepoint",
    "sample_timepoint", "visit_id",
)
_INFLAMMATION_COL_CANDIDATES = (
    "Inflammation_score",  # TAURUS Supp Table 1B (capital I; numeric)
    "inflammation_score", "Inflammation_Score", "inflammation",
    "Inflammation", "inflammation_grade", "infl_score",
    "macroscopic_inflammation", "endoscopic_inflammation",
)
# TAURUS has only 2 cell-type tiers — `major` (broad; 8 classes) and
# `minor` (fine; ~109 classes). No compartment / low / intermediate.
_MAJOR_COL_CANDIDATES = (
    "major", "Major", "cell_type_major", "broad", "cell_type_broad",
)
_MINOR_COL_CANDIDATES = (
    "minor", "Minor", "cell_type_minor", "fine", "cell_type_fine",
    "cell_state", "Cell_State", "annotation", "Annotation",
)

# ---- Canonical broad vocab (single-sourced; same 15 in every loader) ----
# Single-sourced from sibling _broad_vocab module. Loader-local copies
# were a drift risk that 06_concordance would have reported as biology.
# Gate (1) at module load + gate (2) at end-of-load. See DECISIONS (20).
from _broad_vocab import _BROAD_VOCAB
from _qc_policy import EXCLUDE_LINEAGE_AMBIGUOUS_FINE  # QC_STATE_TO_PARENT keys on cell_type_low which TAURUS lacks

# Map from TAURUS `major`-tier labels into the canonical broad vocab
# (single-sourced from _broad_vocab._BROAD_VOCAB — same 15 terms as
# Garrido + Smillie so 06_concordance's string-intersection on
# cell_type_broad aligns across atlases).
#
# 21 major labels enumerated by gate (2) first-run 2026-08-09.
#
# Scope notes:
# - TAURUS has NO "goblet" / "enteroendocrine/tuft" / "epithelial
#   progenitor" bucket at the major tier — those subsets live inside
#   `Non_ileal_epithelium` and only surface at the minor tier. The
#   broad-tier map collapses Non_ileal_epithelium → "colonocyte" per
#   user directive; concordance for goblet/EE/progenitor at broad tier
#   is therefore NOT estimable from TAURUS (fine tier still works).
# - TAURUS has NO "granulocyte" bucket (droplet lysis — structural
#   zero, matches Smillie).
# - All T-lineage majors collapse to "T cell"; `Innate_lymphocytes`
#   goes to "NK/ILC" (the innate-lymphoid bucket, not T cell).
# - All fibroblast/stromal-flavoured majors collapse to "fibroblast";
#   `Pericyte` and `Glial` go to "mural/glia" (the joint mural+glia
#   bucket used by Garrido / Smillie).
# - `Cycling_MNP` goes to "monocyte/macrophage" (MNP is dominantly the
#   macro/mono branch in gut mucosa; cycling DCs would be a small
#   share and don't warrant a separate branch here).
MAJOR_TO_BROAD: dict[str, str] = {
    # B / plasma lineage
    "B":                                 "B cell",
    "Plasma":                            "plasma cell",
    # C3+_ZBTB-blast label was pasted as "C3pos_zbtb...blast" (user
    # abbreviated); listing plausible full spellings so first-run has
    # a good chance of hitting the exact obs value. All map to plasma
    # cell (ZBTB32+ C3+ blasts are pre-plasmablast lineage in gut).
    # If gate (2) fires again, replace with the exact string it prints.
    "C3pos_ZBTB32pos_preplasmablast":    "plasma cell",
    "C3pos_zbtb32pos_preplasmablast":    "plasma cell",
    "C3pos_ZBTB32pos_plasmablast":       "plasma cell",
    "C3pos_zbtb32pos_plasmablast":       "plasma cell",
    "C3pos_ZBTB32pos_Bblast":            "plasma cell",
    "C3pos_zbtb32pos_Bblast":            "plasma cell",
    # T / NK / ILC — all conventional T collapse to "T cell";
    # innate lymphocytes go to the NK/ILC bucket.
    "CD4_T":                             "T cell",
    "CD8_T":                             "T cell",
    "Unconventional_T":                  "T cell",
    "Innate_lymphocytes":                "NK/ILC",
    # Myeloid
    "Mono_macro":                        "monocyte/macrophage",
    "Cycling_MNP":                       "monocyte/macrophage",
    "DC":                                "dendritic cell",
    "Mast":                              "mast cell",
    # Epithelial (see scope note re. Non_ileal_epithelium collapse)
    "Non_ileal_epithelium":              "colonocyte",
    # Vasculature + stroma
    "Endothelium":                       "endothelium",
    "Fibroblast":                        "fibroblast",
    "Epi_fibroblast":                    "fibroblast",
    "LP_fibroblast":                     "fibroblast",
    "Myofibroblast":                     "fibroblast",
    "THY1pos_FAPpos_PDPNpos_fibroblast": "fibroblast",
    "Cycling_stroma":                    "fibroblast",
    "Pericyte":                          "mural/glia",
    "Glial":                             "mural/glia",
}

# Gate (1): every value the map ships must be in the canonical vocab.
# Vacuously true while the map is empty; protects the day we start
# filling it in.
_unmapped_broad = set(MAJOR_TO_BROAD.values()) - _BROAD_VOCAB
if _unmapped_broad:
    raise ValueError(
        f"load_taurus.MAJOR_TO_BROAD ships broad values outside "
        f"_BROAD_VOCAB: {sorted(_unmapped_broad)}. Typo on the value "
        f"side of the map; see canonical_broad_DRAFT.md."
    )
del _unmapped_broad


# ------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------

def _normalize_label(value: object) -> object:
    """Strip + collapse whitespace; preserve case otherwise."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return value
    return re.sub(r"\s+", " ", str(value)).strip()


def _normalize_token(value: object) -> str:
    """Lower-case, strip, collapse whitespace; for substring filter
    matching on region / timepoint / disease."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    s = str(value).strip().lower()
    return re.sub(r"\s+", " ", s)


def _autodetect_column(
    obs: pd.DataFrame, candidates: tuple[str, ...], purpose: str,
) -> str:
    for c in candidates:
        if c in obs.columns:
            return c
    raise KeyError(
        f"TAURUS loader: could not auto-detect {purpose} column. Tried "
        f"{list(candidates)}. obs columns: {sorted(obs.columns)}. Pass "
        f"the explicit name via load(..., {purpose}_col=...)."
    )


def _canonicalize_disease(value: object) -> str | None:
    """Normalize disease/diagnosis values to short tokens: UC, CD, HC."""
    t = _normalize_token(value)
    if not t:
        return None
    if "ulcerative" in t or t == "uc":
        return "UC"
    if "crohn" in t or t == "cd":
        return "CD"
    if "healthy" in t or "control" in t or t in ("hc", "normal", "non-ibd"):
        return "HC"
    return t  # unknown — surface to the gate


def _is_baseline(value: object) -> bool:
    """True iff the normalized value equals a known baseline token.

    Exact match (not substring) because Supp Table 1B uses "Pre"/"Post"
    (3 / 4 chars) which a substring `key in t` would have dropped (no
    member of BASELINE_TIMEPOINT_KEYS is a substring of "pre").
    """
    t = _normalize_token(value)
    return bool(t) and t in BASELINE_TIMEPOINT_KEYS


def _is_colonic(value: object) -> bool:
    """True iff the normalized Ileum_vs_Colon value is Colon or Rectum.

    TAURUS's Ileum_vs_Colon takes exactly {Colon, Rectum, Ileum};
    match is EXACT after lower-case + whitespace-normalize (not
    substring) so a schema drift surfaces as a filter drop, not a
    silent mis-classification.
    """
    return _normalize_token(value) in COLONIC_REGION_VALUES


# ------------------------------------------------------------------------
# Main entry point
# ------------------------------------------------------------------------

def load(
    h5ad_path: str,
    apply_v1_filter: bool = True,
    raw_count_mode: bool = False,
    donor_col: str | None = None,
    disease_col: str | None = None,
    region_filter_col: str | None = None,
    site_col: str | None = None,
    timepoint_col: str | None = None,
    inflammation_col: str | None = None,
    major_col: str | None = None,
    minor_col: str | None = None,
) -> AnnData:
    """Load TAURUS-IBD pooled h5ad and subset to v1 UC × colonic × baseline.

    Parameters
    ----------
    h5ad_path
        Path to ``TAURUS_raw_counts_annotated_final.h5ad`` (md5
        ``c1bd13b92cacb164a401c6c4a4e7912c``; Zenodo
        ``10.5281/zenodo.14007626``).
    apply_v1_filter
        If True (default), apply the UC × colonic × baseline ×
        inflammation>6.5 filter chain and run the donor-structure hard
        invariant (22 UC donors). If False, returns the full atlas with
        only the schema validation — debug/inspection path.
    raw_count_mode
        Must remain False for v1 (DECISIONS 5/7). The TAURUS pooled file
        ships raw counts in X (per the filename); this loader applies
        ``log1p(CP10k)`` and preserves the raw matrix in
        ``layers['counts']``.
    donor_col, disease_col, region_filter_col, site_col, timepoint_col,
    inflammation_col
        Optional explicit obs column names. Auto-detected when not
        provided. ``region_filter_col`` is the coarse column used for
        the colonic filter (TAURUS: ``Ileum_vs_Colon``); ``site_col``
        is the finer column stored as ``obs['region']`` and used for
        cohort validation + covariate structure (TAURUS: ``Site``).
    major_col, minor_col
        Optional explicit obs column names for TAURUS's 2-tier cell-type
        hierarchy (``major``→broad, ``minor``→fine). Auto-detected when
        not provided.

    Returns
    -------
    AnnData
        cells × genes; ``X`` = log1p(CP10k) float; raw counts preserved
        in ``layers['counts']``; obs schema: ``cell_type_fine``,
        ``cell_type_broad``, ``donor``, ``donor_id``, ``disease``,
        ``region`` (=Site), ``timepoint``, ``inflammation_score``,
        ``batch``, ``tissue`` — and ``sample_id`` when the source h5ad
        provides it. var schema depends on the source h5ad's gene
        representation; the final ``ensembl_to_hgnc`` step normalizes
        ``var_names`` to HGNC.
    """
    if raw_count_mode:
        raise ValueError(
            "raw_count_mode=True is not supported for v1 (DECISIONS "
            "5/7): all atlases use uniform log1p(CP10k) input. The "
            "TAURUS pooled file ships raw counts and is normalized on "
            "load."
        )

    h5ad_path = Path(h5ad_path).expanduser()
    if not h5ad_path.exists():
        raise FileNotFoundError(
            f"TAURUS loader: pooled h5ad not found: {h5ad_path}. Download "
            f"{POOLED_FILENAME} from Zenodo {ZENODO_DOI} (md5 "
            f"{POOLED_MD5}). See DECISIONS correction (16)."
        )

    # ---- 1. Read in backed mode so we can filter obs without paying
    #         for the 12.7 GB X matrix until after the subset is known.
    logger.info("Opening TAURUS h5ad in backed mode: %s", h5ad_path)
    adata = ad.read_h5ad(h5ad_path, backed="r")
    logger.info(
        "Backed read: %d cells x %d genes (pre-filter)",
        adata.n_obs, adata.n_vars,
    )

    obs = adata.obs.copy()  # detach from backed file for filter work

    # ---- 2. Auto-detect schema columns ----
    dcol = donor_col or _autodetect_column(obs, _DONOR_COL_CANDIDATES, "donor")
    discol = disease_col or _autodetect_column(obs, _DISEASE_COL_CANDIDATES, "disease")
    rfcol = region_filter_col or _autodetect_column(
        obs, _REGION_FILTER_COL_CANDIDATES, "region_filter"
    )
    scol = site_col or _autodetect_column(obs, _SITE_COL_CANDIDATES, "site")
    tcol = timepoint_col or _autodetect_column(obs, _TIMEPOINT_COL_CANDIDATES, "timepoint")
    icol = inflammation_col or _autodetect_column(
        obs, _INFLAMMATION_COL_CANDIDATES, "inflammation"
    )
    majcol = major_col or _autodetect_column(obs, _MAJOR_COL_CANDIDATES, "major")
    mincol = minor_col or _autodetect_column(obs, _MINOR_COL_CANDIDATES, "minor")
    logger.info(
        "Auto-detected obs columns: donor=%r disease=%r region_filter=%r "
        "site=%r timepoint=%r inflammation=%r major=%r minor=%r",
        dcol, discol, rfcol, scol, tcol, icol, majcol, mincol,
    )

    # Canonicalize disease for the disease-set assertion + filter.
    obs["_disease_canon"] = obs[discol].map(_canonicalize_disease)
    disease_counts = obs["_disease_canon"].value_counts(dropna=False).to_dict()
    logger.info("Disease values (canonicalized): %s", disease_counts)
    expected_diseases = {"UC", "CD", "HC"}
    unknown_diseases = set(disease_counts) - expected_diseases - {None}
    if unknown_diseases:
        raise ValueError(
            f"TAURUS loader: unrecognized disease values after "
            f"canonicalization: {sorted(unknown_diseases)}. Expected only "
            f"{sorted(expected_diseases)}. Extend _canonicalize_disease."
        )

    if not apply_v1_filter:
        logger.warning(
            "apply_v1_filter=False — returning full atlas with only "
            "schema validation; donor invariant and v1 filter chain "
            "skipped."
        )
        # Materialize and skip directly to cell-type schema + normalize.
        return _finalize(
            adata, obs, majcol, mincol,
            dcol, discol, scol, tcol, icol,
            run_donor_assert=False,
        )

    # ---- 3. Four-stage filter chain. Each stage logs n_dropped + n_kept. ----
    n0 = len(obs)

    # Stage A: disease == UC (drop CD + HC).
    keep_disease = obs["_disease_canon"].isin(KEEP_DISEASE)
    n_drop_disease = int((~keep_disease).sum())
    obs = obs[keep_disease].copy()
    logger.info(
        "Filter A (disease in %s): dropped %d, kept %d",
        list(KEEP_DISEASE), n_drop_disease, len(obs),
    )

    # Stage B: colonic region via Ileum_vs_Colon in {Colon, Rectum}.
    keep_region = obs[rfcol].map(_is_colonic).fillna(False).astype(bool)
    n_drop_region = int((~keep_region).sum())
    region_dropped_values = sorted(
        set(obs.loc[~keep_region, rfcol].dropna().astype(str).unique())
    )
    obs = obs[keep_region].copy()
    logger.info(
        "Filter B (%s in %s): dropped %d (values excluded: %s); kept %d",
        rfcol, sorted(COLONIC_REGION_VALUES),
        n_drop_region, region_dropped_values, len(obs),
    )

    # Stage C: baseline timepoint (drop W14 / post-treatment).
    keep_timepoint = obs[tcol].map(_is_baseline).fillna(False).astype(bool)
    n_drop_timepoint = int((~keep_timepoint).sum())
    obs = obs[keep_timepoint].copy()
    logger.info(
        "Filter C (baseline timepoint): dropped %d, kept %d",
        n_drop_timepoint, len(obs),
    )

    # NOTE: inflammation_score > 6.5 is NOT applied here. The paper's
    # inflamed-baseline cutoff is preserved as obs metadata so OPEN_FLAGS
    # F1 can apply one inflamed/non-inflamed/pooled policy uniformly
    # across Smillie + Garrido + TAURUS downstream. See DECISIONS (18)(a)
    # and the module docstring for the rationale.
    infl_summary = pd.to_numeric(obs[icol], errors="coerce")
    logger.info(
        "Carrying inflammation_score through obs (NOT filtering): "
        "n_with_score=%d, n_NaN=%d, range=[%s, %s], n>6.5=%d (the paper's "
        "inflamed-baseline analysis cutoff, applied by F1 downstream).",
        int(infl_summary.notna().sum()), int(infl_summary.isna().sum()),
        f"{infl_summary.min():.2f}" if infl_summary.notna().any() else "—",
        f"{infl_summary.max():.2f}" if infl_summary.notna().any() else "—",
        int((infl_summary > PAPER_BASELINE_INFLAMMATION_MIN).sum()),
    )

    logger.info(
        "Total filter chain: %d -> %d cells (%.1f%% retained)",
        n0, len(obs), 100 * len(obs) / max(n0, 1),
    )
    if len(obs) == 0:
        raise ValueError(
            "TAURUS loader: filter chain produced zero cells. "
            "Inspect the per-stage logs to find which filter eliminated "
            "the cohort; common cause is a column-name auto-detect that "
            "picked up the wrong column. Pass explicit *_col= arguments."
        )

    # ---- 4. Materialize the filtered subset out of backed mode. ----
    kept_index = obs.index
    logger.info("Materializing %d-cell filtered subset out of backed file...", len(kept_index))
    adata_sub = adata[kept_index].to_memory()
    # Replace obs with our augmented one (carries _disease_canon).
    adata_sub.obs = obs

    # ---- 5. Donor + sample invariants from Supp Table 1B dry-run. ----
    # Hard pins derived from the actual Supp Table 1 metadata (DECISIONS
    # 20). All 22 UC patients have at least one Pre colonic sample;
    # total 52 samples across the 22.
    n_uc_donors = int(adata_sub.obs[dcol].nunique())
    breakdown = adata_sub.obs[dcol].value_counts().head(30).to_dict()
    if n_uc_donors == 0 or n_uc_donors > 60:
        raise ValueError(
            f"TAURUS loader: implausible donor count after filter "
            f"(n_uc_donors={n_uc_donors}). Filter chain almost certainly "
            f"misconfigured (e.g., wrong column auto-detected). "
            f"Per-donor cell counts (top 30): {breakdown}."
        )
    if n_uc_donors != EXPECTED_N_UC_DONORS_POST_FILTER:
        # Compute per-donor, per-region diff against Supp Table 1B so the
        # failure points at the actual drift, not just the donor-count
        # delta. v3 of the Zenodo deposit exists because of a metadata
        # fix; an obs-vs-Supp-Table-1B disagreement is meaningful, and we
        # want "donor X missing from region Y", not "expected 22, got 21".
        obs_donors = set(adata_sub.obs[dcol].astype(str).unique())
        exp_donors = set(EXPECTED_UC_COHORT)
        missing_from_obs = sorted(exp_donors - obs_donors)
        unexpected_in_obs = sorted(obs_donors - exp_donors)
        # Per-donor region diff (only for donors that exist on both sides)
        obs_by_donor_region = (
            adata_sub.obs[[dcol, scol]].astype(str)
            .drop_duplicates()
            .groupby(dcol)[scol]
            .agg(lambda s: frozenset(s.dropna()))
            .to_dict()
        )
        region_drifts = []
        for d in sorted(exp_donors & obs_donors):
            exp_r = EXPECTED_UC_COHORT[d]
            obs_r = obs_by_donor_region.get(d, frozenset())
            if exp_r != obs_r:
                missing_r = sorted(exp_r - obs_r)
                extra_r = sorted(obs_r - exp_r)
                region_drifts.append(
                    f"    {d}: expected {sorted(exp_r)}; observed "
                    f"{sorted(obs_r)}; missing {missing_r}; extra {extra_r}"
                )
        region_drift_block = (
            "\n".join(region_drifts) if region_drifts
            else "    (none — observed regions match Supp Table 1B for "
                 "donors present on both sides)"
        )
        raise ValueError(
            f"TAURUS loader: UC cohort mismatch vs Supp Table 1B "
            f"(DECISIONS 20 / 21).\n"
            f"  expected donors: {EXPECTED_N_UC_DONORS_POST_FILTER}, "
            f"observed: {n_uc_donors}\n"
            f"  missing from observed: {missing_from_obs}\n"
            f"  unexpected in observed: {unexpected_in_obs}\n"
            f"  per-donor region drifts:\n{region_drift_block}\n"
            f"  per-donor cell counts (top 30): {breakdown}\n"
            f"v3 of the Zenodo deposit exists *because* of a metadata "
            f"fix vs earlier versions, so an obs-vs-Supp-Table-1B "
            f"disagreement is signal, not noise. Either the h5ad's "
            f"obs columns drifted from the published Supp Table, or "
            f"the filter chain auto-detected a wrong column."
        )
    logger.info(
        "Donor invariant passed: %d UC donors (matches Supp Table 1B).",
        n_uc_donors,
    )

    # Sample count check (soft warning, not hard — the h5ad may have
    # collapsed some sample distinctions post-QC).
    if "sample_id" in adata_sub.obs.columns:
        n_samples = int(adata_sub.obs["sample_id"].nunique())
    elif "Batch" in adata_sub.obs.columns:
        n_samples = int(adata_sub.obs["Batch"].nunique())
    else:
        n_samples = None
    if n_samples is not None and n_samples != EXPECTED_N_UC_SAMPLES_POST_FILTER:
        logger.warning(
            "TAURUS UC sample count %d != expected %d (Supp Table 1B). "
            "Soft tripwire; could be a benign QC drop. Hard donor gate "
            "above has already passed.",
            n_samples, EXPECTED_N_UC_SAMPLES_POST_FILTER,
        )

    if EXPECTED_N_CELLS_HINT is not None and adata_sub.n_obs != EXPECTED_N_CELLS_HINT:
        logger.warning(
            "Cell count %d != expected %d (tripwire only; the donor "
            "invariant above is the hard gate).",
            adata_sub.n_obs, EXPECTED_N_CELLS_HINT,
        )

    return _finalize(
        adata_sub, adata_sub.obs, majcol, mincol,
        dcol, discol, scol, tcol, icol,
        run_donor_assert=True,
    )


def _finalize(
    adata: AnnData,
    obs: pd.DataFrame,
    majcol: str, mincol: str,
    dcol: str, discol: str, scol: str, tcol: str, icol: str,
    run_donor_assert: bool,
) -> AnnData:
    """Build standard obs schema, validate canonical vocab, normalize,
    HGNC-remap. Shared between the apply_v1_filter=True and =False paths.

    TAURUS has only 2 cell-type tiers (major → broad, minor → fine); the
    cell_type_compartment / _low / _intermediate obs columns are omitted.
    """
    # ---- 6. Build the standard obs schema. ----
    new_obs = pd.DataFrame(index=obs.index)
    new_obs["cell_type_fine"] = obs[mincol].astype(str).map(_normalize_label).astype("category")
    new_obs["donor_id"]   = obs[dcol].astype(str).astype("category")
    new_obs["donor"]      = new_obs["donor_id"]
    new_obs["disease"]    = obs.get("_disease_canon", obs[discol].map(_canonicalize_disease)).astype("category")
    new_obs["region"]     = obs[scol].astype(str).map(_normalize_label).astype("category")
    new_obs["timepoint"]  = obs[tcol].astype(str).map(_normalize_label).astype("category")
    new_obs["inflammation_score"] = pd.to_numeric(obs[icol], errors="coerce")
    new_obs["batch"]      = obs[dcol].astype(str).astype("category")  # per-donor batches
    new_obs["tissue"]     = "colonic mucosa"
    # Carry through biopsy-level sample_id when the source h5ad provides
    # it — run_taurus_load.py prefers it over region for the covariate
    # file's sample dummy.
    if "sample_id" in obs.columns:
        new_obs["sample_id"] = obs["sample_id"].astype(str).astype("category")

    # ---- 7. Map major -> canonical broad; assert vocab membership (gate 2). ----
    major = obs[majcol].astype(str).map(_normalize_label)
    unmapped = sorted(set(major.unique()) - set(MAJOR_TO_BROAD))
    if unmapped:
        raise KeyError(
            f"TAURUS loader: {len(unmapped)} major-tier labels have no "
            f"MAJOR_TO_BROAD entry. Extend the map "
            f"(load_taurus.MAJOR_TO_BROAD). Unmapped labels (full list): "
            f"{unmapped}"
        )
    broad = major.map(MAJOR_TO_BROAD)
    emitted = set(broad.dropna().unique())
    unrecognized = emitted - _BROAD_VOCAB
    if unrecognized:
        raise ValueError(
            f"TAURUS loader: emitted cell_type_broad values "
            f"{sorted(unrecognized)} are not in the canonical vocab. "
            f"Fix MAJOR_TO_BROAD value side; see canonical_broad_DRAFT.md."
        )
    new_obs["cell_type_broad"] = broad.astype("category")

    n_fine  = int(new_obs["cell_type_fine"].nunique())
    n_broad = int(new_obs["cell_type_broad"].nunique())
    logger.info(
        "Tier cardinalities: major/broad=%d, minor/fine=%d",
        n_broad, n_fine,
    )
    if not (5 <= n_broad <= 15):
        logger.warning(
            "Broad-tier cardinality %d outside the v1 5-15 target for "
            "TAURUS (paper reports 8 major); review MAJOR_TO_BROAD grouping.",
            n_broad,
        )

    adata.obs = new_obs

    # ---- 7b. Cross-atlas QC policy (DECISIONS 22).
    # EXCLUDE_LINEAGE_AMBIGUOUS_FINE keys on cell_type_fine — applies as
    # for Garrido/Smillie. QC_STATE_TO_PARENT keys on cell_type_low which
    # TAURUS does not have; it is skipped here. Revisit if a
    # minor-tier equivalent surfaces.
    fine_str = adata.obs["cell_type_fine"].astype(str)
    excl_mask = fine_str.isin(EXCLUDE_LINEAGE_AMBIGUOUS_FINE)
    n_excl = int(excl_mask.sum())
    if n_excl:
        excl_labels = sorted(set(fine_str[excl_mask]))
        logger.info(
            "EXCLUDE_LINEAGE_AMBIGUOUS_FINE: dropping %d cells across %d "
            "fine labels %s (DECISIONS 22; revisit in marker-QC).",
            n_excl, len(excl_labels), excl_labels,
        )
        adata = adata[(~excl_mask).values].copy()

    if run_donor_assert:
        per_donor = adata.obs["donor"].value_counts()
        logger.info("UC donors (n=%d) post-filter:", per_donor.size)
        for donor, count in per_donor.items():
            logger.info("  %s: %d cells", donor, count)

    # ---- 8. Normalize raw counts -> log1p(CP10k); preserve raw in layer. ----
    if not np.issubdtype(adata.X.dtype, np.floating):
        adata.X = adata.X.astype(np.float32)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    logger.info("Normalized X to log1p(CP10k); raw counts in layers['counts'].")

    # ---- 9. Ensembl -> HGNC via pinned remap (correction 11). ----
    adata = ensembl_to_hgnc(adata)

    logger.info(
        "TAURUS load complete: %d cells x %d genes, %d donors.",
        adata.n_obs, adata.n_vars, int(adata.obs["donor"].nunique()),
    )
    return adata
