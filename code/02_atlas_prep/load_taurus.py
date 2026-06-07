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

**Annotation hierarchy** (Methods, Extended Data Fig. 1b): 4 levels —
``compartment`` → ``low`` → ``intermediate`` → ``cell_state``. The
``cell_state`` tier is the 109-class finest output (used as
``obs['cell_type_fine']``). The ``low`` tier is the source for
``obs['cell_type_broad']`` via ``LOW_TO_BROAD`` mapping into the
canonical 15-term vocab. All four hierarchy levels are preserved in obs
for downstream use.

**Filter chain (v1)** — three stages, NOT four:

1. Disease == UC (drop CD + HC).
2. Region in colonic set (drop terminal ileum; keep ascending /
   descending / sigmoid / rectum and any other non-ileal label).
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

- Donor-structure hard invariant: 22 UC donors after the full filter.
- Cell-count tripwire (soft): expect ~30-50k cells but TBD on first run;
  log only.
- Canonical-vocab assertions (gates 1 + 2): same two-gate pattern as
  the Garrido + Smillie loaders. ``LOW_TO_BROAD`` ships **empty** in
  this v0 because the ``low``-tier label set is not in the paper's
  Methods text — gate 2 will fail loud on first run with the actual
  labels listed, and the map then gets filled in (one commit per
  Muskaan biology pass).
- Counts pipeline: ``log1p(CP10k)`` on load per DECISIONS (5/7); raw
  integer counts preserved in ``layers['counts']``;
  ``raw_count_mode=True`` unsupported.

**Open before first run** (DECISIONS 16 + load gate will catch):

- The exact obs column names — Methods doesn't name them. This loader
  uses defensive auto-detect against candidate name lists; KeyError
  raises if no candidate matches, with the actual obs.columns dumped
  for triage.
- The ``low``-tier label set → ``LOW_TO_BROAD`` map.
- Exact per-donor cell counts after subset (Supp Table 1 has these;
  load-gate will validate when the values are populated below).

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

# Region values that count as colonic. Match is substring-based after
# whitespace/case normalization (so "ascending colon" matches; so does
# "ascending_colon"). Terminal ileum and any other ileal region drop.
COLONIC_REGION_KEYS: tuple[str, ...] = (
    "ascending", "descending", "transverse", "sigmoid", "rectum", "colon",
)
ILEAL_REGION_KEYS: tuple[str, ...] = ("ileum", "ileal")

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
_REGION_COL_CANDIDATES = (
    "Site",  # TAURUS Supp Table 1B
    "region", "Region", "tissue", "Tissue",
    "anatomical_region", "site", "location", "Location",
    "biopsy_site", "tissue_region",
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
_HIERARCHY_LEVELS = ("compartment", "low", "intermediate", "cell_state")
_HIERARCHY_CANDIDATES = {
    "compartment": ("compartment", "Compartment", "lineage", "level_1",
                    "cell_type_compartment"),
    "low":         ("low", "Low", "cell_type_low", "level_2",
                    "broad_celltype", "cell_class"),
    "intermediate":("intermediate", "Intermediate", "cell_type_intermediate",
                    "level_3", "cell_type", "Celltype"),
    "cell_state":  ("cell_state", "Cell_State", "level_4", "state",
                    "Annotation", "annotation", "cellstate"),
}

# ---- Canonical broad vocab (single-sourced; same 15 in every loader) ----
# Single-sourced from sibling _broad_vocab module. Loader-local copies
# were a drift risk that 06_concordance would have reported as biology.
# Gate (1) at module load + gate (2) at end-of-load. See DECISIONS (20).
from _broad_vocab import _BROAD_VOCAB

# Map from TAURUS `low`-tier labels into the canonical broad vocab.
# Ships EMPTY in this v0 — the `low`-tier label set is not in the paper
# Methods text. Gate (2) at end-of-load will fail loud on first run with
# every unmapped label listed; populate below (one commit per Muskaan
# biology pass) until the gate passes. The Garrido and Smillie loaders
# document analogous maps (FINE_TO_BROAD); they're the template.
LOW_TO_BROAD: dict[str, str] = {
    # TODO(taurus-first-run): populate from the actual `low`-tier label
    # set after the first end-to-end load. Until then, gate (2) raises.
}

# Gate (1): every value the map ships must be in the canonical vocab.
# Vacuously true while the map is empty; protects the day we start
# filling it in.
_unmapped_broad = set(LOW_TO_BROAD.values()) - _BROAD_VOCAB
if _unmapped_broad:
    raise ValueError(
        f"load_taurus.LOW_TO_BROAD ships broad values outside "
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
    t = _normalize_token(value)
    if not t:
        return False
    if any(key in t for key in ILEAL_REGION_KEYS):
        return False
    return any(key in t for key in COLONIC_REGION_KEYS)


# ------------------------------------------------------------------------
# Main entry point
# ------------------------------------------------------------------------

def load(
    h5ad_path: str,
    apply_v1_filter: bool = True,
    raw_count_mode: bool = False,
    donor_col: str | None = None,
    disease_col: str | None = None,
    region_col: str | None = None,
    timepoint_col: str | None = None,
    inflammation_col: str | None = None,
    fine_col: str | None = None,
    low_col: str | None = None,
    intermediate_col: str | None = None,
    compartment_col: str | None = None,
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
    donor_col, disease_col, region_col, timepoint_col, inflammation_col
        Optional explicit obs column names. Auto-detected when not
        provided.
    fine_col, low_col, intermediate_col, compartment_col
        Optional explicit obs column names for the 4-level cell-type
        hierarchy. Auto-detected when not provided.

    Returns
    -------
    AnnData
        cells × genes; ``X`` = log1p(CP10k) float; raw counts preserved
        in ``layers['counts']``; obs schema: ``cell_type_fine``,
        ``cell_type_intermediate``, ``cell_type_low``,
        ``cell_type_compartment``, ``cell_type_broad``, ``donor``,
        ``donor_id``, ``disease``, ``region``, ``timepoint``,
        ``inflammation_score``, ``batch``, ``tissue``. var schema
        depends on the source h5ad's gene representation; the final
        ``ensembl_to_hgnc`` step normalizes ``var_names`` to HGNC.
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
    rcol = region_col or _autodetect_column(obs, _REGION_COL_CANDIDATES, "region")
    tcol = timepoint_col or _autodetect_column(obs, _TIMEPOINT_COL_CANDIDATES, "timepoint")
    icol = inflammation_col or _autodetect_column(
        obs, _INFLAMMATION_COL_CANDIDATES, "inflammation"
    )
    hierarchy_cols = {
        "compartment":  compartment_col  or _autodetect_column(obs, _HIERARCHY_CANDIDATES["compartment"],  "compartment"),
        "low":          low_col          or _autodetect_column(obs, _HIERARCHY_CANDIDATES["low"],          "low"),
        "intermediate": intermediate_col or _autodetect_column(obs, _HIERARCHY_CANDIDATES["intermediate"], "intermediate"),
        "cell_state":   fine_col         or _autodetect_column(obs, _HIERARCHY_CANDIDATES["cell_state"],   "cell_state"),
    }
    logger.info(
        "Auto-detected obs columns: donor=%r disease=%r region=%r "
        "timepoint=%r inflammation=%r hierarchy=%s",
        dcol, discol, rcol, tcol, icol, hierarchy_cols,
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
            adata, obs, hierarchy_cols, dcol, discol, rcol, tcol, icol,
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

    # Stage B: colonic region (drop terminal ileum + any non-colonic).
    keep_region = obs[rcol].map(_is_colonic).fillna(False).astype(bool)
    n_drop_region = int((~keep_region).sum())
    region_dropped_values = sorted(
        set(obs.loc[~keep_region, rcol].dropna().astype(str).unique())
    )
    obs = obs[keep_region].copy()
    logger.info(
        "Filter B (colonic region only): dropped %d (regions excluded: %s); kept %d",
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
        raise ValueError(
            f"TAURUS loader: UC donor count {n_uc_donors} != expected "
            f"{EXPECTED_N_UC_DONORS_POST_FILTER} (Supp Table 1B, "
            f"DECISIONS 20). The Supp Table is the empirical ground "
            f"truth for the v1 subset; a mismatch means metadata drift "
            f"in the h5ad obs vs the paper, OR a filter is wrong. "
            f"Per-donor cell counts (top 30): {breakdown}."
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
        adata_sub, adata_sub.obs, hierarchy_cols, dcol, discol, rcol, tcol, icol,
        run_donor_assert=True,
    )


def _finalize(
    adata: AnnData,
    obs: pd.DataFrame,
    hierarchy_cols: dict[str, str],
    dcol: str, discol: str, rcol: str, tcol: str, icol: str,
    run_donor_assert: bool,
) -> AnnData:
    """Build standard obs schema, validate canonical vocab, normalize,
    HGNC-remap. Shared between the apply_v1_filter=True and =False paths.
    """
    # ---- 6. Build the standard obs schema. ----
    new_obs = pd.DataFrame(index=obs.index)
    new_obs["cell_type_compartment"]  = obs[hierarchy_cols["compartment"]].astype(str).map(_normalize_label).astype("category")
    new_obs["cell_type_low"]          = obs[hierarchy_cols["low"]].astype(str).map(_normalize_label).astype("category")
    new_obs["cell_type_intermediate"] = obs[hierarchy_cols["intermediate"]].astype(str).map(_normalize_label).astype("category")
    new_obs["cell_type_fine"]         = obs[hierarchy_cols["cell_state"]].astype(str).map(_normalize_label).astype("category")
    new_obs["donor_id"]   = obs[dcol].astype(str).astype("category")
    new_obs["donor"]      = new_obs["donor_id"]
    new_obs["disease"]    = obs.get("_disease_canon", obs[discol].map(_canonicalize_disease)).astype("category")
    new_obs["region"]     = obs[rcol].astype(str).map(_normalize_label).astype("category")
    new_obs["timepoint"]  = obs[tcol].astype(str).map(_normalize_label).astype("category")
    new_obs["inflammation_score"] = pd.to_numeric(obs[icol], errors="coerce")
    new_obs["batch"]      = obs[dcol].astype(str).astype("category")  # per-donor batches
    new_obs["tissue"]     = "colonic mucosa"

    # ---- 7. Map low -> canonical broad; assert vocab membership (gate 2). ----
    low = new_obs["cell_type_low"].astype(str)
    unmapped = sorted(set(low.unique()) - set(LOW_TO_BROAD))
    if unmapped:
        raise KeyError(
            f"TAURUS loader: {len(unmapped)} low-tier labels have no "
            f"LOW_TO_BROAD entry. Extend the map "
            f"(load_taurus.LOW_TO_BROAD). Unmapped labels (full list): "
            f"{unmapped}"
        )
    broad = low.map(LOW_TO_BROAD)
    emitted = set(broad.dropna().unique())
    unrecognized = emitted - _BROAD_VOCAB
    if unrecognized:
        raise ValueError(
            f"TAURUS loader: emitted cell_type_broad values "
            f"{sorted(unrecognized)} are not in the canonical vocab. "
            f"Fix LOW_TO_BROAD value side; see canonical_broad_DRAFT.md."
        )
    new_obs["cell_type_broad"] = broad.astype("category")

    n_fine  = int(new_obs["cell_type_fine"].nunique())
    n_low   = int(new_obs["cell_type_low"].nunique())
    n_int   = int(new_obs["cell_type_intermediate"].nunique())
    n_comp  = int(new_obs["cell_type_compartment"].nunique())
    n_broad = int(new_obs["cell_type_broad"].nunique())
    logger.info(
        "Tier cardinalities: compartment=%d, low=%d, intermediate=%d, "
        "cell_state(fine)=%d, broad=%d",
        n_comp, n_low, n_int, n_fine, n_broad,
    )
    if not (10 <= n_broad <= 15):
        logger.warning(
            "Broad-tier cardinality %d outside the v1 10-15 target; "
            "review LOW_TO_BROAD grouping.", n_broad,
        )

    adata.obs = new_obs

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
