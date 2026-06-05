"""Loader for Garrido-Trigo 2023 (Atlas 2, UC subset).

.. warning::
   **SUPERSEDED MATRIX SOURCE (2026-06-04, DECISIONS correction 9).**
   This loader builds against the CELLxGENE deposit, whose ``obs.index``
   is synthetic (``cell1, cell2, ...``) and whose original 10X barcodes
   were stripped at deposit time. That breaks every barcode-join
   strategy below: there is no shared key on the obs side to match
   against ``GSE214695_cell_annotation.csv``. v1 is moving to a loader
   that pulls the matrix from ``GSE214695_RAW.tar`` per-GSM files, where
   barcodes are intact and the CSV's ``Unnamed: 0`` (``SC_xxx_<barcode>``)
   is the unique, deterministic join key. The RAW.tar path also requires
   ``log1p(CP10k)`` normalization on load to match HCA / Pan-GI (the same
   treatment Mennillo receives) — ``--flag-raw-count False`` stays
   uniform per correction 5/7. Rewrite scheduled alongside the Smillie
   SCP259 download. This file is retained as reference; ``load()`` will
   fail loudly at the barcode-join step on the CELLxGENE matrix.

Sources (pre-correction-9, retained for reference)
--------------------------------------------------
- **Matrix:** CELLxGENE deposit ``b1a62801-f509-45f8-b55f-533fbb7e7800.h5ad``
  (log-normalized X; var_names are Ensembl IDs; HGNC symbols in
  ``var['feature_name']``).
- **Labels:** GEO supplementary file ``GSE214695_cell_annotation.csv``
  (91-label fine annotation from the Salas lab; not present in the
  CELLxGENE deposit).

The CELLxGENE-only path is deprecated: it shipped a 5-CL-lineage broad
label with no fine tier, blocking fine-tier cross-atlas concordance for
Garrido-Trigo. The GEO supplementary annotation restores the full
broad + fine tier. See DECISIONS.md correction (4/7), the correction
reversing it (8), and correction 9 superseding the matrix source.

References: ``code/02_atlas_prep/atlas_schemas.md``;
DECISIONS.md corrections 2026-05-20 (2/7), (4/7), (5/7), 2026-06-03 (8),
and 2026-06-04 (9).
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import anndata as ad
from anndata import AnnData
import numpy as np
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

# Expected counts after UC + healthy subset filter
# (DECISIONS.md correction 2026-05-20 (2/7)).
EXPECTED_UC_SUBSET_N_CELLS = 30_068
EXPECTED_UC_SUBSET_N_DONORS = 12

# Candidate column names for auto-detection in the GEO annotation CSV.
# The exact schema of GSE214695_cell_annotation.csv is not pinned; auto-detect
# rather than hard-code so a schema drift surfaces as a loud error.
_BARCODE_COL_CANDIDATES = (
    # "Unnamed: 0" is first because in GSE214695_cell_annotation.csv it
    # holds the unique SC_xxx_<barcode> composite, while "cell_id" is the
    # bare 10X barcode that collides across samples (280 duplicates).
    # _load_annotation_csv prefers candidates whose values are unique.
    "Unnamed: 0",
    "cell_id", "cell", "Cell", "barcode", "Barcode", "cellID", "CellID",
    "index",
)
_ANNOTATION_COL_CANDIDATES = (
    "annotation", "Annotation", "cluster", "Cluster", "Population",
    "cell_type_fine", "celltype", "cell_type", "subset", "label",
)
_DONOR_COL_CANDIDATES = (
    "donor", "Donor", "donor_id", "sample", "Sample", "sample_id",
    "orig.ident", "patient", "Patient",
)


# Ribhi = ribosomal-high cross-lineage transcriptional state, not a lineage.
# Verified empirically: top-20 markers of every Ribhi cluster are dominated
# by RPL*/RPS* genes (plus lineage markers in the 4 stromal/myeloid Ribhis).
# Collapse each Ribhi label to its parent before any tier logic.
RIBHI_TO_PARENT = {
    "B cell Ribhi":     "B cell",
    "Epithelium Ribhi": "epithelial",
    "M0_Ribhi":         "M0",
    "Ribhi T cells 1":  "T",
    "Ribhi T cells 2":  "T",
    "S1 Ribhi":         "S1",
    "Fibroblasts Ribhi":"fibroblast",
    "Mast Ribhi":       "mast",
    "DCs CCL22_Ribhi":  "DCs CCL22",
}

# Fine -> broad roll-up. Built to target the v1 ~10-15 broad-tier budget.
# Keys are the 91 published Garrido-Trigo fine labels *plus* the generic
# parent labels introduced by RIBHI_TO_PARENT ("epithelial", "T",
# "fibroblast", "mast"). Whitespace is normalized before lookup
# (see _normalize_label).
FINE_TO_BROAD = {
    # --- Epithelium ---
    "BEST4 OTOP2":            "colonocyte",
    "Colonocyte 1":           "colonocyte",
    "Colonocyte 2":           "colonocyte",
    "Inflammatory colonocyte":"colonocyte",
    "Laminin colonocytes":    "colonocyte",
    "PLCG2 colonocytes":      "colonocyte",
    "Cycling TA":             "epithelial progenitor",
    "Secretory progenitor":   "epithelial progenitor",
    "Goblet":                 "goblet",
    "Mature goblet":          "goblet",
    "Paneth-like":            "goblet",
    "Enteroendocrine":        "enteroendocrine/tuft",
    "Tuft cells":             "enteroendocrine/tuft",
    "epithelial":             "epithelial progenitor",  # Ribhi parent

    # --- Stroma ---
    "S1":                    "fibroblast",
    "S1.2":                  "fibroblast",
    "S2a":                   "fibroblast",
    "S2b":                   "fibroblast",
    "S3":                    "fibroblast",
    "IER fibroblasts":       "fibroblast",
    "Inflammatory fibroblasts":"fibroblast",
    "MT fibroblasts":        "fibroblast",
    "Myofibroblasts":        "fibroblast",
    "FRCs":                  "fibroblast",
    "fibroblast":            "fibroblast",  # Ribhi parent
    "Endothelium":           "endothelium",
    "Activated endothelium": "endothelium",
    "Lymphatic endothelium": "endothelium",
    "Perycites":             "mural/glia",
    "Glia":                  "mural/glia",

    # --- T cells / NK / ILC ---
    "CD4 ANXA1":      "T cell",
    "CD4 naive":      "T cell",  # whitespace/unicode-normalized form
    "CD8 CTL":        "T cell",
    "CD8 CTL TRM":    "T cell",
    "CD8 FGFBP2":     "T cell",
    "Cycling T cells":"T cell",
    "DN EOMES":       "T cell",
    "DN TNF":         "T cell",
    "MT T cells":     "T cell",
    "S1PR1 T cells":  "T cell",
    "T cells CCL20":  "T cell",
    "ThF":            "T cell",
    "Tregs":          "T cell",
    "gd IEL":         "T cell",
    "T":              "T cell",   # Ribhi parent
    "MAIT":           "T cell",
    "NK":             "NK/ILC",
    "ILC3":           "NK/ILC",

    # --- B cells / plasma cells ---
    "B cell":         "B cell",
    "Memory B cell":  "B cell",
    "Naive B cell":   "B cell",   # whitespace/unicode-normalized form
    "GC B cell":      "B cell",
    "Cycling cells":  "B cell",   # B/plasma cycling pool
    "Cycling cells 2":"B cell",
    "Cycling cells 3":"B cell",
    "PC IER":                     "plasma cell",
    "PC immediate early response":"plasma cell",  # GEO long-form synonym
    "PC IGLL5":                   "plasma cell",
    "PC IgA 1":                   "plasma cell",
    "PC IgA 2":                   "plasma cell",
    "PC IgA 3":                   "plasma cell",
    "PC IgA 4":                   "plasma cell",
    "PC IgA IgM":                 "plasma cell",
    "PC IgA Lambda 1":            "plasma cell",
    "PC IgA heat shock 1":        "plasma cell",
    "PC IgA heat shock 2":        "plasma cell",
    "PC IgG 1":                   "plasma cell",
    "PC IgG 2":                   "plasma cell",
    "Plasmablast IgA Lambda 2":   "plasma cell",
    "Plasmablast IgG":            "plasma cell",
    "Plasmablast IgG Lambda":     "plasma cell",

    # --- Myeloid: monocyte/macrophage ---
    "M0":                   "monocyte/macrophage",
    "M1 ACOD1":             "monocyte/macrophage",
    "M1 CXCL5":             "monocyte/macrophage",
    "M2":                   "monocyte/macrophage",
    "M2.2":                 "monocyte/macrophage",
    "IDA macrophage":       "monocyte/macrophage",
    "Inflammatory monocytes":"monocyte/macrophage",

    # --- Myeloid: dendritic / mast / granulocyte / cycling ---
    "DCs CCL22":      "dendritic cell",
    "DCs CD1c":       "dendritic cell",
    "Mast 1":         "mast cell",
    "Mast 2":         "mast cell",
    "mast":           "mast cell",  # Ribhi parent
    "Neutrophil 1":   "granulocyte",
    "Neutrophil 2":   "granulocyte",
    "Neutrophil 3":   "granulocyte",
    "Eosinophils":    "granulocyte",
    "Cycling myeloid":"monocyte/macrophage",
}


def _normalize_label(value: object) -> object:
    """Strip + collapse internal whitespace; normalize curly to ASCII.

    The GEO CSV contains at least one label with a literal double space
    ('PC immediate early response'), and Salas-lab labels use Latin-1
    characters that can round-trip as mojibake (e.g., 'Na\xefve B cell').
    Normalizing on load prevents silent join breakage.
    """
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return value
    s = str(value)
    # Normalize the few non-ASCII characters seen in the marker xlsx so the
    # FINE_TO_BROAD map can use plain ASCII keys.
    s = s.replace("ï", "i").replace("\xef", "i")  # Naïve -> Naive
    s = s.replace("é", "e")                       # café -> cafe
    s = re.sub(r"\s+", " ", s).strip()
    return s


def _autodetect_column(
    df: pd.DataFrame,
    candidates: tuple[str, ...],
    purpose: str,
    prefer_unique: bool = False,
) -> str:
    """Return the first candidate column that exists in ``df``.

    If ``prefer_unique`` is True, scan all matching candidates and return
    the first one whose values are unique. Falls back to the first match
    if none are unique. Used for the barcode column, where picking the
    bare-barcode ``cell_id`` over the unique ``Unnamed: 0`` composite
    silently breaks the join.
    """
    matches = [c for c in candidates if c in df.columns]
    if not matches:
        raise KeyError(
            f"GSE214695 annotation CSV: could not auto-detect {purpose} "
            f"column. Tried {list(candidates)}. Got columns: "
            f"{list(df.columns)}. Pass the explicit column name via the "
            f"loader's `{purpose}_col` argument."
        )
    if prefer_unique:
        for c in matches:
            if not df[c].duplicated().any():
                return c
    return matches[0]


def _load_annotation_csv(
    annotation_csv_path: str | Path,
    barcode_col: str | None,
    annotation_col: str | None,
    donor_col: str | None,
) -> tuple[pd.DataFrame, str, str, str | None]:
    """Read GSE214695_cell_annotation.csv and standardize key columns.

    Returns
    -------
    (df, barcode_col, annotation_col, donor_col)
        ``df`` indexed by the normalized barcode column, with the
        annotation column whitespace-normalized in place. ``donor_col``
        may be ``None`` if not present in the file.
    """
    logger.info("Reading GEO annotation CSV: %s", annotation_csv_path)
    ann = pd.read_csv(annotation_csv_path)

    bcol = barcode_col or _autodetect_column(
        ann, _BARCODE_COL_CANDIDATES, "barcode", prefer_unique=True
    )
    acol = annotation_col or _autodetect_column(
        ann, _ANNOTATION_COL_CANDIDATES, "annotation"
    )
    dcol = donor_col
    if dcol is None:
        for c in _DONOR_COL_CANDIDATES:
            if c in ann.columns:
                dcol = c
                break

    ann[acol] = ann[acol].map(_normalize_label)
    ann[bcol] = ann[bcol].astype(str).str.strip()
    if dcol is not None:
        ann[dcol] = ann[dcol].astype(str).str.strip()

    if ann[bcol].duplicated().any():
        n_dup = int(ann[bcol].duplicated().sum())
        raise ValueError(
            f"GSE214695 annotation CSV: {n_dup} duplicate barcodes in "
            f"column {bcol!r}; cannot use as join key."
        )

    logger.info(
        "Annotation CSV: %d rows, columns barcode=%r annotation=%r donor=%r",
        len(ann), bcol, acol, dcol,
    )
    return ann, bcol, acol, dcol


def _try_join_keys(
    obs: pd.DataFrame, ann: pd.DataFrame, bcol: str, dcol: str | None,
) -> tuple[pd.Series, str]:
    """Return (annotation-indexed-by-obs.index, strategy_name).

    Strategies tried, in order:
      1. obs.index (raw .h5ad barcode) == ann[bcol].
      2. obs.index after stripping a trailing "-1" suffix == ann[bcol].
      3. obs.index after adding a "-1" suffix == ann[bcol].
      4. ``donor_id + "_" + obs.index`` == ann[bcol].
      5. ``ann[dcol] + "_" + ann[bcol]`` (rebuilt from CSV side) == obs.index.
      6. obs.index after stripping a leading "SC_" prefix == ann[bcol].
      7. ``"SC_" + donor_id + "_" + obs.index`` == ann[bcol].

    Strategies 6/7 bracket the Salas-lab ``SC_xxx + sample`` cell_id
    format used by the per-GSM matrices.

    The first strategy that hits every obs row (zero orphans) wins. If
    none do, raises ValueError with the best partial match for triage.
    """
    obs_index = obs.index.astype(str)
    candidates: list[tuple[str, pd.Series]] = []

    candidates.append(("raw obs.index", obs_index))
    candidates.append((
        "obs.index without -1 suffix",
        obs_index.str.replace(r"-1$", "", regex=True),
    ))
    candidates.append((
        "obs.index with -1 suffix",
        obs_index.where(obs_index.str.endswith("-1"), obs_index + "-1"),
    ))
    if "donor_id" in obs.columns:
        donor = obs["donor_id"].astype(str).to_numpy()
        idx_np = obs_index.to_numpy()
        candidates.append((
            "donor_id + '_' + obs.index",
            pd.Index(donor + "_" + idx_np),
        ))
        candidates.append((
            "'SC_' + donor_id + '_' + obs.index",
            pd.Index("SC_" + donor + "_" + idx_np),
        ))
    candidates.append((
        "obs.index without leading SC_",
        obs_index.str.replace(r"^SC_", "", regex=True),
    ))
    # Strategy 5 inverts: build the CSV-side composite and try to match obs.
    if dcol is not None:
        composite = (ann[dcol].astype(str) + "_" + ann[bcol].astype(str))
        if composite.is_unique:
            composite_ann = ann.set_index(composite)
            hit = obs_index.isin(composite_ann.index)
            if hit.all():
                return composite_ann.loc[obs_index].reset_index(drop=True), \
                       "ann[donor]_ann[barcode] == obs.index"

    ann_indexed = ann.set_index(bcol)
    best_hits = -1
    best_strategy = ""
    for name, keys in candidates:
        hit = keys.isin(ann_indexed.index)
        n_hit = int(hit.sum())
        if n_hit == len(obs):
            joined = ann_indexed.loc[keys.values].reset_index(drop=True)
            return joined, name
        if n_hit > best_hits:
            best_hits = n_hit
            best_strategy = name

    raise ValueError(
        f"Garrido-Trigo barcode join failed: no strategy matched all "
        f"{len(obs)} cells. Best partial match was {best_strategy!r} with "
        f"{best_hits}/{len(obs)} hits. Inspect the CSV barcode format vs "
        f"adata.obs_names (example obs: {obs_index[:3].tolist()}, "
        f"example CSV: {ann[bcol].head(3).tolist()}) and pass an explicit "
        f"`barcode_col` / `donor_col` to load()."
    )


def load(
    h5ad_path: str,
    annotation_csv_path: str | None = None,
    apply_v1_filter: bool = True,
    raw_count_mode: bool = False,
    barcode_col: str | None = None,
    annotation_col: str | None = None,
    donor_col: str | None = None,
) -> AnnData:
    """Load Garrido-Trigo with full broad + fine tier from the GEO annotation.

    Parameters
    ----------
    h5ad_path
        Path to the CELLxGENE .h5ad (log-normalized matrix + obs).
    annotation_csv_path
        Path to ``GSE214695_cell_annotation.csv`` from GEO (the Salas-lab
        91-cluster fine annotation). Required for full broad+fine output.
        If ``None``, the loader falls back to the CELLxGENE 5-CL labels
        only (broad-tier degraded; emits a warning).
    apply_v1_filter
        If True (default), keep only ``disease in {normal, ulcerative colitis}``.
    raw_count_mode
        Must remain False for v1 (DECISIONS.md correction 5/7).
    barcode_col, annotation_col, donor_col
        Optional explicit column names for the GEO CSV. Auto-detected when
        not provided.
    """
    if raw_count_mode:
        raise ValueError(
            "raw_count_mode=True is not supported for v1 "
            "(see DECISIONS.md correction 2026-05-20 (5/7))."
        )

    logger.error(
        "Garrido-Trigo loader: CELLxGENE matrix path is superseded by "
        "DECISIONS correction 2026-06-04 (9). The CELLxGENE deposit's "
        "obs.index is synthetic (cell1, cell2, ...) — barcode join cannot "
        "succeed against the GEO CSV. Use the RAW.tar loader once it "
        "lands; do not trust any output from this code path."
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

    if annotation_csv_path is None:
        logger.warning(
            "No GEO annotation CSV supplied; falling back to CELLxGENE "
            "5-CL broad labels with NO fine tier (degraded mode). Pass "
            "annotation_csv_path=GSE214695_cell_annotation.csv for full "
            "broad+fine."
        )
        obs["cell_type_broad"] = obs["cell_type"].astype("category")
        obs["cell_type_fine"] = pd.Series(pd.NA, index=obs.index, dtype="object")
    else:
        ann, bcol, acol, dcol = _load_annotation_csv(
            annotation_csv_path, barcode_col, annotation_col, donor_col,
        )
        joined, strategy = _try_join_keys(obs, ann, bcol, dcol)
        joined.index = obs.index
        logger.info("Barcode join succeeded via: %s", strategy)

        # Cross-check: CELLxGENE obs['disease'] must agree with the GEO sample
        # prefix (HC*/UC*) for every cell. If they disagree, filter-before-join
        # is silently trusting one source's definition of "UC" over the other.
        if dcol is not None:
            geo_prefix = (
                joined[dcol].astype(str)
                .str.extract(r"^([A-Za-z]+)", expand=False)
                .str.upper()
            )
            prefix_to_disease = {"HC": "normal", "UC": "ulcerative colitis"}
            unknown = sorted(set(geo_prefix.dropna()) - set(prefix_to_disease))
            if unknown:
                raise ValueError(
                    f"Garrido-Trigo loader: unrecognized GEO sample prefixes "
                    f"{unknown} in column {dcol!r}; expected HC*/UC* only "
                    f"after the HC+UC filter. Did the filter-before-join "
                    f"step leak CD donors, or has the sample-naming scheme "
                    f"drifted?"
                )
            expected = geo_prefix.map(prefix_to_disease)
            actual = obs["disease"].astype(str)
            mismatch = expected.values != actual.values
            if mismatch.any():
                n_mis = int(mismatch.sum())
                ex_idx = obs.index[mismatch][:5].tolist()
                ex_actual = actual.values[mismatch][:5].tolist()
                ex_expected = expected.values[mismatch][:5].tolist()
                raise ValueError(
                    f"Garrido-Trigo loader: {n_mis} cells disagree on disease "
                    f"between CELLxGENE obs['disease'] and the GEO sample "
                    f"prefix in {dcol!r}. The HC+UC filter and the annotation "
                    f"source are using different definitions of who's UC. "
                    f"Examples (obs.index / actual / expected): "
                    f"{list(zip(ex_idx, ex_actual, ex_expected))}."
                )
        else:
            logger.warning(
                "GEO CSV has no donor/sample column; skipping the "
                "disease/sample-prefix cross-check."
            )

        fine_raw = joined[acol].map(_normalize_label)
        # Collapse Ribhi clusters into their parent before any tier logic.
        fine_collapsed = fine_raw.replace(RIBHI_TO_PARENT)
        n_ribhi_collapsed = int((fine_raw != fine_collapsed).sum())
        logger.info(
            "Collapsed %d Ribhi cells into parent fine clusters", n_ribhi_collapsed
        )

        # Completeness: every joined cell must carry an annotation. The join
        # itself only guarantees barcodes matched; the matched CSV row could
        # still have NaN in the annotation column. The unmapped check below
        # explicitly drops NaN, so without this gate, NA annotations slip
        # through silently and the cell is later scored with no fine label.
        n_missing_fine = int(fine_collapsed.isna().sum())
        if n_missing_fine:
            raise ValueError(
                f"Garrido-Trigo loader: {n_missing_fine}/{len(fine_collapsed)} "
                f"cells have NaN fine annotation after the barcode join. The "
                f"join hit every barcode but at least one matched CSV row had "
                f"no value in column {acol!r}. Inspect the GEO CSV — partial "
                f"annotation breaks every downstream tier."
            )

        unmapped = sorted(
            set(fine_collapsed.dropna().unique()) - set(FINE_TO_BROAD.keys())
        )
        if unmapped:
            raise KeyError(
                f"Garrido-Trigo loader: {len(unmapped)} fine labels have no "
                f"entry in FINE_TO_BROAD: {unmapped}. Extend the map "
                f"(load_garrido_trigo.FINE_TO_BROAD)."
            )

        broad = fine_collapsed.map(FINE_TO_BROAD)
        n_broad = broad.dropna().nunique()
        n_fine = fine_collapsed.dropna().nunique()
        logger.info(
            "Tier cardinalities: fine=%d (post-Ribhi-collapse), broad=%d",
            n_fine, n_broad,
        )
        if not (10 <= n_broad <= 15):
            logger.warning(
                "Broad-tier cardinality %d is outside the v1 10-15 target; "
                "review FINE_TO_BROAD grouping.", n_broad,
            )

        obs["cell_type_fine"] = fine_collapsed.astype("category")
        obs["cell_type_broad"] = broad.astype("category")

    obs["donor"] = obs["donor_id"].astype("category")
    if "tissue" not in obs.columns:
        obs["tissue"] = "colonic mucosa"
    obs["batch"] = pd.Series(pd.NA, index=obs.index, dtype="object")
    adata.obs = obs

    if apply_v1_filter:
        # Donor structure is a hard invariant: 12 donors, 6 HC + 6 UC. This
        # is fixed by the GEO study design and does not drift with re-pulls
        # or QC nudges; a violation means the filter or annotation is wrong.
        donors_by_disease = (
            adata.obs[["donor", "disease"]]
            .drop_duplicates()
            .groupby("disease", observed=True)
            .size()
        )
        n_donors = int(adata.obs["donor"].nunique())
        n_hc = int(donors_by_disease.get("normal", 0))
        n_uc = int(donors_by_disease.get("ulcerative colitis", 0))
        if (
            n_donors != EXPECTED_UC_SUBSET_N_DONORS
            or n_hc != 6
            or n_uc != 6
        ):
            raise ValueError(
                f"Garrido-Trigo loader: donor-structure invariant violated. "
                f"Got n_donors={n_donors} ({n_hc} HC + {n_uc} UC); expected "
                f"{EXPECTED_UC_SUBSET_N_DONORS} (6 HC + 6 UC). Per-disease "
                f"donor breakdown: {dict(donors_by_disease)}."
            )

        # Cell count is a derived intersection (GEO ∩ CELLxGENE post-QC) and
        # can drift for benign reasons (re-pull, QC nudge). Tripwire only —
        # the hard gates above (completeness, donor structure, no orphans,
        # no unmapped fine labels) are what mean "the join is broken."
        n_cells = adata.n_obs
        if n_cells != EXPECTED_UC_SUBSET_N_CELLS:
            logger.warning(
                "UC-subset cell count %d != expected %d "
                "(DECISIONS.md correction 2/7). Tripwire only; investigate "
                "if this drifts unexpectedly, but the hard gates have "
                "already passed.", n_cells, EXPECTED_UC_SUBSET_N_CELLS,
            )

    logger.info("Post-filter cell count: %d", adata.n_obs)
    per_donor = adata.obs["donor"].value_counts()
    logger.info("Donors (n=%d):", per_donor.size)
    for donor, count in per_donor.items():
        logger.info("  %s: %d cells", donor, count)

    adata = ensembl_to_hgnc(adata)
    return adata
