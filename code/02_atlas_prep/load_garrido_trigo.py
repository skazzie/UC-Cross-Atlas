"""Loader for Garrido-Trigo 2023 (Atlas 2, UC subset).

Sources
-------
- **Matrix:** CELLxGENE deposit ``b1a62801-f509-45f8-b55f-533fbb7e7800.h5ad``
  (log-normalized X; var_names are Ensembl IDs; HGNC symbols in
  ``var['feature_name']``).
- **Labels:** GEO supplementary file ``GSE214695_cell_annotation.csv``
  (91-label fine annotation from the Salas lab; not present in the
  CELLxGENE deposit).

The CELLxGENE-only path is deprecated: it shipped a 5-CL-lineage broad
label with no fine tier, blocking fine-tier cross-atlas concordance for
Garrido-Trigo. The GEO supplementary annotation restores the full
broad + fine tier. See DECISIONS.md correction (4/7) and the correction
reversing it.

References: ``code/02_atlas_prep/atlas_schemas.md``;
DECISIONS.md corrections 2026-05-20 (2/7), (4/7), (5/7), and the
correction reversing (4/7).
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
    "cell_id", "cell", "Cell", "barcode", "Barcode", "cellID", "CellID",
    "index", "Unnamed: 0",
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
    df: pd.DataFrame, candidates: tuple[str, ...], purpose: str
) -> str:
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(
        f"GSE214695 annotation CSV: could not auto-detect {purpose} column. "
        f"Tried {list(candidates)}. Got columns: {list(df.columns)}. "
        f"Pass the explicit column name via the loader's "
        f"`{purpose}_col` argument."
    )


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

    bcol = barcode_col or _autodetect_column(ann, _BARCODE_COL_CANDIDATES, "barcode")
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
        donor = obs["donor_id"].astype(str)
        candidates.append((
            "donor_id + '_' + obs.index",
            (donor + "_" + obs_index).astype(str),
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

        fine_raw = joined[acol].map(_normalize_label)
        # Collapse Ribhi clusters into their parent before any tier logic.
        fine_collapsed = fine_raw.replace(RIBHI_TO_PARENT)
        n_ribhi_collapsed = int((fine_raw != fine_collapsed).sum())
        logger.info(
            "Collapsed %d Ribhi cells into parent fine clusters", n_ribhi_collapsed
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
        n_cells = adata.n_obs
        n_donors = int(adata.obs["donor"].nunique())
        if n_cells != EXPECTED_UC_SUBSET_N_CELLS:
            logger.warning(
                "UC-subset cell count %d != expected %d "
                "(DECISIONS.md correction 2/7). Investigate before trusting "
                "the join.", n_cells, EXPECTED_UC_SUBSET_N_CELLS,
            )
        if n_donors != EXPECTED_UC_SUBSET_N_DONORS:
            logger.warning(
                "UC-subset donor count %d != expected %d "
                "(DECISIONS.md correction 2/7).",
                n_donors, EXPECTED_UC_SUBSET_N_DONORS,
            )

    logger.info("Post-filter cell count: %d", adata.n_obs)
    per_donor = adata.obs["donor"].value_counts()
    logger.info("Donors (n=%d):", per_donor.size)
    for donor, count in per_donor.items():
        logger.info("  %s: %d cells", donor, count)

    adata = ensembl_to_hgnc(adata)
    return adata
