"""Loader for Garrido-Trigo 2023 (Atlas 2, UC subset) — RAW.tar matrix path.

Per DECISIONS correction 2026-06-04 (9), the matrix source is the GEO
supplementary archive ``GSE214695_RAW.tar`` (real 10X barcodes
preserved), NOT the CELLxGENE deposit ``b1a62801-...h5ad``, whose
``obs.index`` was synthetic and stripped of the original barcodes —
making any external-annotation barcode join undoable.

Sources
-------
- **Matrix:** ``GSE214695_RAW.tar`` from GEO. 18 per-GSM 10X v2
  triplets (``GSM{nnnnnnn}_{HC|UC|CD}-{n}_{matrix.mtx,barcodes.tsv,features.tsv}.gz``).
  The 12 HC/UC GSMs are the v1 set; CD-1..CD-6 are skipped at the
  glob step. Matrices are gene-sorted (genes × cells, 33,538 features
  × 737,280 cell barcodes — the full 10X v2 whitelist; the loader
  transposes and inner-joins against the annotation CSV to drop the
  ~99.8% empty droplets).
- **Annotation:** ``GSE214695_cell_annotation.csv`` from GEO (Salas-lab
  91-cluster fine annotation, with ``sample`` + ``cell_id`` columns).

Join
----
Composite-on-composite, deterministic:
- **RAW side**, per cell: ``f"{sample}_{barcode}"`` (e.g.
  ``"HC1_AAACCTGCAAGTCTGT-1"``).
- **CSV side**, per row:  ``f"{sample}_{cell_id}"`` (same shape).

The seven barcode-join strategies from the pre-correction-9 loader are
gone — when we control matrix assembly, the composite is unique on both
sides by construction. Inner-join keeps only CSV-annotated cells, which
drops empty droplets and lands at ~30,068 HC+UC cells (DECISIONS 2/7).

Three sample-naming conventions are reconciled on load:

================  =========  ===========================
form              example    where seen
================  =========  ===========================
dash              ``HC-1``   RAW.tar filenames
no separator      ``HC1``    CSV ``sample`` column
underscore        ``HC_1``   prior CELLxGENE ``donor_id``
================  =========  ===========================

The loader normalizes RAW's dashed form to the CSV's no-separator form
for the join key, and emits the underscore form in ``obs['donor_id']``
to preserve continuity with the pre-correction-9 obs schema.

Counts and genes
----------------
- ``X`` is **raw integer counts** (MatrixMarket ``coordinate integer``).
  The loader applies ``log1p(CP10k)`` so the input state matches the
  other atlases (DECISIONS 5/7); raw counts preserved in
  ``layers['counts']``. ``raw_count_mode=True`` is unsupported.
- ``var_names`` = Ensembl IDs (column 1 of ``features.tsv.gz``);
  ``var['feature_name']`` = HGNC symbol (column 2). The final
  ``ensembl_to_hgnc`` step (pinned per correction 11) takes the
  ``feature_name`` path — same as the prior CELLxGENE deposit.

Donor structure (hard invariant)
--------------------------------
12 donors, 6 HC + 6 UC, fixed by study design. The loader hard-asserts
this; the 30,068 cell count remains a soft tripwire per the correction
9 framing.

References: ``code/02_atlas_prep/atlas_schemas.md``;
DECISIONS corrections 2026-05-20 (2/7), (4/7), (5/7); 2026-06-03 (8);
2026-06-04 (9) [this loader implements (9)]; 2026-06-04 (11) [HGNC pin].
"""

from __future__ import annotations

import gzip
import logging
import re
import tarfile
from pathlib import Path

import anndata as ad
from anndata import AnnData
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io

from hgnc_remap import ensembl_to_hgnc

logger = logging.getLogger(__name__)

# CD samples are skipped at the glob step for v1 (DECISIONS 2/7).
V1_DONOR_PREFIXES: tuple[str, ...] = ("HC", "UC")

# Hard invariants (DECISIONS 2/7).
EXPECTED_UC_SUBSET_N_CELLS = 30_068  # soft tripwire
EXPECTED_UC_SUBSET_N_DONORS = 12     # hard
EXPECTED_N_HC = 6
EXPECTED_N_UC = 6

# Disease harmonization from sample prefix.
DISEASE_MAP = {"HC": "normal", "UC": "ulcerative colitis"}

# Tar entry filename pattern, e.g. ``GSM6614348_HC-1_matrix.mtx.gz``.
_TAR_ENTRY_RE = re.compile(
    r"^(?P<gsm>GSM\d+)_(?P<prefix>HC|UC|CD)-(?P<n>\d+)"
    r"_(?P<kind>matrix|barcodes|features)\.(?:mtx|tsv)\.gz$"
)

# Required columns in GSE214695_cell_annotation.csv.
_CSV_REQUIRED_COLS = ("sample", "cell_id", "annotation")

# Canonical broad vocabulary — loader-local until CANONICAL_BROAD is
# formally locked (see code/_shared/canonical_broad_DRAFT.md). The
# loader runs two assertions against this set:
#   (1) At module load: every value in FINE_TO_BROAD must be a member.
#       Catches typos in the map at authoring time (e.g. the live F8
#       preview "Perycites" -> if the value side had been typo'd, this
#       would have fired at import, not at step 06).
#   (2) At end of load(): every distinct value of obs['cell_type_broad']
#       must be a member. Catches anything that slips past (1) — e.g.
#       a stringified NaN landing in the column via an unexpected path.
# Kept loader-local rather than imported from code/_shared/ to avoid
# premature coupling against a vocab that may shift on lock; promote
# to code/_shared/ when CANONICAL_BROAD locks.
_BROAD_VOCAB: frozenset[str] = frozenset({
    "B cell",
    "NK/ILC",
    "T cell",
    "colonocyte",
    "dendritic cell",
    "endothelium",
    "enteroendocrine/tuft",
    "epithelial progenitor",
    "fibroblast",
    "goblet",
    "granulocyte",
    "mast cell",
    "monocyte/macrophage",
    "mural/glia",
    "plasma cell",
})


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

# Fine -> broad roll-up. Targets the v1 ~10-15 broad-tier budget. Keys are
# the 91 published Garrido-Trigo fine labels plus the generic parent labels
# introduced by RIBHI_TO_PARENT (epithelial, T, fibroblast, mast).
# Whitespace is normalized before lookup (see _normalize_label).
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
    "Pericytes":             "mural/glia",
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

# Gate (1): validate the map's value side against the canonical vocab
# at module import. A typo on the value side would otherwise survive
# every per-cell assertion and land as a phantom broad category at
# step 06 — read as discordance, not as the typo it is.
_unmapped_broad = set(FINE_TO_BROAD.values()) - _BROAD_VOCAB
if _unmapped_broad:
    raise ValueError(
        f"load_garrido_trigo.FINE_TO_BROAD ships broad values outside "
        f"_BROAD_VOCAB: {sorted(_unmapped_broad)}. Likely a typo on the "
        f"value side of the map (see canonical_broad_DRAFT.md for the "
        f"15-term vocabulary)."
    )
del _unmapped_broad


def _normalize_label(value: object) -> object:
    """Strip + collapse internal whitespace; normalize curly to ASCII.

    The GEO CSV contains at least one label with a literal double space
    ('PC  immediate early response'), and Salas-lab labels use Latin-1
    characters that can round-trip as mojibake (e.g., 'Na\xefve B cell').
    Normalizing on load prevents silent join breakage.
    """
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return value
    s = str(value)
    s = s.replace("ï", "i").replace("\xef", "i")  # Naïve -> Naive
    s = s.replace("é", "e")                       # café -> cafe
    s = re.sub(r"\s+", " ", s).strip()
    return s


def _group_tar_entries(
    tar: tarfile.TarFile, keep_prefixes: tuple[str, ...]
) -> dict[str, dict[str, str]]:
    """Walk the tar and group entries by GSM_SAMPLE -> {kind: filename}.

    Filters out CD samples (or any prefix not in ``keep_prefixes``) at
    the glob step. Raises if any kept GSM is missing one of the three
    triplet files.
    """
    groups: dict[str, dict[str, str]] = {}
    for name in tar.getnames():
        m = _TAR_ENTRY_RE.match(name)
        if m is None:
            logger.warning("Tar entry does not match expected pattern: %s", name)
            continue
        prefix = m.group("prefix")
        if prefix not in keep_prefixes:
            continue
        key = f"{m.group('gsm')}_{prefix}-{m.group('n')}"
        groups.setdefault(key, {})[m.group("kind")] = name

    for key, files in groups.items():
        missing = {"matrix", "barcodes", "features"} - set(files)
        if missing:
            raise ValueError(
                f"Garrido-Trigo RAW.tar: GSM {key!r} is missing required "
                f"files {sorted(missing)}. Tar may be truncated or the "
                f"deposit shape has changed."
            )
    return groups


def _read_triplet_from_tar(
    tar: tarfile.TarFile, files: dict[str, str]
) -> tuple[scipy.io.matrix, list[str], pd.DataFrame]:
    """Read one GSM's (matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz)
    in-memory from the open tar.

    Returns ``(mat_genes_x_cells, barcodes, features_df)`` where
    ``features_df`` has columns ``ensembl``, ``symbol``, ``feature_type``.
    """
    with tar.extractfile(files["matrix"]) as raw, gzip.GzipFile(fileobj=raw) as gz:
        mat = scipy.io.mmread(gz).tocsr()  # genes x cells

    with tar.extractfile(files["barcodes"]) as raw, gzip.GzipFile(fileobj=raw) as gz:
        barcodes = [line.decode("utf-8").strip() for line in gz if line.strip()]

    with tar.extractfile(files["features"]) as raw, gzip.GzipFile(fileobj=raw) as gz:
        features = pd.read_csv(
            gz, sep="\t", header=None,
            names=["ensembl", "symbol", "feature_type"],
            dtype=str,
        )

    if mat.shape != (len(features), len(barcodes)):
        raise ValueError(
            f"Garrido-Trigo triplet shape mismatch: matrix {mat.shape} "
            f"vs features={len(features)} x barcodes={len(barcodes)}."
        )
    return mat, barcodes, features


def _load_per_gsm(
    tar: tarfile.TarFile, files: dict[str, str], sample_dash: str
) -> AnnData:
    """Build a per-GSM AnnData (cells × genes) with the composite obs key.

    ``sample_dash`` is the RAW.tar form, e.g. ``"HC-1"``. The function
    normalizes this to the no-separator (CSV-compatible) ``"HC1"`` for
    the composite barcode and the obs ``sample`` column, and to the
    underscore (CELLxGENE-deposit) ``"HC_1"`` for ``donor_id``.
    """
    mat, barcodes, features = _read_triplet_from_tar(tar, files)
    sample_no_sep = sample_dash.replace("-", "")    # "HC1"
    sample_underscore = sample_dash.replace("-", "_")  # "HC_1"
    prefix = sample_dash.split("-")[0]              # "HC"

    composite_index = pd.Index(
        [f"{sample_no_sep}_{b}" for b in barcodes],
        name="cell_id_composite",
    )
    obs = pd.DataFrame(index=composite_index)
    obs["sample"] = sample_no_sep
    obs["donor_id"] = sample_underscore
    obs["disease"] = DISEASE_MAP[prefix]

    var = pd.DataFrame(
        {
            "feature_name": features["symbol"].to_numpy(),
            "feature_type": features["feature_type"].to_numpy(),
        },
        index=pd.Index(features["ensembl"], name="ensembl_id"),
    )

    adata = ad.AnnData(X=mat.T.tocsr().astype(np.float32), obs=obs, var=var)
    return adata


def _load_annotation_csv(
    annotation_csv_path: Path, keep_prefixes: tuple[str, ...],
) -> pd.DataFrame:
    """Read the GEO annotation CSV, filter to the v1 cohort, build the
    composite join key, validate.

    Returns the DataFrame reindexed by ``f"{sample}_{cell_id}"`` (the
    CSV-side composite that matches the RAW-side ``obs.index`` built in
    ``_load_per_gsm``).

    The filter to ``keep_prefixes`` happens **before** the duplicate-key
    check. The CSV has known CD-only duplicates — same cell barcode
    appears twice for one CD3 donor with conflicting fine annotations
    (Salas-lab authoring artifact, traced to inconsistent whitespace in
    the ``Unnamed: 0`` composite). Filtering first means the gate fires
    only on duplicates within the cohort we actually load.
    """
    logger.info("Reading GEO annotation CSV: %s", annotation_csv_path)
    ann = pd.read_csv(annotation_csv_path)
    missing = set(_CSV_REQUIRED_COLS) - set(ann.columns)
    if missing:
        raise KeyError(
            f"Garrido-Trigo annotation CSV: missing required columns "
            f"{sorted(missing)}. Got: {list(ann.columns)}."
        )
    ann["sample"] = ann["sample"].astype(str).str.strip()
    ann["cell_id"] = ann["cell_id"].astype(str).str.strip()
    ann["annotation"] = ann["annotation"].map(_normalize_label)

    n_pre = len(ann)
    sample_prefix = ann["sample"].str.extract(r"^([A-Za-z]+)", expand=False)
    keep_mask = sample_prefix.isin(keep_prefixes)
    ann = ann[keep_mask].copy()
    logger.info(
        "Annotation CSV: %d rows total; %d kept after sample-prefix filter "
        "(keep=%s)", n_pre, len(ann), list(keep_prefixes),
    )

    composite = ann["sample"] + "_" + ann["cell_id"]
    if composite.duplicated().any():
        n_dup = int(composite.duplicated().sum())
        dup_examples = composite[composite.duplicated(keep=False)].head(6).tolist()
        raise ValueError(
            f"Garrido-Trigo annotation CSV: {n_dup} duplicate "
            f"(sample, cell_id) composites in the kept cohort; cannot "
            f"use as join key. Examples: {dup_examples}. The CSV has "
            f"known CD-only duplicates from a Salas-lab authoring bug; "
            f"if this fires on HC/UC, that's a new defect to inspect."
        )
    ann = ann.set_index(pd.Index(composite, name="cell_id_composite"))
    return ann


def load(
    raw_tar_path: str,
    annotation_csv_path: str,
    apply_v1_filter: bool = True,
    raw_count_mode: bool = False,
) -> AnnData:
    """Load Garrido-Trigo from GSE214695_RAW.tar + GSE214695_cell_annotation.csv.

    Parameters
    ----------
    raw_tar_path
        Path to ``GSE214695_RAW.tar`` (the GEO supplementary archive).
        18 per-GSM 10X triplets; the loader reads in-memory via
        ``tarfile`` + ``gzip`` (no extraction).
    annotation_csv_path
        Path to ``GSE214695_cell_annotation.csv`` (the Salas-lab
        91-cluster fine annotation).
    apply_v1_filter
        If True (default), keep only HC + UC GSMs (skip CD at glob time)
        and run the donor-structure invariant (6 HC + 6 UC). If False,
        loads all 18 GSMs and skips the invariant — debug/inspection
        path only.
    raw_count_mode
        Must remain False for v1 (DECISIONS 5/7). The RAW.tar matrices
        are raw counts; the loader normalizes them to log1p(CP10k);
        ``raw_count_mode=True`` (leave X as raw counts) is unsupported.

    Returns
    -------
    AnnData
        cells × genes; ``X`` = log1p(CP10k) float; raw counts in
        ``layers['counts']``. obs schema: ``cell_type_fine``,
        ``cell_type_broad``, ``donor``, ``donor_id``, ``disease``,
        ``sample``, ``batch``, ``tissue``. var schema: index = Ensembl
        ID (input); ``var['feature_name']`` = HGNC symbol;
        ``var['feature_type']``. The final ``ensembl_to_hgnc`` step
        replaces ``var_names`` with HGNC symbols.
    """
    if raw_count_mode:
        raise ValueError(
            "raw_count_mode=True is not supported for v1 (DECISIONS "
            "correction 5/7): all atlases use uniform log1p(CP10k) "
            "input. The Garrido-Trigo RAW.tar matrices are raw counts "
            "and are normalized on load."
        )

    raw_tar_path = Path(raw_tar_path).expanduser()
    annotation_csv_path = Path(annotation_csv_path).expanduser()
    for p in (raw_tar_path, annotation_csv_path):
        if not p.exists():
            raise FileNotFoundError(f"Garrido-Trigo loader: input not found: {p}")

    keep = V1_DONOR_PREFIXES if apply_v1_filter else ("HC", "UC", "CD")
    expected_gsm_count = (
        EXPECTED_N_HC + EXPECTED_N_UC if apply_v1_filter else 18
    )

    # ---- 1. Enumerate tar; group entries by GSM_SAMPLE; verify triplets ----
    logger.info("Opening tar: %s", raw_tar_path)
    with tarfile.open(raw_tar_path, "r:") as tar:
        groups = _group_tar_entries(tar, keep_prefixes=keep)
        if len(groups) != expected_gsm_count:
            raise ValueError(
                f"Garrido-Trigo RAW.tar: expected {expected_gsm_count} GSM "
                f"groups for prefixes {keep}; got {len(groups)} "
                f"({sorted(groups)})."
            )

        # ---- 2. Load each GSM as a per-sample AnnData (composite-keyed) ----
        per_gsm: list[AnnData] = []
        for key in sorted(groups):
            files = groups[key]
            sample_dash = key.split("_", 1)[1]  # "HC-1"
            adata_g = _load_per_gsm(tar, files, sample_dash)
            logger.info(
                "  %s: %d barcodes x %d features",
                key, adata_g.n_obs, adata_g.n_vars,
            )
            per_gsm.append(adata_g)

    # ---- 3. Concat. Gene sets should be identical (single CellRanger ref). ----
    gene_sets = [set(a.var_names) for a in per_gsm]
    shared = set.intersection(*gene_sets)
    union = set.union(*gene_sets)
    if len(shared) != len(union):
        logger.warning(
            "Gene sets differ across GSMs: shared=%d union=%d. Outer-joining; "
            "investigate if this is unexpected (single CellRanger reference "
            "should produce identical gene lists across the deposit).",
            len(shared), len(union),
        )
    adata = ad.concat(
        per_gsm, axis=0, join="outer", fill_value=0, merge="first",
    )
    if not adata.obs_names.is_unique:
        n_dup = int(adata.obs_names.duplicated().sum())
        raise ValueError(
            f"Garrido-Trigo loader: {n_dup} duplicate composite obs keys "
            f"after concat. Sample prefix may not have made the keys unique."
        )
    logger.info(
        "Concatenated %d GSMs: %d cell-barcodes x %d genes (pre-CSV-join)",
        len(per_gsm), adata.n_obs, adata.n_vars,
    )

    # ---- 4. Inner-join against the annotation CSV's composite key. ----
    ann = _load_annotation_csv(annotation_csv_path, keep_prefixes=keep)

    in_csv = adata.obs_names.isin(ann.index)
    n_pre = adata.n_obs
    n_drop = int((~in_csv).sum())
    adata = adata[in_csv].copy()
    logger.info(
        "CSV inner-join: dropped %d empty droplets / unannotated barcodes; "
        "kept %d cells", n_drop, adata.n_obs,
    )
    if adata.n_obs == 0:
        raise ValueError(
            "Garrido-Trigo loader: inner join produced zero cells. The "
            "composite keys did not match — inspect sample-name normalization."
        )

    # Pull the aligned annotation rows onto the kept cells.
    aligned = ann.loc[adata.obs_names]

    # ---- 5. Disease / sample-prefix agreement. The CSV's `sample` column
    # must agree with the per-GSM `disease` we set from the filename prefix.
    # If they disagree, the filename convention and the annotation source
    # are using different definitions of who's HC vs UC.
    geo_prefix = (
        aligned["sample"].astype(str)
        .str.extract(r"^([A-Za-z]+)", expand=False).str.upper()
    )
    unknown = sorted(set(geo_prefix.dropna()) - set(DISEASE_MAP))
    if unknown:
        raise ValueError(
            f"Garrido-Trigo loader: unrecognized GEO sample prefixes "
            f"{unknown} in CSV; expected only {sorted(DISEASE_MAP)} for v1. "
            f"Did filter-before-join leak CD samples, or has the naming "
            f"scheme drifted?"
        )
    expected_disease = geo_prefix.map(DISEASE_MAP)
    actual_disease = adata.obs["disease"].astype(str)
    mismatch = expected_disease.values != actual_disease.values
    if mismatch.any():
        n_mis = int(mismatch.sum())
        ex_idx = adata.obs_names[mismatch][:5].tolist()
        ex_act = actual_disease.values[mismatch][:5].tolist()
        ex_exp = expected_disease.values[mismatch][:5].tolist()
        raise ValueError(
            f"Garrido-Trigo loader: {n_mis} cells disagree on disease "
            f"between the RAW.tar GSM filename prefix and the CSV `sample` "
            f"prefix. The filename convention and the annotation source "
            f"have different definitions of who's UC. Examples "
            f"(obs.index / actual / expected): "
            f"{list(zip(ex_idx, ex_act, ex_exp))}."
        )

    # ---- 6. Fine-label processing: RIBHI collapse, completeness, unmapped. ----
    fine_raw = aligned["annotation"].map(_normalize_label)
    fine_collapsed = fine_raw.replace(RIBHI_TO_PARENT)
    n_ribhi = int((fine_raw != fine_collapsed).sum())
    logger.info("Collapsed %d Ribhi cells into parent fine clusters", n_ribhi)

    n_missing_fine = int(fine_collapsed.isna().sum())
    if n_missing_fine:
        raise ValueError(
            f"Garrido-Trigo loader: {n_missing_fine}/{len(fine_collapsed)} "
            f"cells have NaN fine annotation after the CSV join. Inspect "
            f"the CSV — partial annotation breaks every downstream tier."
        )

    unmapped = sorted(set(fine_collapsed.unique()) - set(FINE_TO_BROAD))
    if unmapped:
        raise KeyError(
            f"Garrido-Trigo loader: {len(unmapped)} fine labels have no "
            f"FINE_TO_BROAD entry: {unmapped}. Extend the map "
            f"(load_garrido_trigo.FINE_TO_BROAD)."
        )

    broad = fine_collapsed.map(FINE_TO_BROAD)

    # Gate (2): every emitted broad label must be in the canonical vocab.
    # The map-side gate (1) catches typos in FINE_TO_BROAD values at
    # import; this gate catches anything that lands a non-canonical
    # string in the broad column via a path other than the map (e.g.
    # stringified NaN). Hard-raise so a typo cannot quietly become a
    # phantom non-matching category at step 06's string intersection.
    emitted = set(broad.dropna().unique())
    unrecognized = emitted - _BROAD_VOCAB
    if unrecognized:
        raise ValueError(
            f"Garrido-Trigo loader: emitted cell_type_broad values "
            f"{sorted(unrecognized)} are not in the canonical vocab. "
            f"Either FINE_TO_BROAD's value side is broken (see gate 1) "
            f"or a non-map path is populating the broad column."
        )

    n_broad = int(broad.nunique())
    n_fine = int(fine_collapsed.nunique())
    logger.info(
        "Tier cardinalities: fine=%d (post-Ribhi-collapse), broad=%d",
        n_fine, n_broad,
    )
    if not (10 <= n_broad <= 15):
        logger.warning(
            "Broad-tier cardinality %d outside the v1 10-15 target; "
            "review FINE_TO_BROAD grouping.", n_broad,
        )

    # ---- 7. Finalize obs schema. ----
    obs = adata.obs
    obs["cell_type_fine"] = fine_collapsed.astype("category")
    obs["cell_type_broad"] = broad.astype("category")
    obs["donor"] = obs["donor_id"].astype("category")
    obs["batch"] = obs["sample"].astype("category")
    obs["tissue"] = "colonic mucosa"
    adata.obs = obs

    # ---- 8. Donor-structure invariant (hard) + cell-count tripwire (soft). ----
    if apply_v1_filter:
        donor_disease = (
            adata.obs[["donor", "disease"]].drop_duplicates()
            .groupby("disease", observed=True).size()
        )
        n_donors = int(adata.obs["donor"].nunique())
        n_hc = int(donor_disease.get("normal", 0))
        n_uc = int(donor_disease.get("ulcerative colitis", 0))
        if (
            n_donors != EXPECTED_UC_SUBSET_N_DONORS
            or n_hc != EXPECTED_N_HC
            or n_uc != EXPECTED_N_UC
        ):
            raise ValueError(
                f"Garrido-Trigo loader: donor-structure invariant violated. "
                f"Got n_donors={n_donors} ({n_hc} HC + {n_uc} UC); expected "
                f"{EXPECTED_UC_SUBSET_N_DONORS} "
                f"({EXPECTED_N_HC} HC + {EXPECTED_N_UC} UC). Per-disease "
                f"donor breakdown: {dict(donor_disease)}."
            )

        n_cells = adata.n_obs
        if n_cells != EXPECTED_UC_SUBSET_N_CELLS:
            logger.warning(
                "UC-subset cell count %d != expected %d (DECISIONS 2/7). "
                "Tripwire only — derived intersection, drifts on benign "
                "re-pulls; the hard gates above have already passed.",
                n_cells, EXPECTED_UC_SUBSET_N_CELLS,
            )

    logger.info("Post-join cell count: %d", adata.n_obs)
    per_donor = adata.obs["donor"].value_counts()
    for donor, count in per_donor.items():
        logger.info("  %s: %d cells", donor, count)

    # ---- 9. Raw counts -> log1p(CP10k); preserve raw in layer. ----
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    logger.info("Normalized X to log1p(CP10k); raw counts in layers['counts'].")

    # ---- 10. Ensembl -> HGNC via the pinned remap (correction 11). ----
    adata = ensembl_to_hgnc(adata)

    logger.info(
        "Garrido-Trigo load complete: %d cells x %d genes, %d donors.",
        adata.n_obs, adata.n_vars, int(adata.obs["donor"].nunique()),
    )
    return adata
