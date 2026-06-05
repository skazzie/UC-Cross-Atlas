"""Loader for Smillie 2019 UC core atlas (Atlas 1, Single Cell Portal SCP259).

Canonical source is **SCP259** (Broad Single Cell Portal), NOT the CELLxGENE
deposit (e373cf41-...), which is a 34,772-cell healthy-epithelial-only subset
and is unusable for v1. See DECISIONS.md correction 2026-05-20 (7/7).

Source layout (after the SCP bulk download)
--------------------------------------------
``<scp259_dir>/`` contains, under hashed sub-directories:

- **Three compartment matrices**, each a gene-sorted 10X triplet::

      expression/<hash>/gene_sorted-Epi.matrix.mtx   (+ Epi.genes.tsv, Epi.barcodes2.tsv)
      expression/<hash>/gene_sorted-Imm.matrix.mtx   (+ Imm.genes.tsv, Imm.barcodes2.tsv)
      expression/<hash>/gene_sorted-Fib.matrix.mtx   (+ Fib.genes.tsv, Fib.barcodes2.tsv)

  The hash directory names are not stable, so the loader globs by the
  ``gene_sorted-<compartment>.matrix.mtx`` filename rather than hard-coding
  them.  **The matrices are gene-sorted = genes x cells**; the loader
  transposes each to cells x genes before assembling the AnnData.

- **One metadata file**::

      metadata/all.meta2.txt

  Tab-separated, with an SCP-specific second header row (``TYPE``) that must
  be skipped.  Columns (captured 2026-06-04):

  ===========  =======  ===================================================
  column       TYPE     meaning
  ===========  =======  ===================================================
  NAME         group    cell id, ``Subject.Sample.barcode`` — join key,
                        identical to the matrix ``*.barcodes2.tsv`` values
  Cluster      group    fine cell type (51 labels) -> ``cell_type_fine``
  nGene        numeric  genes detected -> ``n_genes``
  nUMI         numeric  UMIs -> ``n_counts``
  Subject      group    donor -> ``donor``
  Health       group    Healthy / Inflamed / Non-inflamed -> ``health``
  Location     group    Epi / Imm / Fib -> ``compartment``
  Sample       group    biopsy id -> ``sample`` (also used as ``batch``)
  ===========  =======  ===================================================

Counts / genes
--------------
- ``X`` ships as **raw integer counts** (MatrixMarket ``coordinate integer``).
  Unlike Garrido/Pan-GI/HCA (which ship CELLxGENE log-normalized X), this
  loader applies ``log1p(CP10k)`` itself so the input is cross-atlas
  comparable.  ``raw_count_mode`` must stay ``False`` for v1 (DECISIONS
  correction 5/7) — uniform log1p(CP10k) across atlases.  Raw counts are
  preserved in ``layers['counts']``.
- ``var_names`` are **HGNC symbols** (single-column ``*.genes.tsv``), not
  Ensembl IDs.  The final ``ensembl_to_hgnc`` step takes its symbol-fallback
  path (no ``feature_name`` column): it dedups duplicate symbols by
  max-expression and filters to the NCBI-authoritative symbol set — the same
  gene universe every other atlas converges to, so no separate symbol<->Ensembl
  mapping is needed.

Donor structure (hard invariant, fixed by study design)
-------------------------------------------------------
30 donors: 12 healthy controls + 18 UC patients.  Every UC patient
contributed BOTH an inflamed and a non-inflamed biopsy (paired design), so
``Sample`` != ``Subject`` for UC donors (36 UC samples, 18 UC donors).
Cells: 110,110 Healthy / 125,119 Inflamed / 130,263 Non-inflamed (~365,492).

OPEN FLAGS touched by this loader (see OPEN_FLAGS.md)
----------------------------------------------------
- **F1 (UC tissue definition).** ``Health`` is 3-state. This loader is
  *agnostic* to the decision: it loads all 30 donors, sets the harmonized
  2-state ``disease`` (Healthy->normal, Inflamed & Non-inflamed->ulcerative
  colitis, matching Garrido's vocab) AND preserves the raw 3-state in
  ``obs['health']``.  Any inflamed-vs-pooled subsetting happens downstream
  off ``health``; no re-load needed.
- **F2 / F7 (QC-state + crosswalk REVIEW rows).** A few FINE_TO_BROAD
  assignments are tentative and flagged inline below with ``# REVIEW``:
  ``M cells`` (microfold — no clean broad home), ``Immature Enterocytes
  1/2`` (progenitor/colonocyte boundary), ``Secretory TA`` (progenitor vs
  secretory lineage), and ``MT-hi`` (mitochondrial-high QC state, the
  Smillie analogue of Garrido's MT labels — extends F2).

References: ``code/02_atlas_prep/atlas_schemas.md``;
DECISIONS.md correction 2026-05-20 (7/7) and correction 5/7.
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
import scipy.io

from hgnc_remap import ensembl_to_hgnc

logger = logging.getLogger(__name__)

COMPARTMENTS: tuple[str, ...] = ("Epi", "Imm", "Fib")

# Hard invariants fixed by the SCP259 study design.
EXPECTED_N_DONORS = 30
EXPECTED_N_HC = 12
EXPECTED_N_UC = 18
EXPECTED_N_CELLS = 365_492  # tripwire only (genes get QC'd, cells should not)

HEALTH_VALUES = ("Healthy", "Inflamed", "Non-inflamed")

# Harmonized 2-state disease, matching Garrido's CELLxGENE vocabulary so the
# cross-atlas concordance compares like with like. The inflamed/non-inflamed
# distinction is NOT lost — it is preserved in obs['health'] (see F1).
HEALTH_TO_DISEASE = {
    "Healthy":      "normal",
    "Inflamed":     "ulcerative colitis",
    "Non-inflamed": "ulcerative colitis",
}

# Fine -> broad roll-up into the SHARED 15-level vocab defined in
# load_garrido_trigo.FINE_TO_BROAD / atlas_schemas.md. Smillie populates 14 of
# the 15 (no granulocyte: this taxonomy has no neutrophils/eosinophils).
# Lines marked "# REVIEW" are judgment calls to confirm with markers
# (OPEN_FLAGS F7); they are mapped to their best-guess parent so the loader
# runs, but they are not yet biology-locked.
FINE_TO_BROAD = {
    # --- Epithelial: absorptive (colonocyte lineage) ---
    "Enterocytes":            "colonocyte",
    "Best4+ Enterocytes":     "colonocyte",
    "Immature Enterocytes 1": "colonocyte",            # REVIEW: progenitor/colonocyte boundary
    "Immature Enterocytes 2": "colonocyte",            # REVIEW: progenitor/colonocyte boundary
    "M cells":                "colonocyte",            # REVIEW: microfold cell — no clean broad home

    # --- Epithelial: secretory goblet ---
    "Goblet":                 "goblet",
    "Immature Goblet":        "goblet",

    # --- Epithelial: enteroendocrine / tuft ---
    "Enteroendocrine":        "enteroendocrine/tuft",
    "Tuft":                   "enteroendocrine/tuft",

    # --- Epithelial: progenitor / proliferating ---
    "Stem":                   "epithelial progenitor",
    "TA 1":                   "epithelial progenitor",
    "TA 2":                   "epithelial progenitor",
    "Cycling TA":             "epithelial progenitor",
    "Secretory TA":           "epithelial progenitor",  # REVIEW: progenitor vs secretory lineage
    "Enterocyte Progenitors": "epithelial progenitor",

    # --- Stroma: fibroblast ---
    "WNT2B+ Fos-lo 1":        "fibroblast",
    "WNT2B+ Fos-lo 2":        "fibroblast",
    "WNT2B+ Fos-hi":          "fibroblast",
    "WNT5B+ 1":               "fibroblast",
    "WNT5B+ 2":               "fibroblast",
    "RSPO3+":                 "fibroblast",
    "Inflammatory Fibroblasts":"fibroblast",
    "Myofibroblasts":         "fibroblast",

    # --- Stroma: endothelium ---
    "Endothelial":            "endothelium",
    "Microvascular":          "endothelium",
    "Post-capillary Venules": "endothelium",

    # --- Stroma: mural / glia ---
    "Pericytes":              "mural/glia",
    "Glia":                   "mural/glia",

    # --- T cells ---
    "CD4+ Memory":            "T cell",
    "CD4+ Activated Fos-hi":  "T cell",
    "CD4+ Activated Fos-lo":  "T cell",
    "CD4+ PD1+":              "T cell",
    "CD8+ LP":                "T cell",
    "CD8+ IELs":              "T cell",
    "CD8+ IL17+":             "T cell",
    "Tregs":                  "T cell",
    "Cycling T":              "T cell",
    "MT-hi":                  "T cell",               # REVIEW/QC: mito-high state (F2); verify lineage/exclude

    # --- NK / ILC ---
    "NKs":                    "NK/ILC",
    "ILCs":                   "NK/ILC",

    # --- B cells ---
    "Follicular":             "B cell",
    "GC":                     "B cell",
    "Cycling B":              "B cell",

    # --- Plasma cells ---
    "Plasma":                 "plasma cell",

    # --- Myeloid: monocyte / macrophage ---
    "Macrophages":            "monocyte/macrophage",
    "Inflammatory Monocytes": "monocyte/macrophage",
    "Cycling Monocytes":      "monocyte/macrophage",

    # --- Myeloid: dendritic ---
    "DC1":                    "dendritic cell",
    "DC2":                    "dendritic cell",

    # --- Myeloid: mast ---
    "CD69+ Mast":             "mast cell",
    "CD69- Mast":             "mast cell",
}


def _normalize_label(value: object) -> object:
    """Strip + collapse internal whitespace (defensive, matches Garrido).

    Smillie's labels look clean, but trailing/duplicated whitespace in a
    metadata export would silently break the FINE_TO_BROAD lookup. ``+`` and
    ``-`` (e.g. ``WNT2B+``, ``CD69-``) are preserved.
    """
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return value
    return re.sub(r"\s+", " ", str(value)).strip()


def _find_compartment_triplet(scp259_dir: Path, compartment: str) -> tuple[Path, Path, Path]:
    """Locate (matrix, genes, barcodes) for one compartment by globbing.

    Globs ``**/gene_sorted-<compartment>.matrix.mtx`` so the unstable hash
    directory names do not need to be hard-coded. Raises loudly if the
    matrix is missing or ambiguous, or if either sidecar file is absent.
    """
    matches = sorted(scp259_dir.glob(f"**/gene_sorted-{compartment}.matrix.mtx"))
    if not matches:
        raise FileNotFoundError(
            f"Smillie loader: no gene_sorted-{compartment}.matrix.mtx under "
            f"{scp259_dir}. Did the SCP259 download land elsewhere?"
        )
    if len(matches) > 1:
        raise ValueError(
            f"Smillie loader: {len(matches)} candidate matrices for "
            f"compartment {compartment!r}: {matches}. Expected exactly one."
        )
    mtx = matches[0]
    genes = mtx.parent / f"{compartment}.genes.tsv"
    barcodes = mtx.parent / f"{compartment}.barcodes2.tsv"
    for sidecar in (genes, barcodes):
        if not sidecar.exists():
            raise FileNotFoundError(
                f"Smillie loader: expected sidecar {sidecar} next to {mtx} "
                f"but it is missing."
            )
    return mtx, genes, barcodes


def _load_compartment(scp259_dir: Path, compartment: str) -> AnnData:
    """Read one gene-sorted 10X triplet and return it as cells x genes."""
    mtx_path, genes_path, barcodes_path = _find_compartment_triplet(
        scp259_dir, compartment
    )
    logger.info("Reading %s compartment matrix: %s", compartment, mtx_path)

    # gene_sorted matrix is genes x cells; mmread is slow + memory-heavy on
    # the ~10^8-nnz compartments — run this on a compute node, not the login
    # node (see module/run notes).
    mat = scipy.io.mmread(str(mtx_path)).tocsr()  # genes x cells
    genes = pd.read_csv(genes_path, header=None)[0].astype(str).tolist()
    barcodes = pd.read_csv(barcodes_path, header=None)[0].astype(str).tolist()

    if mat.shape != (len(genes), len(barcodes)):
        raise ValueError(
            f"Smillie loader: {compartment} matrix shape {mat.shape} does not "
            f"match genes={len(genes)} x barcodes={len(barcodes)}. The triplet "
            f"files are inconsistent."
        )

    adata_c = ad.AnnData(
        X=mat.T.tocsr(),  # -> cells x genes
        obs=pd.DataFrame(index=pd.Index(barcodes, name=None)),
        var=pd.DataFrame(index=pd.Index(genes, name=None)),
    )
    logger.info(
        "  %s: %d cells x %d genes (%d nnz)",
        compartment, adata_c.n_obs, adata_c.n_vars, mat.nnz,
    )
    return adata_c


def load(
    scp259_dir: str,
    apply_v1_filter: bool = True,
    raw_count_mode: bool = False,
) -> AnnData:
    """Load Smillie 2019 (SCP259) with full broad + fine tier.

    Parameters
    ----------
    scp259_dir
        Path to the SCP259 download root (the directory that contains the
        ``expression/`` and ``metadata/`` sub-trees), e.g.
        ``~/uc-cross-atlas-data/atlases/SCP259``.
    apply_v1_filter
        If True (default), validate that all cells are ``normal``/``ulcerative
        colitis`` and keep the whole cohort (all 30 donors are v1-relevant).
        This is the hook where an F1 inflamed/non-inflamed subset would later
        be applied; for now it keeps everything.
    raw_count_mode
        Must remain False for v1 (DECISIONS correction 5/7). The Smillie
        matrices ARE raw counts, so the loader normalizes them to log1p(CP10k);
        ``raw_count_mode=True`` (leave X as raw counts) is not supported.

    Returns
    -------
    AnnData
        cells x genes, ``X`` = log1p(CP10k) float, raw counts in
        ``layers['counts']``, with the standard obs schema:
        ``cell_type_fine``, ``cell_type_broad``, ``donor``, ``disease``,
        ``health``, ``sample``, ``compartment``, ``tissue``, ``batch``,
        ``n_genes``, ``n_counts``.
    """
    if raw_count_mode:
        raise ValueError(
            "raw_count_mode=True is not supported for v1 (DECISIONS correction "
            "5/7): all atlases use uniform log1p(CP10k) input. The Smillie "
            "matrices are raw counts and are normalized on load."
        )

    root = Path(scp259_dir).expanduser()
    if not root.exists():
        raise FileNotFoundError(f"Smillie loader: scp259_dir does not exist: {root}")

    # ---- 1. Load + transpose + concatenate the three compartments ----
    per_compartment = [_load_compartment(root, c) for c in COMPARTMENTS]
    gene_sets = [set(a.var_names) for a in per_compartment]
    shared = set.intersection(*gene_sets)
    union = set.union(*gene_sets)
    logger.info(
        "Compartment gene sets: sizes=%s, shared=%d, union=%d",
        [a.n_vars for a in per_compartment], len(shared), len(union),
    )
    # Outer join unions the gene space (zero-fill where a gene is absent from a
    # compartment's matrix); if the three lists are identical this is a no-op.
    adata = ad.concat(
        per_compartment, axis=0, join="outer", fill_value=0, merge="first",
    )
    adata.obs_names_make_unique()  # barcodes are Subject.Sample.barcode -> already unique
    if not adata.obs_names.is_unique:
        raise ValueError("Smillie loader: cell barcodes are not unique after concat.")
    logger.info("Concatenated atlas: %d cells x %d genes", adata.n_obs, adata.n_vars)

    # ---- 2. Join the metadata (direct: barcodes2.tsv == NAME) ----
    meta_matches = sorted(root.glob("**/all.meta2.txt"))
    if not meta_matches:
        raise FileNotFoundError(f"Smillie loader: metadata/all.meta2.txt not found under {root}")
    meta_path = meta_matches[0]
    logger.info("Reading metadata: %s", meta_path)
    # skiprows=[1] drops the SCP 'TYPE' boilerplate row; header (row 0) kept.
    meta = pd.read_csv(meta_path, sep="\t", skiprows=[1])
    expected_cols = {"NAME", "Cluster", "Subject", "Health", "Location", "Sample"}
    missing_cols = expected_cols - set(meta.columns)
    if missing_cols:
        raise KeyError(
            f"Smillie loader: metadata is missing columns {sorted(missing_cols)}. "
            f"Got: {list(meta.columns)}"
        )
    meta["NAME"] = meta["NAME"].astype(str).str.strip()
    meta = meta.set_index("NAME")
    if meta.index.duplicated().any():
        raise ValueError(
            f"Smillie loader: {int(meta.index.duplicated().sum())} duplicate NAME "
            f"values in metadata; cannot use as a join key."
        )

    aligned = meta.reindex(adata.obs_names)
    # Completeness gate: every matrix cell must have a metadata row.
    orphan = aligned["Cluster"].isna()
    if orphan.any():
        ex = adata.obs_names[orphan][:5].tolist()
        raise ValueError(
            f"Smillie loader: {int(orphan.sum())} cells have no metadata row "
            f"after the NAME join (examples: {ex}). Matrix barcodes and "
            f"all.meta2.txt NAME values disagree."
        )

    # ---- 3. Build the standard obs schema ----
    fine = aligned["Cluster"].map(_normalize_label)
    unmapped = sorted(set(fine.unique()) - set(FINE_TO_BROAD))
    if unmapped:
        raise KeyError(
            f"Smillie loader: {len(unmapped)} fine labels have no FINE_TO_BROAD "
            f"entry: {unmapped}. Extend the map (load_smillie.FINE_TO_BROAD)."
        )
    broad = fine.map(FINE_TO_BROAD)

    health = aligned["Health"].astype(str)
    bad_health = sorted(set(health.unique()) - set(HEALTH_VALUES))
    if bad_health:
        raise ValueError(
            f"Smillie loader: unexpected Health values {bad_health}; expected "
            f"{list(HEALTH_VALUES)}."
        )
    disease = health.map(HEALTH_TO_DISEASE)

    obs = pd.DataFrame(index=adata.obs_names)
    obs["cell_type_fine"] = fine.astype("category")
    obs["cell_type_broad"] = broad.astype("category")
    obs["donor"] = aligned["Subject"].astype(str).astype("category")
    obs["disease"] = disease.astype("category")
    obs["health"] = health.astype("category")          # 3-state preserved for F1
    obs["sample"] = aligned["Sample"].astype(str).astype("category")
    obs["compartment"] = aligned["Location"].astype(str).astype("category")
    obs["batch"] = aligned["Sample"].astype(str).astype("category")
    obs["tissue"] = "colonic mucosa"
    if "nGene" in aligned.columns:
        obs["n_genes"] = pd.to_numeric(aligned["nGene"], errors="coerce")
    if "nUMI" in aligned.columns:
        obs["n_counts"] = pd.to_numeric(aligned["nUMI"], errors="coerce")
    adata.obs = obs

    n_broad = int(obs["cell_type_broad"].nunique())
    n_fine = int(obs["cell_type_fine"].nunique())
    logger.info("Tier cardinalities: fine=%d, broad=%d", n_fine, n_broad)
    if not (10 <= n_broad <= 15):
        logger.warning(
            "Broad-tier cardinality %d outside the v1 10-15 target; review "
            "FINE_TO_BROAD.", n_broad,
        )

    # ---- 4. v1 filter + donor-structure invariant ----
    if apply_v1_filter:
        bad_disease = sorted(set(obs["disease"].unique()) - {"normal", "ulcerative colitis"})
        if bad_disease:
            raise ValueError(
                f"Smillie loader: unexpected disease values {bad_disease} after "
                f"Health mapping; expected normal / ulcerative colitis only."
            )
        # (F1 inflamed/non-inflamed subsetting would hook here. v1 keeps all.)

        donor_disease = (
            obs[["donor", "disease"]].drop_duplicates()
            .groupby("disease", observed=True).size()
        )
        n_donors = int(obs["donor"].nunique())
        n_hc = int(donor_disease.get("normal", 0))
        n_uc = int(donor_disease.get("ulcerative colitis", 0))
        if (n_donors != EXPECTED_N_DONORS or n_hc != EXPECTED_N_HC or n_uc != EXPECTED_N_UC):
            raise ValueError(
                f"Smillie loader: donor-structure invariant violated. Got "
                f"n_donors={n_donors} ({n_hc} HC + {n_uc} UC); expected "
                f"{EXPECTED_N_DONORS} ({EXPECTED_N_HC} HC + {EXPECTED_N_UC} UC)."
            )

    if adata.n_obs != EXPECTED_N_CELLS:
        logger.warning(
            "Cell count %d != expected %d (tripwire only; the donor-structure "
            "and completeness gates are the hard checks).",
            adata.n_obs, EXPECTED_N_CELLS,
        )

    # ---- 5. Normalize raw counts -> log1p(CP10k); keep counts in a layer ----
    if not np.issubdtype(adata.X.dtype, np.floating):
        adata.X = adata.X.astype(np.float32)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    logger.info("Normalized X to log1p(CP10k); raw counts in layers['counts'].")

    # ---- 6. Gene symbols -> deduplicated, NCBI-authoritative HGNC set ----
    # No 'feature_name' column -> ensembl_to_hgnc takes its symbol-fallback
    # path: dedup duplicate symbols by max-expression, drop non-approved.
    adata = ensembl_to_hgnc(adata)

    logger.info(
        "Smillie load complete: %d cells x %d genes, %d donors.",
        adata.n_obs, adata.n_vars, int(adata.obs["donor"].nunique()),
    )
    return adata
