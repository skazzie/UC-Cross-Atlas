# Atlas Schemas (v1)

Canonical reference for atlas obs schemas. Loaders read this implicitly
through their parameters. When adding an atlas, add a section here first
and pin the schema, then write the loader.

Schema captured 2026-05-20.

---

## garrido_trigo

> **Matrix source is changing (DECISIONS correction 9, 2026-06-04).** The
> CELLxGENE deposit's `obs.index` is synthetic (`cell1, cell2, ...`) — the
> 10X barcodes are stripped, so the GEO annotation cannot be joined onto
> it deterministically. v1 moves to `GSE214695_RAW.tar` per-GSM matrices,
> where barcodes are intact and the CSV's `Unnamed: 0` (`SC_xxx_<barcode>`)
> is the unique join key. This section still documents the
> pre-correction-9 schema; it will be rewritten when the RAW.tar loader
> lands (next session, alongside Smillie SCP259). Per-donor cell counts
> (n_cells, n_donors below) carry over unchanged because they're fixed by
> the study design.

- **Matrix source**: CELLxGENE (superseded — moving to GEO `GSE214695_RAW.tar`)
- **Matrix download URL**: https://datasets.cellxgene.cziscience.com/b1a62801-f509-45f8-b55f-533fbb7e7800.h5ad
- **Matrix file size**: 176 MB
- **Annotation source**: GEO GSE214695 supplementary file
  `GSE214695_cell_annotation.csv` (Salas-lab 91-cluster fine labels;
  NOT present in the CELLxGENE deposit — see DECISIONS.md correction
  reversing (4/7), and correction (9) superseding the matrix source)
- **n_cells**: 46,700 (full deposit) / 30,068 (UC subset HC + UC)
- **n_donors**: 18 (full) / 12 (UC subset: HC_1-6 + UC_1-6)
- **Assay**: 10x 3' v3

### Filter chain (UC subset)

```python
adata = adata[adata.obs["disease"].isin(["normal", "ulcerative colitis"])]
```

### Tiers

Both tiers are populated by joining the GEO annotation CSV onto the
CELLxGENE matrix by barcode (auto-detected join key; loader raises a
diagnostic error if any cell is orphaned).

- **fine_tier**: 82 surviving labels (91 published - 9 Ribhi collapsed
  into parent lineage). Stored in `obs['cell_type_fine']`. Hyper-specific
  states (PC IgA heat shock 1/2, M1 ACOD1, etc.) won't map cleanly to
  Smillie/Mennillo; the cross-atlas fine intersection is expected to be
  modest — do not oversell "91-way fine concordance."
- **broad_tier**: 15-level roll-up of the surviving fine labels via
  `FINE_TO_BROAD` in `load_garrido_trigo.py`. Stored in
  `obs['cell_type_broad']`. Levels: colonocyte, goblet,
  epithelial progenitor, enteroendocrine/tuft, fibroblast, endothelium,
  mural/glia, T cell, NK/ILC, B cell, plasma cell, monocyte/macrophage,
  dendritic cell, mast cell, granulocyte.
- **CELLxGENE `obs['cell_type']` (5 CL lineages)** is NOT used for tier
  output; the loader uses the GEO annotation as the source of truth.
  Useful only as a coherence cross-check.

### Ribhi collapse policy

The 9 `*Ribhi` clusters are a ribosomal-high cross-lineage transcriptional
state (verified empirically: RPL*/RPS* dominate top-20 markers in
`garrido_trigo_markers.xlsx`). Collapsed in the loader before any tier
logic, via `RIBHI_TO_PARENT`:

| Ribhi label         | Parent fine label |
|---------------------|-------------------|
| B cell Ribhi        | B cell            |
| Epithelium Ribhi    | epithelial        |
| M0_Ribhi            | M0                |
| Ribhi T cells 1/2   | T                 |
| S1 Ribhi            | S1                |
| Fibroblasts Ribhi   | fibroblast        |
| Mast Ribhi          | mast              |
| DCs CCL22_Ribhi     | DCs CCL22         |

Ribhi cells are never entered as standalone fine clusters in cross-atlas
concordance.

### Label normalization on load

Every label from the GEO CSV is passed through `_normalize_label`:

- Collapse internal whitespace (the CSV contains
  `'PC  immediate early response'` with a literal double space).
- Normalize `ï`/`é` → `i`/`e` so map keys can stay ASCII.
- Strip leading/trailing whitespace.

This prevents silent join breakage when label strings vary between
the marker xlsx and the annotation CSV.

### Covariates for scDRS

- `donor_id`, `sex`, `biospsy_or_surgical_resection_area` (NOTE: typo
  "biospsy" is the actual column name in the deposit; preserve it)

### Counts

- `X`: log-normalized float32 (CELLxGENE)
- `raw`: NOT PRESENT (use `--flag-raw-count False`)
- `var_names`: Ensembl IDs (HGNC symbols in `var['feature_name']`)
- **Do NOT switch `--flag-raw-count True` for Garrido-Trigo** even though
  the GEO RAW.tar matrices have true integer counts; uniform
  `--flag-raw-count False` across atlases preserves cross-atlas input
  comparability (DECISIONS.md correction 5/7).

### Markers reference

Marker gene table at `data/atlases/garrido_trigo_markers.xlsx` (5 sheets,
91 fine clusters total). Used by `sanity_check.py` to validate MAGMA
top-gene patterns against known inflammatory subsets, AND by us in the
session that wrote correction (8) to verify Ribhi = ribosomal-high state.

---

## pangi

- **Source**: CELLxGENE
- **Download URL**: https://datasets.cellxgene.cziscience.com/1dcf15ee-c103-4aaa-8b8c-0fc697fcccc8.h5ad
- **File size**: estimated 15-25 GB (verify with `curl -I` before downloading)
- **n_cells**: 1,596,200 (full atlas) / ~150-200k expected after v1 filter
- **Studies integrated**: 25 source studies including Elmentaite2021
  (398,460 cells), Kong2023 (235,327 cells), Kim2022, Madissoon2019, etc.
  See `obs['study']`.

### Filter chain (v1)

```python
adata = adata[
    adata.obs["disease"].isin(["normal", "ulcerative colitis",
                               "inflammatory bowel disease"])
    & adata.obs["organ_unified"].isin([
        "ascending colon", "caecum", "colon", "descending colon",
        "rectum", "sigmoid colon", "transverse colon"
    ])
    & (adata.obs["sample_type"] != "Organ_donor_resection")
]
```

### Tiers

- **broad_tier_column**: `level_2_annot` (~30 categories: Absorptive,
  B_plasma, Conventional_CD4/CD8, Fibroblast, Macrophage, etc.)
- **fine_tier_column**: `level_3_annot` (~70+ subclusters)
- **cl_ontology_column**: `cell_type` (~75 CL terms; harmonization anchor only)

### Covariates for scDRS

- `donorID_unified` (NOT `donor_id`; the latter is study-prefixed and
  inconsistent across studies)
- `sex`, `assay`, `sample_type`

### Sensitivity columns

- `study`: filter "Elmentaite2021" out for HCA Gut overlap sensitivity
- `donorID_unified`: scan for any Smillie 2019 donor IDs (likely none;
  if any found, produce removal sensitivity)

### Counts

- `X`: log-normalized float32 (expect; verify on download)
- `raw`: expected NOT PRESENT (verify on download)
- `var_names`: Ensembl IDs

---

## hca_gut

- **Source**: CELLxGENE
- **Download URL**: https://datasets.cellxgene.cziscience.com/f34d2b82-9265-4a73-bda4-852933bf2a8d.h5ad
- **File size**: estimated 4-8 GB
- **n_cells**: 428,469 (full) / 30k-70k expected after v1 filter

### Filter chain (v1)

```python
adata = adata[
    adata.obs["Age_group"].isin(["Adult", "Adult_MLN"])
    & adata.obs["tissue"].isin([
        "ascending colon", "caecum", "colon", "descending colon",
        "large intestine", "rectum", "sigmoid colon", "transverse colon"
    ])
]
```

### Tiers

- **broad_tier_column**: `category` (9 lineages: Epithelial, Mesenchymal,
  Myeloid, T cells, Plasma cells, B cells, Endothelial, Neuronal, RBC)
- **fine_tier_column**: `author_cell_type` (~120 cell types; REQUIRES
  roll-up via `cl_rollup.py` to ~30-50 for v1)
- **cl_ontology_column**: `cell_type` (~75 CL terms; harmonization anchor)

### Covariates for scDRS

- `donor_id`, `sex`, `assay`, `batch`, `Fraction`
- `Fraction` is mandatory: 5 sorting strategies (SC, SC-45N, SC-45P,
  SC-EPCAMN, SC-EPCAMP) strongly determine cell-type composition.

### Sensitivity

- Excluding `disease == "Crohn disease"` cells (analogous to Pan-GI's
  Elmentaite2021 sensitivity)

### Counts

- `X`: log-normalized float32 (expect; verify on download)
- `raw`: expected NOT PRESENT (verify on download)
- `var_names`: Ensembl IDs

### Important caveat

Cell-type proportions in this atlas are artificially controlled by the
`Fraction` column (sort strategy) and are NOT biologically interpretable.
scDRS per-cell scoring is unaffected, but DO NOT report cell-type
abundance from HCA Gut in figures.

---

## smillie

Schema captured 2026-06-04 from the SCP259 download at
`~/uc-cross-atlas-data/atlases/SCP259/`. See DECISIONS.md correction
(10) for the schema-capture entry and the (7/7) source-lock backstory.

- **Matrix source**: Single Cell Portal SCP259 (NOT CELLxGENE; that
  deposit is a 34,772-cell healthy-epithelial-only subset and is
  unusable for v1)
- **Download URL**: https://singlecell.broadinstitute.org/single_cell/study/SCP259
  (requires SCP account + email verification + browser consent click on
  the study page; bulk download via SCP CLI or web)
- **On-disk location**: `~/uc-cross-atlas-data/atlases/SCP259/`
  (Hummingbird home — outside the repo tree; ~5 GB)
- **Annotation source**: in-band — `metadata/all.meta2.txt` ships with
  the SCP bundle (no external CSV needed, unlike Garrido)
- **n_cells**: 365,492 (110,110 Healthy + 125,119 Inflamed + 130,263 Non-inflamed)
- **n_donors**: 30 (12 HC + 18 UC, paired — every UC patient
  contributed BOTH an inflamed and a non-inflamed biopsy, so 36 UC
  samples / 18 UC donors)
- **Assay**: 10x 3' v2 / v3 (per Smillie 2019)

### On-disk layout

Three gene-sorted 10X triplets under hashed sub-directories (hash names
NOT stable; the loader globs by filename), plus one metadata file:

```
expression/<hash>/  gene_sorted-Epi.matrix.mtx  Epi.genes.tsv  Epi.barcodes2.tsv
expression/<hash>/  gene_sorted-Imm.matrix.mtx  Imm.genes.tsv  Imm.barcodes2.tsv
expression/<hash>/  gene_sorted-Fib.matrix.mtx  Fib.genes.tsv  Fib.barcodes2.tsv
metadata/all.meta2.txt
cluster/{Epi,Imm,Fib}.tsne.txt   (coordinates; not used by the loader)
```

- Matrices are **gene-sorted = genes × cells** — the loader transposes
  each to cells × genes before assembling the AnnData. Epi header is
  20,028 genes × 123,006 cells × 174.4M nnz; the others are similar
  scale. mmread + transpose + outer-join concat needs ~64–96 GB RAM;
  run on a compute node, not the login node.
- `var_names` are **HGNC symbols**, single-column (e.g. `7SK`, `A1BG`),
  NOT Ensembl. `ensembl_to_hgnc` takes its symbol-fallback path on load.
- `X` is **raw integer counts** (MatrixMarket `coordinate integer`).
  Unlike Garrido/HCA/Pan-GI (which ship CELLxGENE log-normalized X),
  this loader applies `log1p(CP10k)` itself to keep input state uniform
  across atlases (DECISIONS correction 5/7). Raw counts preserved in
  `layers['counts']`.

### Metadata schema (`metadata/all.meta2.txt`)

Tab-separated. **Line 2 is an SCP-specific `TYPE` boilerplate row that
must be skipped on read (`pd.read_csv(..., sep="\t", skiprows=[1])`);
header on line 1 is kept.** Columns:

| column   | TYPE    | role                                                   |
|----------|---------|--------------------------------------------------------|
| NAME     | group   | cell id `Subject.Sample.barcode` — join key, equals `*.barcodes2.tsv` exactly (direct join, no reconstruction) |
| Cluster  | group   | fine cell type (51 labels) → `cell_type_fine`          |
| nGene    | numeric | genes detected → `n_genes`                             |
| nUMI     | numeric | UMIs → `n_counts`                                      |
| Subject  | group   | donor (e.g. `N7`) → `donor`                            |
| Health   | group   | `Healthy` / `Inflamed` / `Non-inflamed` → `health`     |
| Location | group   | `Epi` / `Imm` / `Fib` → `compartment`                  |
| Sample   | group   | biopsy id (e.g. `N7.EpiA`) → `sample` (also `batch`)   |

### Filter chain (v1)

```python
# No subset filter — all 30 donors are v1-relevant.
# F1 (UC tissue definition) hooks here downstream off obs['health'];
# the loader keeps the 3-state health column for that decision.
```

### Tiers

Both tiers come from the `Cluster` column via the joined metadata.

- **fine_tier**: 51 published Smillie labels. Stored in
  `obs['cell_type_fine']`.
- **broad_tier**: 14-level roll-up via `FINE_TO_BROAD` in
  `load_smillie.py`. Levels: colonocyte, goblet,
  enteroendocrine/tuft, epithelial progenitor, fibroblast, endothelium,
  mural/glia, T cell, NK/ILC, B cell, plasma cell, monocyte/macrophage,
  dendritic cell, mast cell. **No granulocyte** — this taxonomy has no
  neutrophils/eosinophils, so Smillie populates 14 of the shared 15-level
  vocab from Garrido. (Cross-atlas concordance handles the missing level
  per-atlas; no broad-vocab change needed.)
- **Disease** (2-state, harmonized to Garrido vocab):
  `Healthy → normal`, `Inflamed | Non-inflamed → ulcerative colitis`.
  The 3-state `Health` is preserved in `obs['health']` so any future F1
  inflamed-vs-pooled subsetting is downstream-only, no re-load.

### REVIEW rows (OPEN_FLAGS F7 — marker-QC at M2 to lock)

Four `FINE_TO_BROAD` rows are best-guess judgment calls, not yet
biology-locked. The loader runs with these provisional assignments
flagged inline; M2 marker-QC moves them into DECISIONS once resolved.

- `M cells` → colonocyte (microfold; no clean broad home)
- `Immature Enterocytes 1/2` → colonocyte (progenitor/colonocyte boundary)
- `Secretory TA` → epithelial progenitor (progenitor vs secretory lineage)
- `MT-hi` → T cell (mitochondrial-high QC state — Smillie analogue of
  Garrido's MT labels; F2 governs whether MT-state labels collapse vs
  exclude uniformly across atlases)

---

## mennillo (DEFERRED to next session)

- **Source**: GEO GSE229072 (verify accession)
- **Paper**: Mennillo 2024, Nat Commun 15:1493
- **Expected n_cells**: ~50-100k (verify on download)
- **Schema**: TBD (capture after download)

Placeholder loader (`load_mennillo.py`) created with TODO markers; will be
filled in after GEO download.
