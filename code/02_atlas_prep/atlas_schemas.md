# Atlas Schemas (v1)

Canonical reference for atlas obs schemas. Loaders read this implicitly
through their parameters. When adding an atlas, add a section here first
and pin the schema, then write the loader.

Schema captured 2026-05-20.

---

## garrido_trigo

- **Source**: CELLxGENE
- **Download URL**: https://datasets.cellxgene.cziscience.com/b1a62801-f509-45f8-b55f-533fbb7e7800.h5ad
- **File size**: 176 MB
- **n_cells**: 46,700 (full deposit) / 30,068 (UC subset HC + UC)
- **n_donors**: 18 (full) / 12 (UC subset: HC_1-6 + UC_1-6)
- **Assay**: 10x 3' v3

### Filter chain (UC subset)

```python
adata = adata[adata.obs["disease"].isin(["normal", "ulcerative colitis"])]
```

### Tiers

- **broad_tier_column**: `cell_type` (5 CL lineages: T cell of anorectum,
  plasma cell, colon epithelial cell, stromal cell of lamina propria of
  colon, myeloid cell)
- **fine_tier_column**: NONE in CELLxGENE deposit. Documented limitation.
  Garrido-Trigo contributes broad-tier only to cross-atlas concordance.

### Covariates for scDRS

- `donor_id`, `sex`, `biospsy_or_surgical_resection_area` (NOTE: typo
  "biospsy" is the actual column name in the deposit; preserve it)

### Counts

- `X`: log-normalized float32
- `raw`: NOT PRESENT (use `--flag-raw-count False`)
- `var_names`: Ensembl IDs (HGNC symbols in `var['feature_name']`)

### Markers reference

Marker gene table at `data/atlases/garrido_trigo_markers.xlsx` (5 sheets,
91 fine clusters total). Used by `sanity_check.py` to validate MAGMA
top-gene patterns against known inflammatory subsets.

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

## smillie (DEFERRED to next session)

- **Source**: Single Cell Portal SCP259 (NOT CELLxGENE; that deposit is
  healthy-epithelial-only and unusable)
- **Download URL**: TBD (requires SCP account + browser consent click)
- **Expected n_cells**: 366,650 across 30 donors (18 UC + 12 HC)
- **Schema**: TBD (capture after download)

Placeholder loader (`load_smillie.py`) created with TODO markers; will be
filled in after SCP259 is on disk.

---

## mennillo (DEFERRED to next session)

- **Source**: GEO GSE229072 (verify accession)
- **Paper**: Mennillo 2024, Nat Commun 15:1493
- **Expected n_cells**: ~50-100k (verify on download)
- **Schema**: TBD (capture after download)

Placeholder loader (`load_mennillo.py`) created with TODO markers; will be
filled in after GEO download.
