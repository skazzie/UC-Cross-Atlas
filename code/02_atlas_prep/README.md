# 02_atlas_prep — atlas standardization + donor-attribution metadata

Download, harmonize, and HGNC-remap the 5 atlases (UC trio + 2 broad
comparators) into a uniform `.h5ad` format. Generate per-atlas covariate
files for scDRS and donor-attribution metadata for Pan-GI donor-overlap
analysis.

## Atlases

- **Smillie 2019** (Single Cell Portal SCP259) — 366,650 colon mucosa
  cells, 18 UC + 12 healthy.
- **Garrido-Trigo 2023 UC subset** (GEO GSE214695 / CELLxGENE) — 6 healthy
  + 6 active UC. (Previously listed as "Kong 2023" — see DECISIONS.md
  correction 2026-05-20 (2/7).)
- **TAURUS-IBD** (Thomas et al. 2024, Zenodo
  `10.5281/zenodo.14007626` (v3, pinned 2026-06-06; md5 of pooled file
  `c1bd13b92cacb164a401c6c4a4e7912c`)) — longitudinal anti-TNF single-cell atlas (UC + CD).
  Subset to **UC donors only** and a **single time-point per donor**
  (pre-treatment baseline preferred). Replaces the previously-planned
  Mennillo 2024 (anti-integrin) per DECISIONS 16. Verify ≥8 UC donors
  after subsetting.
- **HCA Gut Cell Atlas** (cellxgene, Elmentaite 2021) — large-intestine
  subset only. M1 must verify zero donor overlap with the UC trio.
- **Pan-GI** (cellxgene, Oliver 2024) — Extended+ slice (1.6M cells),
  the only slice with all lineages. Filtered to large-intestine cells
  by `load_pangi.py`. Contains Elmentaite2021 (= HCA Gut, 398,460 cells)
  and the Kong 2023 CD atlas (235,327 cells); no Smillie2019. See
  DECISIONS.md correction 2026-05-20 (3/7) for the overlap policy.

## HGNC pin

A dated NCBI `gene_info` snapshot is committed at
`data/reference/gene_info.tsv.gz` (current pin: **2026-05-21**). The
approved symbol set is built from its `Symbol` column **only** —
synonyms are deliberately excluded so deprecated aliases cannot slip
through the membership filter unchanged (if alias resolution is ever
needed, it goes in as an explicit alias→approved remap, not a
membership test).

`hgnc_remap.ensembl_to_hgnc` runs as the last step of every loader and
enforces a **hard canonical-hit survival gate**: ≥95% of the canonical
UC GWAS hits (IL23R, JAK2, TYK2, NKX2-3, ATG16L1) must be present in
`adata.var_names` after dedup + symbol-validity filter. A miss raises
loudly — a missing IBD/UC locus post-remap silently biases every
downstream score. To refresh the pin: download a fresh
`Homo_sapiens.gene_info.gz`, update `GENE_INFO_PIN_DATE` in
`hgnc_remap.py`, replace the committed snapshot in one commit, and
document the date bump in DECISIONS.md.

## Atlas preprocessing policy (Option B)

Take published preprocessing as-is. Do **not** re-run uniform doublet
removal / ambient correction / QC across atlases. Disclose preprocessing
heterogeneity in Methods Limitations. The marker-gene QC step at end of
M2 catches the worst label-mismatch problems.

## Per-atlas covariate file generation

`<atlas>_covariates.tsv` per `code/03_scdrs/README.md` covariate set.

## Donor-attribution metadata (v2-setup decision)

For every atlas, save
`data/atlases/donor_metadata/<atlas>_donor_metadata.csv` with columns:

- `donor_id`
- `originating_study` — atlas name (or upstream study for integrated
  atlases like Pan-GI)
- `tissue` — colon, large intestine, ileum, etc.
- `disease_status` — UC active, UC remission, healthy, etc.
- `assay_protocol` — 10X v2/v3, dropseq, etc.
- `tissue_state` (optional) — inflamed / uninflamed / unspecified.
  Used for Methods Table 1.

Pan-GI is the immediate v1 consumer (donor-overlap exclusion at M5).
Same schema scales to AMP-RA Phase 2 and any future integrated atlas in v2.

## Pan-GI / HCA Gut log-normalized matrix verification (M1)

Verify each atlas exposes a usable **log-normalized** expression matrix
(in `.X` or a named layer). scDRS does NOT require raw counts: its
`preprocess` operates on the log-scale matrix (gene stats are computed in
both log and non-log scale). Raw counts are only relevant when
`flag_raw_count=True`, which applies CP10k + log1p internally; v1 locks
`--flag-raw-count False` uniformly across all five atlases (DECISIONS
correction 5/7), so the log-normalized matrix is the input contract.
Loaders that ship raw counts (Smillie SCP259, TAURUS Zenodo, Garrido-Trigo
RAW.tar once the correction-9 rewrite lands) apply `log1p(CP10k)` on load
and preserve raw counts in `layers['counts']`.

## Output

- `data/atlases/<atlas>.h5ad` — standardized AnnData with HGNC-pinned symbols.
- `data/atlases/<atlas>_covariates.tsv` — scDRS covariate file.
- `data/atlases/donor_metadata/<atlas>_donor_metadata.csv` — donor attribution.
- `data/atlases/<atlas>_sce.rds` — SingleCellExperiment for seismicGWAS.

## v1 loader scripts (session 2026-05-20)

Canonical schema reference: `atlas_schemas.md` (commit with this batch).
DECISIONS.md corrections 2026-05-20 (1/7) through (7/7) document the v1
locked policy.

| Script | Atlas | Status | Filter chain |
|--------|-------|--------|--------------|
| `load_garrido_trigo.py` | Garrido-Trigo 2023 (UC core) | Production | 12 HC/UC GSMs from `GSE214695_RAW.tar` (CD skipped at glob); composite inner-join `{sample}_{barcode}` against `GSE214695_cell_annotation.csv`; `log1p(CP10k)` applied on load (raw in `layers['counts']`); 86 fine (post-Ribhi) → 15 broad; verified 30,068 cells / 12 donors end-to-end. See DECISIONS corrections (9), (12). |
| `load_pangi.py` | Pan-GI Extended+ (broad comparator) | Production | UC + IBD diseases x colon organs x non-organ-donor |
| `load_hca_gut.py` | HCA Gut / Elmentaite 2021 (broad reference) | Production | `Age_group in {Adult, Adult_MLN}` x colon tissues |
| `load_smillie.py` | Smillie 2019 (UC core) | Loader complete; first compute-node run pending | 30 donors (12 HC + 18 UC paired); `Health` 3-state preserved (`obs['health']`); harmonized 2-state `disease` via `HEALTH_TO_DISEASE`; `log1p(CP10k)` applied on load (raw in `layers['counts']`); 51 fine → 14 broad. See DECISIONS correction 10. |
| `load_taurus.py` | TAURUS-IBD (Thomas 2024) (UC core) | Skeleton | Zenodo download deferred; replaces Mennillo per DECISIONS (16); subset to UC + pre-treatment baseline on load |

Each production loader exposes a `load(h5ad_path, apply_v1_filter=True,
raw_count_mode=False)` function returning a standardized AnnData with
obs columns: `cell_type_broad`, `cell_type_fine`, `donor`, `sex`,
`tissue`, `disease`, `assay`, `batch` (plus atlas-specific extras).

Paired sensitivity loaders:

- `load_pangi.load_pangi_no_elmentaite()` — drops Elmentaite2021 cells
  (HCA Gut overlap test).
- `load_pangi.load_pangi_no_smillie()` — drops any cells matching the
  Smillie donor-ID pattern (expected no-op based on inspection).
- `load_hca_gut.load_hca_gut_no_crohn()` — drops Crohn-disease cells.

## Shared utilities

| Script | Purpose |
|--------|---------|
| `hgnc_remap.py` | Convert Ensembl `var_names` to HGNC symbols via `var['feature_name']`; drop duplicates / non-canonical entries. Called as the last step of every loader. |
| `cl_rollup.py` + `cl_rollup_maps.yaml` | Collapse high-cardinality fine annotations (e.g. HCA Gut's ~120 `author_cell_type`) down to ~30-50 categories. Adds `cell_type_fine_rolled` to obs; preserves the original column. Override mappings per atlas in the YAML. |
| `../03_scdrs/aggregate_null_draws.py` | Aggregate scDRS per-cell null-draw z-scores to per-cell-type tensors (input to Brown's-method empirical correlation matrix in M5). Invoked by the SLURM wrapper. |

## Data references

- `data/atlases/garrido_trigo_markers.xlsx` — Salas-lab marker gene table
  (91 fine clusters across 5 compartments), used by sanity_check.py to
  validate MAGMA top-gene patterns and used to verify Ribhi = ribosomal-high
  state (RPL*/RPS* dominate top-20 markers in every Ribhi cluster).
- `GSE214695_cell_annotation.csv` (GEO supplementary, downloaded to
  Hummingbird scratch alongside the matrix) — the per-cell 91-cluster
  fine annotation. Joined onto the CELLxGENE matrix by barcode in
  `load_garrido_trigo.py`. See DECISIONS.md correction reversing (4/7).
