# 02_atlas_prep — atlas standardization + donor-attribution metadata

Download, harmonize, and HGNC-remap the 5 atlases (UC trio + 2 broad
comparators) into a uniform `.h5ad` format. Generate per-atlas covariate
files for scDRS and donor-attribution metadata for Pan-GI donor-overlap
analysis.

## Atlases

- **Smillie 2019** (Single Cell Portal SCP259) — 366,650 colon mucosa
  cells, 18 UC + 12 healthy.
- **Kong 2023 UC subset** (GEO GSE214695) — 6 healthy + 6 active UC.
- **Mennillo 2024** (GEO accession to confirm in M1) — anti-integrin
  therapy time course; **subset to pre-treatment baseline samples only**.
  Verify ≥8 donors after subsetting; if fewer, swap to Garrido-Trigo 2023.
- **HCA Gut Cell Atlas** (cellxgene, Elmentaite 2021) — large-intestine
  subset only. M1 must verify zero donor overlap with the UC trio.
- **Pan-GI** (cellxgene, Oliver 2024) — large-intestine subset only.
  Substantial donor overlap with UC trio (integrated Smillie 2019,
  anchored Kong 2023).

## HGNC pin

Pin a single HGNC release (e.g., HGNC 2024-Q1) and remap every atlas.
Verify ≥95% of canonical UC GWAS hits (IL23R, JAK2, TYK2, NKX2-3,
ATG16L1) appear under the same symbol after harmonization.

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

## Pan-GI / HCA Gut raw counts verification (M1)

Verify `adata.layers['counts']` (or equivalent) contains raw counts for
each. If only normalized/integrated counts are accessible, drop the
comparator to stretch — scDRS requires raw counts.

## Output

- `data/atlases/<atlas>.h5ad` — standardized AnnData with HGNC-pinned symbols.
- `data/atlases/<atlas>_covariates.tsv` — scDRS covariate file.
- `data/atlases/donor_metadata/<atlas>_donor_metadata.csv` — donor attribution.
- `data/atlases/<atlas>_sce.rds` — SingleCellExperiment for seismicGWAS.
