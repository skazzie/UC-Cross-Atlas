# 04_seismic — seismicGWAS regime-1 pipeline

Lai et al. 2025 (*Nature Communications*) seismicGWAS implementation. R-side
package (`devtools::install_github("ylaboratory/seismicGWAS")`). Runs on
each atlas under each GWAS at both granularity tiers.

## Pipeline (per atlas)

```r
library(seismicGWAS)
library(SingleCellExperiment)

sce <- readRDS("data/atlases/<atlas>_sce.rds")

# Once per atlas — GWAS-independent
spec_broad <- calc_specificity(sce, ct_label_col = "cell_type_broad")
spec_fine  <- calc_specificity(sce, ct_label_col = "cell_type_fine")

# Once per (atlas, GWAS, granularity)
res_broad_delange <- get_ct_trait_associations(spec_broad, magma_z_delange)
res_fine_delange  <- get_ct_trait_associations(spec_fine,  magma_z_delange)
res_broad_liu     <- get_ct_trait_associations(spec_broad, magma_z_liu)
res_fine_liu      <- get_ct_trait_associations(spec_fine,  magma_z_liu)
```

## Confounders

`get_ct_trait_associations()` defaults — gene length (log-transformed),
gene-gene LD score, transcript count. **Verify in M1 by inspecting package
source**; if defaults differ, override and apply this explicit set.
Documented in DECISIONS.md.

## Specificity score serialization (v2-setup decision)

Save specificity as long-format `.rds` (or feather) per
(atlas, granularity, cell_type, gene) into `specificity_long/`. Schema:

```r
data.frame(
  atlas       = "smillie",
  granularity = c("broad", "fine"),
  cell_type   = "T_cell",
  gene        = "IL23R",
  specificity = 0.42
)
```

**Stored alongside regression results, NOT embedded in them.** Required
because v2's cross-tissue claim compares UC-T-cell specificity to
RA-T-cell specificity for the same gene set; embedding inside
regression-result objects forces full re-run later.

## Permuted-Z-score null draws for Brown's method

M = 1000 permutations of the gene-Z-score vector per (atlas, GWAS, tier).
Re-run `get_ct_trait_associations()` under each permutation, save
per-permutation cell-type-level test statistics to feather indexed by
`(permutation_idx, cell_type)`. ~30 min per (atlas, GWAS, tier);
~5 hours total for the locked core. Inputs to Brown's empirical correlation
matrix in M5.

## Test-retest gate

ρ ≥ 0.999 (deterministic method given fixed cell-type labels and expression
matrix; deviation indicates a bug).

## Inputs

- `data/atlases/<atlas>_sce.rds` (SingleCellExperiment with logcounts).
  Build from h5ad via `zellkonverter` or `anndata2ri`.
- `data/gwas/<gwas>_gene_z.tsv` — MAGMA gene-level Z-score table from
  `code/01_magma/` (MHC-excluded, autosomes only).

## Output

- `results/seismic/<atlas>_<gwas>_<tier>.tsv` — regression results.
- `code/04_seismic/specificity_long/<atlas>_specificity.feather` — long-format specificity.
- `results/seismic/permutations/<atlas>_<gwas>_<tier>_permnulls.feather` —
  per-permutation null draws for Brown's.
