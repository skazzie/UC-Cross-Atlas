# 05_magma_geneprop — MAGMA gene-property sanity track

Watanabe-style MAGMA gene-property regression as a sanity check on the
upstream pipeline. Validates that MAGMA gene-level Z-scores produce
reasonable cell-type rankings on Smillie at broad tier under de Lange.

**Scope:** Smillie only, broad tier, de Lange only. One sanity run.

## Why this is locked-core sanity

scDRS and seismicGWAS both consume MAGMA gene-Z as upstream input. If
MAGMA gene-property regression on Smillie does not produce a sensible
cell-type ordering for UC, the entire pipeline has a bug — investigate
before proceeding to M3 Gate 1.

## Pipeline

For each broad cell type in Smillie:

1. Compute mean log-normalized expression per gene per cell type.
2. Regress MAGMA gene-Z on the cell-type expression vector with
   confounders (gene length, gene-gene LD, transcript count).
3. Cell types with positive regression coefficients and FDR < 0.05 are
   "MAGMA gene-property-prioritized" for UC.

Reference: Watanabe K et al. 2019, *Nat Commun* 10:3222.

## Expected output

UC-relevant cell types should rank near the top — mature enterocytes,
plasma cells, T cells, monocyte/macrophages, DCs. If the top cell types
are biologically unrelated to UC (e.g., smooth muscle, endothelium), the
upstream pipeline has a bug.

## Output

`results/magma_geneprop/smillie_broad_delange.tsv` — one row per broad
cell type with regression coefficient, p-value, BH-FDR.

## Driver script

`code/05_magma_geneprop/run_magma_geneprop.py`. Smillie / broad / de Lange.
CLI:

```bash
python code/05_magma_geneprop/run_magma_geneprop.py \
    --atlas smillie \
    --h5ad-path data/atlases/smillie.h5ad \
    --magma-z results/magma/delange_gene_z.tsv \
    --tier broad \
    --gwas delange \
    --out results/magma_geneprop/smillie_broad_delange.tsv
```

Runs on the login node (single-pass OLS; no SLURM wrapper needed).
