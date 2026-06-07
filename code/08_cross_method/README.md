# 08_cross_method — scDRS vs seismicGWAS within-atlas concordance

Concordance axis #2 from PLAN.md §"Five concordance axes": within each
(atlas, GWAS, tier) combination, compare scDRS Z-scores to seismicGWAS
regression coefficients on the shared cell-type set.

**Scope:** 3 atlases × 2 GWAS × 2 tiers = 12 comparisons.

## Caveat

Both methods consume the same upstream MAGMA gene-level Z-scores. Cross-method
agreement is partly mechanical; cross-method *disagreement* is the more
informative signal. Methods text states this explicitly.

## Pipeline

For each (atlas, GWAS, tier):

1. Load cell-type-level scDRS Z-scores (mean per-cell Z within cell type)
   from `results/scdrs/<atlas>_<gwas>/cell_type_<tier>/`.
2. Load seismicGWAS regression coefficients from
   `results/seismic/<atlas>_<gwas>_<tier>.tsv`.
3. Compute Spearman ρ on shared cell-type intersection (with min-50-cells
   filter), with bootstrap 95% CIs (1000 iters, percentile, seed=42).
4. Top-k Jaccard on cell types ranked by Z (scDRS) vs by coefficient
   (seismicGWAS).
5. Cohen's κ on FDR-significance agreement (BH-FDR within atlas, with
   κ@FDR<0.01 contingency for ≥80% saturation).

Uses the shared library in `code/06_concordance/metrics.py`.

## Output

`results/cross_method/cross_method_concordance.tsv` — long format with
columns `atlas`, `gwas`, `tier`, `spearman_rho`, `ci_lo`, `ci_hi`,
`jaccard_top5`, `jaccard_top10` (+ `jaccard_top20` at fine tier), `kappa`,
`kappa_threshold` (0.05 or 0.01 contingency), `n_sig_scdrs`,
`n_sig_seismic`, `n_common`.

## Driver script

`code/08_cross_method/run_cross_method.py`. CLI:

```bash
python code/08_cross_method/run_cross_method.py \
    --atlases smillie garrido_trigo taurus \
    --gwas delange liu \
    --tiers broad fine \
    --scdrs-dir results/scdrs \
    --seismic-dir results/seismic \
    --out results/cross_method/cross_method_concordance.tsv
```

Runs on login node (lightweight aggregation; no SLURM wrapper needed).
