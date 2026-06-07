# 06_concordance — cross-atlas concordance metrics

Implements the locked v1 concordance bundle from PLAN.md §"Three statistical
metrics" / DECISIONS.md.

## What's locked

- **Headline metric: Spearman ρ on cell-type-level Z-scores (scDRS) /
  regression coefficients (seismicGWAS) — NOT p-values.** Tied ranks
  broken by average-rank (scipy default).
- **Bootstrap 95% CIs on every reported ρ:** 1000 iterations, percentile
  method, seed = 42 (locked). BCa deferred to revision.
- **Concordance computed on shared cell-type intersection per pair, with
  ≥50 cells in BOTH atlases.** Atlas-specific cell types reported
  separately, not entered into concordance.
- **Top-k Jaccard:** k = 5, 10 at broad tier; k = 5, 10, 20 at fine tier.
- **Cohen's κ on FDR-significance** with **marginal-saturation
  contingency**: if ≥80% of cell types pass FDR < 0.05 in both atlases,
  report κ @ FDR < 0.01 as the headline κ instead.

This module is the shared library used by all five concordance axes:
`07_regime2_meta`, `08_cross_method`, `09_cross_gwas`, `10_broad_atlas_hca`,
`11_broad_atlas_pangi` — and by this directory's own `compute_concordance.py`
for the cross-atlas axis.

## Files

- `metrics.py` — library: `spearman()`, `top_k_jaccard()`,
  `fdr_concordance()`, `bootstrap_spearman_ci()`, `concordance()` (the
  bundle).
- `compute_concordance.py` — CLI driver. Loads per-(atlas, method, GWAS,
  tier) cell-type test-statistic tables; emits the long-format headline
  table.
- `test_metrics.py` — pytest suite covering perfect/partial/disjoint
  agreement, marginal-saturation contingency, min-cell-count filter,
  fine-tier Jaccard k=20, locked-seed determinism.

## Input format

For each (atlas, method, GWAS, tier), one TSV with columns:

```
cell_type   score   pval   qval   n_cells
T cell      0.42    1.2e-08 1.5e-07 12480
B cell      0.18    0.003   0.011   8412
...
```

- `score` is the **larger-is-stronger** statistic for that method:
  - **scDRS:** mean per-cell Z within cell type (the
    `cell_type_z` / `mean_z` aggregate from the group analysis).
  - **seismicGWAS:** regression coefficient.
- `pval` and `qval` are within-atlas p- and q-values (BH-FDR).
- `n_cells` is the cell count per cell type in that atlas, used by the
  min-50-cells filter.

scDRS's `.scdrs_group.cell_type` output uses different column names —
adapt with a one-liner pandas rename or via the `--score-col`,
`--pval-col`, `--qval-col`, `--n-cells-col` flags.

## Run (cross-atlas axis)

```bash
python code/06_concordance/compute_concordance.py \
    --input \
        results/scdrs/smillie_delange_broad/UC.cell_type.tsv:smillie:scdrs:delange:broad \
        results/scdrs/garrido_trigo_delange_broad/UC.cell_type.tsv:garrido_trigo:scdrs:delange:broad \
        results/scdrs/taurus_delange_broad/UC.cell_type.tsv:taurus:scdrs:delange:broad \
        results/scdrs/smillie_delange_fine/UC.cell_type.tsv:smillie:scdrs:delange:fine \
        results/scdrs/garrido_trigo_delange_fine/UC.cell_type.tsv:garrido_trigo:scdrs:delange:fine \
        results/scdrs/taurus_delange_fine/UC.cell_type.tsv:taurus:scdrs:delange:fine \
        results/scdrs/smillie_liu_broad/UC.cell_type.tsv:smillie:scdrs:liu:broad \
        ... \
        results/seismic/smillie_delange_broad.tsv:smillie:seismic:delange:broad \
        ... \
    --out results/concordance/cross_atlas_table.csv
```

Output is one row per (method, GWAS, tier, atlas_a, atlas_b) comparison
with all metrics + bootstrap CI + saturation flag — directly suitable
for the headline figure (3×3 atlas-pair Spearman ρ heatmap, faceted by
method × granularity × GWAS).

## Tests

```bash
pytest code/06_concordance/test_metrics.py -v
```

## Notes

- This module's API is **larger-is-stronger** by default. p-value-style
  ranking is supported via `larger_is_stronger=False` on `top_k_jaccard`,
  but headline metrics in the paper are Z-score-based.
- The bootstrap here resamples cell types only — it does not capture
  cell-level sampling noise. Donor-level uncertainty is reported
  separately as **LOO jackknife ranges** (see PLAN.md §"Donor-LOO
  uncertainty intervals"); LOO code lives next to the scDRS / seismicGWAS
  runners since it has to re-invoke scoring.
- Multiple-testing correction across the secondary comparison battery
  (~60 comparisons across the four axes outside the primary 3×3 panel)
  is applied **downstream**, in the figure-generation step, not here.
  This module emits raw uncorrected p-values; the caller applies BH-FDR
  across the appropriate set.
