<!--
Status (2026-06-13, handoff/laptop-validation): code-complete;
3 E2E tests green on synthetic fixtures (scDRS path, seismic+BH
path, missing-input clean-exit). Awaits HB regime-1 outputs under
both de Lange and Liu (Liu gated on DECISIONS 27(c) ancestry-LD
decision).
-->

# 09_cross_gwas — de Lange vs Liu within-atlas concordance

Concordance axis #3 from PLAN.md §"Five concordance axes": within each
(atlas, method, tier) combination, compare cell-type rankings under
de Lange 2017 vs Liu 2023.

**Scope:** 3 atlases × 2 methods × 2 tiers = 12 comparisons.

## Caveat

de Lange and Liu are both UC GWAS, with partially overlapping disease
cohorts. Cross-GWAS concordance is *expected* to be high — the interesting
case is cell types where they disagree. This is *not* a cross-trait or
cross-disease test (PLAN.md §"What this paper does not claim").

## Pipeline

For each (atlas, method, tier):

1. Load cell-type-level statistics under de Lange and Liu from
   `results/<method>/<atlas>_<gwas>_<tier>.tsv`.
2. Compute Spearman ρ on shared cell-type intersection (with min-50-cells
   filter), with bootstrap 95% CIs (1000 iters, percentile, seed=42).
3. Top-k Jaccard.
4. Cohen's κ on FDR-significance agreement (with κ@FDR<0.01 contingency
   for ≥80% saturation).

Uses the shared library in `code/06_concordance/metrics.py`.

## Output

`results/cross_gwas/cross_gwas_concordance.tsv` — long format with columns
`atlas`, `method`, `tier`, `spearman_rho`, `ci_lo`, `ci_hi`,
`jaccard_top5`, `jaccard_top10` (+ `jaccard_top20` at fine tier), `kappa`,
`kappa_threshold`, `n_sig_delange`, `n_sig_liu`, `n_common`.

## Driver script

`code/09_cross_gwas/run_cross_gwas.py`. CLI:

```bash
python code/09_cross_gwas/run_cross_gwas.py \
    --atlases smillie garrido_trigo taurus \
    --methods scdrs seismic \
    --tiers broad fine \
    --scdrs-dir results/scdrs \
    --seismic-dir results/seismic \
    --out results/cross_gwas/cross_gwas_concordance.tsv
```

Runs on login node.
