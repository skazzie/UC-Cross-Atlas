<!--
Status (2026-06-13, handoff/laptop-validation): code-complete;
8 tests green (math + E2E on synthetic null tensors). Awaits HB
null draws (results/null_draws/*_nulls.npz) + regime-1 TSVs to
produce real combined p-values.
-->

# 07_regime2_meta — Brown's method per-cell-type meta-analysis

Implements regime 2 from PLAN.md / DECISIONS.md: combining cell-type
p-values across the three UC atlases via Brown's method (analytical
Kost-McDermott), with the cross-atlas correlation matrix estimated
empirically from null draws.

**Scope:** de Lange GWAS only, both methods (scDRS + seismicGWAS), both
granularity tiers. Four regime-2 outputs total.

## Inputs

- Per-(atlas, GWAS, tier) cell-type test statistics from regime 1
  (`results/scdrs/...`, `results/seismic/...`).
- Per-null-draw cell-type-level statistics serialized in M3:
  - scDRS Monte Carlo null aggregations (feather, indexed by
    `(null_draw_idx, cell_type)`).
  - seismicGWAS gene-Z permutation aggregations (feather, indexed by
    `(permutation_idx, cell_type)`).

## Pipeline

1. For each (method, tier) under de Lange, build the per-cell-type
   empirical correlation matrix: Pearson correlation of N = 1000 null
   statistic pairs `(T_c^A_i, T_c^B_i)` between atlas pairs.
2. Sanity check: off-diagonal entries should be positive. Investigate any
   cell type with near-zero or negative cross-atlas null correlation.
3. Apply edge-case fallback: cell types with null-statistic SD below the
   5th percentile in any atlas use the median cross-atlas correlation
   instead of their own.
4. Run `EmpiricalBrownsMethod::kostsMethod()` per cell type.
5. Handle fine-tier missing-cell-type cases (3/3, 2/3, 1/3) per
   DECISIONS.md.

## Stretch fallback

If empirical correlation matrices are pathological (negative or near-zero
off-diagonals), activate stretch #1: permutation-based covariance with
~300 additional scDRS runs. Decide at M5 only.

## Output

`results/regime2/<method>_<tier>_delange.tsv` with columns:
`cell_type`, `combined_pval`, `n_atlases_combined`, `correlation_fallback`
(boolean for the median-substitution edge case).

## Driver script

`code/07_regime2_meta/run_brown.py`. One (method, tier, gwas) per invocation.
Implements the Kost-McDermott analytical form of Brown's method inline
using Brown's (1975) polynomial covariance approximation — no R or extra
PyPI dependency. Matches the math in R's `EmpiricalBrownsMethod::kostsMethod()`.

CLI:

```bash
python code/07_regime2_meta/run_brown.py \
    --method scdrs \
    --tier broad \
    --gwas delange \
    --regime1-dir results/scdrs \
    --null-draws-dir results/null_draws \
    --out results/regime2/scdrs_broad_delange.tsv
```

Runs on the login node.
