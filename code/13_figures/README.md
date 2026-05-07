# 13_figures — publication figures

Figure-generation code for the paper. Notebooks live in `notebooks/`; the
final figure scripts that are reproducible CLI entry points live here.

## Headline figures (locked v1)

Per PLAN.md §"Five concordance axes" + M4 plan:

1. **3×3 atlas-pair Spearman ρ heatmap** with bootstrap 95% CIs — faceted
   by method × granularity × GWAS. The single primary analysis is the
   broad-tier × scDRS × de Lange panel.
2. **Cell-type forest plot for broad tier** (de Lange) with donor-LOO
   95% jackknife ranges. Both methods.
3. **Cell-type forest plot for fine tier** (de Lange) with donor-LOO ranges.
4. **Cross-method concordance comparison panel** with bootstrap CIs —
   scDRS vs seismicGWAS within atlas.
5. **Cross-GWAS concordance comparison panel** with bootstrap CIs —
   de Lange vs Liu within atlas.
6. **Regime 1 vs regime 2 comparison** — combined cell-type rankings under
   Brown's method vs single-atlas rankings.
7. **HCA Gut concordance panel** — independent broad-atlas comparator.
8. **Pan-GI dual-analysis panel** — concordance with vs without
   donor-overlap exclusion.

## Style

- ggplot2 + ggsci or matplotlib + seaborn; consistent colour palette across
  panels.
- All ρ values reported with bootstrap 95% CI brackets, not raw point
  estimates.
- All κ values reported alongside marginals (`n_sig_a`, `n_sig_b`).

## Output

`results/figures/figure_<n>_<short_name>.{pdf,png}`
