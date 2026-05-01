# 07_concordance — cross-atlas concordance metrics

Implements the three metrics from spec §2.5 / §2.7:

- **Spearman ρ** between cell-type rank vectors (with bootstrap 95% CI)
- **Top-k Jaccard** for k = 5 and k = 10
- **Cohen's κ** on FDR-significance calls, with marginals reported alongside

`metrics.py` is the library. `compute_concordance.py` is a CLI that ingests
per-(atlas, method) cell-type p-value/q-value tables and emits the headline
table from §2.5.

## Input format

For each (atlas, method) pair, one TSV with these columns:

```
cell_type    pval    qval
T cell       1.2e-08 1.5e-07
B cell       0.003   0.011
...
```

scDRS's `.scdrs_group.cell_type` output uses different column names — adapt
with a one-line pandas rename or pass `--cell-type-col`/`--pval-col`/`--qval-col`.

## Run

```bash
python code/07_concordance/compute_concordance.py \
    --input \
        results/scdrs/smillie/UC.cell_type.tsv:smillie:scdrs \
        results/scdrs/kong/UC.cell_type.tsv:kong:scdrs \
        results/scdrs/mennillo/UC.cell_type.tsv:mennillo:scdrs \
        results/seismic/smillie/UC.cell_type.tsv:smillie:seismic \
        results/seismic/kong/UC.cell_type.tsv:kong:seismic \
        results/seismic/mennillo/UC.cell_type.tsv:mennillo:seismic \
    --out results/concordance/headline_table.csv
```

Output is one row per (method, atlas_a, atlas_b) pair, with all three metrics
plus marginals — directly suitable for the §2.5 headline summary table.

## Tests

```bash
pytest code/07_concordance/test_metrics.py -v
```

13 tests covering perfect agreement, perfect disagreement, partial overlap,
non-shared cell types, the all-significant κ degenerate case, and bootstrap CI
sanity.

## Notes

- Pass **p-values or q-values** (smaller = stronger), not Z-statistics — top-k
  Jaccard uses smallest-first ranking.
- Cohen's κ is unstable when marginals are extreme. The function returns `NaN`
  when every cell type is significant in both atlases (§2.5 edge case); always
  print `n_sig_a` and `n_sig_b` next to κ in tables.
- The bootstrap CI here resamples cell types only — it does not capture
  cell-level sampling noise. For headline numbers, also run cell-level
  bootstraps (resample cells, re-run scDRS, re-compute ρ; ~100 iterations
  per spec §2.7). That code lives next to the scDRS runner since it has to
  re-invoke scoring.
