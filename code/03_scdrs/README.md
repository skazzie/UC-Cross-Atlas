# 03_scdrs — scDRS regime-1 pipeline

Zhang et al. 2022 (*Nature Genetics*) scDRS. Pre-commits to all the
configuration choices in DECISIONS.md.

## Pipeline (per atlas, per GWAS)

```bash
# Once per (atlas, GWAS) — generates per-cell scores (~30-45 min, 15-20 GB RAM)
scdrs compute-score \
    --h5ad-file data/atlases/<atlas>.h5ad \
    --gs-file data/gwas/<gwas>_top1000.gs \
    --cov-file data/atlases/<atlas>_covariates.tsv \
    --flag-filter-data True \
    --flag-raw-count True \
    --n-ctrl 1000 \
    --out-folder results/scdrs/<atlas>_<gwas>/ \
    --random-seed 42

# Group analysis at broad (post-processing, fast)
scdrs perform-downstream \
    --score-file results/scdrs/<atlas>_<gwas>/UC.score.gz \
    --group-analysis cell_type_broad \
    --out-folder results/scdrs/<atlas>_<gwas>_broad/

# Group analysis at fine (same per-cell scores)
scdrs perform-downstream \
    --score-file results/scdrs/<atlas>_<gwas>/UC.score.gz \
    --group-analysis cell_type_fine \
    --out-folder results/scdrs/<atlas>_<gwas>_fine/
```

`compute-score` is the heavy step. Per-cell scores feed both granularity
outputs. **Fine-tier promotion adds essentially zero compute at the scoring
level.**

## Cross-atlas headline metric uses Z-scores, not p-values

scDRS p-values are computed against an atlas-specific Monte Carlo null
distribution; same biology in two atlases produces different p-values
purely from compositional differences. Headline cross-atlas Spearman ρ
uses **cell-type-level Z-scores** (mean per-cell Z within cell type).
Within-atlas FDR-significance counts and Cohen's κ stay on p-values.

## Save per-null-draw cell-type Z-scores

scDRS internally generates ~1000 null gene-set Monte Carlo draws and
computes per-cell null scores. Aggregate to cell-type level (mean per-cell
Z within cell type, per null draw, per atlas) and serialize as feather
files indexed by `(null_draw_idx, cell_type)`. These are the inputs to
Brown's method correlation-matrix estimation in M5.

## Covariate set (locked)

`<atlas>_covariates.tsv` columns:
- `cell_id`
- `log_n_genes` — log10 of detected genes per cell
- `log_n_counts` — log10 of total counts per cell
- `donor_<id>` — one-hot encoding of donor
- `sample_<id>` — one-hot encoding of sample, if multiple samples per donor
- `sex_<F/M/unknown>` — one-hot encoding if metadata exists; document
  omission per atlas in DECISIONS.md if missing

## All-cells policy

`compute-score` and group analysis run on **all cells** in the atlas,
regardless of disease status. Smillie/Garrido-Trigo/Mennillo all include both UC
and reference samples; running on the full mucosa is the right answer to
the biological question and is consistent with Zhang et al. 2022's
standard practice.

## Random seeds

- Headline run: `seed = 42`.
- Test-retest: seeds 1, 2, 3.

## Pan-GI / HCA Gut counts policy

Use raw count layer (`adata.layers['counts']` from cellxgene Discover, or
`adata.raw.X`) with **published cell-type labels** — no re-clustering.

## MHC-included sensitivity run

One supplementary run on Smillie × de Lange with the MHC-included `.gs`
file (built by `code/01_magma/make_scdrs_gs.py --keep-mhc`). Confirms MHC
inclusion shifts antigen-presenting cell rankings as expected.

## Output

- `results/scdrs/<atlas>_<gwas>/UC.full_score.gz` — per-cell scores.
- `results/scdrs/<atlas>_<gwas>_<tier>/UC.scdrs_group.cell_type_<tier>` —
  group analysis.
- `results/scdrs/<atlas>_<gwas>/null_aggregations_<tier>.feather` —
  per-null-draw cell-type Z aggregations for Brown's.
