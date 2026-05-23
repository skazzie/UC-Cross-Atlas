# 10_broad_atlas_hca — HCA Gut comparator (independent broad atlas)

Concordance axis #5a from PLAN.md §"Five concordance axes": broad-atlas
substitution test using HCA Gut Cell Atlas (Elmentaite 2021), restricted
to the large-intestine subset.

**Role:** independent broad comparator. M1 must verify zero donor overlap
with the UC trio (Smillie 2019 / Kong 2023 / Mennillo 2024). If overlap
exists, HCA Gut becomes a second integration-pipeline-robustness comparator
(not an independent comparator) and the framing in §2.4.3 of the plan
changes accordingly.

**Scope:** scDRS + seismicGWAS, 2 GWAS each = 4 runs total.

## Counts policy

Use raw count layer (`adata.layers['counts']` from cellxgene Discover, or
`adata.raw.X`) with **published cell-type labels** — no re-clustering.
Methods text limits results to "broad-comparator robustness check," not
independent validation.

## Caveat

HCA Gut is healthy-only; it lacks UC-specific inflammation-induced cell
states (BEST4+ disease-state enterocytes, inflammatory fibroblasts).
Concordance metrics computed on the **shared cell-type intersection**
between HCA Gut and each UC trio atlas (with min-50-cells threshold).

## Output

`results/broad_atlas/hca_gut_concordance.tsv` — comparison of HCA Gut
cell-type rankings to each UC trio atlas under each method × GWAS.

## v1 filter chain (DECISIONS.md 2026-05-20 (6/7))

Defined in `code/02_atlas_prep/load_hca_gut.py`:

```
Age_group in {Adult, Adult_MLN}
tissue in {ascending colon, caecum, colon, descending colon,
           large intestine, rectum, sigmoid colon, transverse colon}
```

Excludes the fetal, pediatric, and pediatric-IBD subjects (HCA Gut is
used as a v1 *healthy reference*; it does not need UC samples for its
cell-type prioritization role). Expected post-filter: 30k-70k cells.

## scDRS covariates

`assay`, `batch`, `Fraction`, `sex`.

`Fraction` is **mandatory**: HCA Gut has 5 sorting strategies
(`SC`, `SC-45N`, `SC-45P`, `SC-EPCAMN`, `SC-EPCAMP`) that strongly
determine cell-type composition. Omitting this covariate would conflate
sort-induced expression differences with biology.

## v1 sensitivity

**No-Crohn run** — re-run with `disease == "Crohn disease"` cells
removed. Tests whether residual Crohn signal in this broad reference
affects UC cell-type prioritization.
Loader: `code/02_atlas_prep/load_hca_gut.load_hca_gut_no_crohn()`.

## Important caveat

Cell-type proportions in HCA Gut are artificially controlled by the
`Fraction` column (sort strategy) and are NOT biologically interpretable.
scDRS per-cell scoring is unaffected, but DO NOT report cell-type
abundance from HCA Gut in figures.

## Driver script

`code/10_broad_atlas_hca/run_hca_comparison.py`. CLI:

```bash
python code/10_broad_atlas_hca/run_hca_comparison.py \
    --hca-results-base results/hca_gut \
    --hca-no-crohn-base results/hca_gut_no_crohn \
    --uc-atlases smillie garrido_trigo mennillo \
    --gwas delange liu \
    --methods scdrs seismic \
    --tiers broad fine \
    --scdrs-dir results/scdrs \
    --seismic-dir results/seismic \
    --out-headline results/broad_atlas/hca_gut_concordance.tsv \
    --out-no-crohn results/broad_atlas/hca_gut_concordance_no_crohn.tsv
```

Loads pre-existing scDRS/seismicGWAS results — does NOT run the methods
itself. Runs on login node.
