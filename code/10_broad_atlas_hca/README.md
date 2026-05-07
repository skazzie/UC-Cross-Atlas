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
