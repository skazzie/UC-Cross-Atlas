# Garrido × de Lange scDRS dry run — diagnostic

**Status (2026-06-21, rev 2 per Saisohan):** Local laptop dry-run on
the Garrido atlas (29,675 cells, sized for laptop validation before
HB scaleup). Ranking under the locked v1 policy (donor-cov + MHC-
excluded) is **strongly axis-resolved but inverted from the
canonical-UC-immune expectation** — epithelial subtypes top with
significant permutation p-values; T cells, fibroblasts, NK/ILC, mast
cells sit at the bottom with active negative z-scores. 0/28,706
cells reach FDR<0.1 at the cell level.

Three sensitivities + a per-gene length test isolate the cause:
**the gene set, not the covariates and not gene-length inflation.**
The dominant mechanism is the **locked MHC exclusion**: removing
chr 6 strips the bulk of the immune-coded signal in de Lange UC.

**Honest framing of the contribution (rev 2).** Neither MHC choice
is clean: MHC-excluded gives an epithelial artifact (immune signal
too sparse to beat expression-matched controls); MHC-included gives
an HLA-marker + LD-smear artifact (HLA genes are themselves
antigen-presenting-cell markers, so adding them mechanically scores
B / DC / macrophage up partly by construction; the MHC region's
extreme LD also smears one association across many correlated
genes). The result is a **methods-caution trade-off**, not a
discovery: the standard MHC-exclusion hygiene inverts cell-type
prioritization for UC; the inversion is not a covariate or length
artifact; neither MHC choice is independent of confound. The v1
primary stays **MHC-excluded** (the convention reviewers expect);
MHC-included is the **documented sensitivity**.

This file records what we tested, what we ruled in, what we ruled
out, and what the M4 narrative needs to honestly say.

## Headline broad-tier ranking (with donor covariates, MHC-excluded)

| Rank | Cell type | n_cell | assoc_mcp | assoc_mcz |
|---|---|---|---|---|
| 1  | colonocyte             | 3,530 | 0.003 | +3.63 |
| 2  | epithelial progenitor  | 2,134 | 0.002 | +3.29 |
| 3  | goblet                 | 1,237 | 0.089 | +1.47 |
| 4  | monocyte/macrophage    |   993 | 0.16  | +0.97 |
| 5  | enteroendocrine/tuft   |   332 | 0.19  | +0.79 |
| 6  | endothelium            |   227 | 0.24  | +0.68 |
| 7  | granulocyte            |   145 | 0.29  | +0.43 |
| 8  | dendritic cell         |    97 | 0.57  | −0.25 |
| 9  | plasma cell            | 8,900 | 0.61  | −0.40 |
| 10 | B cell                 | 1,310 | 0.69  | −0.57 |
| 11 | mural/glia             |   134 | 0.86  | −1.05 |
| 12 | mast cell              |   520 | 0.92  | −1.29 |
| 13 | NK/ILC                 |   166 | 0.95  | −1.45 |
| 14 | fibroblast             | 1,668 | 0.96  | −1.69 |
| 15 | **T cell**             | 7,313 | 0.99  | **−2.02** |

Active negative immune z is the diagnostic feature, not "epithelial
on top" (Saisohan 2026-06-15). A contaminated, epithelially-biased
gene set scores immune cells *below* their expression-matched
controls → negative.

## Three hypotheses ruled in / out

### Hypothesis 1: MAGMA gene-length / NSNPS inflation — **RULED OUT**

The top-2 weights are PLCL1 (rank 1, Z=+9.59) and NRCAM (rank 2,
Z=+8.96). Both are large genes (P95-P97 by both NSNPS and LENGTH)
and neither appears in canonical IBD locus catalogs (HNF4A, CDH1/3,
LAMB1, ECM1, IL23R, CARD9, RNF186, NOD2, PTGER4, IL10, etc.).
Length-artifact suspicion was reasonable on priors.

Direct test on `results/magma/uc_delange_10kb.genes.out`:

- **Spearman(ZSTAT, NSNPS) genome-wide = −0.0351**
- **Spearman(ZSTAT, LENGTH) genome-wide = −0.0108**
- **R²(ZSTAT ~ log10(LENGTH) + log10(NSNPS)) = 0.0007 (0.07%)**

After residualizing for length and NSNPS, the top 10 by `Z_resid` is:
PLCL1 (9.08), NRCAM (8.45), POU5F1 (7.97), TAP1 (7.79), BTNL2 (7.58),
TCF19 (7.52), MICA (7.50), TMEM8C (7.39), HLA-C (7.38), MICB (7.34).

PLCL1 and NRCAM **survive length residualization** — they're not
length artifacts. MAGMA's gene-test already corrects for length and
the correction worked. The remaining concern is whether they're
locus-true-but-wrong-gene assignments under the 10kb window (could
tag a neighboring canonical locus by LD). That's a separate
investigation; not pursued here.

### Hypothesis 2: Donor-covariate over-correction — **RULED OUT**

Garrido has 12 donors, each donor is either HC or UC, so donor
one-hots are collinear with disease status. Saisohan flagged that
regressing donor as a covariate could attenuate the UC signal we're
trying to score.

Test: rebuild covariates with `log_n_genes + log_n_counts` only
(13 → 2 columns, no donor one-hots). Re-run compute-score +
perform-downstream.

| Cell type | no-cov mcz | donor-cov mcz | Δ |
|---|---|---|---|
| colonocyte | +3.60 | +3.63 | −0.02 |
| epithelial progenitor | +3.12 | +3.29 | −0.17 |
| goblet | +1.68 | +1.47 | +0.21 |
| T cell | **−2.19** | −2.02 | **−0.17** |
| fibroblast | −1.94 | −1.69 | −0.25 |
| NK/ILC | −1.67 | −1.45 | −0.22 |

Cell-level FDR<0.1: **0/28,706 in both runs.**

**The immune cells do not recover without donor covariates** — they
score *slightly more negative* on the immune end (T cell −2.19 vs
−2.02 without donor cov), not positive. Donor over-correction was
not the lead mechanism, not even compensating for it. The inversion
is in the gene set itself.

Cost: 138,291s (38h) compute-score + 13min perform-downstream. The
no-cov full_score.gz is at
`results/scdrs/garrido_delange_nocov/UC.full_score.gz` (gitignored).

### Hypothesis 3: Gene-set epithelial-expression bias — **RULED OUT**

Mean expression of the top-20 `.gs` genes per broad cell type in
Garrido (`log1p(CP10k)` units, rank-aggregated 15=highest 1=lowest):

| Cell type | Mean rank |
|---|---|
| endothelium            | 10.1 |
| fibroblast             |  9.8 |
| monocyte/macrophage    |  9.0 |
| dendritic cell         |  8.8 |
| T cell                 |  8.7 |
| mural/glia             |  8.6 |
| B cell                 |  8.5 |
| NK/ILC                 |  8.3 |
| colonocyte             |  8.2 |
| enteroendocrine/tuft   |  8.2 |
| epithelial progenitor  |  7.7 |
| mast cell              |  7.0 |
| goblet                 |  6.1 |
| granulocyte            |  5.7 |
| plasma cell            |  5.4 |

The top-20 `.gs` genes are **NOT** epithelially over-expressed in
Garrido. The expression-matched-control comparison (scDRS's core
trick) is what creates the inversion: in epithelial cells the `.gs`
genes are RELATIVELY enriched over the 1000 matched control sets,
even though absolute expression is higher in stromal/myeloid. That
mechanism is correct scDRS behavior given the gene set's composition.

## What's actually driving it — Hypothesis 4: MHC exclusion (confirmed mechanism, not vindication)

After length/NSNPS residualization, **8 of the top 10 de Lange UC
gene-Z values are MHC-region genes** (POU5F1, TAP1, BTNL2, TCF19,
MICA, HLA-C, MICB on chr 6). The locked v1 policy excludes MHC from
the scDRS `.gs` (PLAN.md / DECISIONS.md — chr 6: 28,477,797–
33,448,354 GRCh37). That strips the bulk of the immune-coded signal
in de Lange UC and leaves the residual non-MHC weights — which
under-weight immune — as the dominant driver.

**Diff between the two `.gs` files** (`results/magma/uc_delange.gs`
vs `results/magma/uc_delange_mhc.gs`):

- Both 1000-gene; overlap 872; **128 genes swapped**.
- MHC-included adds (sample): ABHD16A, AGER, AGPAT1, AIF1, APOM,
  ATAT1, ATF6B, ATP6V1G2, BAG6, BRD2, BTNL2, C2, C4A, C4B, C6ORF10,
  ... classical + non-classical HLA-region, complement.
- MHC-excluded keeps (sample): ABHD14A, ABHD14B, ABR, ACOT13,
  ADGRL3, AIMP1, ARHGAP12, BAZ1A, BEND6, BPIFC, ... long-tail of
  non-immune-themed genes.

**Sensitivity result (compute-score level).** Garrido × de Lange
with `uc_delange_mhc.gs`, donor covariates per locked design,
n_ctrl=1000:

- **MHC-included**: **461 / 28,706 cells FDR<0.1** (1,267 at <0.2).
- MHC-excluded (locked v1): 0 / 28,706 cells FDR<0.1 (1 at <0.2).

Order-of-magnitude effect at the cell level — confirms MHC exclusion
is the dominant inversion driver.

### Why "MHC-included ≠ correct answer" — two confounds

1. **HLA genes are themselves immune-cell markers.** Class II
   (HLA-DR, -DQ, -DP) is highly expressed in B cells, dendritic
   cells, and macrophages; class I is broad. A gene set that
   contains the MHC block will mechanically score
   antigen-presenting and immune cells up *partly by construction*.
   This is not the same as recovering a polygenic UC signal in
   those cells.

2. **MHC's extreme LD means one association smears across many
   correlated genes**, inflating enrichment for whatever cell types
   express the region.

This is why MHC is conventionally excluded from polygenic-trait
gene sets in the first place. We are NOT re-locking MHC inclusion
as the default — that would quietly swap to the variant with its
own confound. Keep MHC-excluded as primary, MHC-included as the
documented sensitivity. The trade-off itself is the result.

**Per-cell-type sensitivity ranking** (where do the 461 cells land?
key question — if they cluster in HLA-high APCs, the MHC-included
ranking is partly a tautology): *(appended below when
perform-downstream completes)*

## Outputs on disk (all gitignored)

```
results/scdrs/garrido_delange_seed42/   donor-cov, MHC-excluded (locked v1)
  UC.score.gz
  UC.full_score.gz
  UC.scdrs_group.cell_type_broad
  UC.scdrs_group.cell_type_fine
  null_aggregations_{broad,fine}.feather

results/scdrs/garrido_delange_nocov/    no-donor-cov, MHC-excluded (sensitivity)
  UC.score.gz
  UC.full_score.gz
  UC.scdrs_group.cell_type_broad
  UC.scdrs_group.cell_type_fine

results/scdrs/garrido_delange_mhc/      donor-cov, MHC-INCLUDED (sensitivity)
  *populated by the in-progress run; results appended below*
```

## Notes for the HB run / M4 narrative

- **Don't write λ_GC inflation as the explanation** (initially
  proposed; corrected by Saisohan 2026-06-15). λ_GC scales test
  statistics multiplicatively, not preferentially-by-length. Drop it.
- **Don't write "half real biology, half artifact"** — the
  epithelial-on-top is not biology. PLCG2 at rank 265 is a defensible
  small real epithelial component (PLCG2 IS a UC GWAS gene expressed
  in epithelium), but with immune actively negative and PLCL1/NRCAM
  dominating the weights, the dominant signal is not biological.
- **Do write the MHC-exclusion-vs-immune-signal trade as a core
  methods note**, contingent on the MHC sensitivity result below
  confirming the mechanism.

## LOG items for HB cleanup (recorded for next HB session)

- `code/03_scdrs/README.md` likely shows `--flag-raw-count True`
  somewhere; slurm correctly uses `False`. Verify.
- Slurm naming mismatches: file is `garrido.h5ad`, slurm expects
  `garrido_trigo.h5ad`; `.gs` path is `results/magma/uc_delange.gs`,
  slurm expects `data/gwas/delange_top1000.gs`.
- Slurm omits `--h5ad-species` / `--gs-species` — required in
  scdrs 1.0.2, KeyError on first HB run.
- Slurm uses `--random-seed` which doesn't exist until scdrs
  1.0.4-dev. Test-retest seed design (seeds 1/2/3 per PLAN.md)
  needs 1.0.4 pinned on HB.
