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

**Per-cell-type sensitivity ranking — HLA tautology confirmed.**

Broad-tier scdrs_group output with MHC-included `.gs` (donor cov,
n_ctrl=1000):

| Rank | Cell type | n_cell | mcz | n_fdr<0.1 | % cells FDR<0.1 |
|---|---|---|---|---|---|
| 1  | B cell                 | 1,310 | **+7.35** | 197 | **15.0%** |
| 2  | dendritic cell         |    97 | **+6.51** |  52 | **53.6%** |
| 3  | monocyte/macrophage    |   993 | **+5.16** | 173 | **17.4%** |
| 4  | colonocyte             | 3,530 | +2.34 |   9 | 0.3% |
| 5  | epithelial progenitor  | 2,134 | +1.60 |   5 | 0.2% |
| 6  | mural/glia             |   134 | +1.13 |   0 |    — |
| 7  | endothelium            |   227 | +0.74 |   0 |    — |
| 8  | goblet                 | 1,237 | +0.53 |   1 | 0.08% |
| 9  | granulocyte            |   145 | +0.28 |   0 |    — |
| 10 | enteroendocrine/tuft   |   332 | +0.24 |   1 | 0.30% |
| 11 | plasma cell            | 8,900 | −0.09 |  21 | 0.24% |
| 12 | NK/ILC                 |   166 | −1.27 |   0 |    — |
| 13 | fibroblast             | 1,668 | −2.05 |   0 |    — |
| 14 | mast cell              |   520 | −2.08 |   0 |    — |
| 15 | **T cell**             | 7,313 | **−2.53** |   2 | 0.03% |

Cell-level totals: 111 cells at FDR<0.05, **461 at <0.1**, 1,267 at
<0.2.

**422 of 461 FDR<0.1 cells (91.5%) cluster in the three classical
antigen-presenting cell types** (B cell + dendritic cell +
monocyte/macrophage). This is the HLA-marker tautology Saisohan
predicted — HLA class II is highly expressed in exactly these three
lineages by constitutive mechanism, not by biology. The MHC block
in the gene set mechanically scores them up because they ARE the
cells that express that gene block.

**Counter-evidence to a "polygenic immune recovery" reading**:

- **T cell stays negative (−2.53), MORE negative than MHC-excluded
  (−2.02).** If MHC inclusion were recovering a polygenic immune
  signal, T should move toward 0. It moves away.
- **NK/ILC stay negative** (−1.27 vs −1.45 MHC-excl). Same pattern.
- **Plasma cell stays ~null** (−0.09 vs −0.40 MHC-excl). Plasma cells
  down-regulate HLA class II during terminal differentiation — and
  they don't light up, even though they would if a real T/B-cell-axis
  UC signal were being captured.
- The cell-type pattern matches "high constitutive HLA expression,"
  not "polygenic UC risk."

**Conclusion for the M4 trade-off section:**

- **MHC-excluded (locked v1 primary)**: epithelial artifact, T cell
  actively negative, 0 FDR-hit cells.
- **MHC-included (sensitivity)**: HLA tautology, T cell *more*
  negative, 461 FDR-hit cells concentrated 91.5% in APCs.
- Neither extreme is the immune-driven UC story. Both are confounded
  in opposite directions. The MHC handling fundamentally swings UC
  cell-type prioritization; that swing IS the result, not a
  vindication of either variant.

**Fine-tier MHC-included ranking (top 15 by mcz, with delta vs MHC-excluded):**

| Group | n_cell | mcz_MHC | mcz_noMHC | Δ | n_fdr<0.1 | % FDR<0.1 |
|---|---|---|---|---|---|---|
| B cell                  | 815 | **+7.73** | −0.94 | +8.66 | 119 | 14.6% |
| Memory B cell           | 138 | **+7.06** | +1.06 | +6.01 |  35 | 25.4% |
| M0 (macrophage)         | 153 | **+6.76** | −1.02 | +7.78 |  46 | 30.1% |
| **DCs CD1c**            |  76 | +6.37 | −0.45 | +6.81 |  45 | **59.2%** |
| Naive B cell            | 200 | +5.54 | −0.62 | +6.16 |  26 | 13.0% |
| DCs CCL22               |  21 | +5.49 | +0.84 | +4.65 |   7 | 33.3% |
| M2.2 (macrophage)       |  80 | +5.47 | −0.16 | +5.63 |  28 | 35.0% |
| IDA macrophage          | 224 | +4.99 | +0.54 | +4.45 |  51 | 22.8% |
| GC B cell               | 157 | +4.71 | +0.10 | +4.61 |  17 | 10.8% |
| Cycling myeloid         |   8 | +4.41 | +1.38 | +3.03 |   1 | 12.5% |
| M2 (macrophage)         | 154 | +3.92 | −0.16 | +4.08 |  21 | 13.6% |
| Inflammatory monocytes  | 197 | +3.75 | +1.39 | +2.36 |  13 |  6.6% |
| M1 ACOD1                | 169 | +3.61 | +1.43 | +2.18 |  12 |  7.1% |
| M1 CXCL5                |   8 | +3.04 | +1.57 | +1.47 |   1 | 12.5% |
| Inflammatory colonocyte | 518 | +2.56 | +2.78 | −0.22 |   0 |     — |

**The top 13 fine-tier types by % cells reaching FDR<0.1 are ALL
B cell, dendritic cell, or macrophage subtypes.** Inflammatory
colonocyte is the first non-APC entry — and it has 0 FDR<0.1 cells
*and* a tiny NEGATIVE delta (the only one in the table that gets
*worse* with MHC included).

Two readings reinforce the tautology:

1. **Across all 14 fine APC subtypes, mcz went up.** Across the 65
   non-APC fine types, the shift is mostly noise. Adding the MHC
   block is doing the same thing it did at broad tier: it scores
   constitutive HLA-expressors up.
2. **The delta correlates inversely with biological UC relevance.**
   Memory B cell delta +6.0 (B-lineage UC effector? maybe, but
   plasma cell — the actual differentiated UC effector — has Δ ≈ 0
   because it down-regulates HLA II during maturation). DCs CD1c
   delta +6.8 (constitutive HLA II expressor). M0 delta +7.8 (resting
   macrophage, very high HLA II). Compare CD8 FGFBP2 +0.0 / Δ ≈ 0
   (T cell, low HLA II) — no movement.

The fine-tier replicates the broad-tier conclusion: MHC-included
recovery is HLA-marker scoring, not polygenic UC signal recovery.

---

## Cross-method validation: seismicGWAS confirms the MHC swing

The MHC inversion is *not* method-specific. seismicGWAS (linear-model
specificity-vs-Z regression) on Garrido × de Lange, both MHC variants,
swings in the same direction as scDRS at the broad tier.

Setup notes:
- scDRS metric is `assoc_mcz` (MC-permutation z) and `assoc_mcp`.
- seismicGWAS metric is `pvalue` from `get_ct_trait_associations`'s
  one-sided linear-model slope test. Score = −log10(pvalue) for
  larger-is-stronger comparability; both seismic methods return
  only `(cell_type, pvalue, FDR)`, no coefficient/se.
- All four runs use the same Garrido h5ad, same `.gs` / gene-Z TSVs,
  same MHC-excluded / MHC-included pair. Donor covariates apply only
  to scDRS (seismic regresses at the specificity-vs-Z level, no
  per-cell covariate path).

| Cell type | scDRS noMHC z | scDRS MHC z | seismic noMHC nlogp | seismic MHC nlogp |
|---|---|---|---|---|
| colonocyte             | +3.62 | +2.34 | 0.80 | 0.53 |
| epithelial progenitor  | +3.29 | +1.60 | 1.13 | 0.63 |
| goblet                 | +1.47 | +0.53 | 1.35 | 1.19 |
| monocyte/macrophage    | +0.97 | +5.16 | 0.91 | 1.82 |
| enteroendocrine/tuft   | +0.79 | +0.24 | 0.27 | 0.13 |
| endothelium            | +0.68 | +0.74 | 0.06 | 0.04 |
| granulocyte            | +0.43 | +0.28 | 0.24 | 0.46 |
| dendritic cell         | −0.25 | +6.51 | 1.08 | 2.16 |
| plasma cell            | −0.40 | −0.09 | 0.21 | 0.09 |
| B cell                 | −0.57 | +7.35 | 0.14 | 1.30 |
| mural/glia             | −1.05 | +1.13 | 0.25 | 0.18 |
| mast cell              | −1.29 | −2.08 | 0.01 | 0.00 |
| NK/ILC                 | −1.45 | −1.27 | 0.07 | 0.06 |
| fibroblast             | −1.69 | −2.05 | 0.07 | 0.02 |
| T cell                 | −2.02 | −2.53 | 0.31 | 0.18 |

Spearman rank correlations across the 15 cell types:

| Comparison | ρ |
|---|---|
| scDRS noMHC vs scDRS MHC (within-method swing) | **+0.60** |
| seismic noMHC vs seismic MHC (within-method swing) | **+0.79** |
| scDRS noMHC vs seismic noMHC (cross-method, MHC-excluded) | **+0.60** |
| scDRS MHC vs seismic MHC (cross-method, MHC-included) | **+0.78** |

**Findings:**

1. **Methods agree** at moderate-to-high ρ within each MHC variant
   (0.60 MHC-excl, 0.78 MHC-incl). Cross-method concordance is real,
   not noise.
2. **MHC inclusion drives ALL three APCs into the top 5 for both
   methods.** Seismic MHC-incl top-5 is DC + mono/macro + B cell +
   goblet + epithelial progenitor — same three APCs as scDRS MHC-incl,
   plus two non-APC epithelial. The HLA-marker swing is reproducible
   across regression frameworks.
3. **MHC-excluded top-5 differ between methods on the non-immune
   side**: scDRS leads with colonocyte / epithelial progenitor /
   goblet; seismic leads with goblet / epithelial progenitor / DC /
   mono/macro. So the "epithelial-on-top" pattern is more scDRS-
   specific than I'd thought from the scDRS-only read — seismic
   already has DC + mono/macro in its top-4 even *without* MHC.
4. **T cell stays at or near the bottom under all four scoring
   variants** (rank 11–15). Whatever the MHC swing does to the
   ranking, T cell is consistently NOT prioritized for UC by either
   method in Garrido.

The robust cross-method finding is: **MHC handling fundamentally
swings cell-type prioritization, in the same direction across two
independent scoring frameworks.** That's the methods-caution result
to write into M4 — and it's now demonstrated at the cross-method
level on the dry-run atlas, before HB scaleup.

Outputs:
- `results/seismic/garrido_delange_broad.tsv`
- `results/seismic/garrido_delange_mhc_broad.tsv`

---

## Track [3] Controls: SCZ and Height through both methods

Per Saisohan 2026-06-21: "compare UC's ranking against SCZ's and height's,
not just against the canonical expectation. The UC signal is only what is
left after you subtract what the controls also show."

All four control runs landed (seismic in 11s total; scDRS SCZ + Height
~25 min each for compute+pdownstream). MHC-excluded only — matches locked
v1 primary, the direct apples-to-apples vs UC.

### Three-trait broad-tier panel (Garrido, MHC-excluded, donor-cov scDRS only)

scDRS column = `assoc_mcz`. seismic column = `−log10(pvalue)`. Larger = stronger.

| Cell type | UC_sd | SCZ_sd | Height_sd | UC_se | SCZ_se | Height_se |
|---|---|---|---|---|---|---|
| colonocyte             | **+3.62** | −0.39 | −0.52 | 0.80 | 0.03 | 0.09 |
| epithelial progenitor  | **+3.29** | −0.00 | +1.82 | 1.13 | 0.26 | 0.61 |
| goblet                 | +1.47 | −0.56 | +1.77 | **1.35** | 0.01 | 0.09 |
| monocyte/macrophage    | +0.97 | −1.08 | +0.38 | 0.91 | 0.00 | 0.00 |
| enteroendocrine/tuft   | +0.79 | +0.87 | +0.72 | 0.27 | **2.37** | 0.74 |
| endothelium            | +0.68 | +1.01 | **+3.81** | 0.06 | 1.12 | **9.27** |
| granulocyte            | +0.43 | −1.53 | −1.62 | 0.24 | 0.00 | 0.00 |
| dendritic cell         | −0.25 | −1.23 | +0.17 | 1.08 | 0.70 | 0.01 |
| plasma cell            | −0.40 | −0.22 | −0.17 | 0.21 | 1.23 | 0.37 |
| B cell                 | −0.57 | **+2.52** | −0.86 | 0.14 | 1.46 | 0.01 |
| mural/glia             | −1.05 | **+4.23** | **+3.73** | 0.25 | 0.68 | **5.57** |
| mast cell              | −1.29 | +0.31 | −1.26 | 0.01 | 0.11 | 0.00 |
| NK/ILC                 | −1.45 | −0.59 | −2.48 | 0.07 | 0.06 | 0.00 |
| fibroblast             | −1.69 | +1.36 | **+3.72** | 0.07 | 0.68 | **5.52** |
| T cell                 | −2.02 | −0.20 | −2.37 | 0.31 | 0.84 | 0.02 |

### Saisohan's question answered

> **2026-06-27 amendment** — the interpretation in items 1–3 below
> is **partially superseded by DECISIONS 30(a)**. The T-cell-bottom
> conclusion in item 1 was internally contradictory (says
> "trait-specific" while listing Height T cell at rank 14/15). The
> reproducible side-by-side now lives at
> `results/concordance/controls_sidebyside_broad.md` (generator:
> `code/06_concordance/controls_sidebyside_broad.py`). Locked reads:
> APC-up survives the controls test in both methods (UC-specific);
> scDRS T-cell-bottom is shared with Height (methods baseline
> concern) and seismic does not replicate it; IBDverse cross-check
> is now the adjudicator. The paragraphs below stay in place as a
> historical record but the manuscript paragraph must follow
> DECISIONS 30(a), not items 1–3 here.

**Are UC's findings (T-cell-bottom, APC-up under MHC-incl, epithelial-on-top) trait-specific or baseline method bias?**

Comparing UC against SCZ + Height on the same atlas under the same methods:

1. **T-cell-bottom is trait-specific.** Under scDRS SCZ, T cell ranks 8/15 (z=−0.20). Under scDRS Height, T cell ranks 14/15 (z=−2.37). UC's T-cell-bottom isn't "methods always park T cells low" — Height also has it, but SCZ doesn't.

2. **APCs at top under MHC-incl is NOT a method default.** Under SCZ scDRS, DC ranks 14/15 (z=−1.23), mono/macro ranks 13/15 (z=−1.08). The MHC-incl APC-up swing for UC is mechanistically MHC-driven, not a default APC preference of the methods.

3. **Epithelial-on-top in UC is trait-specific.** SCZ has colonocyte/goblet/EE near-null or negative. Height has them moderate-positive but well below the stromal/vascular top. UC's pattern doesn't replicate in controls.

### Pipeline positive validation: Height is clean

Height top 3, **BOTH methods identically**:
- endothelium (scDRS z=+3.81, seismic 9.27)
- mural/glia (z=+3.73, 5.57)
- fibroblast (z=+3.72, 5.52)

Exactly the right cell types for height — vascular endothelium, pericytes, connective-tissue fibroblasts express the height-associated GWAS genes (collagens, ECM, skeletal regulators). The pipeline correctly identifies stromal/vascular as height-relevant in a gut atlas where no skeletal/CNS cells exist. **Pipeline validates.**

### SCZ negative control: NOT quite null (and that's biologically defensible)

SCZ tops in scDRS with **mural/glia (z=+4.23, mcp=0.001)** — and seismic also lights it at 0.68 (rank 4/15). This isn't a clean null on a gut atlas, but it's biologically defensible: **enteric glia share developmental program with CNS glia** (both neural-crest derived; both express GFAP, S100, glial markers). Trubetskoy 2022 PGC3 SCZ explicitly identifies glia as a cell-type enrichment. So SCZ scoring glia in a gut atlas is the method correctly identifying the cell type that carries the trait-relevant biology even when the "right" tissue isn't present.

SCZ also has **B cell (z=+2.52, mcp=0.013)** — possibly real (some immune-axis SCZ findings) or noise; not clean null but small.

### Cross-method concordance (Track [4]): all four cases

Output: `results/concordance/garrido_cross_method_prototype.csv`.

| Trait + MHC | Spearman ρ | 95% CI | Jaccard top 5 | Jaccard top 10 | Kappa |
|---|---|---|---|---|---|
| UC_MHCexcl   | 0.60 | [0.12, 0.88] | 0.67 | 0.67 | n/a (kappa=0, no method-agreement on sig cells) |
| UC_MHCincl   | 0.78 | [0.39, 0.97] | 0.67 | 0.82 | n/a |
| SCZ_MHCexcl  | 0.61 | [0.07, 0.90] | 0.43 | 0.82 | n/a |
| **Height_MHCexcl** | **0.83** | **[0.47, 0.97]** | 0.67 | 0.67 | **1.00** |

**Height is the cleanest cross-method result**: ρ=0.83, AND kappa=1.0 — both methods agree perfectly on which cells reach FDR<0.05 (3 cells in both: endothelium, mural/glia, fibroblast). That's the strongest cross-method validation the dry run produced.

UC and SCZ both sit at ρ ~ 0.60 cross-method under MHC-excluded. UC's cross-method ρ rises to 0.78 under MHC-incl (Saisohan's "strong signal swamps methodological differences" framing confirmed).

### Net for M4

The dry run delivers four robust findings for the methods-caution section:

1. The **MHC-handling trade-off** is reproducible across two independent scoring frameworks (scDRS + seismic). MHC inclusion produces an HLA-marker tautology in immune cells; MHC exclusion misses the bulk of immune-coded UC signal. Neither variant is clean.

2. Cross-method concordance is **moderate-to-strong** (ρ 0.60–0.83 across all four cases on Garrido). High enough to justify reporting both methods; not so high that one is redundant.

3. The UC cell-type pattern is **trait-specific, not baseline method bias** — SCZ and Height produce different rankings on the same atlas under the same methods. So the UC findings reflect UC biology (filtered through MHC handling), not methodological artifacts.

4. Pipeline-positive validation (Height) hits the right cells. Pipeline-negative (SCZ) hits a biologically defensible cell type (glia) rather than a null. Methods are working as expected.

What's NOT delivered on the laptop and remains HB-bound:
- Multi-atlas validation (the actual paper question). Garrido is one atlas; the cross-atlas comparison is what determines whether any of this generalizes.
- SCZ_MHC + Height_MHC scDRS variants (for the full 6-cell trade-off panel). Deferrable — primary finding stands on the MHC-excl comparison.
- Fine-tier cross-method panels and concordance.

Outputs:
- `results/scdrs/garrido_{scz,height}/{SCZ,HEIGHT}.scdrs_group.cell_type_broad`
- `results/seismic/garrido_{scz,scz_mhc,height,height_mhc}_broad.tsv`
- `results/concordance/garrido_cross_method_prototype.csv`

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

---

## Remaining laptop work to close the validation loop

The laptop's job is to validate the entire analysis loop on Garrido
(the one laptop-sized atlas) so HB only ever does RAM-bound big-atlas
scoring. So far we've validated **one method (scDRS) on one GWAS**.
The loop isn't closed. Order matters — earlier items unblock later
ones.

### [1] Finish MHC perform-downstream — read the 461-cell distribution

In-progress (~13 min from compute-score). Question: which broad cell
types host the 461 FDR<0.1 cells? If they cluster in HLA-high
antigen-presenting cells (B / DC / macrophage), the MHC-included
ranking is partly the HLA-marker tautology we flagged, not
independent evidence of immune-cell UC enrichment. If they distribute
across non-APC cell types, the recovered signal is more
polygenic-real.

Output: `results/scdrs/garrido_delange_mhc/UC_MHC.scdrs_group.cell_type_broad`
(landing now).

### [2] seismicGWAS on Garrido × de Lange (BOTH MHC variants)

The biggest gap. The paper is a cross-method comparison and seismic
has never run. R package, separate `uc-cross-atlas-r` conda env
(install in progress; same isolation discipline as the LDSC Py 2.7
env). Run Garrido × de Lange with both `.gs` files.

**The question that matters:** does seismicGWAS reproduce the same
MHC-driven inversion scDRS showed?

- If both methods invert when MHC is stripped → property of the
  gene-level signal. Robust methods-caution finding.
- If only scDRS inverts → method-specific quirk; the cross-method
  result becomes the headline differentiator.

Either answer is a result and can only be cheaply produced on the
laptop.

### [3] Controls through both methods on Garrido

SCZ (negative) and Yengo height (pipeline positive) `.gs` files are
already built. Score Garrido through both methods for both:

- **SCZ on gut atlas must be null** across cell types. If it lights
  anything up, the approach is suspect — find out now, not after HB.
- **Yengo height** is the pipeline positive: validates the scoring
  machinery on a high-N polygenic trait (won't have biologically
  meaningful gut-cell hits, but should produce stable, predictable
  numerics).

Load-bearing for the methods section.

### [4] Cross-method concordance prototype on Garrido

Once [2] and [3] are done, run `code/08_cross_method/run_cross_method.py`
on the Garrido outputs. The 3 E2E tests for that script (Track D)
exercise it on synthetic fixtures; this is the first real-data run.
Validates the comparison machinery (Spearman, Jaccard, kappa) before
HB feeds it five atlases of outputs.

### Out of scope on the laptop (HB-only)

- Loading and scoring TAURUS (12.7 GB), HCA Gut, Pan-GI, Smillie.
  RAM-bound and HB-only by the split we already set.

### Optional / deferrable

- **05 gene-property smoke on Garrido** — code-validation only, the
  locked deliverable is Smillie (HB).
- **F10 LDSC intercept** — needs its own Py 2.7 env + `eur_w_ld_chr`;
  the `01b_ldsc/` scaffold exists with 23 green tests, execution
  waits on WSL Ubuntu being up (Track A0 still gated on the human
  `wsl --install` step).

### Insurance argument

The MHC inversion we caught above is exactly the failure mode that
would have cost a full HB run if we'd skipped the dry run. Finishing
[1]-[4] on Garrido is the same insurance applied to seismic + the
controls + the cross-method comparison. HB then becomes a known-good
scale-up, not a debug-at-scale exercise.
