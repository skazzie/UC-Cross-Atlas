# UC Cross-Atlas Reproducibility v1 — Operational Summary

A 12-month, two-author dry-lab project at UCSC. This file is a tight
operational summary; the canonical plan is the PDF in `docs/`:

- [`docs/uc-cross-atlas-v1-plan.pdf`](docs/uc-cross-atlas-v1-plan.pdf)
- [`docs/uc-cross-atlas-pre-project-curriculum.pdf`](docs/uc-cross-atlas-pre-project-curriculum.pdf)

All locked methodological commitments live in [`DECISIONS.md`](DECISIONS.md)
(append-only). When this summary and the PDF disagree, the PDF wins.

---

## What we are doing

Comparing two operational regimes for GWAS-driven cell-type prioritization
on ulcerative colitis colon mucosa, across three independent atlases, two
GWAS, two methods, and two cell-type granularity tiers.

- **Regime 1 (single-atlas):** run scDRS or seismicGWAS on whichever atlas
  you have.
- **Regime 2 (per-atlas-then-meta-analyze):** run each method per atlas,
  then combine cell-type p-values across atlases via Brown's method
  (correlation-aware Fisher's combination).

We additionally compare to two broad multi-tissue gut atlases — HCA Gut
(independent of our trio) and Pan-GI (which integrated Smillie 2019 and
anchored Kong 2023, providing an integration-pipeline-robustness check
rather than independent replication).

## What this paper does *not* claim

- Not "regimes generalize across diseases" — UC-only.
- Not "regime 3 (integration) is gold standard" — deferred to v2.
- Not "concordance measures whether the same biology is identified" —
  scDRS group-level p-values use atlas-specific reference distributions, so
  some apparent disagreement is structural.
- Not "scDRS and seismicGWAS are statistically independent" — both consume
  the same MAGMA gene-Z. Cross-method *disagreement* is the informative
  signal.
- Not "Pan-GI is independent broad-atlas replication" — 2/3 of our trio's
  donors are inside Pan-GI. HCA Gut is the only fully-independent broad
  comparator (subject to M1 donor-overlap verification).
- Not "Liu 2023 cross-GWAS robustness implies cross-disease generalization."
  Liu 2023 is the same disease in a partially overlapping cohort.

## Locked v1 scope

- **Disease:** ulcerative colitis only.
- **GWAS:** de Lange 2017 UC (primary, GCST004131) + Liu 2023 multi-ancestry
  UC arm (cross-GWAS sensitivity). Trubetskoy 2022 schizophrenia GWAS as
  negative control on Smillie at broad tier.
- **Atlases:** Smillie 2019, Kong 2023 UC, Mennillo 2024 (pre-treatment-only
  baseline; verify ≥8 donors at M1, else fall back to Garrido-Trigo 2023).
  Plus HCA Gut (independent comparator) and Pan-GI (integration-pipeline
  comparator with donor-overlap dual analysis).
- **Methods:** scDRS + seismicGWAS. *scPagwas is NOT in v1.*
- **Regimes:** Regime 1 (single-atlas) + Regime 2 (Brown's analytical
  Kost-McDermott with empirical correlation matrix from null draws).
  Regime 3 (scANVI integration) is deferred to v2.
- **Granularity:** broad (~10–15 cell types) AND fine (~30–50 cell states),
  both primary.
- **MHC region** (chr 6: 28,477,797–33,448,354 GRCh37) **excluded** from the
  scDRS top-1000 gene set and the seismicGWAS gene-Z-score table.
  Sensitivity analysis with MHC retained: one supplementary scDRS run on
  Smillie × de Lange.
- **Autosomes only** (chr 1–22). X chromosome excluded.
- **LD reference:** 1000G EUR for both UC GWAS (acknowledged approximation
  for multi-ancestry Liu 2023).
- λ_GC reported in DECISIONS.md after M1; if > 1.1, flag for revision
  response.

## Five concordance axes

1. Cross-atlas (within method, GWAS, tier) — 24 ranking comparisons.
2. Cross-method (within atlas, GWAS, tier) — 12 comparisons.
3. Cross-GWAS (within atlas, method, tier) — 12 comparisons.
4. Regime 1 vs regime 2 (within method, GWAS, tier) — 4 comparisons (regime
   2 on de Lange only).
5. Tissue-matched vs broad atlas — HCA Gut (independent) and Pan-GI
   (donor-overlap dual analysis).

## Three statistical metrics

- **Headline metric:** Spearman ρ on cell-type-level Z-scores (scDRS) /
  regression coefficients (seismicGWAS) — *not* p-values. Computed via
  scipy/R defaults (average-rank tie-breaking) on the shared cell-type
  intersection between the two atlases. Each cell type in the intersection
  must have ≥50 cells in both atlases.
- **Top-k Jaccard:** k = 5, 10 at broad tier; k = 5, 10, 20 at fine tier.
- **Cohen's κ on FDR-significance:** within-atlas BH-FDR. Marginal-saturation
  contingency: report κ at FDR < 0.01 as headline if ≥80% of cell types
  pass FDR < 0.05 in both atlases.

**Bootstrap 95% CIs** on every reported Spearman ρ — 1000 iterations,
percentile method, resampling over cell types within shared intersection,
seed = 42.

**Multiple-testing correction:** the single primary analysis is the 3×3
atlas-pair Spearman ρ heatmap at broad tier under de Lange via scDRS,
reported uncorrected. All ~60 secondary comparisons receive
Benjamini-Hochberg FDR < 0.05 across the full battery.

**Headline 3×3 heatmap caveat:** the matrix has only three lower-triangle
entries (Smillie–Kong, Smillie–Mennillo, Kong–Mennillo). Pattern claims
across these three values are qualitative descriptions of n = 3
observations, not statistical claims. Methods text says so explicitly.

## Brown's method (regime 2)

Analytical Kost-McDermott via `EmpiricalBrownsMethod::kostsMethod()`. The
cross-atlas correlation matrix is **estimated empirically** from
per-null-draw cell-type-level statistics (scDRS Monte Carlo null draws +
seismicGWAS gene-Z permutations, both serialized at M3). Per-cell-type
correlation matrix is the Pearson correlation across N = 1000 null draws
between atlas pairs. Off-diagonal entries should be positive; investigate
any cell type with near-zero or negative cross-atlas null correlation
before combining.

**Edge-case fallback for low-variance null statistics:** if a cell type's
null-statistic SD in any atlas falls below the 5th percentile of
null-statistic SDs across all cell types in that (method, tier, GWAS)
combination, replace its correlation matrix with the median cross-atlas
correlation across well-behaved cell types in that combination. Affected
cell types flagged in supplementary tables.

**Heterogeneous per-atlas N:** Smillie 30, Kong 12, Mennillo ~10–15 donors.
Use unweighted Brown's via `kostsMethod()` and document the
equal-weighting as defensible because per-atlas N ratios are < 3×.
Stouffer-weighted Brown's deferred to revision if reviewers push back.

**Cell types missing from some atlases (fine tier):**
- 3/3 atlases — full 3×3 empirical correlation matrix.
- 2/3 — combine with `n_atlases = 2` flag, 2×2 correlation submatrix.
- 1/3 — do not combine; report regime-1 with `n_atlases = 1`. Excluded from
  regime-2 ranking.

## Sanity scaffolding (locked, do not drop)

- Test-retest baselines: 3 seeds × 3 atlases × 2 methods, de Lange only.
  Pass gates: scDRS ρ ≥ 0.9; seismicGWAS ρ ≥ 0.999 (deterministic).
- scDRS positive control: Tabula Muris × Yengo 2022 height GWAS.
- scDRS negative control: Trubetskoy 2022 schizophrenia × Smillie at broad
  tier — no colon cell type should achieve FDR < 0.05.
- MAGMA gene-property sanity track on Smillie at broad tier under de Lange.
- Donor-LOO uncertainty intervals on broad-tier headline metrics under de
  Lange for both methods (~180 LOO runs). Liu = point estimates only. For
  Kong (12 donors, 6 UC + 6 healthy), enforce ≥5 cases per group after LOO.

## scDRS configuration

- **Covariates** (`<atlas>_covariates.tsv`): `log_n_genes`, `log_n_counts`,
  donor (one-hot), sample (one-hot if multiple per donor), sex (one-hot if
  metadata exists; document omission per atlas in DECISIONS.md).
- **All-cells policy:** `compute-score` and group analysis run on all cells
  in each atlas regardless of disease status.
- **Top-1000 MAGMA gene cutoff** for `.gs` (scDRS default), MHC-excluded.
- **Headline random seed:** `seed = 42` for Monte Carlo null sampling and
  bootstrap resampling. Test-retest uses seeds 1, 2, 3.
- **Min cell-count threshold:** cell types with < 50 cells in any atlas in
  a comparison are excluded from concordance metrics for that comparison
  (reported in supplementary tables).

## seismicGWAS configuration

- Confounders in `get_ct_trait_associations()`: gene length (log), gene-gene
  LD score, transcript count (package defaults; verify in M1 and override
  if defaults differ).
- Specificity scores serialized as long-format `.rds`/feather per (atlas,
  granularity, cell_type, gene), stored alongside regression results — not
  embedded in them. v2 cross-tissue claim depends on this.
- Gene-Z permutation null draws: M = 1000 per (atlas, GWAS, tier),
  serialized for Brown's empirical correlation matrix.

## v2 trajectory architectural decisions (M0)

If the PI commits to v2 (UC + CD + RA cross-tissue with regime 3, 2027–2028),
v1 carries three lightweight setup decisions:

1. **CL-ontology-aware harmonization (M2, ~1 day).** Crosswalk includes
   `cl_term` and `cl_parent_chain` columns. Immune-cell CL assignments use
   parent terms admitting both gut and synovium daughters.
2. **seismicGWAS specificity score serialization (M3, ~half day).**
   Long-format files stored alongside regression results.
3. **Generalized donor-attribution metadata (M1, near-zero cost beyond
   Pan-GI).** `donor_metadata.csv` per atlas with columns `donor_id`,
   `originating_study`, `tissue`, `disease_status`, `assay_protocol`.

If PI does *not* commit to v2 at M0, drop all three. v1 stands alone.

## Stretch ladder (compact)

1. **Brown's-method permutation-based covariance** — conditional fallback
   if empirical Kost-McDermott matrix is pathological (negative or
   near-zero off-diagonals). ~300 additional scDRS runs.
2. **scANVI integration as regime 3** — deferred to v2 by default.
3. **LDL × Tabula Sapiens hepatocytes** complementary positive control.
   Closes the human-only-code-path gap. ~1 day.

Stretch N+1 only starts when stretch N is fully clean and written up.

## Compute budget

- Locked v1: ~112 scDRS runs + ~109 seismicGWAS runs ≈ **69–142 node-hours**.
- With all stretches: ≈ **234–472 node-hours**.
- Pan-GI (~1.1M cells) may need ~30 GB memory and a higher-memory node.
- **Descope order if compute is short:** stretches → drop seismic donor-LOO
  → reduce scDRS donor-LOO. *Do not* drop scDRS broad-tier donor-LOO under
  de Lange.

## Gating decisions

- **Gate 1 (M3):** regime 1 working for both methods, both GWAS, both tiers.
- **Gate 2 (M5):** locked core clean.
- **Gate 3 (M7):** writing on track.

## Three explicit fork points for descope decisions

- **M3 fork** — drop fine-tier from primary, drop seismicGWAS, drop Liu, or
  drop the v2-setup decisions if any are blocking.
- **M5 fork** — drop Pan-GI to stretch if donor-overlap analysis is
  unclear; defer regime-2 cleanup to M6.
- **M7 fork** — drop all stretches; M7–M11 is pure writing.

**Do NOT drop:** test-retest baselines, scDRS positive control on Tabula
Muris, donor-LOO on broad-tier scDRS headline under de Lange,
DECISIONS.md updates.

## Deliverable

bioRxiv preprint targeting *Bioinformatics Advances* (primary) or *NAR
Genomics and Bioinformatics* (secondary). Submit M11–M12. First decision
M14–M15. Published M17–M19.

## Repository structure

```
uc-cross-atlas-v1/
├── DECISIONS.md                # locked methodological commitments (append-only)
├── PLAN.md                     # this file
├── README.md
├── docs/
│   ├── uc-cross-atlas-v1-plan.pdf
│   └── uc-cross-atlas-pre-project-curriculum.pdf
├── data/
│   ├── gwas/                   # de Lange + Liu 2023 + Trubetskoy
│   ├── atlases/                # Smillie, Kong-UC, Mennillo, HCA Gut, Pan-GI
│   │   └── donor_metadata/     # per-atlas donor_metadata.csv (v2-setup pattern)
│   └── reference/              # 1000G EUR, CL ontology, HGNC pin
├── code/
│   ├── 01_magma/               # ×2 GWAS + negative control
│   ├── 02_atlas_prep/          # download, harmonize, HGNC remap, donor-attribution
│   ├── 03_scdrs/               # regime 1 (×2 GWAS, ×2 tiers)
│   ├── 04_seismic/             # regime 1 (×2 GWAS, ×2 tiers)
│   │   └── specificity_long/   # serialized specificity for v2 reuse
│   ├── 05_magma_geneprop/      # Watanabe-style sanity track
│   ├── 06_concordance/         # 5 axes
│   ├── 07_regime2_meta/        # Brown's method (analytical, de Lange only)
│   ├── 08_cross_method/        # scDRS vs seismicGWAS within-atlas
│   ├── 09_cross_gwas/          # de Lange vs Liu within-atlas
│   ├── 10_broad_atlas_hca/     # HCA Gut comparator (independent)
│   ├── 11_broad_atlas_pangi/   # Pan-GI comparator (with/without donor-overlap)
│   ├── 12_positive_control/    # Tabula Muris × height
│   └── 13_figures/
├── results/
├── notebooks/
└── manuscript/
```
