# DECISIONS.md

Append-only log of analytical commitments. Never edit existing entries —
only append, with date and rationale, when a commitment changes. The
canonical plan is `docs/uc-cross-atlas-v1-plan.pdf`.

---

## YYYY-MM-DD (M0): Locked v1 scope

> Replace `YYYY-MM-DD` with the actual M0 lock-in date when this commit is
> finalized at the project kickoff. All `[X]` placeholders must be filled
> in at M0 (Hummingbird allocation, Cell Ontology release date, etc.).

**Project:** Cross-atlas reproducibility of GWAS-driven cell-type
prioritization in UC, comparing regime 1 (single-atlas) and regime 2
(per-atlas-then-meta-analyze via Brown's method).

**Disease:** Ulcerative colitis only. Cross-disease + cross-tissue
(UC + CD + RA) generalization is v2's job.

### GWAS (locked core)

- **Primary:** de Lange et al. 2017 UC, GWAS Catalog GCST004131.
- **Cross-GWAS sensitivity:** Liu et al. 2023, UC arm only. Per-SNP N
  column verified in M1 (or fixed-N approximation documented if absent).
- **Negative control:** Trubetskoy 2022 schizophrenia GWAS — run on Smillie
  at broad tier under the same MAGMA pipeline. No colon cell type should
  achieve FDR < 0.05.
- **LD reference:** 1000G EUR for both UC GWAS (acknowledged in Methods as
  approximate for multi-ancestry Liu 2023).
- **MHC region** (chr 6: 28,477,797–33,448,354 GRCh37) **excluded** from
  the scDRS top-1000 gene set and the seismicGWAS gene-Z-score table.
  Sensitivity analysis with MHC retained: one supplementary scDRS run on
  Smillie × de Lange.
- **Autosomes only** (chr 1–22). X chromosome excluded.
- **λ_GC verified ≤ 1.1** for both UC GWAS in M1 (or flagged for revision
  response if higher).
- **Regime 2 on de Lange only.**

### Atlases (locked core)

- **UC trio:** Smillie 2019 (SCP259), Kong 2023 UC (GSE214695), Mennillo
  2024 (GEO accession to confirm in M1; pre-treatment baseline samples
  only — verify ≥8 donors after subsetting, else fall back to
  Garrido-Trigo 2023).
- **Independent broad comparator:** HCA Gut Cell Atlas / Elmentaite 2021.
  M1 must verify zero donor overlap with the UC trio.
- **Integration-pipeline-robustness comparator:** Pan-GI / Oliver 2024
  (with/without donor-overlap analysis). Pan-GI integrated Smillie 2019 and
  anchored Kong 2023; 2/3 of trio donors are inside Pan-GI.
- **Fallbacks:** Garrido-Trigo 2023, Boland 2020, Devlin/Zhao 2023.
- **Atlas preprocessing:** take published preprocessing as-is (Option B).
  Doublet/ambient/QC differences across atlases are disclosed in Methods
  Limitations. Marker-gene QC at end of M2 catches the worst label-mismatch
  problems.
- **Inflamed/uninflamed sample composition** documented per atlas in donor
  metadata; reported in Methods Table 1.

### Methods

- **scDRS** and **seismicGWAS** in locked core.
- **scPagwas NOT in v1** (compute infeasible).

### Regimes

- **Regime 1** (single-atlas) + **Regime 2** (Brown's method analytical
  Kost-McDermott via `EmpiricalBrownsMethod::kostsMethod()` with empirical
  cross-atlas correlation matrix from null draws).
- Regime 3 (scANVI integration) deferred to v2.

### Granularity

- Broad (~10–15 cell types) AND fine (~30–50 cell states), **both primary**.

### Cell-type harmonization

- **Vocabulary:** gut mucosa.
- **Strategy:** scOntoMatch + manual CSV crosswalk + marker-gene
  confirmation, applied across all 5 atlases.
- **Sub-clustering depth at fine tier:** roll up to lowest-resolution atlas
  in any pairwise comparison. Higher-resolution sub-states retained in
  supplementary tables, not used in headline metrics.
- **Marker-gene QC required at end of M2.** For each broad cell type,
  compute mean expression of canonical markers per atlas. Flag any cell
  type where one atlas's mean marker expression deviates by > 2 SD from
  the others; investigate and re-map before locking the crosswalk.
- **Cell Ontology release pinned** in this file at M2 (record OWL file
  release date downloaded for harmonization work). Currently: `[X]`.
- **Hard deadlines:** broad-tier end of M2 (including marker-gene QC);
  fine-tier mid-M3.

### scDRS configuration

- **Covariates** in `compute-score`: `log_n_genes`, `log_n_counts`, donor
  (one-hot), sample (one-hot if multiple per donor), sex (one-hot if
  donor-level metadata exists). Sex omission documented per atlas if
  metadata is missing.
- **All-cells policy:** `compute-score` and group analysis run on all cells
  in each atlas regardless of disease status.
- **Top-1000 MAGMA gene cutoff** for `.gs` file (scDRS default),
  MHC-excluded. Sensitivity analysis with MHC retained: one supplementary
  scDRS run on Smillie × de Lange.
- **Headline random seed: `seed = 42`** for Monte Carlo null sampling and
  bootstrap resampling. Test-retest uses seeds 1, 2, 3.
- **Min cell-count threshold:** cell types with < 50 cells in any atlas in
  a comparison are excluded from concordance metrics for that comparison.

### Statistical metrics

- **Cross-atlas headline metric: Spearman ρ on cell-type-level Z-scores
  (scDRS) and regression coefficients (seismicGWAS) — NOT p-values.**
  scDRS p-values are computed against an atlas-specific Monte Carlo null
  distribution; same biology in two atlases produces different p-values
  purely from compositional differences.
- Cohen's κ on FDR-significance with marginals; top-5/10 Jaccard at broad,
  top-5/10/20 Jaccard at fine.
- **Tied-rank handling for Spearman ρ: average-rank tie-breaking
  (`scipy.stats.spearmanr` default; R `cor(method="spearman")` default).**
- **Concordance computed on shared cell-type intersection per pair, with
  minimum cell-count threshold (≥50 cells in both atlases per cell type).**
  Atlas-specific cell types reported separately, not entered into
  concordance metrics in either direction.
- **"Detected" operationally defined:** harmonized label present in atlas
  crosswalk AND atlas has ≥50 cells assigned to that label.
- **Bootstrap 95% CIs on every reported Spearman ρ:** 1000 iterations,
  resampling over cell types within shared intersection, percentile method,
  seed = 42. BCa deferred to revision if reviewers request.
- Five concordance axes: cross-atlas, cross-method, cross-GWAS,
  regime-1-vs-2, broad-comparator (HCA Gut + Pan-GI).
- **Headline 3×3 heatmap is summarizing three pairwise ρ values.** Pattern
  claims across these three values are qualitative descriptions of n = 3
  observations, not statistical claims. Methods text says so explicitly.
- **Donor-LOO uncertainty intervals** on broad-tier headline metrics under
  de Lange for both methods. Reported as "95% LOO jackknife range," not
  "95% CI." Liu sensitivity reported as point estimates only. Kong (12
  donors, 6 UC + 6 healthy) enforces a minimum-after-LOO threshold of
  ≥5 cases per group; LOO iterations violating this are dropped.

### Multiple-testing strategy

- **Single primary analysis:** 3×3 atlas-pair Spearman ρ at broad tier
  under de Lange via scDRS — reported uncorrected.
- **All other secondary comparisons (~60 across the four other axes):**
  Benjamini-Hochberg FDR < 0.05 across the full battery for any claim of
  statistical significance. Uncorrected p-values reported in supplementary
  tables for transparency.

### Cohen's κ marginal-imbalance contingency

If ≥80% of cell types pass FDR < 0.05 in both atlases, **report κ at
FDR < 0.01** as the headline κ instead. Report the FDR < 0.05 κ and
marginals in supplementary tables.

### Pan-GI donor-overlap policy

Report Pan-GI concordance with each UC atlas in two flavors:
- **With donor overlap retained:** sanity check (mostly the same data).
- **With donor-overlap exclusion** (subset Pan-GI to donors not from
  Smillie 2019 or Kong 2023, using donor-attribution metadata from M1):
  **this is the headline Pan-GI panel.**

Donor-overlap-excluded Pan-GI is still an scVI-integrated multi-study
atlas — just without our trio's donors. It is **not** "Pan-GI without
integration." Methods text states this explicitly.

### Pan-GI and HCA Gut counts policy

Use the **raw count layer** (`adata.layers['counts']` from cellxgene
Discover, or `adata.raw.X` if available) with **published cell-type
labels** — no re-clustering. Verify accessibility during M1 atlas
standardization; if only normalized/integrated counts are accessible, drop
the comparator to stretch.

### Brown's method (regime 2) implementation

- **Cross-atlas correlation matrix estimated empirically** from
  per-null-draw cell-type-level statistics (scDRS Monte Carlo null draws +
  seismicGWAS gene-Z permutations, both serialized in M3).
- **Per-cell-type correlation matrix:** Pearson correlation across N = 1000
  null draws between atlas pairs; the (A, B) entry for cell type c is the
  Pearson correlation of N null statistic pairs (T_c^A_i, T_c^B_i).
  Kost-McDermott combination performed per cell type using that
  cell-type-specific correlation matrix.
- **Sanity check:** off-diagonal entries should be positive (atlases share
  GWAS inputs and tissue biology). If any cell type has near-zero or
  negative cross-atlas null correlation, investigate before combining
  (likely harmonization mismatch or insufficient sample size).
- **Edge-case fallback:** for any cell type where null-statistic SD in any
  atlas falls below the 5th percentile of null-statistic SDs across all
  cell types in that (method, tier, GWAS) combination, replace the
  cell-type-specific correlation matrix with the **median cross-atlas
  correlation across well-behaved cell types** (those above the
  5th-percentile threshold) in that combination. Affected cell types
  flagged in supplementary tables.
- **Heterogeneous per-atlas sample sizes** (Smillie 30 / Kong 12 / Mennillo
  ~10–15 donors): use **unweighted Brown's** via `kostsMethod()`. Document
  equal-weighting as defensible because per-atlas N ratios are < 3×.
  Stouffer-weighted Brown's deferred to revision (~1 day of custom code) if
  reviewers push back.
- **Stretch #1 (permutation Brown's) is a conditional fallback** activated
  only if the empirical correlation matrix shows pathological entries
  (negative or near-zero off-diagonals where shared MAGMA inputs and
  shared etiology guarantee positive).
- **Cell types missing from some atlases (fine tier):**
  - 3/3 — full 3×3 empirical correlation matrix.
  - 2/3 — combine with `n_atlases = 2` flag, 2×2 correlation submatrix.
  - 1/3 — do not combine; report regime-1 with `n_atlases = 1`. Excluded
    from regime-2 ranking.

### seismicGWAS configuration

- **Confounders in `get_ct_trait_associations()`:** gene length
  (log-transformed), gene-gene LD score, transcript count. Verified in M1
  by inspecting package source; if package defaults differ, override and
  apply this explicit set.
- **Specificity score serialization (v2-setup decision).** Saved as
  long-format `.rds` (or feather) per (atlas, granularity, cell_type, gene)
  with columns `atlas`, `cell_type`, `granularity`, `gene`,
  `specificity_score`. Stored alongside regression results, not embedded.
- **Permuted-Z-score null draws for Brown's method:** M = 1000 permutations
  of the gene-Z-score vector per (atlas, GWAS, tier); per-permutation
  cell-type-level test statistics saved to feather.
- **Test-retest gate: ρ ≥ 0.999** (deterministic method given fixed
  cell-type labels and expression matrix; deviation indicates a bug).

### Sanity scaffolding (locked core, do not drop)

- Test-retest baseline: 3 seeds × 3 atlases × 2 methods, de Lange only.
  Pass: scDRS ρ ≥ 0.9; seismicGWAS ρ ≥ 0.999.
- scDRS **positive** control on Tabula Muris × Yengo 2022 height GWAS.
- **scDRS negative control on Smillie × Trubetskoy 2022 schizophrenia GWAS
  at broad tier.** Caveat: schizophrenia hits include MHC and
  complement-pathway genes; "absence of enrichment" rather than "exact
  null" is the operational expectation. MHC excluded per the MHC policy
  above.
- MAGMA gene-property sanity track on Smillie at broad tier under de Lange.

### Compute budget

- Locked v1: ~112 scDRS + ~109 seismicGWAS = ~69–142 node-hours.
- With all stretches: ~442 scDRS + ~109 seismicGWAS = ~234–472 node-hours.
- Hummingbird allocation as of `[date]`: `[X] node-hours`.
- Pan-GI may require ~30 GB memory.
- **Descope order:** stretches → seismic donor-LOO → reduce scDRS LOO.
  Do not drop scDRS broad-tier donor-LOO under de Lange.

### Stretch ladder (priority order)

1. Brown's-method permutation-based covariance (conditional fallback if
   empirical Kost-McDermott matrix is pathological).
2. scANVI integration as regime 3 (deferred to v2 by default).
3. LDL × Tabula Sapiens hepatocytes complementary positive control.

### Gating decisions

- Gate 1 (M3): regime 1 working for both methods, both GWAS, both tiers.
- Gate 2 (M5): locked core clean.
- Gate 3 (M7): writing on track.

### Authorship and division of labor

- Documented in README. First author = primary manuscript drafter; co-first
  if contributions are genuinely indistinguishable. Discuss with PI mentor
  before work starts. Revisit at M4.

---

## v2 trajectory architectural decisions (commit at M0 if PI is committed to v2)

**v2 scope:** UC + CD + RA cross-tissue generalization with regime 3 (scANVI
integration) as primary methodological contribution. 18–24 month follow-up
project after v1 ships.

**v1 architectural decisions to enable v2 cheaply:**

1. **CL-ontology-aware harmonization (M2, ~1 day).**
   - Crosswalk CSV includes `cl_term` and `cl_parent_chain` columns for
     every immune cell type.
   - Immune-cell CL assignments use parent terms that admit both gut and
     synovium daughters (e.g., gut "T cell" CL:0000084 →
     CL:0000542 "lymphocyte").
   - Walk parent chain until tissue-agnostic; document chain in CSV.
   - Time-boxed to 1 day during M2 harmonization. If running over, take
     crosswalk as-is; v2 retrofits later. v1 figures unaffected.

2. **seismicGWAS specificity score serialization (M3, ~half day).**
   - Save specificity scores per (atlas, granularity, cell_type, gene) as
     long-format `.rds` (or feather), in `code/04_seismic/specificity_long/`.
   - Schema: `atlas`, `cell_type`, `granularity`, `gene`, `specificity_score`.
   - Stored alongside regression results, not embedded in them.
   - Required because v2 cross-tissue claim compares UC-T-cell specificity
     to RA-T-cell specificity for the same gene set. Embedding inside
     regression-result objects forces full re-run later.

3. **Generalized donor-attribution metadata (M1, near-zero cost beyond Pan-GI).**
   - Per-atlas `data/atlases/donor_metadata/<atlas>_donor_metadata.csv`
     with columns: `donor_id`, `originating_study`, `tissue`,
     `disease_status`, `assay_protocol` (+ optional `tissue_state`
     ∈ {inflamed, uninflamed, unspecified}).
   - Pan-GI is the immediate v1 consumer (donor-overlap exclusion at M5).
   - Same schema scales to AMP-RA Phase 2 and any future integrated atlas
     in v2.

**These three decisions are droppable** under v1 schedule pressure without
compromising v1. Drop priority: CL parent chain first (least painful
retrofit) → generalized donor-attribution second (Pan-GI metadata still
required regardless) → specificity score serialization last (most painful
retrofit).

**v2 decisions deferred to v2 planning** (typically M9–M12 of v1 when
bandwidth permits):

- Specific RA atlases (likely AMP-RA Phase 1 + Phase 2; verify access at v2 M0).
- Specific RA GWAS (Okada 2014 vs Ishigaki 2022 multi-ancestry).
- CD atlas additions (Kong 2023 CD arm; Martin 2019; possibly Elmentaite 2020 CD subset).
- Three-tier harmonization architecture full design (cross-tissue /
  within-tissue broad / within-tissue fine).
- Regime 3 implementation choices (scANVI parameter tuning, integration QC
  metric thresholds).
- Whether to add a third author for v2.

**PI commitment status as of M0:** `[explicitly recorded — yes / no / conditional]`.
If yes, the three v1 architectural decisions above are committed. If no,
they are dropped and v1 stands alone.

---

## [Append future entries here, dated, with rationale]
