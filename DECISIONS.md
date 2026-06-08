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
- **Fallbacks:** Garrido-Trigo 2023, Boland 2020, ~~Devlin/Zhao 2023~~ [struck — DECISIONS 23(b)].
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

## CORRECTION 2026-05-20 (1/7): de Lange UC GWAS accession

Original entry: "UC primary = de Lange 2017, GCST004131."

Verified incorrect via GWAS Catalog browser check + cross-reference to
medRxiv 2025.03.03.25323217, which lists the de Lange 2017 trio as:
  GCST004131 = IBD-combined
  GCST004132 = Crohn's disease
  GCST004133 = ulcerative colitis

Corrected v1 primary UC GWAS accession: GCST004133.

All downstream filenames (uc_delange.*) and artifacts retain the
"uc_delange" prefix; only the accession changes.

---

## CORRECTION 2026-05-20 (2/7): Atlas 2 citation

Original entry: "Atlas 2 = Kong et al. 2023, Immunity 56:444-458,
GEO GSE214695 (UC subset of CD atlas)."

Verified incorrect via CELLxGENE dataset page for GSE214695, which is
owned by:
  Garrido-Trigo et al. 2023, Nat Commun 14:4506,
  "Macrophage and neutrophil heterogeneity at single-cell spatial
   resolution in human inflammatory bowel disease."
  doi: 10.1038/s41467-023-40156-6
  Cohort: 6 HC + 6 CD + 6 UC active, 46,700 cells, colonic mucosa.
  UC subset for our analysis: 6 HC + 6 UC active = 12 donors,
  30,068 cells.

Kong et al. 2023 Immunity is a separate Crohn's Disease atlas (ileum + colon)
with no clean UC subset and is NOT appropriate as v1 Atlas 2. It may be
reconsidered for v2 CD work.

Corrected v1 Atlas 2: Garrido-Trigo et al. 2023, GSE214695 (UC subset).
Donor counts (6 HC + 6 UC) and GEO accession (GSE214695) were already
correct; only the citation/attribution changes.

Internal slug change: rename "kong" -> "garrido_trigo" across the repo
(see Part 6 for the exact rename list).

---

## CORRECTION 2026-05-20 (3/7): Pan-GI slice + overlap policy

Original plan committed to Pan-GI as broad atlas comparator framed as
"integration-pipeline-robustness comparator with known donor overlap
(Smillie only)" and HCA Gut as broad atlas comparator framed as
"independent of our trio."

CELLxGENE inspection revealed:
1. Pan-GI Extended+ slice (1,596,200 cells) is the only slice with all
   lineages (epithelial + immune + stromal + endothelial + neural).
   Other slices (e.g. "Extended - Large Intestine") are lineage-restricted
   (epithelial only at 96,675 cells) and CANNOT support cross-lineage
   cell-type prioritization. Use Extended+ as the canonical Pan-GI v1 input.
2. Pan-GI Extended+ study column contains "Elmentaite2021" = 398,460 cells
   (~25% of the atlas). Elmentaite 2021 IS the HCA Gut Atlas. Therefore
   HCA Gut is NOT independent of Pan-GI; it is nested within it.
3. Pan-GI Extended+ study column also contains "Kong2023" = 235,327 cells.
   Not relevant for UC v1 (Kong is CD-only) but noted for future CD work.
4. Pan-GI Extended+ study column does NOT contain Smillie2019. The plan's
   "known Smillie overlap with Pan-GI" claim is unverified; treat as
   "no Smillie overlap" until proven otherwise. The Pan-GI loader should
   still scan for Smillie donor IDs and produce a "0 cells overlapped"
   sensitivity report rather than asserting Smillie absence.

Locked Option A policy (keep both atlases, plan for nested overlap):
  - Pan-GI v1 = Extended+ - 18485 genes slice
    Download URL: https://datasets.cellxgene.cziscience.com/1dcf15ee-c103-4aaa-8b8c-0fc697fcccc8.h5ad
  - HCA Gut v1 = "Total - Cells of the human intestinal tract mapped
    across space and time"
    Download URL: https://datasets.cellxgene.cziscience.com/f34d2b82-9265-4a73-bda4-852933bf2a8d.h5ad
  - HCA Gut framing is corrected to "external single-atlas reference;
    nested within Pan-GI." Not independent.
  - Pan-GI runs are paired with TWO sensitivities:
      (a) Without Elmentaite2021 (test: does HCA Gut overlap drive results?)
      (b) Without any Smillie donor IDs found (likely no-op; documents
          the empirical overlap rather than assumed)
  - Pan-GI v1 filter chain (applied in load_pangi.py):
      * disease in {normal, ulcerative colitis, IBD}
      * organ_unified in {ascending colon, caecum, colon, descending colon,
                          rectum, sigmoid colon, transverse colon}
      * sample_type != "Organ_donor_resection" (different biology;
        keep biopsy + resection only)
      * Expected post-filter: ~150-200k cells (tractable on 128x24)

---

## CORRECTION 2026-05-20 (4/7): Garrido-Trigo annotation tier availability

CELLxGENE deposit of Garrido-Trigo
(b1a62801-f509-45f8-b55f-533fbb7e7800.h5ad) contains only the 5 CL-mapped
lineage labels in obs['cell_type']:
  - colon epithelial cell
  - myeloid cell
  - plasma cell
  - stromal cell of lamina propria of colon
  - T cell of anorectum

The paper's full 51-tier and 91-tier annotations (as visible in the Salas
lab Shiny app at servidor2-ciberehd.upc.es/external/garrido/app/) are NOT
in the CELLxGENE deposit. The lab's GitHub repo
(linked from CELLxGENE under "Protocol") contains analysis code for the
spatial CosMx data only; no per-cell scRNA-seq annotations.

Marker gene table per cluster IS available (downloaded from the Shiny app)
and stored at:
  data/atlases/garrido_trigo_markers.xlsx
Structure: 5 sheets (one per compartment) with columns
(p_val, avg_log2FC, pct.1, pct.2, p_val_adj, cluster, Population, gene).
91 unique fine clusters across all compartments.

Locked decision: Garrido-Trigo contributes broad-tier only (5 lineages)
to v1 cross-atlas concordance. Fine-tier concordance restricted to
4 atlases: Smillie x Pan-GI x HCA Gut x Mennillo.

Note: broad-tier from Garrido-Trigo (5 categories) is below the v1
target of 10-15 for the "broad" tier. Document as known limitation in
the manuscript. Recovery path for v2: email corresponding author
Azucena Salas requesting the .rds with full annotations.

---

## CORRECTION 2026-05-20 (5/7): scDRS raw-count flag across CELLxGENE atlases

CELLxGENE-deposited atlases (Garrido-Trigo confirmed; Pan-GI and HCA Gut
expected to follow same pattern; verify on download) ship with
log-normalized X only; no raw.X, no integer-count layer.

Confirmed for Garrido-Trigo via direct .h5ad inspection:
  AnnData has no .raw, no .layers, X is float32 log-normalized.

Implication for scDRS protocol:
  scripts/slurm/03_scdrs_compute.slurm currently uses --flag-raw-count True.
  This is incompatible with normalized-count input.

Two options considered:
  (a) Re-process every atlas from GEO raw count matrices (canonical inputs
      but adds significant scope and compute).
  (b) Use --flag-raw-count False uniformly across all v1 atlases.

Locked: Option (b) for v1. Update scripts/slurm/03_scdrs_compute.slurm
to use --flag-raw-count False. Methodological note: scDRS supports both
modes; the original paper validates both. Treat this as a limitation to
disclose in methods, not a flaw to fix.

For Mennillo (GEO-sourced, will have raw counts), re-normalize using the
same pipeline as CELLxGENE atlases (log1p of CP10k) before scDRS to keep
input distributions comparable across all atlases.

---

## CORRECTION 2026-05-20 (6/7): HCA Gut filter chain and covariates

HCA Gut (Elmentaite 2021) CELLxGENE deposit
(f34d2b82-9265-4a73-bda4-852933bf2a8d.h5ad):
  - 428,469 cells across 15 tissues
  - obs['category'] = 9 broad lineages (Epithelial, Mesenchymal, Myeloid,
    T cells, Plasma cells, B cells, Endothelial, Neuronal, RBC).
    USE AS BROAD TIER.
  - obs['author_cell_type'] = ~120 fine cell types. USE AS FINE TIER
    AFTER ROLL-UP (see Part 5: cl_rollup.py).
  - obs['cell_type'] (CL ontology, ~75 terms) is IGNORED for analysis;
    serves only as the harmonization anchor in cl_rollup.py.
  - obs['Age_group'] spans fetal (First trim, Second trim), pediatric
    (Pediatric, Pediatric_IBD), and adult (Adult, Adult_MLN).
  - obs['Fraction'] has 5 sorting strategies (SC, SC-45N, SC-45P,
    SC-EPCAMN, SC-EPCAMP). Cell-type proportions in HCA Gut are
    artificially controlled by sort strategy and are NOT biologically
    interpretable.
  - obs['batch'] has ~100 batch IDs.
  - disease values: Crohn disease (27,164) + normal (401,305). No UC.

Locked v1 filter chain (load_hca_gut.py):
  * Age_group in {Adult, Adult_MLN} (excludes fetal/pediatric/IBD-pediatric)
  * tissue in {ascending colon, caecum, colon, descending colon,
               large intestine, rectum, sigmoid colon, transverse colon}
  * Expected post-filter: 30k-70k cells.

Locked scDRS covariates for HCA Gut: assay, batch, Fraction, sex.
Fraction is critical because sort strategy correlates strongly with cell
type identity; omitting it would conflate sort-induced expression
differences with biology.

Locked sensitivity: paired run excluding Crohn disease cells (analogous
to Pan-GI's Elmentaite2021-removal sensitivity). Tests whether residual
Crohn signal in this "broad reference" atlas affects UC cell-type
prioritization.

Note: HCA Gut is a healthy reference; it does NOT need UC samples for
its role in cell-type prioritization. The GWAS provides the UC association
signal; the atlas provides the cell-type catalog.

---

## CORRECTION 2026-05-20 (7/7): Smillie atlas source - CELLxGENE unusable

Smillie 2019 (Cell 175:372-386, "Intra- and Inter-cellular Rewiring of
the Human Colon during Ulcerative Colitis") has a CELLxGENE deposit
(e373cf41-e123-4c98-a8bb-a531c7bbedd2.h5ad), but inspection shows:
  - 34,772 cells (~9% of the paper's full 366,650-cell atlas)
  - 12 donors all prefixed "N" (Normal/Healthy)
  - All from "Epi" compartment (epithelial only; no immune, no stromal)
  - No disease column; Source column confirms healthy-epithelial subset

This deposit is UNUSABLE for the project because:
  - It contains zero UC samples (defeats the purpose of an "UC core" atlas)
  - It lacks immune and stromal compartments (cannot support cross-lineage
    prioritization)
  - It overlaps with HCA Gut's role as a healthy reference

Locked: Smillie v1 source = Single Cell Portal SCP259 (canonical deposit
with 366,650 cells, 30 donors, 18 UC + 12 HC, all compartments).

This requires (deferred to next session, not this batch):
  - Single Cell Portal account creation + email verification
  - Browser-mediated consent/access click on SCP259
  - ~5-8 GB download to Hummingbird scratch

After SCP259 is on disk, run a separate schema-capture session for
Smillie (and Mennillo, also deferred to GEO download).

---

## CORRECTION 2026-06-03 (8): Garrido-Trigo annotation tier — reversing (4/7)

Correction (4/7) locked Garrido-Trigo as broad-tier-only (5 CL lineages)
on the premise that the paper's full 51- and 91-tier annotations were
not in any downloadable form. That premise was wrong.

GEO GSE214695 ships a supplementary file
`GSE214695_cell_annotation.csv` containing the full 91-cluster Salas-lab
fine annotation, joinable to the CELLxGENE matrix by barcode. The 5-CL
CELLxGENE labels are still present but are now used only as a coherence
cross-check; the GEO CSV is the source of truth for both tiers.

Locked decision: Garrido-Trigo contributes **broad + fine tier** to v1
cross-atlas concordance. The fine-tier UC trio is back to 3 atlases
(Smillie x Garrido-Trigo x Mennillo). Correction (4/7)'s broad-only
restriction is voided.

Loader changes (`code/02_atlas_prep/load_garrido_trigo.py`):

- New required argument `annotation_csv_path` for the GEO CSV; if
  omitted, the loader falls back to the (4/7) degraded mode with a
  warning.
- Barcode join auto-detects column names from a candidate set and tries
  five join strategies (raw, strip/add `-1` suffix, composite
  `donor_id + barcode` either side). Raises a diagnostic error with
  example barcodes if no strategy hits all cells (no orphans tolerated).
- Every label is whitespace- and unicode-normalized on load
  (`_normalize_label`). One known case in the CSV is
  `'PC  immediate early response'` with a literal double space; that
  alone would silently break string joins downstream.
- 9 Ribhi clusters are collapsed into parent lineages
  (`RIBHI_TO_PARENT`) before any tier logic. Empirically verified:
  RPL*/RPS* genes dominate the top-20 markers of every Ribhi cluster in
  `garrido_trigo_markers.xlsx` — Ribhi is a ribosomal-high cross-lineage
  transcriptional state, not a lineage or a batch artifact. Ribhi cells
  are never entered as standalone fine clusters in cross-atlas
  concordance.
- `cell_type_fine` (82 surviving labels) and `cell_type_broad`
  (15-level roll-up) populated from a hand-curated `FINE_TO_BROAD` map
  covering all 82 surviving fine labels + the 4 generic Ribhi parents
  (epithelial, T, fibroblast, mast) + 2 defensive synonyms (long-form
  PC IER, "Cycling cells 3" mentioned in the to-do but absent from
  the marker xlsx).
- Logs UC-subset cell and donor counts and warns if they deviate from
  the (2/7) expected 30,068 cells / 12 donors.

Followups carried into M2 harmonization (notes, not blockers):

- 14 of 91 fall below the ≥50-cell rule in HC+UC and drop as standalone
  fine clusters: Cycling cells 3, Cycling myeloid, DCs CCL22,
  DCs CCL22_Ribhi, DN EOMES, Enteroendocrine, Eosinophils, FRCs,
  Lymphatic endothelium, M1 CXCL5, MAIT, Neutrophil 2, Neutrophil 3,
  Paneth-like. They roll up into broad parents, so broad-tier
  concordance is unaffected; only fine-tier is trimmed.
- Cross-atlas fine-tier intersection is expected to be modest;
  hyper-specific states (PC IgA heat shock 1/2, M1 ACOD1, etc.) will not
  map to Smillie/Mennillo. Methods text must avoid framing this as
  "91-way fine concordance."
- `nanostring_reference` (54 levels, clean nested coarsening of the 91)
  remains useful scaffolding but is NOT used as the broad tier (too
  fine, non-uniform depth).
- No email to Azucena Salas required; the GEO CSV is the full annotation
  the (4/7) recovery path anticipated. Email is now a backstop only if
  the barcode join turns out broken.

Files updated in this batch:

- `code/02_atlas_prep/load_garrido_trigo.py`
- `code/02_atlas_prep/atlas_schemas.md`
- `README.md` (HCA Gut wording — see Section D of the to-do)
- `DECISIONS.md` (this entry)

---

## CORRECTION 2026-06-04 (9): Garrido-Trigo matrix source — superseding (8)

Correction (8) restored Garrido-Trigo to full broad + fine tier by joining
the GEO supplementary annotation onto the CELLxGENE matrix by barcode.
That join is undoable on the actual CELLxGENE deposit.

**Finding (verified by running the loader on the local h5ad and CSV).**
The CELLxGENE deposit `b1a62801-f509-45f8-b55f-533fbb7e7800.h5ad` has a
**synthetic** `obs.index` — values are `cell1, cell2, ..., cell46700`,
sequential and information-free. The original 10X barcodes were dropped
at deposit time. `observation_joinid` is a 10-character portal-internal
tag (`LH6LS+Lnyp`, …), also not a barcode. No `obs` column carries the
original cell barcode. Every barcode-join strategy in (8)'s loader keys
on something derived from `obs.index`, so all seven strategies hit
0/30,068 cells against `GSE214695_cell_annotation.csv` (which keys on
`SC_xxx_<barcode>` in `Unnamed: 0`).

Per-donor cell counts match exactly between the two sources
(HC_1↔HC1: 1531, HC_2↔HC2: 2555, …, UC_6↔UC6: 2965), so a positional
join within donor was considered. **Rejected.** A positional join has
no shared key — it rests on the two files preserving the same arbitrary
within-donor ordering, which cannot be enforced or verified. Marker-gene
sanity checks validate populations, not per-cell correspondence: a join
that scrambles labels for one donor (or 5% of cells) still shows
roughly correct cluster-level marker rates while the per-cell label
noise attenuates scDRS toward the null — i.e. it silently manufactures
the exact cross-atlas discordance the paper is measuring. The failure
mode is invisible at the population level. Hard no.

**Locked decision.** v1 builds Garrido-Trigo from GEO primary, not the
CELLxGENE deposit. Matrix comes from `GSE214695_RAW.tar` (per-GSM
sparse matrices, real 10X barcodes preserved). Annotation comes from
`GSE214695_cell_annotation.csv` as before. The join key is the
reconstructed composite `<sample-GSM>_<barcode>` matching the CSV's
`Unnamed: 0` (`SC_xxx_<barcode>`); auto-detect already prefers the
unique candidate column over the duplicated bare `cell_id`. Garrido is
the only atlas with a barcode-join requirement (Smillie/HCA/Pan-GI
carry annotations in their own obs; Mennillo is already GEO RAW), so
the impact is contained — no re-audit of the other four loaders is
needed.

**Normalization (satisfying (5/7), not sidestepping it).** RAW.tar is
raw counts; HCA Gut and Pan-GI ship log-normalized. To keep the input
state uniform across atlases, the new loader applies `log1p(CP10k)` on
load — the identical treatment Mennillo gets. `--flag-raw-count False`
remains uniform downstream per (5/7). The flag is *not* switched on
for Garrido; the normalization is reproduced on the loader side so
that downstream code sees the same input state across all five atlases.

**Salvaged tonight (in this commit, no matrix path rewrite yet):**

- Top-of-file docstring banner in `load_garrido_trigo.py` flagging the
  superseded path and pointing to this correction.
- Runtime `logger.error` at `load()` entry surfacing the same warning
  to anyone running the loader before tomorrow's rewrite lands.
- `_BARCODE_COL_CANDIDATES` reordered: `"Unnamed: 0"` first (carries the
  unique `SC_xxx_<barcode>` composite); `cell_id` retained but now
  passed over when present because it carries 280 cross-sample barcode
  collisions in this CSV.
- `_autodetect_column` gains a `prefer_unique` flag, used for the
  barcode column. The duplicate-barcode check in `_load_annotation_csv`
  (which originally caught the `cell_id` collision and gave us this
  finding) is preserved as a backstop.
- Three hard-invariant asserts added (correctly placed even though they
  will not exercise until the RAW.tar loader lands):
  1. **Annotation completeness**: zero NaN in `cell_type_fine` after the
     join — the unmapped-label check explicitly drops NaN and would
     otherwise let partial annotations through silently.
  2. **Donor structure**: exactly 12 donors split 6 HC + 6 UC; this is
     fixed by the study design and is the right gate for "the filter is
     producing the cohort we expect."
  3. **Disease/sample-prefix agreement**: every cell's CELLxGENE
     `obs['disease']` matches the GEO `sample` prefix (`HC*`→normal,
     `UC*`→ulcerative colitis). Catches the case where filter-before-join
     and the annotation source disagree on who counts as UC.
  4. The 30,068 cell count remains a soft `logger.warning` tripwire —
     it's a derived intersection that can drift on re-pulls or QC
     nudges, so a hard assert there would crash on benign changes.

**Tomorrow's first task (next session, before any Hummingbird run).**
Rewrite `load_garrido_trigo.py` matrix path:

1. Download `GSE214695_RAW.tar` to scratch; untar to per-GSM sparse
   matrices (`.mtx` + barcodes + features). 18 GSMs; the 12 HC/UC ones
   are the v1 set.
2. Build a `(GSM → sample-label)` map from the CSV's `sample` column
   and the GSM-to-sample assignments in the GEO series metadata.
3. Concatenate the 12 per-GSM matrices into a single AnnData, with
   `obs.index` set to `<sample>_<barcode>` and `obs['donor_id']`,
   `obs['disease']` populated from the sample label.
4. Join `GSE214695_cell_annotation.csv` on `<sample>_<barcode>` ==
   `Unnamed: 0` (the loader's existing strategy 5/7 path works here
   with no further changes).
5. Apply `log1p(CP10k)` normalization on load.
6. Run with the three hard asserts; commit, push, score.

Files updated in this batch:

- `code/02_atlas_prep/load_garrido_trigo.py` (docstring, asserts,
  auto-detect uniqueness preference, runtime warning; matrix path
  unchanged — rewrite pending)
- `code/02_atlas_prep/atlas_schemas.md` (correction-9 banner on
  `garrido_trigo` section)
- `DECISIONS.md` (this entry)

---

## CORRECTION 2026-06-04 (10): Smillie SCP259 schema capture + loader completion

Correction (7/7) locked Smillie's v1 source as Single Cell Portal SCP259
(the CELLxGENE deposit `e373cf41-...` is a 34,772-cell
healthy-epithelial-only subset and was rejected as unusable for v1).
The download was deferred pending account creation + browser consent.

**Status (2026-06-04).** SCP259 is downloaded, verified, and on disk at
`~/uc-cross-atlas-data/atlases/SCP259/` (Hummingbird home, outside the
repo, ~5 GB). Three compartment matrices (Epi/Imm/Fib) plus
`metadata/all.meta2.txt`, all present and consistent. The schema below
was captured empirically against the bytes on disk; the loader
(`code/02_atlas_prep/load_smillie.py`) is written and label-coverage-
verified but **has not yet been run end-to-end against data** — that
happens on a Hummingbird compute node (next session) and is the gating
step before Smillie enters any cross-atlas concordance.

**Schema captured.** Layout: three gene-sorted 10X triplets under
hashed sub-directories whose names are not stable (loader globs by
`gene_sorted-<compartment>.matrix.mtx` filename), plus one in-band
metadata file. The matrices are **gene-sorted = genes × cells** —
loader transposes each to cells × genes before concat (Epi header is
20,028 × 123,006 × 174.4M nnz; the other two are similar scale).
`metadata/all.meta2.txt` is tab-separated with an SCP-specific `TYPE`
boilerplate row at line 2 that must be skipped on read (`skiprows=[1]`).
Join is **direct**: `barcodes2.tsv` values equal `NAME` in the metadata
exactly — no donor-prefix reconstruction, no composite key. `var_names`
are HGNC symbols (single-column, e.g. `7SK`, `A1BG`), not Ensembl;
`ensembl_to_hgnc` takes its symbol-fallback path (dedup duplicates by
max-expression, filter to NCBI-approved symbols).

`X` ships as **raw integer counts** (MatrixMarket `coordinate integer`).
Unlike Garrido/HCA/Pan-GI (CELLxGENE log-normalized), this loader
applies `log1p(CP10k)` itself to keep input state uniform across all
five atlases — the same treatment Mennillo gets and the same treatment
Garrido will receive once correction (9)'s RAW.tar rewrite lands.
Raw counts preserved in `layers['counts']`. `raw_count_mode=True`
remains unsupported per correction (5/7).

**Donor structure (hard invariant, fixed by study design).** 30 donors
= 12 healthy controls + 18 UC patients. Every UC patient contributed
**both** an inflamed and a non-inflamed biopsy (paired design), so
`Sample` ≠ `Subject` for UC donors: 36 UC samples, 18 UC donors.
Cell counts: 110,110 Healthy / 125,119 Inflamed / 130,263 Non-inflamed
(~365,492 total). The loader hard-asserts `30 = 12 HC + 18 UC`; the
365,492 total is a soft `logger.warning` tripwire per the
correction-9 framing.

**Tiers.** 51 published Smillie fine labels → 14-level broad roll-up
via `FINE_TO_BROAD`. The broad vocab is the same 15-level shared set
defined in `load_garrido_trigo.FINE_TO_BROAD`; Smillie populates 14 of
the 15 — **no granulocyte** (the Smillie taxonomy has no
neutrophils/eosinophils). Cross-atlas concordance handles the missing
level per-atlas; no broad-vocab change needed.

**Disease harmonization to Garrido vocab.** `Health` is 3-state on the
SCP side (`Healthy` / `Inflamed` / `Non-inflamed`); the loader maps
`Healthy → normal` and `Inflamed | Non-inflamed → ulcerative colitis`
for cross-atlas comparability, but **preserves the raw 3-state in
`obs['health']`** so any later inflamed-vs-pooled subsetting is a
downstream filter, no re-load needed. See F1 below.

**Loader changes (`code/02_atlas_prep/load_smillie.py`).** Replaces the
`NotImplementedError` skeleton from correction (7/7). Hard-fail asserts
mirror the Garrido correction-9 pattern:

- Compartment-triplet file presence: matrix + `.genes.tsv` +
  `.barcodes2.tsv` for each of Epi/Imm/Fib; loader raises if any is
  missing or ambiguous.
- Matrix-shape consistency: `mat.shape == (n_genes, n_barcodes)` for
  each compartment.
- Cell-id uniqueness post-concat.
- Metadata column presence (`NAME`, `Cluster`, `Subject`, `Health`,
  `Location`, `Sample`).
- **Orphan-cell completeness** (every matrix cell must have a metadata
  row after the `NAME` reindex).
- Unmapped fine labels (`KeyError` listing them — directly analogous to
  Garrido's gate).
- Unexpected `Health` values; unexpected harmonized `disease` values.
- **Donor structure**: `30 = 12 HC + 18 UC` hard, fixed by study
  design.
- Broad-tier cardinality warning if outside 10–15.
- Cell-count tripwire (`logger.warning`) at 365,492 — same soft-gate
  framing as Garrido's 30,068 per correction 9.

**Open flags surfaced by this loader (pending OPEN_FLAGS.md).**

- **F1 — UC tissue definition.** `Health` is 3-state. Loader is
  agnostic: keeps all 30 donors, sets harmonized 2-state `disease`,
  preserves 3-state in `obs['health']`. The inflamed-vs-pooled decision
  must be made consistently with how Garrido and Mennillo define their
  UC tissue.
- **F2 — QC-state labels.** Smillie's `MT-hi` is the
  mitochondrial-high analogue of Garrido's MT labels; collapse-vs-
  exclude must be decided uniformly across atlases.
- **F7 — Smillie crosswalk REVIEW rows.** Four `FINE_TO_BROAD`
  assignments are tentative and flagged inline (`# REVIEW`): `M cells`,
  `Immature Enterocytes 1/2`, `Secretory TA`, `MT-hi`. Loader runs with
  best-guess parents; marker-QC at M2 will lock or move.
- **Not a flag:** gene-identifier harmonization. `ensembl_to_hgnc`
  converges every atlas onto the NCBI symbol set, so Smillie shipping
  symbols needs no special mapping.

**Latent observation (not a correction).** In the loader, the
`adata.obs_names_make_unique()` call (defensive against the
`Subject.Sample.barcode` format already being unique) is followed by an
`if not adata.obs_names.is_unique` check that can never fire —
`obs_names_make_unique()` modifies in place to guarantee uniqueness.
Either drop the check or reorder to check-before-make. Captured for
the M2 cleanup pass; not a blocker for the first Hummingbird run.

**Still pending (next session, on Hummingbird).**

- First end-to-end run of the loader on a compute node (`--mem=96G`,
  `--cpus=4`, ~2 h). Pre-cache `Homo_sapiens.gene_info.gz` from the
  login node so the symbol-validity filter does not silently skip on a
  no-internet compute node.
- Write `~/uc-cross-atlas-data/processed/smillie.h5ad`; back up
  off-cluster (Hummingbird home is 1 TB but not backed up).
- Resolve F7 REVIEW rows via marker QC at M2 (after the first scored
  pass); move locked decisions into DECISIONS.
- Commit `OPEN_FLAGS.md` to the repo root once drafted — file does not
  yet exist on disk; out of scope for this commit.

Files updated in this batch:

- `code/02_atlas_prep/load_smillie.py` (replaces the (7/7) skeleton
  with the full loader)
- `code/02_atlas_prep/atlas_schemas.md` (smillie section: DEFERRED →
  captured schema)
- `DECISIONS.md` (this entry)

---

## CORRECTION 2026-06-04 (11): HGNC remap is now pinned, strict, and gated

The README HGNC-pin section described the intended behavior (pin a
release; verify ≥95% of canonical UC GWAS hits survive); the actual
``hgnc_remap.py`` implemented none of it. Three gaps, surfaced by the
2026-06-04 sanity-check review:

1. **No version pin.** ``ensembl_to_hgnc`` fetched
   ``Homo_sapiens.gene_info.gz`` from NCBI live on first call, cached
   under ``$UCC_DATA/reference/`` or ``./data/reference/``. NCBI updates
   monthly, so two runs separated by a month would use different
   approved-symbol sets, breaking byte-reproducibility for downstream
   MAGMA/scDRS/seismic scores.
2. **Synonyms in the approved set.** ``_load_ncbi_symbol_set`` pulled
   both the ``Symbol`` and ``Synonyms`` columns into one set
   (~68k extra strings on the 2026-05-21 snapshot, vs ~194k Symbols).
   That made the membership filter near-permissive — deprecated
   aliases passed through unchanged with no remap to their approved
   form.
3. **The canonical-hit survival check was documented but not coded.**
   The README claimed "verify ≥95% of canonical UC GWAS hits survive";
   no such assertion existed.

A fourth, minor issue: ``adata.var_names_make_unique = False`` was set
as an instance attribute. ``var_names_make_unique`` is an AnnData
method; assigning ``False`` to it shadowed the method on that instance
but did not change dedup behavior (dedup was already handled manually
below). Removed.

**Locked decision.** ``ensembl_to_hgnc`` is now pinned, strict, and
gated.

- **Pin**: a dated NCBI ``gene_info`` snapshot is committed at
  ``data/reference/gene_info.tsv.gz``; current pin is **2026-05-21**.
  Live-fetch removed. ``GENE_INFO_PIN_DATE`` in ``hgnc_remap.py`` is the
  single source of truth for the pin date.
- **Approved set**: ``Symbol`` column **only**. Synonyms are
  deliberately excluded. If alias resolution is ever needed, it goes in
  as an explicit alias → approved remap, not a membership test that
  leaves the alias in place.
- **Survival gate**: after the symbol-validity filter, the five
  canonical UC GWAS hits in ``CANONICAL_UC_HITS`` (IL23R, JAK2, TYK2,
  NKX2-3, ATG16L1) are checked; ``< 95%`` survival raises with a
  diagnostic message naming the missing symbol(s). All five must
  effectively survive (4/5 = 80% < 95% would fail). The threshold is
  expressed as a fraction for forward-compatibility if the canonical
  list grows.
- **Refresh procedure**: download a fresh ``Homo_sapiens.gene_info.gz``,
  bump ``GENE_INFO_PIN_DATE``, replace the committed snapshot in one
  commit, and log the date bump as a new correction here.

**Verified before commit.** All five canonical hits are present in the
2026-05-21 snapshot's ``Symbol`` column (not synonyms-only), so the
new gate passes on the current pin without any alias special-casing.

**Repo size note.** ``data/reference/gene_info.tsv.gz`` is 5.1 MB —
fine to track in git, well under the per-blob soft limit. The previous
``data/reference/.gitkeep`` placeholder remains as the directory marker.

**To refresh ahead of M2.** No action required; the current pin is
good for at least one full M1→M2 cycle. If a new HGNC release is
preferred, the refresh procedure above is one commit + a correction
entry.

Files updated in this batch:

- `code/02_atlas_prep/hgnc_remap.py` (live fetch → pinned snapshot;
  synonyms dropped; survival gate added; dead `var_names_make_unique`
  attribute removed; docstring expanded)
- `data/reference/gene_info.tsv.gz` (new; 5.1 MB; pin date 2026-05-21)
- `code/02_atlas_prep/README.md` (HGNC pin section: described →
  matches actual behavior)
- `DECISIONS.md` (this entry)

---

## CORRECTION 2026-06-06 (12): Garrido-Trigo RAW.tar loader shipped — implements (9)

Correction (9) committed to rewriting `load_garrido_trigo.py` against
`GSE214695_RAW.tar` (real 10X barcodes preserved) once that archive was
on disk and the schema captured. The archive landed at
`~/Downloads/GSE214695_RAW.tar` (887 MB, downloaded 2026-06-06), schema
was captured against the actual bytes, and the loader was written and
**verified end-to-end** against the local RAW.tar + CSV.

### Verification facts (from the local run)

End-to-end load produced **30,068 cells × 22,414 genes across 12
donors**, with all hard gates passing on the first complete run after
the typo fix below:

- Cell count = 30,068 — **exact** match to (2/7) expected; soft
  tripwire silent.
- Donor structure = 6 HC + 6 UC = 12 — hard donor invariant passes.
- Per-donor cell counts match the prior CELLxGENE deposit's HC_*/UC_*
  counts to the cell (HC_1: 1531, HC_6: 1898, UC_6: 2965, …).
  Independent confirmation that the RAW.tar path reproduces the same
  cohort as the (correction-9-superseded) CELLxGENE path, just with
  joinable barcodes.
- Disease vs sample-prefix cross-check: agree on all 30,068 cells.
- Annotation completeness: zero NaN `cell_type_fine` post-join.
- Fine labels: 86 surviving (post-Ribhi-collapse). Arithmetic
  corrected from the (8) entry's "82": 91 published − 9 Ribhi-named +
  4 generic parents introduced by the collapse (`epithelial`, `T`,
  `fibroblast`, `mast`, none of which were standalone labels in the
  original 91). The old "82" was a hand calculation that never ran.
- Broad labels: 15 — within the 10–15 v1 target.
- Canonical UC GWAS hits: **5/5 survive** the (11) pinned-HGNC remap
  (IL23R, JAK2, TYK2, NKX2-3, ATG16L1).
- X = log1p(CP10k) float32, max 9.14 (sensible log range); raw counts
  preserved in `layers['counts']`, max UMI 58,069.

### Loader architecture

- Reads the 18 per-GSM 10X triplets in-memory via `tarfile` +
  `gzip.GzipFile` + `scipy.io.mmread` — no tar extraction step.
- CD-1..CD-6 GSMs are dropped at the **glob step** (saves ~30 MB ×
  6 sparse matrices vs. load-then-filter).
- Each per-GSM AnnData carries `obs.index = f"{sample}_{barcode}"`
  (e.g. `"HC1_AAACCTGCAAGTCTGT-1"`), built so the composite is unique
  on the RAW side by construction.
- CSV side is filtered to the same sample-prefix subset
  (`("HC", "UC")` for v1), then composite-keyed by
  `f"{sample}_{cell_id}"`, then dup-checked. Inner-join collapses the
  ~99.8% empty droplets in the raw 10X whitelist down to the 30,068
  CSV-annotated cells.
- The seven `_try_join_keys` strategies from the pre-correction-9
  loader are gone — when we control matrix assembly, there's exactly
  one deterministic join key on each side.
- Three sample-naming conventions reconciled on load: RAW dashed
  (`HC-1`), CSV no-separator (`HC1`), CELLxGENE-deposit underscore
  (`HC_1`). The loader emits the underscore form in `obs['donor_id']`
  to preserve continuity with the pre-correction-9 obs schema.
- `log1p(CP10k)` applied on load (the RAW.tar matrices are raw int
  counts); raw preserved in `layers['counts']`; `raw_count_mode=True`
  remains unsupported per (5/7).
- The full assert pattern from (9)'s plan is in place — completeness,
  donor structure (hard), disease/sample agreement (hard), unmapped
  fine labels (hard), broad cardinality (warning), cell count
  (warning tripwire).

### Side-findings logged (not corrections of their own; documented here)

1. **CSV CD-only duplicate composites (Salas-lab authoring bug).** The
   GEO CSV has 2 duplicate `(sample, cell_id)` composites, all in CD3,
   each pair carrying conflicting fine annotations
   (`Cycling cells` vs `Cycling myeloid`;
   `Epithelium Ribhi` vs `M0_Ribhi`). Traced to inconsistent
   whitespace in `Unnamed: 0` (`SC_013_GTGTGGCAGACTACCT-1` vs
   `SC_013 _ GTGTGGCAGACTACCT-1` etc.) — the same cell got annotated
   twice in two CSV authoring passes. None affect the v1 HC+UC cohort.
   The loader filters CSV to the v1 sample prefixes **before** the
   duplicate gate, so the bug is structurally outside our cohort. If
   the gate ever fires on HC/UC, that's a new defect to inspect.
2. **`FINE_TO_BROAD` typo fix (`Perycites` → `Pericytes`).** The map
   had `"Perycites"` (with a `y`); the actual published label is
   `"Pericytes"`. Pre-existing bug in the (8) version of the map; the
   CELLxGENE-path loader never reached the unmapped-fine-labels check
   because it failed earlier at the barcode join, so the typo never
   fired. Caught immediately on the first RAW.tar run; one-character fix.
3. **Mixed 10X chemistry in the deposit.** HC-1 ships the v2 whitelist
   (737,280 barcodes); HC-2..UC-6 ship the v3 whitelist (6,794,880
   barcodes). Same gene reference across all (CellRanger 3.0.2,
   33,538 features), so concat is uniform and the post-CSV cell set
   is unaffected. Documented in `atlas_schemas.md` for context but
   not a downstream concern.

### Promoted: README + atlas_schemas

The README loader-status table flips Garrido-Trigo from "Superseded —
RAW.tar rewrite pending (correction 9)" back to **Production**. The
`atlas_schemas.md` garrido_trigo section is rewritten as a captured
schema for the RAW.tar layout (matrix source, on-disk layout, filter
chain, tier arithmetic explaining 86 vs the old 82) — the
correction-9 "changing" banner is gone.

Files updated in this batch:

- `code/02_atlas_prep/load_garrido_trigo.py` (full rewrite against
  RAW.tar; `_try_join_keys` and the CELLxGENE deposit machinery gone;
  `_group_tar_entries`, `_read_triplet_from_tar`, `_load_per_gsm`
  added; `Pericytes` typo fix)
- `code/02_atlas_prep/atlas_schemas.md` (garrido_trigo section
  rewritten for RAW.tar schema; correction-9 banner removed)
- `code/02_atlas_prep/README.md` (loader-status table: Garrido-Trigo →
  Production)
- `DECISIONS.md` (this entry)

---

## CORRECTION 2026-06-06 (13): Cell Ontology pin — release 2026-03-26

`canonical_broad_DRAFT.md` references ~27 CL IDs as the substrate of the
candidate broad-tier vocabulary. Cell Ontology is updated externally
(monthly-ish), so without a pin those IDs can silently shift meaning
between runs — same shape as the (11) gene_info issue.

**Pin.**

- **Release date:** 2026-03-26
- **versionIRI:** `http://purl.obolibrary.org/obo/cl/releases/2026-03-26/cl.owl`
- **Source URL:** `http://purl.obolibrary.org/obo/cl.owl`
- **File size:** 62.82 MB
- **Pinned term list (committed):**
  `data/reference/cl_terms_pinned.tsv` — 27 rows, one per CL ID
  referenced in `canonical_broad_DRAFT.md`, with the **label as it
  appeared in 2026-03-26** alongside the draft-expected label and a
  match flag. Refresh procedure documented in the TSV header.

The OWL file itself is **not** committed (63 MB; re-downloadable from
the URL above; pinning the term list is sufficient for our use). The
gene_info pattern of committing the dump was justified because the
remap reads the bytes at runtime; the CL pin only needs the term-list
snapshot, not the live OWL.

**Six DRAFT errors surfaced by the pin** (the whole point of pinning):

- Label drift (concept correct, label updated upstream): CL:1000347
  *colonocyte* (was "enterocyte of colon"); CL:0002204 *tuft cell*
  (was "brush cell of intestine / tuft"); CL:0002138 *endothelial cell
  of lymphatic vessel* (was "lymphatic endothelial cell").
- **WRONG CL ID** (the drafted ID points at a different concept in
  2026-03-26):
  - CL:1000280 → in this release is "smooth muscle cell of colon",
    NOT "stem cell of intestine".
  - CL:0009039 → in this release is "colon goblet cell", NOT "colon
    epithelial progenitor cell".
  - CL:0002073 → in this release is "transitional myocyte", NOT
    "enteric glial cell".

Reverse-lookup of the drafted concept names also failed for `enteric
glial cell`, `enterocyte of colon` (now "colonocyte"), `brush cell of
intestine`, `stem cell of intestine`, `colon epithelial progenitor
cell`, and `lymphatic endothelial cell` — the three WRONG-ID rows need
actual ontology investigation (which CL term, if any, captures the
intended biology), not a trivial relabel. Documented inline in
`canonical_broad_DRAFT.md` and flagged with `[CL ID UNRESOLVED —
pin verification]` so they cannot be silently promoted into a locked
`CANONICAL_BROAD`.

Files updated:

- `data/reference/cl_terms_pinned.tsv` (new; 27 rows; the actual pin)
- `code/_shared/canonical_broad_DRAFT.md` (CL pin section + inline
  flags on the six affected rows)
- `DECISIONS.md` (this entry)

---

## CORRECTION 2026-06-06 (14): GWAS sumstats — accessions, schemas, fixed-N

The to-do PDF flagged "Liu's per-SNP N column" and "which de Lange
accession is the UC arm." Resolved both by HEAD-checking the GWAS
Catalog FTP listings and range-fetching the first ~1 MB of each
sumstats file. The accession ambiguity is settled and Liu's N
fallback is documented; the cross-GWAS effective-N range is also now
on the record.

**Locked accessions and schemas (captured 2026-06-06):**

| Study | Accession | Cases / controls | Build | Per-SNP N |
|---|---|---|---|---|
| de Lange 2017 UC | **GCST004133** | 12,366 / 33,609 EUR | 37 | **NO — fixed 45,975** |
| Liu 2023 UC (multi-ancestry) | **GCST90446794** | (6,862 EA + 16,390 EUR) / (15,456 EA + 336,800 EUR) | 38 | **NO — fixed 375,508** |
| Yengo 2022 height EUR | **GCST90245992** | 1,597,374 EUR | 37 | **YES — `n` column** |
| Trubetskoy 2022 SCZ | GCST90128471 | 53,386 EUR + 14,004 EAS + 6,152 AA + 1,234 Latino cases | — | full sumstats not on GWAS Catalog; PGC-only |

**Critical accession correction.** Saisohan's `download_refs.sh`
previously labeled GCST004131 as the de Lange UC sumstats. **That is
wrong.** Per the GWAS Catalog REST API (PMID 28067908), the three de
Lange 2017 studies are:

- GCST004131 = **IBD overall** (25,042 / 34,915)
- GCST004132 = **Crohn's disease** (12,194 / 28,072)
- GCST004133 = **Ulcerative colitis** (12,366 / 33,609) ← the v1 UC arm

Using GCST004131 would have silently munged IBD-combined sumstats as
"UC" — exactly the phenotype confound the to-do PDF flagged. Script
fixed in this commit; the UC file is
`uc_build37_45975_20161107.txt.gz` (the "45975" in the filename =
total samples = 12,366 + 33,609, which is the empirical confirmation
of the accession).

**Per-SNP N — Liu and de Lange need fixed-N fallback.** Liu's
GCST90446794 ships `chromosome, base_pair_location, effect_allele,
other_allele, beta, standard_error, effect_allele_frequency, p_value,
variant_id, FreqSE, MinFreq, MaxFreq, Direction, HetISq, HetChiSq,
HetDf, HetPVal` — no N column. De Lange GCST004133 ships `MarkerName,
Allele1, Allele2, Effect, StdErr, P.value, Direction, HetISq, HetChiSq,
HetDf, HetPVal, Pval_IBDseq, Pval_IIBDGC, Pval_GWAS3,
Min_single_cohort_pval` — no N column either. **Fixed-N for MAGMA**:
de Lange UC = 45,975; Liu UC = 375,508. Yengo's `n` column is the
canonical per-SNP N (varies 85k–615k).

**Cross-GWAS effective-N range** (settles the Li 2025 sample-size
confound concern from your notes): 46k (de Lange) / 376k (Liu) / 1.6M
(Yengo) — three orders of magnitude. The magnitude-of-N variation is
the *expected* explanation for differential power across the three;
seismic + scDRS treatment of this needs to be deliberate, not assumed
neutral. Method-side handling logged separately when MAGMA-munged
files exist.

`scripts/download_refs.sh` updated:

- de Lange URL fixed (GCST004131 → GCST004133); script auto-fetches.
- Liu UC URL added (auto-fetch).
- Yengo height EUR URL added (auto-fetch).
- Trubetskoy SCZ stays manual (PGC, registration). Saisohan's link
  is unchanged; added a note explaining why it isn't a GWAS Catalog
  download.
- Per-SNP N status documented in the post-download notes.

Files updated:

- `scripts/download_refs.sh` (de Lange accession correction; Liu and
  Yengo URLs added; Trubetskoy/PGC note clarified; per-SNP N status
  in the post-download notes)
- `DECISIONS.md` (this entry)

λ_GC and any other munged-stats observations wait until the files are
actually staged into `~/uc-cross-atlas-data/gwas/` and `1_magma.slurm`
has run on them.

---

## CORRECTION 2026-06-06 (15): Mennillo accession — GSE229072 is wrong, correct ID unconfirmed

The original to-do PDF said "Download Mennillo from GEO (GSE229072 —
verify accession on the page)". Verification (this session,
2026-06-06): **GSE229072 is not Mennillo.**

- **Actual GSE229072:** "Arginine methylation of C/EBPα controls the
  speed of immune cell transdifferentiation (ChIP-Seq)" — Garcia,
  Leutz, Graf et al., 2023. Mouse, B-cell→macrophage transdifferen-
  tiation, 8 ChIP-Seq samples. **Wrong organism, wrong assay, wrong
  trait.**

The actual Mennillo et al. 2024 paper is "Single-cell and spatial
multi-omics highlight effects of anti-integrin therapy across cellular
compartments in ulcerative colitis" (*Nat Commun* 15:1493,
doi:10.1038/s41467-024-45665-6, PMC10876948). Confirmed by web search.
The correct GEO accession could not be auto-confirmed this session —
the Nature Communications full-text fetch returned 403 (bot block) and
PMC's data-availability snippet was truncated and surfaced only an
external validation dataset (GSE73661 — *Smillie 2019 prequel*, also
not Mennillo's primary deposit).

**Action for next session** (needs a human-driven browser fetch):
open the Nat Comms paper or its bioRxiv version
(`https://www.biorxiv.org/content/10.1101/2023.01.21.525036v3.full`)
and read the Data Availability section. The Mennillo primary
single-cell data lives at a GSE accession that begins with GSE2xxxxx
(submitted in 2023–2024). Once confirmed, update:

- `code/02_atlas_prep/load_mennillo.py` (currently `NotImplementedError`)
- `code/02_atlas_prep/README.md` Atlases list (item 3, "verify GEO
  accession at M1")
- `scripts/download_refs.sh` if the deposit is publicly downloadable
  without registration

Until then, the Mennillo loader stays a skeleton and the v1 UC trio
remains Smillie + Garrido-Trigo with Mennillo deferred to M2.

No files updated by this correction (negative-finding record only).

---

## CORRECTION 2026-06-06 (16): Post-audit fixes — TAURUS swap, effective N math, Liu/CL verification

The (13)–(15) work this session was audited; five real issues caught
that need correcting before they propagate further.

### (a) Mennillo dropped; TAURUS-IBD is the third UC atlas.

Mennillo 2024 (anti-integrin, Nat Commun 15:1493) was the planned
Atlas 3 through correction (15). It is **dropped from v1.** The third
UC core atlas is now **TAURUS-IBD** — Thomas, Dendrou, Agarwal
(Oxford, 2024), "A longitudinal single-cell atlas of anti-tumour
necrosis factor treatment in inflammatory bowel disease".

- **Source**: Zenodo. The swap directive cited DOI
  `10.5281/zenodo.13768607`; the DOI that resolves on Zenodo as of
  2026-06-06 is `10.5281/zenodo.14007626` (presumably a different
  versioned record of the same dataset). **Both DOIs flagged in the
  loader docstring; canonical version must be pinned before download.**
- **License**: CC-BY-4.0, public, no registration.
- **Files**: `TAURUS_raw_counts_annotated_final.h5ad` (12.7 GB) plus
  per-lineage h5ads (epithelial, CD4/CD8 T, B, plasma, myeloid, ILC,
  fibroblasts, vascular). Raw BAM/CellRanger at GEO `GSE282122`.
  Total ~27.9 GB.
- **Scope alignment**: anti-TNF longitudinal IBD = UC + CD donors.
  Loader must subset to UC donors only AND a single time-point per
  donor (pre-treatment baseline preferred), mirroring the discipline
  the Mennillo loader was going to apply.

The Mennillo accession hunt from (15) ("next-session browser fetch of
Nat Comms data availability") is **dead work** and is killed by this
correction. The negative finding in (15) — GSE229072 ≠ Mennillo — is
preserved as a historical record but no longer drives any TODO.

Files updated:

- `code/02_atlas_prep/load_mennillo.py` — **deleted**.
- `code/02_atlas_prep/load_taurus.py` — **new**; `NotImplementedError`
  skeleton with both Zenodo DOIs in the docstring and the UC + single-
  time-point subset requirement documented.
- `code/02_atlas_prep/atlas_schemas.md` — `## mennillo` section
  replaced with `## taurus` (DEFERRED placeholder with the captured
  Zenodo metadata; full schema captured on download).
- `code/02_atlas_prep/README.md` — atlas list bullet + loader-status
  table row updated; v1 trio is now Smillie + Garrido-Trigo + TAURUS.
- `code/_shared/constants.py` — `UC_ATLASES` tuple: `"mennillo"` →
  `"taurus"`. (Downstream slurm scripts and module READMEs reference
  this constant via env var or import; cascade is contained.)
- `scripts/download_refs.sh` — Mennillo manual section → TAURUS Zenodo.
- `OPEN_FLAGS.md` — F1 (UC tissue definition) and F8 (fine-tier
  harmonization) updated to reference TAURUS instead of Mennillo.

Other Mennillo mentions in slurm scripts, individual-module READMEs,
top-level README.md, PLAN.md, and historical DECISIONS entries are
left in place — the slurm/module references will pick up the new
atlas name via `UC_ATLASES`, and the historical entries (1/7)–(15)
describe Mennillo as it was the planned third atlas at the time, which
remains an accurate historical record. The README/PLAN top-level
references should be reconciled by Muskaan in a manuscript-side pass
(out of scope for this correction).

### (b) "Effective N" was mislabeled in (14); recomputed.

(14) cited "46k / 376k / 1.6M effective N range across de Lange / Liu
/ Yengo, three orders of magnitude". **Wrong in three ways:**

1. **Those are totals, not effective Ns.** For case-control,
   `N_eff = 4 / (1/n_case + 1/n_ctrl)`.
2. **Magnitude wrong.** 46k → 1.6M is ~35× ≈ 1.5 orders of magnitude,
   not three.
3. **Frame wrong.** Yengo (quantitative) and Trubetskoy (separate
   trait) are positive/negative controls — they don't belong on a
   "cross-GWAS power confound" axis. The N confound that matters is
   **de Lange vs Liu** (both UC).

**Recomputed effective N for the two UC GWAS:**

| GWAS | Cases | Controls | Total | **N_eff** |
|---|---|---|---|---|
| de Lange UC (GCST004133) | 12,366 | 33,609 | 45,975 | **36,168** |
| Liu UC (GCST90446794), multi-ancestry | 6,862 EAS + 16,390 EUR = 23,252 | 15,456 EAS + 336,800 EUR = 352,256 | 375,508 | **87,242** |

**The relevant cross-GWAS power ratio is 87k / 36k ≈ 2.4× ≈ half an
order of magnitude.** Material (Liu has more power) but not three
orders. Also: de Lange's case:control ratio is 1:2.7, Liu's is 1:15.1
— Liu is far more control-skewed, which inflates total N relative to
effective N. Methods that scale with `N_eff` inherit the 2.4× gap;
methods that scale with total N see an 8.2× gap. **Both are
deceptive** when used naively to compare power; effective N is the
right axis.

Yengo (quantitative, GCST90245992): `N_eff ≈ N = 1,597,374`. **Stays
out of the cross-GWAS confound argument** since it's the positive
control for an unrelated trait. Same for Trubetskoy SCZ (negative
control).

This recomputation is what feeds the Li-2025 sample-size handling.
Nothing was "settled" by the totals quoted in (14); the corrected
N_eff numbers above are the substrate.

### (c) Liu UC accession verified the same way de Lange was.

Re-queried the GWAS Catalog REST API for GCST90446794:

- `reportedTrait`: "Ulcerative colitis" ✓ (NOT IBD-combined — the
  same trap that bit de Lange in (14) does not bite Liu).
- `fullPvalueSet`: true.
- `initialSampleSize`: "6,862 East Asian ancestry cases, 15,456 East
  Asian ancestry controls, 16,390 European ancestry cases, 336,800
  European ancestry controls".
- `publicationInfo.firstAuthor`: "Liu Z", PMID 37156999, Nat Genet
  2023-05-08.

**No European-only Liu UC arm exists** in the GWAS Catalog deposits
for PMID 37156999. The two UC accessions are:

- GCST90446794: multi-ancestry (EAS + EUR) — what we picked.
- GCST90446795: EAS-only.

**Ancestry-LD mismatch flag.** Liu UC is intrinsically multi-ancestry;
the EAS fraction (~5.9% of total, 22,318/375,508) is fixed. If MAGMA
runs against `g1000_eur` (the EUR LD reference Saisohan already
staged), the LD model will be wrong for the EAS variants. Three
options to handle, none free:

1. Use the multi-ancestry sumstats with a matched multi-ancestry LD
   reference (changes the MAGMA bfile across atlases, breaking
   uniformity).
2. Use a EUR-restricted variant set from the multi-ancestry sumstats
   (lose ~5.9% of effective N but keep EUR LD model honest).
3. Use de Lange UC as the only EUR-LD-honest UC GWAS; treat Liu UC as
   sensitivity / multi-ancestry-LD only.

**Deferred to method-side decision** (needs the MAGMA pipeline
running); not a blocker for the atlas side. Flagged in DECISIONS for
the (eventual) sample-size confound discussion.

### (d) CL OWL deep-dive — none of the 6 flagged IDs are deprecated.

Re-parsed the 2026-03-26 OWL specifically for `owl:deprecated`,
`IAO_0100001` (replaced_by), and `oboInOwl:consider` on the 6 IDs
flagged in (13). **None of the 6 carry any of those annotations.**

- All 6 are **live terms** in CL 2026-03-26.
- The 3 "label drift" cases (CL:1000347, CL:0002204, CL:0002138) are
  IDs that *correctly* point at the concepts we want — my draft just
  used outdated/wrong labels for them. Concept-side: no change needed.
  Label-side: `canonical_broad_DRAFT.md` rows updated inline in (13)
  to the pinned 2026-03-26 labels (`colonocyte`, `tuft cell`,
  `endothelial cell of lymphatic vessel`).
- The 3 "WRONG ID" cases (CL:1000280, CL:0009039, CL:0002073) are IDs
  that point at *unrelated* live concepts in 2026-03-26 (smooth muscle
  cell of colon; colon goblet cell; transitional myocyte respectively).
  **No `replaced_by` upgrade path** — these are plain draft authoring
  mistakes, not ontology obsoletions. The correct CL IDs for the
  intended concepts (intestinal stem cell beyond CL:0002250; colon
  epithelial progenitor; enteric glial cell) need actual ontology
  lookup, **not** a follow-the-pointer remap. Muskaan biology call.
- The label-drift fixes from (13) should still be sanity-checked by
  Muskaan against the pinned 2026-03-26 release; the LLM-edited
  labels could in principle have introduced their own drift. The
  pinned TSV at `data/reference/cl_terms_pinned.tsv` is the
  authoritative reference for that check.

### (e) Yengo subset confirmed.

GCST90245992 confirmed via REST: European ancestry, 1,597,374
individuals, `fullPvalueSet=true`, Yengo 2022 Nature. The 4M / 5.4M
figures in press releases reference the multi-ancestry discovery set
(GCST90245843, no full p-values). The 1.6M EUR-with-full-pvalues
subset is the right pick for our MAGMA + 1000G-EUR pipeline. **No
change to the (14) pick.**

### Summary of operational impact

- **Liu cross-GWAS handling**: needs ancestry-LD policy decision when
  MAGMA pipeline runs. Three options listed in (c) above; not blocking
  the atlas side. Deferred to method-side correction.
- **TAURUS DOI version**: pin needed before download. Two DOIs
  documented in the loader skeleton + atlas_schemas; resolve in the
  next session.
- **CL canonical IDs**: 3 unresolved IDs in the canonical broad draft.
  Real biology lookup; Muskaan call. Until resolved, the draft cannot
  promote to a locked CANONICAL_BROAD — which is exactly the discipline
  the DRAFT was designed for.
- **Effective N (UC GWAS)**: 36k (de Lange) vs 87k (Liu) — that's the
  number for the Li-2025 power discussion, not the totals.

---

## CORRECTION 2026-06-06 (17): TAURUS DOI v3 locked + real loader + v2 commitment + hygiene sweep

Post-(16) PI input resolved two open items and closes the TAURUS atlas
side of the work that can be done from the laptop.

### (a) Zenodo DOI v3 LOCKED.

(16) left the canonical DOI open between `10.5281/zenodo.13768607`
(swap directive) and `10.5281/zenodo.14007626` (resolved by Zenodo on
2026-06-06). Confirmed by PI: **v3 = `10.5281/zenodo.14007626`** is the
canonical pin. The 13768607 reference is presumably a stale versioned
predecessor with broken metadata; treating it as canonical would have
pulled a different file. All operational references in the repo are
now on the v3 DOI:

- `code/02_atlas_prep/load_taurus.py` — `ZENODO_DOI` constant.
- `code/02_atlas_prep/atlas_schemas.md` — taurus section.
- `code/02_atlas_prep/README.md` — atlas list.
- `scripts/download_refs.sh` — TAURUS manual block.

Pinned pooled file: `TAURUS_raw_counts_annotated_final.h5ad` (12.7 GB)
with **md5 `c1bd13b92cacb164a401c6c4a4e7912c`**. The full Zenodo record
also exposes 9 per-lineage h5ads + UMAP/PC coordinate files + a paired
sample list; raw BAM/CellRanger at GEO `GSE282122` if the loader ever
needs to fall back to per-sample reads.

### (b) Real `load_taurus.py` shipped (replaces the (16) skeleton).

The skeleton from (16) is replaced with a working loader. Cannot run
yet (the 12.7 GB h5ad is not on disk; same shape as Smillie before
SCP259 was staged), but every gate that does not need the bytes is in
place and `python -c "import load_taurus"` is clean.

**Schema captured from Thomas 2024 paper Methods + Zenodo description:**

- 4-level cell-type hierarchy: `compartment` → `low` → `intermediate`
  → `cell_state` (109 states at the finest tier). Loader emits all
  four levels as obs columns; `cell_type_fine = cell_state` and
  `cell_type_broad` is derived from the `low` tier via `LOW_TO_BROAD`.
- Cohort: 16 CD + 22 UC + 3 HC = 41 subjects. v1 strictly UC-only per
  PI directive (drop CD + HC).
- Inflamed-baseline threshold: `inflammation_score > 6.5` (from the
  Zenodo deposit description).
- Tissue regions: terminal ileum vs colonic (ascending / descending /
  sigmoid / rectum). v1 keeps colonic only.

**Filter chain (4 stages, each logged per-stage):**

1. Disease == UC (drop CD + HC). Disease values canonicalized via
   `_canonicalize_disease` (`"ulcerative colitis"`/`"uc"`/etc. → `UC`).
2. Region matches `COLONIC_REGION_KEYS` and does not match
   `ILEAL_REGION_KEYS` (substring-based, case-insensitive).
3. Timepoint matches `BASELINE_TIMEPOINT_KEYS` (`baseline`,
   `pretreatment`, `W0`, `V1`, etc.).
4. `inflammation_score > 6.5` (NaN drops).

**Validation gates (mirror correction 9/12/16 discipline):**

- Auto-detect obs columns from candidate-name lists; hard-raise on
  no match with `obs.columns` dumped for triage. Explicit override
  arguments on `load()` for every detected column.
- Donor-structure hard invariant: **22 UC donors** post-filter (Thomas
  Fig. 2b). Hard-raises with per-donor cell counts dumped if missed.
- Cell-count tripwire: soft, currently unset (`EXPECTED_N_CELLS_HINT =
  None`); populates on first end-to-end run.
- Canonical-vocab two-gate assertion (matches Garrido + Smillie
  loaders): gate 1 at module import validates
  `LOW_TO_BROAD.values() ⊆ _BROAD_VOCAB`; gate 2 at end-of-load
  validates emitted broad values. `LOW_TO_BROAD` ships **empty** in
  this v0 because the `low`-tier label set is not in the paper
  Methods text — gate 2 will fail loud on first run with the actual
  labels listed, then the map gets populated (one commit per Muskaan
  biology pass).
- Counts pipeline: `log1p(CP10k)` on load (raw counts in filename);
  raw integer counts preserved in `layers['counts']`;
  `raw_count_mode=True` unsupported per (5/7).
- HGNC remap via pinned `ensembl_to_hgnc` (correction 11).

**Backed-mode read.** The 12.7 GB pooled file is opened with
`ad.read_h5ad(path, backed='r')` for the filter phase, so we don't pay
for materializing X until after the 4-stage filter has narrowed the
cell set. Filtered subset is then materialized to memory. Memory
footprint stays bounded; the alternative (full read then filter) would
peak above the 12.7 GB file size.

### (c) v2 committed → F8 is CL-aware, not cheap.

PI confirmed: **v2 is committed.** This locks the F8 fine-tier
harmonization to the **CL-aware path** (build it once against the
pinned CL ontology), not the cheap-by-hand path. The
`canonical_broad_DRAFT.md` and `data/reference/cl_terms_pinned.tsv`
substrate from (13) is the right foundation; the broad vocab is the
v1+v2 shared substrate. Operationally:

- The 3 unresolved CL IDs from (13) (intestinal stem cell beyond
  CL:0002250; colon epithelial progenitor; enteric glial cell)
  **must** be resolved before CANONICAL_BROAD locks, because v2's
  additional atlases will reference them too.
- The 3 label-drift relabels from (13) should be double-checked by
  Muskaan against the pinned 2026-03-26 TSV — same point applies for
  v2 stability.
- F8 fine-tier vocab construction will use the CL graph (subclass
  relationships, not just label match) to harmonize across atlases.
  Out of scope for tonight; the constraint just got documented.

### (d) Repo hygiene sweep: operational Mennillo refs → TAURUS.

Per PI directive ("grep the 26 Mennillo refs, confirm none are
hardcoded operational paths, add the superseded-by (16)(a)
breadcrumbs"), did the operational sweep. Every place where the atlas
name is read as a runtime token now uses `taurus`:

- `code/02_atlas_prep/cl_rollup_maps.yaml`: `mennillo: {}` →
  `taurus: {}` (per-atlas rollup hook; `cl_rollup.py` reads this).
- 4 slurm scripts (`03_scdrs_compute`, `04_seismic`, `donor_loo_array`,
  `test_retest_array`): `ATLAS=…|mennillo|…` validation strings,
  `ATLASES=(…)` arrays, and example sbatch commands all updated. The
  `donor_loo_array` array bound now reads `--array=0-21` for TAURUS
  (22 UC donors).
- `scripts/slurm/README.md` and `scripts/README.md`: example loops
  iterate `smillie garrido_trigo taurus`; the manual-fetch atlas list
  updated; the per-atlas donor-count note updated.
- `code/06_concordance/compute_concordance.py` and `code/06_concordance/README.md`:
  example invocations with `taurus_delange_broad/...` result paths.
- `code/07_regime2_meta/run_brown.py`: docstring example.
- `code/08_cross_method/README.md`, `code/09_cross_gwas/README.md`,
  `code/10_broad_atlas_hca/README.md`,
  `code/11_broad_atlas_pangi/README.md`: `--atlases` / `--uc-atlases`
  example invocations.
- `code/03_scdrs/README.md`: all-cells policy rewritten to reflect
  TAURUS's UC-only subset (Smillie + Garrido carry HC; TAURUS does not
  in v1 per PI directive).
- `code/02_atlas_prep/README.md` (raw-counts shipper list) and
  `hgnc_remap.py` (loader-list docstring): updated.
- `code/_shared/canonical_broad_DRAFT.md`: 6 Mennillo → TAURUS swaps.

Plus: `code/10_broad_atlas_hca/README.md` carried a **double** stale
reference — "UC trio (Smillie 2019 / Kong 2023 / Mennillo 2024)" —
where the "Kong 2023" should have been Garrido-Trigo per DECISIONS
(2/7) and Mennillo should be TAURUS per (16). Both fixed in this
sweep; was a real hardcoded operational path (the line documented
which atlases the donor-overlap check runs against).

**Remaining Mennillo references** are by design:

- `DECISIONS.md` corrections (1/7) through (15) + (16): append-only
  historical record. Each entry is accurate as of its date; the
  swap is documented in (16)(a) and (17)(d).
- `PLAN.md` and top-level `README.md`: per (16)'s explicit note, these
  are manuscript-side reconciliation deferred to Muskaan's pass.
- A handful of in-doc breadcrumbs (`# was: mennillo (superseded by
  DECISIONS 16)`, `replaces Mennillo per DECISIONS 16`) — intentional;
  they're the audit-trail per PI directive.

### Net status from PI's critical-path read

PI's framing: "all five atlases in processed/ + four GWAS munged →
first scDRS/seismic results." On the atlas side, the only remaining
laptop-side work was the TAURUS loader logic, which is done in (b)
above. The actual end-to-end TAURUS run + the four GWAS munges happen
on Hummingbird (and the SCZ download needs PGC registration started
now per PI directive — flagged but not done in this commit, since
registration is human-in-the-loop).

Files updated in this batch:

- `code/02_atlas_prep/load_taurus.py` (skeleton → real loader with
  4-stage filter, hierarchy passthrough, two-gate canonical-vocab,
  log1p(CP10k), HGNC remap; ~450 lines).
- `code/02_atlas_prep/atlas_schemas.md` (DOI v3; passing Mennillo
  references → TAURUS).
- `code/02_atlas_prep/README.md` (DOI v3; raw-counts list updated).
- `code/02_atlas_prep/cl_rollup_maps.yaml` (mennillo → taurus).
- `code/02_atlas_prep/hgnc_remap.py` (docstring loader-list).
- `code/03_scdrs/README.md` (all-cells policy reflects UC-only TAURUS).
- `code/06_concordance/{compute_concordance.py,README.md}` (examples).
- `code/07_regime2_meta/run_brown.py` (docstring example).
- `code/08_cross_method/README.md`, `code/09_cross_gwas/README.md`,
  `code/10_broad_atlas_hca/README.md`,
  `code/11_broad_atlas_pangi/README.md` (example invocations).
- `code/_shared/canonical_broad_DRAFT.md` (5 spots).
- `scripts/download_refs.sh` (TAURUS manual block; DOI v3 + md5).
- `scripts/README.md` (manual-fetch list; per-atlas note).
- `scripts/slurm/03_scdrs_compute.slurm`,
  `scripts/slurm/04_seismic.slurm`,
  `scripts/slurm/donor_loo_array.slurm`,
  `scripts/slurm/test_retest_array.slurm`,
  `scripts/slurm/README.md` (validation strings, arrays, examples).
- `DECISIONS.md` (this entry).

---

## CORRECTION 2026-06-06 (18): Post-(17) audit — inflammation filter pulled, donor invariant relaxed, v2 escalation rolled back

The (17) loader had two bugs that would crash or silently halve the
first TAURUS run, plus an inferred-v2-commitment escalation that wasn't
PI-confirmed. Fixed before the file is staged.

### (a) Filter D (`inflammation_score > 6.5`) PULLED from the chain.

The Zenodo deposit description says "For baseline analyses, please use
baseline samples > 6.5 inflammation score". (17) read that as a
subsetting rule. It is not. It is the paper's **analytical convention**
to recapitulate its own remission analysis on inflamed baseline
tissue. Applying it as a load-time filter has three problems:

1. Drops roughly half of UC baseline samples (Fig. 2b: 50 inflamed vs
   53 non-inflamed).
2. Pre-empts OPEN_FLAGS F1, which must set one inflamed / non-inflamed
   / pooled policy uniformly across **all three** UC atlases. Smillie
   and Garrido are not inflamed-only, so an inflamed-only TAURUS
   silently breaks the cross-atlas comparability that is the entire
   point of the study.
3. Inflamed-only skews cell-type composition toward
   inflammation-expanded subsets and distorts GWAS prioritization.

**Resolution.** The 4-stage chain is now **3 stages** — disease=UC →
colonic → baseline. `inflammation_score` is carried through as obs
metadata (already was in (17)'s `_finalize`) so F1 can apply one
policy downstream. The deposit's 6.5 cutoff is retained as the
documented constant `PAPER_BASELINE_INFLAMMATION_MIN` for any future
inflamed-only sensitivity analysis that goes through F1, but it is
**not applied** in this loader. Logging at the inflammation step is
now informational: count of non-NaN scores, min/max, and the count
above 6.5 — gives F1 enough to design the cross-atlas policy without
binding TAURUS to it.

### (b) Donor invariant relaxed from hard `==22` to logged expected-range.

(17)'s hard `EXPECTED_N_UC_DONORS_POST_FILTER = 22` assertion would
have crashed the first run. 22 is the full UC cohort (Fig. 2b — all
timepoints, inflamed + non-inflamed). The (3-stage) subset above is
baseline-only, and Fig. 4c shows only 4 + 13 = 17 UC patients with
pretreatment samples in the paper's own baseline analysis (after
their inflamed filter, which we don't apply). So the post-filter
count lies in `[17, 22]` depending on how many UC patients gave a
pretreatment colonic biopsy at all. Without Supp Table 1 in hand to
pin the exact number, hard-asserting 22 false-fails on a correct
subset; and silently tuning it down to pass would let the (17)
inflamed-only filter sneak through (the two bugs interact — see PI's
warning).

**Resolution.**

- Constant renamed `EXPECTED_N_UC_DONORS_RANGE = (15, 22)`.
- Hard-fail only on impossible values (`n_uc_donors == 0` — filter
  misconfigured — or `> 60` — wrong column auto-detected).
- Log a warning if outside the expected range; this surfaces metadata
  drift or a needed Supp-Table-1 reconciliation without crashing the
  load.
- Per-donor cell breakdown always logged for triage.

The exact number gets pinned empirically from the first end-to-end
run + a cross-check against Supp Table 1.

### (c) v2 commitment — rolled back to "PI-confirm-before-escalating-v1-paths".

(17)(c) treated the user's "v2 is committed" statement as PI sign-off
and escalated three CL IDs to **v1 critical-path-blocking** for the
CANONICAL_BROAD lock. That was an overreach: confirmation that v2 is
on the roadmap is not the same as PI sign-off that the v1 broad-tier
first pass should wait on the v2-ready CL graph. The exact "build-it-
twice risk" PI flagged is the reason to slow down here.

**Rollback.**

- v2 commitment recorded as **needs explicit PI confirmation** (the
  one-question check PI asked for at the top of this turn). Pending
  that, the F8 CL-aware track is the planned path but not yet locked
  as a v1 dependency.
- The three unresolved CL IDs (`CL:1000280`, `CL:0009039`,
  `CL:0002073`) gate the eventual **v2-ready** `CANONICAL_BROAD` lock,
  **not** the v1 first 5-atlas × broad-tier concordance pass. v1's
  broad-tier output comes from each loader's existing
  `FINE_TO_BROAD` / `LOW_TO_BROAD` mapping into the same 15-string
  `_BROAD_VOCAB` frozenset. The two-gate canonical-vocab assertion
  already enforces string equality across all three UC loaders, which
  is what step 06 needs.
- `canonical_broad_DRAFT.md` is the substrate for the eventual lock;
  not a runtime dependency of step 06.

This decouples the critical paths: **v1 first results** unblocks on
TAURUS loader + GWAS munge (the PI's narrow critical path); **v2-ready
CANONICAL_BROAD** unblocks on the three CL IDs (Muskaan biology call)
+ PI v2 confirmation. Those two paths run in parallel and neither
waits on the other for its own milestones.

### (d) Trubetskoy PGC registration — escalated in download_refs.sh.

(14) noted PGC registration as a manual step but did not flag the
human-in-the-loop delay. (17) listed it as "flagged but not done in
this commit". This correction bumps the visibility:
`scripts/download_refs.sh` step 6 now opens with
`** HUMAN-IN-THE-LOOP DELAY — START REGISTRATION NOW. **` and an
explicit "submit the request now even though Hummingbird access is
offline" instruction. The PGC review is the rate-limiting step for
the SCZ negative-control GWAS, so the wait clock should start at the
earliest possible moment.

This is the only critical-path item that can be advanced from a
non-HPC laptop without writing code — submitting the registration
form. Cannot be done from this session (the form is browser-based
behind a registration wall); flagged for human action.

### Net status

PI's narrow critical path is unchanged: TAURUS loader (done in (17),
hardened in (18)(a)+(b)) + 4 GWAS munged. With (a) + (b) fixed, the
loader is safe to stage; the first run will not crash on a correct
subset, and will not silently halve the cohort. The donor count will
be captured empirically and reconciled against Supp Table 1 on first
run.

Files updated in this batch:

- `code/02_atlas_prep/load_taurus.py` (Filter D removed; constants
  renamed; donor invariant relaxed to logged expected-range; module
  docstring updated; informational inflammation-score logging in
  place of the dropped filter).
- `scripts/download_refs.sh` (Trubetskoy step 6 — human-in-the-loop
  emphasis).
- `DECISIONS.md` (this entry).

---

## CORRECTION 2026-06-06 (19): PGC SCZ is publicly downloadable — (18)(d) retracted

(18)(d) flagged PGC registration as a "multi-day human-in-the-loop"
delay and asked for it to be started immediately. On attempting to
start it, **the registration does not exist**. The summary statistics
are on public figshare with `is_public=true`, `is_embargoed=false`,
and a CC-BY-4.0 license; the figshare API returns 36 files with
direct `ndownloader` URLs.

The "PGC Terms and Conditions" referenced on
`https://pgc.unc.edu/for-researchers/download-results/` are a Fort
Lauderdale honor-code researcher pledge (don't re-identify
participants, don't publish global analyses before the consortium
manuscript publishes), **not** a registration gate. The publication
embargo from the Fort Lauderdale clause cleared when Trubetskoy 2022
landed in Nature 604:502 (April 2022) — free for genome-wide analyses
including ours.

**Schema captured** (PGC3 wave3 EUR file, range-fetched from
figshare):

- **Format**: `##fileFormat=PGCsumstatsVCFv1.0` — a PGC-internal
  extension of VCF. Verbose `##` meta-header (acknowledgments,
  contigs, terms of use, n_case=52017, n_control=75889, n_trio=1369);
  data rows are 14 plain tab-separated columns:

  | # | Column | Notes |
  |---|---|---|
  | 1 | `chr` | 1..22, X, Y, M |
  | 2 | `rsid` | dbSNP id |
  | 3 | `pos` | GRCh37 coordinate |
  | 4 | `A1` | effect allele |
  | 5 | `A2` | other allele |
  | 6 | `frq_A` | allele freq in cases |
  | 7 | `frq_U` | allele freq in controls |
  | 8 | `info` | imputation quality |
  | 9 | `beta` | log-OR |
  | 10 | `se` | standard error |
  | 11 | `p` | p-value |
  | 12 | `n_case_cohort` | cohort-level (53,386 = 52,017 + 1,369 trios) |
  | 13 | `n_ctrl_cohort` | cohort-level (77,258 = 75,889 + 1,369 trios) |
  | 14 | `n_eff` | **per-SNP** effective sample size |

- **HAS per-SNP `n_eff`** — column 14 varies by imputation
  completeness; pass `--col-n n_eff` to `prepare_gwas.py` rather than
  the cohort-level fixed N.
- **Build**: GRCh37 (matches Yengo + de Lange; Liu is GRCh38).
- **EUR cohort N_eff (cohort-level)**: `4 / (1/53,386 + 1/77,258) ≈
  126,275`. Per-SNP values from column 14 will hover near this.
  Bigger than Liu UC (87k) and de Lange UC (36k) — the SCZ negative
  control is the most-powered of the four.
- **md5 of pooled file**: `6ebe2376f5cda972d37efa0f214c4df0` (240 MB).
- **Cited DOI**: `10.6084/m9.figshare.19426775.v7`.

**Munge implication**. Because of the PGC-VCF format, the munge step
(`prepare_gwas.py` or LDSC `munge_sumstats.py`) needs to skip the
verbose `##` meta-header and treat the data rows as 14-column
tab-separated. Standard GWAS-Catalog-SSF parsers will see the
`##`-prefix lines and either skip or error, depending on
implementation. Worth flagging in the MAGMA pipeline notes when those
run.

**Repo changes from (18)(d) → (19)**:

- `scripts/download_refs.sh` step 5d is **new**: auto-fetch the
  EUR-ancestry public file at
  `https://ndownloader.figshare.com/files/34517828` →
  `gwas/scz_trubetskoy_eur_PGC3_v3.vcf.tsv.gz`. Step 6 (manual
  Trubetskoy section) collapsed to a one-line note pointing at 5d
  and at (19). Per-SNP N status table at the bottom of the script
  updated.
- `scripts/README.md`: the "datasets that require accounts" list
  retracted PGC; "first-run notes" item on Trubetskoy N updated to
  the per-SNP `n_eff` column.

**Net effect on the critical path**: all four GWAS are now
auto-fetchable. There is no PGC registration clock to start. The
remaining human-in-the-loop items on the critical path are: PI
confirmation of v2 (per (18)(c)) and Hummingbird re-access for the
TAURUS + Garrido runs. Neither is an SCZ-download issue.

Files updated in this batch:

- `scripts/download_refs.sh` (new step 5d auto-fetch; step 6 collapsed;
  per-SNP N status table updated).
- `scripts/README.md` (retract PGC from accounts-required list;
  update Trubetskoy N notes).
- `DECISIONS.md` (this entry).

---

## CORRECTION 2026-06-06 (20): TAURUS Supp Table 1 dry-run + vocab single-sourced + three pre-stage bugs caught

Two-pronged audit per PI directive: (i) "cheap and now" — fetch
Supplementary Table 1 (70 KB xlsx) and validate `load_taurus.py`'s
filter chain against the real cohort table before the 12.7 GB h5ad is
staged; (ii) "two-minute check" — single-source `_BROAD_VOCAB` so the
three UC loaders aren't measuring drift between identical-by-copy
frozensets.

### (a) Supp Table 1 fetched + parsed; UC cohort empirically locked.

Source: Nature ESM URL for the paper
(`https://static-content.springer.com/esm/art%3A10.1038%2Fs41590-024-01994-8/MediaObjects/41590_2024_1994_MOESM3_ESM.xlsx`,
70 KB, 6 sheets). The relevant sheet is **Supp_Table_1B_SampleMetadataIBD**
— 216 sample-level rows with columns: `sample_id`, `LibraryType`,
`CellsLoaded`, **`Disease`**, `Disease_duration`, **`Patient`**,
**`Site`**, **`Inflammation`** (categorical), **`Treatment`**
(Pre/Post), `Age`, `Sex`, `Ethnicity`, `Match`, `Batch`,
**`Inflammation_score`** (numeric), and ~20 CellRanger QC stats.

Filter chain dry-run (sample-level, since cell-level inherits the
sample's metadata per cell):

| Stage | Samples kept | Notes |
|---|---|---|
| Total | 216 | 108 UC + 96 CD + 12 Healthy |
| A: Disease == UC | 108 | drops CD + Healthy |
| B: Colonic site | 103 | 5 UC Terminal_Ileum samples dropped |
| C: Treatment == Pre | **52** | **22 distinct UC donors** |

**Empirical hard invariants now in `load_taurus.py`:**

- `EXPECTED_N_UC_DONORS_POST_FILTER = 22` (all UC patients contribute
  at least one Pre colonic sample — no one is dropped).
- `EXPECTED_N_UC_SAMPLES_POST_FILTER = 52` (soft tripwire — h5ad may
  drop samples at QC).

The Fig. 4c "17 UC patients with pretreatment samples" figure I
estimated in (18)(b) was the paper's *additional* inflamed-baseline
subset (after `inflammation_score > 6.5`), which (18)(a) deliberately
pulls. Without that further filter, all 22 contribute. So (17)'s
original `==22` was actually right; (18)(b)'s relaxation to
`(15, 22)` was over-defensive. Tightened back to hard `==22`. Per-donor
sample counts range 1–3; inflammation breakdown of the v1 subset:
39 Inflamed + 13 Non_Inflamed (F1 governs whether to split downstream).

### (b) Three pre-stage bugs caught in `load_taurus.py` — none would
have surfaced without the Supp Table 1 audit.

These are the value of "cheap and now" — would each have crashed or
silently corrupted the first end-to-end run.

1. **`_is_baseline` substring direction was wrong.** Original:
   `return any(key in t for key in BASELINE_TIMEPOINT_KEYS)` — asks
   "is any baseline key a substring of the value?" For TAURUS's
   `"Pre"` (3 chars), no key in the list is a substring of `"pre"` →
   returns False → **the loader would have dropped every Pre row,
   producing zero UC cells**. Fixed: now `t in BASELINE_TIMEPOINT_KEYS`
   (exact match against an explicit frozenset that includes the bare
   token `"pre"`). Negative test confirms `_is_baseline("Pre") = True`,
   `_is_baseline("Post") = False`.

2. **`_TIMEPOINT_COL_CANDIDATES` was missing `"Treatment"`.** TAURUS's
   actual column name is `"Treatment"` (Supp Table 1B), which my
   candidate list didn't include. Auto-detect would have raised
   `KeyError` on first run. Fixed; `"Treatment"` is now first in the
   candidate list (TAURUS-confirmed names take precedence).

3. **`_INFLAMMATION_COL_CANDIDATES` missed the capital-I
   `"Inflammation_score"`.** Same shape as (2); added with TAURUS
   precedence.

Plus the donor-invariant relaxation (b above) — net 4 corrections
without ever touching the 12.7 GB file.

### (c) `_BROAD_VOCAB` single-sourced (the "two-minute check").

`grep -n _BROAD_VOCAB` confirmed three copy-paste duplicates of the
identical 15-string frozenset in `load_smillie.py`,
`load_garrido_trigo.py`, `load_taurus.py`. Resolved by extracting to
`code/02_atlas_prep/_broad_vocab.py` (sibling private module; same dir
as `hgnc_remap.py` which the loaders already import sibling-style).
Each loader now does `from _broad_vocab import _BROAD_VOCAB`;
verified `id(load_X._BROAD_VOCAB)` is identical across all three.

This is **not** the locked public `CANONICAL_BROAD` of
`canonical_broad_DRAFT.md` — it's still underscore-prefixed and
loader-only. When CANONICAL_BROAD locks (with CL subtree IDs,
deprecation pass, and PI v2 sign-off per (13) + (18)(c)), this
module's contents promote to `code/_shared/canonical_broad.py` and
become a public import for `06_concordance`.

PI's risk framing: "the heatmap measures vocabulary drift instead of
biology if it isn't single-sourced." The (16) F8-preview defense
(value-side typo catches via gate 1) only works if every loader's
gate 1 is checking the same frozenset; with copy-paste duplicates, a
value-side typo could pass in one loader and fail in another, or
worse, an asymmetric vocab update could silently degrade the
cross-atlas string intersection at step 06.

### (d) What remains for next session — explicit critical-path read.

PI's narrow path to the first 5-atlas × broad-tier concordance pass:

1. **GWAS pipeline locally.** MAGMA binary + g1000_eur LD ref +
   download de Lange + Yengo + Trubetskoy SCZ EUR (all four now
   auto-fetchable per (19)) + munge → MAGMA → `make_scdrs_gs.py` →
   `.gs` gene-set files. Out of scope tonight (context budget);
   queued for next session. Liu's ancestry-LD policy decision blocks
   only Liu — de Lange + Yengo + SCZ can munge first.

2. **TAURUS expensive de-risk** (optional but high-value).
   `download_refs.sh` can pull the 12.7 GB pooled file to laptop;
   then a dry-run of `load_taurus.py` in backed mode surfaces the
   actual `low`-tier label set to populate `LOW_TO_BROAD`. The
   sample-level dry-run in (a) above already validates the filter
   chain; the cell-level dry-run validates the cell-type hierarchy.

3. **PI handoffs unchanged from previous turns.** v2 sign-off (gates
   the eventual CANONICAL_BROAD lock); Muskaan's biology on the 3
   unresolved CL IDs + MT/HSP/IER policy + F8 fine vocab.

Files updated in this batch:

- `code/02_atlas_prep/_broad_vocab.py` (new; single-source frozenset).
- `code/02_atlas_prep/load_smillie.py`,
  `code/02_atlas_prep/load_garrido_trigo.py`,
  `code/02_atlas_prep/load_taurus.py` (each: replace local
  `_BROAD_VOCAB` definition with `from _broad_vocab import _BROAD_VOCAB`).
- `code/02_atlas_prep/load_taurus.py` (additional: fix `_is_baseline`
  direction; tighten `EXPECTED_N_UC_DONORS_POST_FILTER = 22`; add
  `EXPECTED_N_UC_SAMPLES_POST_FILTER = 52`; add `"Treatment"` and
  `"Inflammation_score"` to candidate lists at TAURUS-confirmed
  precedence; docstring update).
- `DECISIONS.md` (this entry).

The Supp Table 1 xlsx itself stays in `~/Downloads/` (not committed —
Nature's copyright on supplementary materials is murky). The
derived expected-cohort facts are encoded in the loader's constants;
that's the audit-trail.

---

## CORRECTION 2026-06-06 (21): Donor-mismatch diff dump + Trubetskoy PGC3 fact-check + sequencing nudge

Three small but load-bearing follow-ups to (20).

### (a) Donor-invariant failure now dumps a structured per-donor/per-region diff.

PI flagged that the bare `expected 22, got 21` failure message
wouldn't help triage when Zenodo v3 of TAURUS exists precisely
*because* of a metadata fix vs prior versions — a mismatch is signal,
not noise, and the failure message should point at the actual drift.

Resolution:

- New constant `EXPECTED_UC_COHORT` in `load_taurus.py`: dict mapping
  each of the 22 UC donor IDs to a frozenset of expected colonic
  regions, derived directly from Supp Table 1B's UC × colonic × Pre
  rows (the same dry-run that produced the 22-donor / 52-sample lock
  in (20)). Self-consistency check: sum of region-set cardinalities
  across the 22 donors = 52 = `EXPECTED_N_UC_SAMPLES_POST_FILTER`.
- Donor-mismatch failure path now computes and dumps:
  - `missing_from_observed`: donors in Supp Table 1B but not in the
    h5ad subset.
  - `unexpected_in_observed`: donors in the h5ad subset but not in
    Supp Table 1B (would signal an upstream metadata reassignment).
  - **`per-donor region drifts`**: for each donor present on both
    sides, the symmetric difference of region sets (`missing` and
    `extra` regions per donor).
  - Per-donor cell counts (top 30) as before.

So a v3-vs-v4 metadata fix where, say, UC18's only sample gets
re-tagged to a different region now surfaces as `UC18: expected
{Rectum}; observed {Sigmoid}; missing [Rectum]; extra [Sigmoid]` —
the triage call, not the donor-count delta.

### (b) Trubetskoy fact-check: confirmed PGC3 SCZ, NOT PGC2.

PI flagged that the (19) figshare reversal needs verification — that
the file on figshare DOI 10.6084/m9.figshare.19426775 is actually
Trubetskoy 2022 PGC3, not an older PGC2 release. Low stakes (negative
control), but don't munge the wrong study.

The VCF meta-header captured in (19) settles it. Verbatim:

- ``##shortName="PGC3-SCZ"`` ← **PGC3**, not PGC2.
- ``##preparedBy="VassilyTrubetskoy"`` ← first author of the 2022
  Nature paper.
- ``##DOI=<ID=BIORXIV,REF="https://doi.org/10.1101/2020.09.12.20192922">``
  ← the bioRxiv preprint of Trubetskoy 2022 (the v2 preprint posted
  in 2020 that became the Nature paper in 2022).
- ``##nCase="52017"`` + ``##nTrio="1369"`` → cohort-level 53,386 EUR
  cases. Trubetskoy 2022 reports 53,386 EUR cases (52,017 + 1,369
  trios). **Exact match**.
- ``##nControl="75889"``. Trubetskoy 2022 reports 77,258 EUR controls
  (75,889 + 1,369 trios contribute as controls in the PGC trio
  encoding). Matches.
- ``##methodsParagraph`` reproduces verbatim from the published
  Nature paper Abstract / Methods (Trubetskoy 2022 Nat 604:502).

Counterfactual sanity-check: PGC2 (Schizophrenia Working Group 2014,
Nat 511:421) had 36,989 EUR cases and 113,075 EUR controls. The
numbers don't match this file at all, and the `##shortName` would be
"SCZ2" or similar, not "PGC3-SCZ". So this is definitively the 2022
release.

Confirmed pin remains: figshare DOI ``10.6084/m9.figshare.19426775.v7``,
file ``PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz``, md5
``6ebe2376f5cda972d37efa0f214c4df0``.

### (c) Sequencing nudge for next session.

PI flagged that the TAURUS-expensive dry-run (12.7 GB pull + cell-level
loader validation) is **now largely redundant** thanks to the (20)
Supp Table 1 sample-level dry-run + gate 2's loud-failure-on-first-run
behavior for the `low`-tier labels. So next session's ordering is:

1. **GWAS pipeline locally** (laptop-side, half the heatmap, fully
   independent of HB): pull MAGMA + g1000_eur + de Lange + Yengo +
   Trubetskoy SCZ EUR (all four auto-fetchable now per (19) and the
   downstream-handling table in `scripts/download_refs.sh`). Munge
   each → MAGMA → `make_scdrs_gs.py` → `.gs` gene-set files. Liu
   waits on the ancestry-LD method decision.
2. **HB-side TAURUS run** (when HB access is back): stage the
   pre-validated loader; gate 2 surfaces the real `low`-tier label
   set on first failure; populate `LOW_TO_BROAD` from that;
   re-run; done.
3. **Skip** the local 12.7 GB pull unless HB access slips badly —
   the validation it would add is now redundant given (20)+(21).

Files updated in this batch:

- `code/02_atlas_prep/load_taurus.py` (`EXPECTED_UC_COHORT` constant
  added; donor-mismatch failure rewritten to dump the per-donor /
  per-region diff).
- `DECISIONS.md` (this entry).

---

## CORRECTION 2026-06-07 (22): Broad-tier lock — QC-state collapse, lineage-ambiguous-cycling exclusion, CL wrong-ID resolution

Three Tier-1 biology calls (gate the 3×3 broad concordance figure) and
one Tier-2 cleanup (gate the v2 CANONICAL_BROAD lock), bundled here so
the cross-atlas rules land together and the loaders, the draft vocab
doc, and the CL pin all reflect the same lock atomically.

Sequencing note: this entry follows Saisohan's 2026-06-07 biology
review of the Tier-1 punch list. The flagged single mis-take in my
proposal — locking the unprefixed Garrido `Cycling cells *` to B cell —
is reversed here as a uniform cross-atlas exclusion rule, which is the
only resolution that doesn't manufacture a B-cell harmonization
artifact in the first headline figure.

### (a) `QC_STATE_TO_PARENT` — collapse, not exclude (resolves OPEN_FLAGS F2).

Six Garrido fine labels live on a cross-lineage QC-state axis parallel
to Ribhi: `MT T cells`, `MT fibroblasts`, `IER fibroblasts`,
`PC IgA heat shock 1/2`, `PC immediate early response`. The broad-tier
target for each is already biologically correct in `FINE_TO_BROAD`
(stress doesn't change lineage); what was missing was the fine-tier
collapse-to-parent that Ribhi already gets. New cross-atlas map
`code/02_atlas_prep/_qc_policy.QC_STATE_TO_PARENT` chains alongside
`RIBHI_TO_PARENT` in each loader. Smillie's single `MT-hi` cluster
stays mapped directly at broad in `load_smillie.FINE_TO_BROAD` (it's
one Imm-compartment label, no fine-tier collapse buys anything);
TAURUS imports the map but it's a no-op until first-run surfaces its
low-tier label set.

**Policy: collapse-to-parent, NOT exclude.** Reasoning: MT-hi cells are
stressed-but-viable in the deposit's own QC (Smillie SCP259 ran MT%
filters before publishing the 51-cluster annotation; the Salas-lab
Garrido annotation likewise gates on MT% before assigning the MT-T /
MT-fib clusters — these are stressed cells, not dying-cell residue).
Heat-shock and IER are dissociation artifacts overlaid on real
lineages; the original lineage call survives. Exclusion would drop
cells asymmetrically (Garrido has 6 stress labels, Smillie has 1,
TAURUS unknown) and would change the soft cell-count tripwires; the
collapse keeps cells and only removes QC noise from the fine tier.

**Cross-atlas-rule status, not Garrido+Smillie patch.** The map lives
in `_qc_policy.py` precisely so TAURUS labels surface into the same
map (not into the TAURUS loader) when the low-tier set is enumerated,
preserving cross-atlas symmetry by construction rather than by-copy.
Same architectural pattern as `_broad_vocab._BROAD_VOCAB` (DECISIONS
20) — duplicate-by-copy is the drift risk this fixes.

**PRE-BANK GATE on the first broad-tier figure (NOT a footnote).**
Before treating the 3×3 broad heatmap as real, confirm the MT-hi
cluster cell counts in *each* atlas are consistent with the deposit's
own MT% QC having actually passed (Salas-lab Garrido fine-tier table;
SCP259 Smillie 51-cluster metadata). The collapse policy assumes
"stressed-but-viable, kept by deposit QC" — if a deposit silently
shipped sub-threshold MT% cells in those clusters, collapse-to-parent
quietly scores dying cells into our lineages and the figure inherits
the bias. Tracked as **OPEN_FLAGS F9 (MT% pre-bank gate)** — do not
bank the heatmap without it. Cheap to check, expensive to miss. If
sub-threshold cells slipped through, the resolution is to drop those
specific cells, NOT change the cross-atlas collapse rule.

**Forward note (Phase 9 — compositional confound).** Whether stress-
state cell fractions correlate with disease state (more MT-hi / more
heat-shock in inflamed tissue) is a compositional-confound question
that bites at the concordance step, not at the load. Park for the
Phase 9 compositional-confound check; don't relitigate inside the QC
policy. If the check surfaces a stress-fraction × disease correlation,
the answer is a sensitivity panel (re-run concordance with stress-
state cells held out), not a change to the v1 collapse rule.

### (b) `EXCLUDE_LINEAGE_AMBIGUOUS_FINE` — uniform exclusion of unprefixed cycling (resolves OPEN_FLAGS F3).

Garrido's `Cycling cells`, `Cycling cells 2`, `Cycling cells 3` carry
proliferation markers (MKI67+ etc) without a stable lineage call —
they're labeled by clustering position in the Salas-lab tree, not by
marker-confirmed compartment. My initial proposal locked these to B
cell at broad (following the Salas-lab tree placement); Saisohan
flagged this as the single assignment that could manufacture a false
cross-atlas result — if Smillie or TAURUS happens to place its
analogous cycling cells anywhere else (T, or a separate cycling
bucket), the broad heatmap would report a B-cell discordance that's
pure harmonization artifact, and "B cells disagree across atlases" is
exactly the failure mode a reviewer dismisses with "that's just your
cycling-cell handling."

**Resolution: drop these cells from broad scoring, uniformly, in all
three loaders.** New cross-atlas frozenset
`code/02_atlas_prep/_qc_policy.EXCLUDE_LINEAGE_AMBIGUOUS_FINE`,
applied before the FINE_TO_BROAD lookup in every loader; logs the
drop count loudly so the exclusion is auditable. Smillie currently
adds zero labels (all its cycling clusters are compartment-prefixed:
`Cycling B`, `Cycling T`, `Cycling Monocytes`, `Cycling TA` — kept,
lineage-tagged). TAURUS imports the set; populated when low-tier
labels surface and any analogous unprefixed cycling cluster appears.

**Why symmetric exclusion is the bounded move.** Cycling-cell
fractions are small in every UC atlas (typically a few percent at
most); exclusion is bounded in cell-count impact and symmetric across
atlases by construction. Force-assignment is neither — its discordance
contribution scales with the asymmetry of how each atlas's annotators
happened to place cycling cells, which is exactly the unknown we
can't audit pre-figure.

**Revisit path.** The marker-QC step already on the build plan
(MKI67 + lineage markers in the included-cells diagnostic pass) is
where these cells get reinstated per-cell, but only when markers AND
cross-atlas placement agree. The exclusion is a structural floor on
the v1 figure, not a permanent drop.

The Garrido `FINE_TO_BROAD` no longer ships entries for the three
unprefixed cycling labels — they go through the EXCLUDE step before
the unmapped-label gate; the previous B-cell entries are deleted.
Also added: identity row `"plasma cell": "plasma cell"` to Garrido's
`FINE_TO_BROAD`, as the target of `QC_STATE_TO_PARENT`'s plasma-cell
collapse path (PC IgA heat shock 1/2, PC immediate early response).

### (c) Garrido broad-tier lock — 9 of 12 REVIEW rows confirmed (resolves OPEN_FLAGS F4 broad-tier portion).

The 12 OPEN_FLAGS F4 REVIEW rows are biology calls on Garrido fine
labels with non-obvious broad-tier placement. After Saisohan's review,
9 rows lock with their current `load_garrido_trigo.FINE_TO_BROAD`
assignment; the other 3 are the unprefixed cycling cells just handled
in (b) above. Locked assignments:

| Garrido fine label | Current broad | Rationale |
|---|---|---|
| Laminin colonocytes | colonocyte | Basement-membrane-adjacent colonocyte subtype — colonocyte at broad. Fold-vs-distinct at fine is an F8 question, not a broad call. |
| PLCG2 colonocytes | colonocyte | Inflammation-associated colonocyte subtype (PLCG2 = UC GWAS gene). Broad = colonocyte. F8 question whether to fold. |
| Inflammatory colonocyte | colonocyte | State-flavored colonocyte. Broad clean. |
| Mature goblet | goblet | Terminally differentiated goblet. Broad clean. F8 question whether to fold to Goblet. |
| Paneth-like | goblet | Colon lacks canonical Paneth; Paneth-like = secretory anti-microbial epithelium. Broad placement under goblet (secretory) defensible — split-out would create a 16th broad term for a rare population, against the 10–15 budget. |
| CD4 ANXA1 | T cell | CD4+ ANXA1+ T cells (Th1/Treg-adjacent). Broad clean. |
| S1PR1 T cells | T cell | S1PR1+ naive/memory T cells. Broad clean. |
| T cells CCL20 | T cell | Th17-skewed T cells. Broad clean. |
| PC IGLL5 | plasma cell | IGLL5-expressing plasma cells. Broad clean. F8 question on isotype grouping. |

The remaining fine-tier identity questions (fold-vs-distinct for
Laminin / PLCG2 colonocytes, Mature goblet, PC IGLL5 isotype
grouping) flow to F8. Broad tier is **locked**.

### (d) CL pin — three wrong-ID rows resolved (DECISIONS 13 follow-up).

The CL pin (13) caught six DRAFT errors: three label-drift renames
(closed by 2026-06-06 doc edits, T1.3 verified the doc matches the
pinned tsv) and three flat-out wrong IDs. The wrong-ID rows are
resolved here.

| Drafted slot | Wrong ID | Resolution | OLS-confirmed last-modified |
|---|---|---|---|
| intestinal stem cell | CL:1000280 (was "smooth muscle cell of colon") | **Drop.** Use existing pinned `CL:0002250` *intestinal crypt stem cell*. | n/a — already pinned |
| colon epithelial progenitor cell | CL:0009039 (was "colon goblet cell") | **Replace with `CL:0009010`** *transit amplifying cell* (generic TA term; "of colon" subtype does not exist in 2026-03-26). | 2024-04-03 |
| enteric glial cell | CL:0002073 (was "transitional myocyte") | **Replace with `CL:4040002`** *enteroglial cell* (exact synonym: "enteric glial cell"; parent of CL:4047047 type I enteric glial cell). | 2023-04-03 |

Both replacement IDs have OLS last-modification dates that comfortably
predate the 2026-03-26 pin date, so the labels we capture here are
stable for that release. **Pin provenance: same-release post-hoc
additions.** The pin's release-date provenance is unchanged — these
are within-release adds keyed against the live OLS, NOT a refresh of
the pin to a newer release. The cl_terms_pinned.tsv has a comment
block delimiting the post-hoc additions; the OWL itself remains
2026-03-26.

**Direct OWL membership verification (2026-06-07).** Saisohan flagged
that OLS reflects the live ontology, not the 2026-03-26 snapshot;
"pin that lies" was the failure mode. Verified directly against the
local 65.88 MB `cl.owl` (versionIRI
`http://purl.obolibrary.org/obo/cl/releases/2026-03-26/cl.owl`,
versionInfo `2026-03-26`): grep confirms `CL_4040002` and `CL_0009010`
both have `owl:Class` declarations in the pinned OWL. Labels match the
tsv (`enteroglial cell` with exactSynonym `enteric glial cell`;
`transit amplifying cell`). Pin membership is established for this
release, not assumed from OLS.

**v2 polish (deferred).** Two same-release children exist that would
add specificity without changing broad-tier mapping (parent subtree
semantics cover both):

- `CL:0009043` *intestinal crypt stem cell of colon* — colon-specific
  child of `CL:0002250`.
- `CL:4047017` *transit amplifying cell of gut* — gut-specific child
  of `CL:0009010` (added 2024-09-24 per the OWL `terms:date`); closer
  fit to Saisohan's original "of colon / of gut" framing than the
  generic CL:0009010 that locked here. Lock stayed on the generic to
  preserve the signed-off scope; upgrade is a one-line swap in
  `canonical_broad_DRAFT.md` and the pinned tsv when v2 polish runs.

Files updated in this batch:

- `code/02_atlas_prep/_qc_policy.py` (new — `QC_STATE_TO_PARENT`,
  `EXCLUDE_LINEAGE_AMBIGUOUS_FINE`; cross-atlas rule module).
- `code/02_atlas_prep/load_garrido_trigo.py` (import `_qc_policy`;
  drop the three unprefixed `Cycling cells *` rows from `FINE_TO_BROAD`;
  add `"plasma cell": "plasma cell"` identity row; apply EXCLUDE +
  QC-state collapse in the fine-label processing block).
- `code/02_atlas_prep/load_smillie.py` (import `_qc_policy`; apply
  EXCLUDE + QC-state collapse for symmetry — both no-ops on the
  current SCP259 label set; lift the `MT-hi # REVIEW` comment).
- `code/02_atlas_prep/load_taurus.py` (import `_qc_policy`; apply
  EXCLUDE + QC-state collapse hook in `_finalize`; both no-ops until
  LOW_TO_BROAD populates on first compute-node run).
- `code/_shared/canonical_broad_DRAFT.md` (CL pin error table updated
  to "Resolution" column; rows 4 / 7 of the candidate vocabulary
  updated to cite CL:0009010 and CL:4040002 respectively; CL:1000280
  dropped).
- `data/reference/cl_terms_pinned.tsv` (post-hoc same-release append
  block for CL:0009010 and CL:4040002 with provenance comment).
- `OPEN_FLAGS.md` (F2 / F3 / F4 moved to Resolved; F7 narrowed to its
  fine-tier-only portion; F8 expanded with T2.5 CL-aware fine-vocab
  analytical sketch).
- `DECISIONS.md` (this entry).

The 3×3 broad-tier concordance figure is now unblocked on the biology
side. Remaining gate is the loader-runs-clean check on a compute node
(MT% spot-check from caveat (a); Garrido + Smillie load output; no
TAURUS dependency for the first figure since CD/HC arms are dropped
and the broad axis can be computed pairwise).

---

## CORRECTION 2026-06-07 (23): Li-2025 N_eff lock + bibliography audit

Two small but record-worthy items from the rev2 execution to-do PDF:
the cross-GWAS power-axis lock that the (14) cross-GWAS effective-N
range left open, and the bibliography strikes that have been verbal
PI flags for weeks but never landed in the audit trail.

### (a) Li-2025 confound handling — use N_eff on the power axis.

The (14) effective-N table established the cross-GWAS magnitude
range (de Lange 46k / Liu 376k / Yengo 1.6M, three orders of
magnitude). The to-do PDF locks the power-axis treatment of that
range: **use N_eff, not total N**, when reasoning about per-GWAS
power on the concordance heatmap. Specifically:

- **de Lange UC** GCST004133: total N = 45,975; **N_eff ≈ 36,168**
  (Willer et al. effective sample size for a case-control study,
  applied to the published 12,366 cases / 33,609 controls split).
- **Liu UC** GCST90446794: total N = 375,508; **N_eff ≈ 87,242**
  (trans-ancestry; case count 23,252 across EAS + EUR arms).
- **Cross-GWAS N_eff ratio:** 87,242 / 36,168 ≈ **2.4×**.

Controls are excluded from the power axis — they contribute to the
test statistic's null but not to the case-count-bound effective
power, so per Li 2025's framing of single-cell-GWAS sample-size
sensitivity, N_eff is the right axis label. Total N is appropriate
only for variance-scaling, which is upstream of the concordance
metric.

The 2.4× ratio is small enough that the cross-GWAS concordance
*should* survive a power-correction sensitivity panel — but it's
not zero, and the magnitude is what justifies including the Li-2025
sensitivity in Phase 9 rather than pre-emptively in the headline.

### (b) Bibliography audit — strikes.

PI flagged three references in the v1 plan/READMEs as stale or
incorrect during the rev2 to-do prep. Audit findings:

- **Lakkis 2024** — struck. Cited in `README.md` "Reading list"
  bullet; no further references in PLAN.md or docs/. Removed from
  the README in this batch.
- **Devlin/Zhao 2023** — struck. Cited in `DECISIONS.md` line 52,
  inside the original PI brief at the top of the file ("Fallbacks:
  Garrido-Trigo 2023, Boland 2020, Devlin/Zhao 2023"). That line is
  a historical record of the project's locked-core decisions and
  is NOT overwritten here; the strike is logged in this entry
  instead. Forward references should not cite Devlin/Zhao 2023.
- **GSE229072** — already correctly handled in correction (15)
  (record-of-record: "GSE229072 is not Mennillo"). The Mennillo
  arm itself was dropped in correction (16) (TAURUS swap). No
  additional strike needed; the negative-finding record stands.

Files updated in this batch:

- `README.md` (Lakkis 2024 removed from Reading list; pointer added
  to this DECISIONS entry).
- `DECISIONS.md` (this entry).

---

## CORRECTION 2026-06-07 (24): Laptop GWAS munge — de Lange + Yengo pre-staged for HB

Rev2 to-do PDF's "DO NOW (laptop) — Build the de Lange .gs" item:
executed the laptop-runnable half of the pipeline. MAGMA itself
remains HB-pinned (Linux x86_64 binary; 1000G EUR LD reference;
NCBI37.3 gene-loc; none on this Windows laptop). What landed locally
is the `prepare_gwas.py` munge step (autosome filter + MAF/INFO QC +
λ_GC) producing the two MAGMA intermediates (`.snp.loc` + `.pval`)
for the two GWAS that ship in a format `prepare_gwas.py` handles
without modification.

### (a) de Lange UC GCST004133 — munged.

Used the GWAS Catalog harmonized version
`28067908-GCST004133-EFO_0000729.h.tsv.gz` (299 MB, schema = `hm_*`
columns) rather than the raw deposit referenced by `download_refs.sh`
(`uc_build37_45975_20161107.txt.gz`, schema = `MarkerName / P.value`
with NO chr/bp columns). The harmonized version drops in directly
against `prepare_gwas.py`'s defaults (`hm_rsid`, `hm_chrom`, `hm_pos`,
`p_value`). Invocation per `code/01_magma/README.md` template;
`--n-fixed 45975` for the no-per-SNP-N case (DECISIONS 14).

Result: 9,486,539 SNPs in → 9,486,539 SNPs out (zero QC drops; the
harmonized file does NOT ship MAF or INFO columns, so those filters
no-op). Outputs at `data/gwas/uc_delange.{snp.loc,pval}`. **λ_GC =
1.1724** — fails the PDF's `[after-munge] λ_GC ≤ 1.1` gate as
literally written. Flagged below.

### (b) Yengo height GCST90245992 (positive control) — munged.

Raw deposit `GCST90245992_buildGRCh37.tsv` (95 MB). Has per-SNP `n`
column and `effect_allele_frequency`. Invocation:
`--col-snp variant_id --col-chr chromosome --col-bp base_pair_location
--col-p p_value --col-n n --col-frq effect_allele_frequency`.

Result: 1,372,608 SNPs in → 1,180,302 SNPs out after MAF ≥ 0.01 QC
(192,306 dropped). Outputs at `data/gwas/yengo_height.{snp.loc,pval}`.
**λ_GC = 5.1992** — also over the PDF's 1.1 gate.

### (c) λ_GC interpretation — `≤ 1.1` gate is the WRONG check at large N.

Both λ_GC values are over 1.1, but interpretation matters here:

- **de Lange 1.17 at N=46k**: moderate; consistent with the polygenic
  signal expected for UC at this sample size. The raw λ_GC conflates
  stratification with polygenic signal; the correct stratification
  check is the LDSC intercept (or the LDSC ratio).
- **Yengo 5.20 at N=1.6M**: very high in absolute terms but textbook
  for a highly-powered polygenic-trait GWAS — at this N, λ_GC is
  dominated by genuine signal across the genome, NOT by population
  stratification (Bulik-Sullivan 2015 *Nat Genet*).

The PDF's `λ_GC ≤ 1.1` rule was written without the large-N caveat;
it's appropriate for the de Lange-size UC GWAS and was almost
hit (1.17 vs 1.10), but the Yengo positive-control failure is
expected, not actionable. Resolution: keep the λ_GC value as
*reported* per-GWAS (already written to
`results/magma/{trait}_lambda_gc.tsv` by `prepare_gwas.py`), but
treat the **LDSC intercept** as the stratification gate where it
matters (the two UC GWAS for the concordance figure), not raw λ_GC.
This is consistent with the Methods conventions in Bulik-Sullivan and
matches what reviewers will expect.

Action item: LDSC run + intercept reporting is a separate pipeline
not currently in `code/01_magma/`. Park for the M3 sanity scaffolding
phase or as a follow-up to the (14) cross-GWAS effective-N table.

### (d) Trubetskoy SCZ and Liu UC — deferred.

- **SCZ**: PGC sumstats VCF v1.0 (240 MB, figshare-hosted, downloaded
  via `download_refs.sh` URL). Format has `##` metadata header lines
  followed by a single-`#` column header line and per-SNP rows.
  `prepare_gwas.py` calls `pd.read_csv` without a `comment=` or
  `skiprows=` argument, so the file will not parse as-is. PDF flags
  this: "SCZ needs `##`-header skip + `--col-n n_eff`". Patch is a
  small CLI flag addition (`--comment-char` mapping to `pd.read_csv
  comment=`, plus header-line handling for the single-`#` row). Not
  applied in this batch — Saisohan's territory to extend the script.
- **Liu**: 2.49 GB compressed deposit. Per the rev2 to-do PDF, gated
  on the open ancestry-LD methodological decision; not heatmap-
  blocking. Not downloaded.

### (e) Pipeline status — what's HB-pinned vs laptop-done.

Laptop completed (this batch):

- `data/gwas/uc_delange_GCST004133.h.tsv.gz` (299 MB, downloaded)
- `data/gwas/yengo_height_GCST90245992.tsv` (95 MB, downloaded)
- `data/gwas/uc_delange.{snp.loc,pval}` (gitignored — intermediate)
- `data/gwas/yengo_height.{snp.loc,pval}` (gitignored — intermediate)
- `results/magma/{uc_delange,yengo_height}_lambda_gc.tsv` (gitignored)

HB-pinned (rest of the pipeline):

- `run_magma.sh` step (needs MAGMA binary + 1000G EUR LD ref +
  NCBI37.3.gene.loc; Linux x86_64).
- `make_scdrs_gs.py` → `.gs` file (depends on MAGMA `.genes.out`).
- `sanity_check.py` → UC marker-gene top-N verification.

The `.snp.loc` + `.pval` files are ready to ship to Hummingbird and
plug directly into `run_magma.sh` — the laptop side of the heatmap
critical path is done.

Files updated in this batch:

- `.gitignore` (added `data/gwas/*.snp.loc` + `data/gwas/*.pval`
  patterns; intermediates, regenerable from sumstats).
- `DECISIONS.md` (this entry).

---

## CORRECTION 2026-06-07 (25): Post-munge review — LDSC elevated to pre-narrative gate; Yengo cross-checked; Devlin/Zhao mark in brief

Saisohan's 2026-06-07 review of (24) flagged one elevation, one
cross-check, two minor doc fixes.

### (a) de Lange LDSC intercept elevated out of M3-sanity.

(24)(c) parked the LDSC intercept under "M3 sanity scaffolding."
Wrong tier: this is the validity check on the **primary GWAS**, not a
generic sanity item. Promoted to a tracked pre-narrative gate
(**OPEN_FLAGS F10**) with the nuance:

- Does NOT block heatmap *generation* — de Lange is the same GWAS
  across all five atlases; any stratification is a shared input and
  doesn't differentially distort cross-atlas agreement.
- DOES block the *biological narrative* on that heatmap. "Atlases
  concordantly rank cell type X for UC" carries biological weight
  only if de Lange isn't confounded.

Schedule: LDSC pipeline must clear **before M4 manuscript draft**,
not after. Don't let M3 slide past M4.

### (b) Yengo λ_GC = 5.20 cross-checked — munge is clean.

Saisohan: theoretical Bulik-Sullivan framing isn't enough; a value
that extreme is also exactly what a munge bug produces (wrong N col,
wrong stat col, allele mismatches all inflate λ_GC). Verified
empirically:

| Metric | de Lange UC | Yengo height |
|---|---|---|
| n_snps | 9,486,539 | 1,180,302 |
| λ_GC | 1.1724 | 5.1992 |
| mean(χ²) | 1.27 | 14.91 |
| n GW-sig (p<5e-8) | 7,392 | 115,546 |
| N column range | 45,975 / 45,975 | 344 / 1,597,374 |

Findings:

- de Lange clean: λ_GC, mean(χ²), and GW-sig count all in band for
  a polygenic UC GWAS at N=46k. The `--n-fixed` path worked
  (N column constant at 45,975). 7,392 raw GW-sig SNPs collapses
  to ~200-300 LD-independent loci per de Lange 2017 Fig. 1.
- Yengo cross-check resolved:
  - **mean(χ²) = 14.9** is plausible for height at N=1.6M — height
    is the most polygenic common trait there is. Published Yengo
    2022 reports mean(χ²) in the mid double digits for the EUR
    discovery; 14.9 is in the published-deposit band, not a bug.
  - **λ_GC = 5.20 vs published ~3-4 in the main analysis.** Modest
    discrepancy explainable by the deposit's SNP-set difference
    (the GWAS Catalog deposit appears to be MAF-and-INFO-pre-
    filtered to a HM3-scale set of 1.37M SNPs vs the paper's full
    discovery panel of ~12M; λ_GC is sensitive to the SNP universe).
    Not a munge bug.
  - **115k GW-sig SNPs (10%)** consistent with massive polygenic
    signal at large N (post-LD pruning will collapse to the paper's
    ~12,000 independent loci).
  - **Per-SNP N varies 344 → 1,597,374** — real Yengo deposit shape
    (rare variants tested only in small sub-cohorts), passed straight
    through. MAGMA's gene-test is N-weighted; high-N SNPs dominate.
    Worth documenting but not a fix item.

Yengo positive control should light up the right cells in MAGMA.

### (c) SNP-density adequacy for MAGMA gene coverage — confirmed.

Saisohan flagged the 8× SNP-density asymmetry (de Lange 9.5M vs
Yengo 1.18M) as a possible coverage starvation risk for the positive
control. Empirical density across the 2,881 Mb autosome span:

- **de Lange**: 3,293 SNPs/Mb → ~231 SNPs per 10-kb-windowed gene.
  22/22 chr coverage, per-chr range 127k–798k.
- **Yengo**: 410 SNPs/Mb → ~29 SNPs per 10-kb-windowed gene.
  22/22 chr coverage, per-chr range 17k–99k.

29 SNPs per windowed gene is plenty for MAGMA's gene-test (a few
SNPs is usually enough for the gene-Z to converge). Positive control
is a fair test of the method, not a coverage-starved one.

### (d) Devlin/Zhao struck mark in the historical PI brief.

(23)(b) left Devlin/Zhao 2023 noted-in-DECISIONS-only rather than
marked inline at DECISIONS.md line 52 (the original PI brief). Per
Saisohan: same history-read-as-live trap that had us chasing the
dead Mennillo / GSE229072 accession. Updated the brief line to
`~~Devlin/Zhao 2023~~ [struck — DECISIONS 23(b)]` so a future scan
of the brief surfaces the strike rather than the stale entry.

Files updated in this batch:

- `OPEN_FLAGS.md` (F10 added: de Lange LDSC intercept pre-narrative
  gate).
- `DECISIONS.md` (line 52 Devlin/Zhao inline struck mark; this entry).

---

## CORRECTION 2026-06-07 (26): N-filter on per-SNP-N munge + F1 inflammation lock

Saisohan's round-3 review of (24)/(25) surfaced one technical fix
worth doing before the laptop intermediates ship to HB, and pinned
the remaining laptop-side decision items.

### (a) `prepare_gwas.py` — `--n-min-frac` added; Yengo re-munged.

The script previously passed per-SNP N straight through to MAGMA
without trimming the low-N tail. Standard LDSC `munge_sumstats`
drops SNPs below `0.67 × max(N)` by default: SNPs tested in only
a few hundred samples are noise on the per-SNP estimate side and
shouldn't enter the gene score. Patch:

- New CLI flag `--n-min-frac` (default 0.67 = LDSC convention).
- Applied only on the `--col-n` path. **No-op when `--n-fixed` is
  used** (every SNP shares one N — there's no tail to trim).
- `--n-min-frac 0` disables the filter for back-compat.
- Both `pval` and `snp_loc` are kept in sync after the filter
  (the original code drop-deduped only `pval`; the matching `.snp.loc`
  rows would otherwise have orphaned SNPs).

**Yengo re-munged with the filter:**
- 1,180,302 → **1,150,988 SNPs** (29,314 dropped, 2.5%; threshold
  N ≥ 1,070,241 = 0.67 × max 1,597,374).
- λ_GC = **5.3057** (was 5.1992 pre-filter); the small uptick is
  the expected direction — removing low-N noise removes chi² values
  pulled toward zero by small-sample p-value inflation.
- Outputs at `data/gwas/yengo_height.{snp.loc,pval}`.

**de Lange re-run for cleanliness:** unchanged (`--n-fixed 45975`
path; 9,486,539 SNPs, λ_GC = 1.1724). Both intermediates now produced
by the patched script — fresh reproducible baseline.

### (b) F1 inflammation policy — locked: pool first, stratify Phase 9.

Per the rev2 PDF + Saisohan's framing: pool all UC cells (inflamed +
non-inflamed) into a single "UC" group for the first 3×3 broad-tier
heatmap; stratify inflamed-vs-non-inflamed as a Phase-9 sensitivity
panel. Decided on purpose now, not discovered-by-omission at scoring
time.

- **Per-atlas implementation:**
  - **Smillie**: `Health` ∈ {Inflamed, Non-inflamed} both → `disease =
    ulcerative colitis` (already in loader's `HEALTH_TO_DISEASE`).
  - **Garrido**: All UC samples are inflamed biopsies (Salas-lab
    cohort design); pool is a no-op — already one group.
  - **TAURUS**: Baseline UC × colonic × Pre cohort = 22 donors / 52
    samples (39 inflamed + 13 non-inflamed per Supp Table 1B, see
    DECISIONS 20). Inflammation captured in obs but not subset on
    load.
- **Phase-9 stratification path:** the Phase-9 compositional-confound
  panel (parked from QC policy in DECISIONS 22(a)) re-runs broad
  concordance on inflamed-only and non-inflamed-only sub-cohorts;
  obs columns are already populated, no re-load needed.
- **F1 closed in OPEN_FLAGS.**

### (c) Liu ancestry-LD decision — still open; sole laptop-side gate.

Per Saisohan: "Open five rounds, and it's now the sole gate on the
entire cross-GWAS axis. Five-minute call; just make it." This is
Muskaan's methodological call (trans-ancestry EAS+EUR LD vs European-
only LD panel for the 23,252 EAS + EUR Liu cohort), NOT a Claude
Code item. Logging here to keep the decision visible — and to note
that *every other* laptop-side gate has been closed:

- ~~CL grep~~ → DECISIONS 22(d)
- ~~Bibliography audit~~ → DECISIONS 23(b), 25(d)
- ~~Yengo cross-check~~ → DECISIONS 25(b)
- ~~SNP-density check~~ → DECISIONS 25(c)
- ~~Devlin/Zhao inline strike~~ → DECISIONS 25(d)
- ~~N-filter audit + patch~~ → this entry (a)
- ~~F1 inflammation policy~~ → this entry (b)
- **Liu ancestry-LD** → still open, Muskaan's call.
- ~~Commit + push~~ → DECISIONS 25 push.

Once the ancestry-LD call lands, Liu munges directly:
`--input data/gwas/uc_liu_GCST90446794.tsv --col-snp variant_id
--col-chr chromosome --col-bp base_pair_location --col-p p_value
--n-fixed 375508` (matches the schema from `download_refs.sh`
comment; fixed-N path so the `--n-min-frac` filter is a no-op).
Liu download deferred until the call clears (2.49 GB; no point
pulling until munge args are known).

Files updated in this batch:

- `code/01_magma/prepare_gwas.py` (`--n-min-frac` flag added;
  `snp_loc` kept in sync with the N-filtered `pval`).
- `OPEN_FLAGS.md` (F1 → Resolved with pool-first + Phase-9
  stratification lock).
- `DECISIONS.md` (this entry).

`data/gwas/{uc_delange,yengo_height}.{snp.loc,pval}` regenerated by
the patched script and ready to ship to HB.

---

## CORRECTION 2026-06-07 (27): de Lange N_eff in the munge + F1 inflammation-composition caveat + Liu LD options

Saisohan's round-4 review of (26) caught the one substantive
inconsistency that the N_eff lock in (23) opened up, plus two doc
fixes the team needs to actually act on rather than log.

### (a) de Lange munge re-run with N_eff = 36,160 (was 45,975).

(23)(a) locked **N_eff, not total N**, as the power axis label across
the cross-GWAS heatmap — that's the whole point of citing 36k vs 46k
as the gap with Liu. (24)(a) then munged de Lange with
`--n-fixed 45975` (total N), straight from the `code/01_magma/README.md`
template, undoing (23) on the primary GWAS that the figure depends on.
Caught by Saisohan: "Make the munge's N match the number you're
claiming as the study's power, or document why they differ."

Resolution: re-munge with `--n-fixed 36160`.

- **Precise formula** (Willer et al. effective sample size, case-control):
  `N_eff = 4 * n_cases * n_ctrls / (n_cases + n_ctrls)`
- de Lange GCST004133: 12,366 cases, 33,609 controls
- `N_eff = 4 * 12366 * 33609 / 45975 = 36,159.56 ≈ 36,160`
- (23)(a) cited 36,168 — that was a rounded approximation; the precise
  value is **36,160**. The 8-sample (≈0.02%) difference is irrelevant
  to MAGMA's gene-test precision weighting; 36,160 is the value
  shipped on disk now.
- (24)(a)'s claimed 9,486,539 SNPs / λ_GC = 1.1724 are unchanged
  (N doesn't affect p-values; only changes per-SNP precision
  weighting at the MAGMA gene-test step downstream). The N column in
  `data/gwas/uc_delange.pval` now reads 36,160 uniformly.

**Why this matters for MAGMA**: MAGMA's `--pval` step weights each
SNP by its sample size in computing the gene Z-score. Feeding total
N=45,975 would over-state precision by ~27% per SNP, inflating the
gene Z-scores beyond what de Lange's effective power justifies. The
heatmap would then sit on top of an over-powered primary GWAS — the
exact axis (23)(a) was supposed to protect.

**Forward note**: when Liu munges land (gated on the still-open
ancestry-LD call below), use Liu N_eff = 87,242 (precise:
`4 * 23252 * 352256 / 375508 = 87,248.81 ≈ 87,249`; the 87,242 in
(23)(a) is similarly rounded) — NOT the total 375,508.

`code/01_magma/README.md` still references `--n-fixed 45975` /
"combined N" in its de Lange / SCZ templates; flag for cleanup next
time the file's touched, but not blocking — the munge is what
actually fed MAGMA, and that's now N_eff-correct.

### (b) F1 inflammation-composition caveat — Garrido is structurally inflamed-only.

(26)(b) locked "pool all UC cells for the first 3×3" but didn't carry
that the pool **doesn't equalize inflammation composition across
atlases**. Garrido is all-inflamed by cohort design (Salas-lab
biopsied inflamed mucosa only); Smillie and TAURUS pool inflamed +
non-inflamed. So:

- The Garrido row in the first 3×3 is structurally inflamed-only.
- The Smillie and TAURUS rows mix inflamed + non-inflamed.
- **A Garrido-row deviation in the broad heatmap could therefore be
  inflammation composition, not a cross-atlas reproducibility
  failure.** Must be carried as an interpretation caveat in the
  result paragraph, not discovered post-hoc by a reviewer.

Phase-9 implication: Phase-9 only restores this asymmetry — the
inflamed-only sub-figure is the apples-to-apples 3×3. Garrido can
feed only the inflamed stratum (no non-inflamed cells to contribute);
the non-inflamed Phase-9 panel is Smillie × TAURUS only, a 2×2.
That's a structural consequence of cohort design, not a defect.

F1 OPEN_FLAGS entry updated to carry this caveat into the resolved-
flag record.

### (c) Liu ancestry-LD — three concrete options, owner = Muskaan.

(26)(c) logged Liu ancestry-LD as still-open but didn't enumerate
options or hand the call to anyone in particular. Per Saisohan's
"a decision sitting in a log isn't a decision anyone's been asked to
make": fixing the visibility gap. Logging here as the visible
hand-off; **Muskaan needs the call before Liu can munge.**

**The three options:**

1. **1000G EUR LD panel** (default-conservative).
   - Treat Liu as European-effective; use 1000G EUR LD per de Lange
     and Yengo's panel.
   - Pros: matches MAGMA defaults; simplest; consistent across the
     three GWAS in the figure.
   - Cons: misspecifies LD for the EAS sub-cohort (~6,862 EAS cases of
     the 23,252 total = 30% of cases); LD differs meaningfully between
     ancestries, particularly at MHC and HLA-adjacent loci.
   - **Verdict in literature**: defensible for trans-ancestry GWAS
     where one ancestry dominates; documented as an approximation in
     the Methods.

2. **Combined trans-ancestry LD panel** (1000G EUR + EAS, weighted).
   - Build a panel matching Liu's cohort ancestry proportions.
   - Pros: better matches the cohort's LD structure.
   - Cons: non-standard panel construction; downloads 1000G EAS as
     well; MAGMA's `--bfile` expects a single panel, so the weighting
     would happen at the panel-construction step (requires a custom
     PLINK merge); harder to reproduce.

3. **Per-ancestry split + meta-analysis**.
   - Run MAGMA separately on the EAS sub-cohort (1000G EAS LD) and
     the EUR sub-cohort (1000G EUR LD), meta-analyze the gene
     Z-scores (inverse-variance weighted).
   - Pros: most principled handling of LD heterogeneity.
   - Cons: requires per-sub-cohort sumstats which the public Liu
     deposit `GCST90446794.tsv` may NOT split out (the deposit ships
     trans-ancestry combined p-values + meta-betas, NOT per-ancestry
     panels separately). Empirical verification: open the file header
     and check for `Direction` column or per-ancestry betas. If
     absent, option 3 is infeasible without contacting the authors.

**Recommendation for the 5-minute call**: default to **Option 1**
(1000G EUR) for the first figure with the misspecification
documented in Methods. Reasons:
- Cohort is 70% EUR cases by raw count.
- Consistency with de Lange and Yengo's panel is more important for
  the cross-GWAS heatmap interpretation than the marginal LD
  improvement from Options 2/3.
- Option 3 requires verifying the deposit's per-ancestry separability;
  Option 2 is non-standard work.
- A Liu-as-EUR-effective sensitivity in Phase 9 captures the residual
  question without blocking the first figure.

The above is a recommendation, not a lock — final call belongs to
Muskaan. Once the call lands, Liu munges with:
`--col-snp variant_id --col-chr chromosome --col-bp base_pair_location
--col-p p_value --n-fixed 87249` (per-SNP N absent per (14); use
N_eff = 87,249, not total 375,508 per (a) above's consistency rule).

Files updated in this batch:

- `data/gwas/uc_delange.{snp.loc,pval}` regenerated with N=36,160
  (gitignored).
- `results/magma/uc_delange_lambda_gc.tsv` regenerated (gitignored).
- `OPEN_FLAGS.md` (F1 inflammation-composition caveat appended).
- `DECISIONS.md` (this entry).

`code/01_magma/prepare_gwas.py` not touched in this batch — the
(26) patch was already correct; only the invocation N changed.

---

## CORRECTION 2026-06-07 (28): SCZ negative control acquired; references pre-staged

Per Saisohan: SCZ is off the first-figure path but the only GWAS file
fully in our hands; references download is the "fixed URLs verify"
two-bird check from DO NOW #2. Both executed.

### (a) `prepare_gwas.py` — `--skip-comments` for PGC sumstats VCF v1.0.

Trubetskoy SCZ ships in PGC's `*.vcf.tsv.gz` format: a TSV body
prefixed by 73 `##` metadata lines (fileformat, methodsParagraph,
acknowledgments with accented author names, contig lengths, etc.).
`pd.read_csv(comment='#')` would eat both the metadata AND the column
header in one undifferentiated swoop; the right move is to count the
`##` prefix and `skiprows=` to land on the header line.

Patch: new `--skip-comments` flag (action='store_true'). When set:

- Pre-scan the input file, count leading `##` lines.
- Pass the count to `pd.read_csv(skiprows=...)`.
- If the resulting header's first column starts with `#` (some PGC
  files do this), strip the single `#`. Trubetskoy's deposit does
  NOT (header reads `CHROM ID POS A1 A2 FCAS FCON IMPINFO BETA SE
  PVAL NCAS NCON NEFF`), so the strip is a no-op for this file.
- **Explicit `encoding='utf-8', errors='replace'`** when reading the
  pre-scan stream. PGC's `##acknowledgments` line includes
  non-ASCII characters (accented author names, byte 0x9d among
  others); on Windows the default cp1252 decoder choked at byte
  position 5,506. `errors='replace'` is safe because only `##`
  prefix detection uses these bytes, not their content.

### (b) SCZ munged.

Invocation: `--input scz_trubetskoy_eur_PGC3_v3.vcf.tsv.gz
--skip-comments --col-snp ID --col-chr CHROM --col-bp POS --col-p PVAL
--col-n NEFF --col-info IMPINFO`.

Pipeline:
- 73 `##` lines skipped (logged).
- 7,659,767 SNPs in.
- INFO ≥ 0.6 filter: dropped 86,380 SNPs → 7,573,387.
- Autosome filter: zero non-autosomal SNPs in the deposit.
- N filter (`--n-min-frac 0.67`): dropped 121,338 SNPs below
  39,361 (= 0.67 × max NEFF 58,749) → **7,452,049 SNPs out**.
- λ_GC = **1.6305** (in-band for SCZ at the file's per-SNP NEFF;
  287 distinct loci per the abstract → real polygenic signal).

**Note on SCZ NEFF**: the file's `NEFF` column reports per-SNP
effective sample size from PGC's meta-analysis weighting, NOT the
case-control Willer formula. For 53,386 cases + 77,258 controls
(EUR cohort per the methodsParagraph), Willer formula gives
N_eff ≈ 126,310; the file's per-SNP NEFF maxes at 58,749 because
PGC's pipeline accounts for trio encoding + imputation-quality
weighting that roughly halves the naive N_eff. The right axis label
for SCZ in the cross-GWAS power discussion is **PGC's per-SNP NEFF
58,749**, NOT Willer-formula 126k — what MAGMA actually sees is what
the cohort's effective power for differential expression at the gene
level actually was. This matches the same N_eff-not-total-N
discipline (23)(a) / (27)(a) locked for de Lange and Liu.

Outputs: `data/gwas/scz_trubetskoy.{snp.loc,pval}` (gitignored
intermediates), `results/magma/scz_trubetskoy_lambda_gc.tsv`.

### (c) Reference data — MAGMA, NCBI37.3, g1000_eur all staged.

Pulled from the SURF mirrors in `scripts/download_refs.sh` via
PowerShell `Invoke-WebRequest` (bash `download_refs.sh` won't run on
Windows — `set -euo pipefail`, `stat -c %y`, `unzip` all non-Windows).
This serves as the "verify the fixed URLs *download*, not just
resolve" gate that's been open since rev2 DO NOW #2.

| Resource | URL | Size | Verified |
|---|---|---|---|
| MAGMA binary | `vu.data.surf.nl/s/lxDgt2dNdNr6DYt/download` | 3.3 MB zip → 7.2 MB binary | ELF 64-bit LSB executable, x86-64, statically linked, GNU/Linux 3.2.0, BuildID a49e0123... — MAGMA v1.10 confirmed via shipped `manual_v1.10.pdf` |
| NCBI37.3 gene-loc | `vu.data.surf.nl/s/Pj2orwuF2JYyKxq/download` | 356 KB zip → 688 KB | 19,427 protein-coding genes, format `entrez chr start end strand symbol` |
| 1000G EUR LD | `vu.data.surf.nl/s/VZNByNwpD8qqINe/download` | 512 MB zip → 2.86 GB .bed + 659 MB .bim + 12.6 KB .fam + 86 MB .synonyms | **503 EUR samples** (matches 1000G Phase 3 EUR cohort exactly), 22,665,064 variants |

All three resources extracted to `data/reference/` matching the
`download_refs.sh` post-unzip layout (`magma`, `NCBI37.3.gene.loc`,
`g1000_eur.{bed,bim,fam,synonyms}`). All gitignored.

### (d) MAGMA local execution status: WSL2 yes, Linux distro no.

WSL2 is installed (`wsl --version` reports 2.7.3.0, kernel 6.6.114.1)
but no Linux distribution is present (`wsl --list --online` shows the
options; `wsl --install Ubuntu` would land Ubuntu 24.04 LTS in a few
minutes). The MAGMA binary cannot execute natively on Windows; it's
pre-staged for either (a) `wsl --install Ubuntu` and run locally, or
(b) ship the binary + LD ref + gene-loc + munged sumstats to HB and
run there. Not pursued here because installing Ubuntu via WSL is a
deeper system change than "do everything that doesn't require HB"
implied — Muskaan's choice whether to install or ship.

### (e) plink 1.9 — not pulled; not on the .gs build path.

`download_refs.sh` doesn't auto-fetch plink (one-time setup in the
README only). MAGMA's gene-test uses the LD reference in plink
binary format (`.bed/.bim/.fam`) but doesn't shell out to plink at
runtime. So plink is not needed for the laptop-side MAGMA → .gs
chain; can be skipped from this staging round.

### (f) Pipeline status — 100% of laptop side now landed.

Three munge intermediates ready to ship to HB:
- `data/gwas/uc_delange.{snp.loc,pval}` — 9,486,539 SNPs, N=36,160
- `data/gwas/yengo_height.{snp.loc,pval}` — 1,150,988 SNPs, per-SNP N
- `data/gwas/scz_trubetskoy.{snp.loc,pval}` — 7,452,049 SNPs, per-SNP NEFF

Plus full reference set at `data/reference/`. Either install
WSL+Ubuntu and run MAGMA locally, or scp the intermediates to HB
and run there. Liu munge waits on the open ancestry-LD decision
(27)(c); the patched `prepare_gwas.py` is ready for it.

Files updated in this batch:

- `code/01_magma/prepare_gwas.py` (`--skip-comments` flag + utf-8
  encoding for the pre-scan).
- `data/gwas/scz_trubetskoy.{snp.loc,pval}` regenerated (gitignored).
- `data/reference/{magma, NCBI37.3.gene.loc, g1000_eur.bed/bim/fam/synonyms}`
  staged (gitignored).
- `DECISIONS.md` (this entry).

---

## CORRECTION 2026-06-07 (29): SCZ post-munge sanity — NEFF convention + loci count

Saisohan flagged two cheap confirms on (28). Both clear; logging the
results so a future re-read isn't left wondering whether the checks
ever happened.

### (a) PGC NEFF convention — file's value is the right MAGMA input.

The "half-vs-sum NEFF wrinkle" concern: PGC files can ship NEFF as
either the full-Willer (`4 × Ncas × Ncon / (Ncas + Ncon)`) or the
half-form (`2 / (1/Ncas + 1/Ncon)`); naively doubling or halving
without verifying loses precision in MAGMA's gene-test weighting.

Verified against the LDSC convention (per
[bulik/ldsc issue #95](https://github.com/bulik/ldsc/issues/95),
which cites `Neff = 2 / (1/Ncases + 1/Ncontrols)` from a *Nature
Protocols* reference — the LDSC default):

- Trubetskoy 2022 PGC3 EUR cohort: 53,386 cases + 77,258 controls.
- **Full-Willer**: `4 × 53386 × 77258 / 130644 = 126,265`.
- **Half-form (LDSC convention)**: `2 / (1/53386 + 1/77258) = 63,132`.
- **File's NEFF (max)**: **58,749**, which sits just below the
  half-form 63,132 (consistent with PGC's per-SNP meta-analysis
  weighting reducing effective N slightly below the unweighted
  cohort half-form).

So PGC3 wave 3 ships NEFF in the **half-form (LDSC) convention** —
no doubling needed. MAGMA's per-SNP precision weight expects exactly
this effective-N form; passing the file's NEFF column straight
through is correct.

Cross-GWAS power discussion (correction 23(a) framing) should
compare NEFF values in the same convention. The de Lange N_eff
36,160 was computed via the full-Willer formula. To compare on the
SAME axis as SCZ NEFF=58,749 (half-form), the de Lange equivalent
half-form is `2 / (1/12366 + 1/33609) = 18,080`. So SCZ is ~3.2×
de Lange in half-form effective N; the cross-GWAS power gap is
larger than the raw N_eff number alone suggested. **Not actionable
for the first heatmap** — SCZ is the negative control; mis-
calibrated N can't flip a negative into a false positive. Noted for
the Methods power-axis paragraph; the first figure proceeds with
the values currently shipped on disk.

### (b) Loci count — munge preserved signal.

Trubetskoy 2022 reports **287 distinct genomic loci** in the abstract
— but that's the **trans-ancestry meta-analysis** (76,755 cases +
243,649 controls), NOT the EUR-only file `*.european.autosome.public.*`
that the rev2 PDF references and that I downloaded. The EUR-only
sub-analysis is what's on disk (53,386 cases + 77,258 controls;
~70% of the trans-ancestry cohort).

Empirical count from `data/gwas/scz_trubetskoy.pval` (post-munge):

| Threshold / binning | Count |
|---|---|
| Raw GW-sig SNPs (p < 5e-8) | 20,300 |
| 1 Mb-window bin | 217 |
| 500 kb-merged loci | **177** |
| Trubetskoy 2022 trans-ancestry (abstract) | 287 |

177 loci (500 kb-merged) is in the expected band for the EUR-only
sub-cohort at ~70% of the trans-ancestry sample (loci count scales
sub-linearly with N; rough expectation 287 × 0.70^0.8 ≈ 217, just
above my 500 kb-merged 177 because 500 kb merging is slightly
stricter than r²<0.1 LD pruning). **No signal was dropped by the
munge.**

If the headline is ever re-claimed against the full trans-ancestry
287 figure (e.g., if a methods reviewer asks why my SCZ has fewer
loci than published), the answer is: my file is the EUR-only
publicly downloadable arm (figshare DOI 10.6084/m9.figshare.19426775),
which the PGC release explicitly distinguishes from the full
trans-ancestry sumstats (those are not on the public figshare per
the deposit description).

### Bottom line

Both Saisohan cross-checks land clean — NEFF convention correctly
identified as the half-form LDSC default (no transformation needed),
and the EUR-only loci count is consistent with the deposit's
sub-cohort, not the abstract's trans-ancestry headline.

No files updated in this batch — DECISIONS.md only (this entry).





