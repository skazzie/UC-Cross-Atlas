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
