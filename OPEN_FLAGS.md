# Open Flags — deferred decisions to commit to later

Running list of judgment calls and unresolved questions surfaced during the build but
not yet decided. The companion to DECISIONS.md: that file is *locked* decisions, this is
the *still-open* side. When a flag is resolved, log it in DECISIONS.md and strike it here.

_Last updated: June 7, 2026._

---

## ~~F1 · UC tissue definition — inflamed vs non-inflamed~~ — RESOLVED 2026-06-07 (DECISIONS 26)
Locked: **pool all UC cells for the first 3×3 broad heatmap** (inflamed
+ non-inflamed combined as one "UC" group). Inflamed-vs-non-inflamed
stratification is a **Phase-9 sensitivity panel**, not a v1-figure
gate. Applies uniformly across Smillie / Garrido / TAURUS — TAURUS's
baseline-pretreatment timepoint includes both inflamed and
non-inflamed biopsies per the rev2 PDF cohort breakdown (39 inflamed
+ 13 non-inflamed of the 52 baseline samples). Loaders already capture
the inflammation state in obs (`health` / `inflammation_score`), so
the Phase-9 split is a downstream subset, no re-load needed.

## ~~F2 · Garrido QC-state labels (MT / heat-shock / IER)~~ — RESOLVED 2026-06-07 (DECISIONS 22)
Resolved as a cross-atlas `QC_STATE_TO_PARENT` collapse policy in
`code/02_atlas_prep/_qc_policy.py`. Collapse-to-parent (not exclude); applied
identically by Garrido / Smillie / TAURUS loaders. Smillie's `MT-hi` stays
mapped at broad without a fine-tier collapse (single Imm-compartment cluster).
Phase-9 compositional-confound check parked: whether stress-state fractions
correlate with disease state is sensitivity-panel territory, not a v1 policy
question.

## ~~F3 · Garrido lineage-ambiguous cycling labels~~ — RESOLVED 2026-06-07 (DECISIONS 22)
Resolved as a cross-atlas `EXCLUDE_LINEAGE_AMBIGUOUS_FINE` rule in
`code/02_atlas_prep/_qc_policy.py`. Unprefixed `Cycling cells` / `Cycling
cells 2` / `Cycling cells 3` (Garrido) drop before the broad-tier mapping;
all loaders import the same set so the exclusion is symmetric. Force-
assignment was the failure mode (would have created a B-cell discordance
that's pure harmonization artifact); reinstatement is the marker-QC step
(MKI67 + lineage markers + cross-atlas placement agree). Smillie's cycling
clusters are all compartment-prefixed — kept, no labels match the
exclusion set today; TAURUS hook in place pre-first-run.

## ~~F4 · Garrido crosswalk REVIEW rows (12)~~ — RESOLVED 2026-06-07 (DECISIONS 22)
9 of 12 REVIEW rows locked at broad with their `load_garrido_trigo.FINE_TO_BROAD`
assignment; 3 unprefixed cycling rows handled via F3 exclusion above. Per-row
broad-tier rationale documented in DECISIONS (22)(c). Fine-tier fold-vs-distinct
questions (Laminin / PLCG2 colonocytes, Mature goblet, PC IGLL5 isotype
grouping) flow to F8 (fine vocab).

## F5 · Broad-tier vocabulary size
Draft harmonized broad tier is 15 categories (top of the 10–15 target). Optional: fold
Pericyte → Fibroblast and reconsider Glia to land ~13. Interacts with what the other four
atlases can resolve, so settle once their native labels are in.
- **Bites at:** lock broad tier (end M2).

## F6 · Persistent storage + backup (standing risk)
Hummingbird scratch auto-purges after 2 weeks; home (1 TB) and lab folders are **not
backed up**. Decide where persistent project data lives (`~/uc-cross-atlas-data/` vs a
lab/group folder) and back up irreplaceable outputs (processed atlases, results) off
the cluster. Reference data and processed atlases must not be left on scratch.
- **Bites at:** ongoing. Reference data (`g1000_eur`) is currently sitting on scratch.

## F7 · Smillie crosswalk REVIEW rows — fine-tier only after 2026-06-07
Broad-tier portion locked 2026-06-07 (DECISIONS 22): `MT-hi` resolved
via the QC policy (kept mapped to T cell at broad; SCP259's own MT% QC
gate is the viability check); `M cells` / `Immature Enterocytes 1/2` /
`Secretory TA` confirmed at broad with their current
`load_smillie.FINE_TO_BROAD` assignments — same logic as Garrido's F4
9-of-12 lock (the broad call holds even when the fine-tier identity
is ambiguous). What remains:
- `M cells` — microfold cell, no clean fine-tier home outside colonocyte
- `Immature Enterocytes 1/2` — progenitor/colonocyte fine boundary
- `Secretory TA` — progenitor vs secretory fine lineage
These three are fine-tier identity questions, not broad. **Bites at:**
F8 fine vocab (T2.5).

## F8 · Fine-tier cross-atlas harmonization layer (blocks the fine concordance axis)
`06_concordance` compares atlases on the *string intersection* of cell-type
labels (`set(scores_a) & set(scores_b)`). Broad tier works because
CANONICAL_BROAD forces identical strings; fine tier has no equivalent — loaders
emit each atlas's *native* fine labels (Smillie 51 `Cluster`, Garrido
91-collapsed, TAURUS 109 cell states pending) and no canonical-fine
vocab/crosswalk exists. Native labels won't string-match across studies (even
"Inflammatory Fibroblasts" vs "Inflammatory fibroblasts"), so the fine
intersection is near-empty and the README's k=20 fine-tier Jaccard is
undeliverable as wired.
Need: a fine-harmonization layer (canonical fine vocab, or per-atlas crosswalk
to shared names). Achievable shared fine set = subtypes biologically alignable
AND ≥ MIN_CELLS_PER_TYPE (50) in all three atlases ≈ a few dozen well-defined
types; study-specific inflammatory states excluded/caveated. Cannot be confirmed
until the TAURUS loader exists (3-way intersection bounded by worst atlas).
- **Bites at:** fine-tier concordance (step 06); gates the fine-tier novelty
  axis (the gap Li et al. 2025 flag as future work). Scope before banking fine
  tier as a headline result.

### F8 · T2.5 analytical sketch — CL-aware candidate fine vocab (2-of-3 atlases known)
Saisohan confirmed v2 is committed → fine vocab is CL-aware. Each fine bucket
anchors on a CL subtree, not just a string. The candidate below uses the two
atlases whose fine labels are on disk (Smillie 51, Garrido 91-collapsed) plus
the post-QC-policy filters (DECISIONS 22: exclude unprefixed cycling; collapse
MT/heat-shock/IER to parent). TAURUS column is structural placeholder — fills
on first compute-node run when LOW_TO_BROAD surfaces and the cell_state (109)
labels are enumerated against this candidate set.

Selection rule: a fine bucket lands here only if (i) it has a CL subtree
anchor in the 2026-03-26 pin, (ii) Smillie + Garrido each have ≥1 native
label that maps in (≥50 cells per atlas check happens at first compute-node
run, not pre-emptively), and (iii) the biology is shared across UC atlases
(study-specific inflammatory states tagged as caveated, not used as the
primary harmonization key).

Tentative 14-bucket candidate (refine after TAURUS surfaces):

| # | Fine bucket | CL anchor | Smillie candidate(s) | Garrido candidate(s) | Notes |
|---|---|---|---|---|---|
| 1 | absorptive enterocyte | CL:1000347 colonocyte | Enterocytes, Best4+ Enterocytes | Colonocyte 1/2, BEST4 OTOP2 | Fold inflammatory / Laminin / PLCG2 subtypes into this bucket per DECISIONS 22(c) fine-tier carry-over |
| 2 | immature enterocyte / TA-adjacent | CL:0009010 transit amplifying cell | Immature Enterocytes 1/2, Enterocyte Progenitors | Cycling TA, Secretory progenitor | F7 boundary call lands here |
| 3 | crypt stem | CL:0002250 intestinal crypt stem cell | Stem | (none resolvable after Ribhi collapse) | Garrido under-represents at fine; Smillie carries this bucket |
| 4 | goblet | CL:0000160 goblet cell | Goblet, Immature Goblet | Goblet, Mature goblet | Fold Mature goblet here per F4 fine call |
| 5 | Paneth-like (secretory antimicrobial) | (no CL anchor — secretory parent only) | (none) | Paneth-like | Garrido-only; flag as study-specific, not used as harmonization key |
| 6 | enteroendocrine | CL:0000164 enteroendocrine cell | Enteroendocrine | Enteroendocrine | Clean 1:1 |
| 7 | tuft | CL:0002204 tuft cell | Tuft | Tuft cells | Clean 1:1; split-vs-keep with #6 above is the open biology call (canonical_broad_DRAFT.md "Open biology call" section) |
| 8 | fibroblast (stromal) | CL:0000057 fibroblast | WNT2B+ / WNT5B+ / RSPO3+ family | S1 / S1.2 / S2a/b / S3 | High cardinality on both sides; fine-tier subtype mapping (WNT2B+ Fos-hi vs S1 etc) needs marker overlap, not just name match |
| 9 | inflammatory fibroblast | CL:0000057 (state subtype) | Inflammatory Fibroblasts | Inflammatory fibroblasts | One of the few inflammatory states alignable across atlases by name AND markers |
| 10 | endothelial | CL:0000115 endothelial cell | Endothelial, Microvascular | Endothelium, Activated endothelium | Vascular fold |
| 11 | lymphatic endothelial | CL:0002138 endothelial cell of lymphatic vessel | (Post-capillary Venules — boundary; revisit) | Lymphatic endothelium | Smillie boundary case; mark as Garrido-strong |
| 12 | T cell (broad lymphoid placeholder — refine post-TAURUS) | CL:0000084 T cell | CD4+ Memory, CD8+ LP, etc. | CD4 ANXA1, CD8 CTL, etc. | High cardinality; fine subtype alignment (CD4 vs CD8, naive vs memory) needs marker-based mapping, not naming |
| 13 | plasma cell (isotype-collapsed) | CL:0000786 plasma cell | Plasma | PC IgA 1-4, PC IgG 1-2, PC IGLL5, Plasmablast variants | Garrido has high isotype-resolution; Smillie has one label. F4 fine-tier isotype-grouping call lands here. |
| 14 | enteroglial | CL:4040002 enteroglial cell | Glia | Glia | Clean 1:1; small N both sides; ≥50 gate may fail in one or both |

Not in v1 fine vocab (defer or exclude):
- Granulocytes (Smillie structural zero; preserved at broad as exclusion-doc'd).
- Unprefixed cycling (DECISIONS 22 exclusion — see F3).
- Lineage-prefixed cycling (Cycling T / Cycling B / Cycling Monocytes /
  Cycling TA) — small N each, fold into their lineage bucket at fine for v1.
- Study-specific inflammatory states beyond #9 (e.g., M1 ACOD1 / CXCL5 in
  Garrido) — caveat as study-specific, exclude from cross-atlas fine
  intersection, present as per-atlas in supplementary.

Validation path (not on the critical path for the first figure):
1. Wait for TAURUS first-run; gate 2 surfaces low + cell_state labels.
2. Walk TAURUS's 109 cell states against this 14-bucket candidate; add/adjust
   buckets per coverage.
3. Per-bucket per-atlas cell-count gate (≥50); document structural zeros
   explicitly.
4. Promote to `code/_shared/canonical_fine_DRAFT.md` and have Saisohan
   red-line before the fine-tier figure runs.

This is analytical sketch only — no loader retarget, no constant promotion,
no `CANONICAL_FINE` import path created until TAURUS has populated and
Saisohan signs off on the bucket list. **Does NOT block** the broad-tier
figure (DECISIONS 22 unblocks that side).

## F10 · de Lange LDSC intercept — pre-narrative gate (NOT M3-sanity)
**Validity check on the primary GWAS, elevated tier vs generic sanity
scaffolding.** Per Saisohan 2026-06-07: λ_GC conflates polygenic
signal with stratification confound (Bulik-Sullivan 2015 *Nat Genet*);
the LDSC intercept is what separates them. de Lange UC GCST004133
munged with λ_GC = 1.1724 — slightly over the rule-of-thumb 1.10 but
in-band for a well-powered polygenic UC GWAS. **The intercept is what
gates interpretation, not generation.**

Nuance — what this gate does and doesn't block:
- **Does NOT block** the 3×3 broad-tier concordance heatmap. de Lange
  is the same GWAS across all five atlases; any stratification is a
  shared input and doesn't differentially distort cross-atlas
  agreement. Compute the heatmap, look at the structure.
- **DOES block** the *biological narrative* on the heatmap.
  "Atlases concordantly rank cell type X for UC" only carries
  biological weight if de Lange isn't confounded. Without the
  intercept clear, the heatmap is a methods result, not a biology
  result.

Resolution path: LDSC pipeline. Not currently in `code/01_magma/` —
needs setup (LDSC python package + EUR LD scores) before the
intercept can land. **Do this before the M4 narrative write-up,
not after.** Liu UC + Yengo height should get the same treatment
once the pipeline is in place (Liu for the cross-GWAS sensitivity
narrative; Yengo as positive-control calibration).

- **Bites at:** the heatmap-to-narrative transition (M4). Pre-narrative
  gate, NOT pre-figure gate. Don't let M3 slide past M4 — if the
  intercept hasn't cleared before the manuscript draft, the result
  paragraph cannot be written without an explicit caveat.

## F9 · MT% pre-bank gate (must run before the first broad heatmap)
**Gate, not footnote.** Before treating the first 3×3 broad-tier
concordance figure as real, eyeball MT-hi cell counts + MT%
distribution against each deposit's own QC threshold, and confirm the
QC_STATE_TO_PARENT collapse-to-parent is scoring "stressed-but-viable"
cells — not dying cells the deposit silently shipped sub-threshold.
Cheap to check, expensive to miss: if a deposit shipped sub-threshold
cells in the MT clusters, collapse-to-parent quietly scores dying cells
into your lineages and the figure inherits the bias.

Per-atlas check:
- **Smillie SCP259**: `MT-hi` cluster cell count (~442 per the published
  metadata) and the MT% distribution within that cluster. SCP259's own
  QC gate is the viability anchor — confirm cluster MT% is below it.
- **Garrido GSE214695**: `MT T cells` and `MT fibroblasts` per-cluster
  cell counts and MT% distributions. Salas-lab gates these on MT%
  before annotation — confirm by spot-check.
- **TAURUS**: same check on whatever MT-flagged low-tier labels surface
  on first compute-node run.

Resolution path: a single notebook cell (per atlas) computing
`(obs['cell_type_fine'] == 'MT-hi').sum()` plus a histogram of MT% for
those cells against the deposit's published threshold. Log the
threshold and the spot-check result alongside the broad-tier figure as
a pre-bank certification. If sub-threshold cells DID slip through,
the collapse policy needs revisiting (likely outcome: drop those
specific cells, NOT change the cross-atlas rule).

- **Bites at:** broad-tier figure (step 06). PRE-BANK gate — do not
  bank the heatmap until this check passes.

> Note: gene-identifier harmonization is **not** a separate flag. Smillie ships
> HGNC symbols and Garrido/Pan-GI/HCA ship Ensembl, but every loader's final
> `ensembl_to_hgnc` step converges all atlases onto the NCBI-authoritative
> symbol set — that's the shared gene space. No symbol↔Ensembl map needed.

---

### Resolved (moved to DECISIONS.md)
- **F1** — UC tissue inflamed/non-inflamed/pooled → DECISIONS 26(b).
  Locked: pool for first heatmap; Phase-9 stratification sensitivity.
- **F2** — Garrido QC-state labels (MT / heat-shock / IER) → DECISIONS 22(a).
  Cross-atlas `QC_STATE_TO_PARENT` collapse policy in `_qc_policy.py`.
- **F3** — Garrido lineage-ambiguous cycling labels → DECISIONS 22(b).
  Cross-atlas `EXCLUDE_LINEAGE_AMBIGUOUS_FINE` rule in `_qc_policy.py`.
  **Resolved for the broad-tier figure only** — cells excluded so the
  harmonization-artifact failure mode is structurally prevented. The
  underlying lineage question (which lineage are these cycling cells?)
  is NOT closed; it is parked in the marker-QC step (MKI67 + lineage
  markers + cross-atlas placement agree). Reinstate cells per-cell, not
  per-label, when that pass runs.
- **F4** — Garrido crosswalk REVIEW rows (12) → DECISIONS 22(c) + 22(b).
  9 of 12 broad-tier-locked; 3 unprefixed cycling handled by F3 rule.
  Fine-tier identity questions flow to F8.
