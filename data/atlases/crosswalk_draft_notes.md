# Crosswalk draft notes — companion to `celltype_crosswalk.tsv`

## Status: DRAFT (paper-derived, NOT obs-validated)

The TAURUS, HCA Gut, and Pan-GI rows in `celltype_crosswalk.tsv` are
**paper-derived drafts** as of 2026-06-27. **They have not been
validated against the actual `.obs` of the deposited h5ads** — what's
in a figure legend and what's in `.obs` routinely differ (the rule
set in Phase 3 was: map from deposit labels, not paper labels;
papers are the pre-work, deposits are the truth). The drafts exist
to front-load the boring part and surface granularity mismatches
before HB scoring runs at the wrong granularity, NOT to replace the
first-HB-load reconciliation.

**Anyone reading this in a future session: do NOT treat these 406
TSV rows as a done crosswalk.** ~177 native fine labels across the
three atlases are still flagged `<pending-HB-unseen>` in the TSV;
the unseen-count drops as obs reconciliation lands per atlas. The
drafted-but-seen labels still require obs verification because
spelling variants, renames, and silent post-publication edits all
happen between paper and deposit.

The build-script status output reports `DRAFT (paper-supplement;
verify against obs)` for any atlas whose rows include
`paper-supplement` sources, distinct from `complete (loader-extracted)`
for Garrido and Smillie.

---

## Source

Generated 2026-06-27 alongside `scripts/build_celltype_crosswalk.py`
draft-block additions. Covers TAURUS, HCA Gut, Pan-GI — the three atlases
whose obs lives on Hummingbird only and whose crosswalk rows are drafted
from the published papers / supplementary tables / portal landing pages
rather than from h5ad obs.

The TSV holds the machine-readable mapping. This file holds the
human-readable narrative: source URLs, granularity flags, uncertain rows,
and per-bucket structural-zero rationale.

**Hard rule (per DECISIONS 31 + Saisohan 2026-06-24 review):** every
drafted row is marked `source = paper-supplement` and notes begin with
`VERIFY ON LOAD:`. The h5ad obs is authoritative if it disagrees with the
paper. A draft-to-verify is strictly safer than a blank-to-fill, and it
catches granularity mismatches BEFORE the cluster operator scores at the
wrong granularity — which is the whole point of doing the crosswalk early.

---

## TAURUS-IBD (Thomas et al., Nat Immunol 25:2152-2165, 2024)

**Sources:**
- Paper full-text: https://pmc.ncbi.nlm.nih.gov/articles/PMC11519010/
- Zenodo v3 deposit (h5ad on HB): https://zenodo.org/records/14007626
- DOI: 10.1038/s41590-024-01994-8

**Coverage:** ~45 of 109 paper-reported `cell_state` labels drafted
(~58% unseen at draft time). The full 109-class catalog lives in
Supplementary Table 1 (xlsx); MOESM endpoints were unreachable from the
laptop fetch. Recommend dumping `adata.obs['cell_state'].unique()` on HB
as ground truth and extending `_TAURUS_DRAFT_FINE` accordingly.

**Compartment count (paper text):** 9 compartments — Epithelial (ileal),
Epithelial (colonic), Stromal, CD4 T, CD8/innate-T/NK/ILC, B cell, Plasma
cell, Myeloid, and a 9th that did not surface verbatim in fetched text
(likely Mast/granulocyte split or GC-B-as-own). **Verify on HB.**

**Granularity flags:**
- 14 T-cell-compartment `cell_state` labels collapse into 1 canonical
  broad bucket. F8 fine vocab has only `t_cell_broad` for the T-cell
  branch — none of Th/Tfh/Tph/Treg/CD8-GZMK*/MAIT/gd are individually
  cross-atlas alignable. Expected per F8 T2.5 intent.
- 9 myeloid `cell_state` labels collapse to `monocyte/macrophage` +
  `dendritic cell` — no separate `inflammatory_macrophage` fine bucket
  exists in F8 vocab; granularity loss.

**Structural zero (marked in TSV):**
- `granulocyte` — biopsy + droplet 10x lyses granulocytes (same caveat
  as Smillie).

**Uncertain rows (marked as `paper-supplement-uncertain`):**
- `mast cell` — no mast label surfaced in ED Fig 1; may sit inside the
  myeloid compartment under a label not extracted.
- `mural/glia` — no enteric glia label surfaced. The "Stromal"
  compartment in ED Fig 1h appears to contain only Fibroblast / THY1+
  FAP+ PDPN+ / Pericyte / Vascular. Enteric glia may be structurally
  absent (mucosal biopsies do not sample submucosa / muscularis).
- `lymphatic_endothelial` fine bucket — "Vascular" appears to lump
  LEC + BEC; lymphatic subset likely not recoverable from
  `cell_state` alone.

**Naming caveats (per row):**
- `Enterocyte` → `canonical_broad = colonocyte` because the canonical
  vocab reuses `colonocyte` as the absorptive bucket (DECISIONS 22,
  canonical_broad_DRAFT.md row 1). SI absorptive in TAURUS is not
  literally "colon-ocyte" — collapse is intentional.
- `Vascular` → `endothelium` — confirm marker gene set on HB; if it
  bundles pericyte-adjacent, may need reassignment.

---

## HCA Gut (Elmentaite et al., Nature 597:250, 2021)

**Sources:**
- Paper full-text: https://pmc.ncbi.nlm.nih.gov/articles/PMC8426186/
- Nature DOI: https://www.nature.com/articles/s41586-021-03852-1
- CELLxGENE deposit (h5ad on HB):
  `f34d2b82-9265-4a73-bda4-852933bf2a8d`
  (portal facets were unreachable from laptop fetch)

**Coverage:** ~90 of ~120 paper-reported `author_cell_type` labels
drafted (~25% unseen at draft time). PMC excerpt was thin on the myeloid
compartment — recommend dumping
`adata.obs['author_cell_type'].unique()` on HB and extending
`_HCA_DRAFT_FINE` for any unseen myeloid subtypes.

**Region + age filter (CRITICAL):** HCA Gut spans fetal → adult and 11
anatomical regions. ~15-20 native labels are developmental-only
(Mesoderm 1/2, Pro-B, Pre-B, CLP, NEUROG3+ progenitors, mLTo, Branch A/B
enteric neurons) — `load_hca_gut.py` already filters to `Age_group ∈
{Adult, Adult_MLN}` and colonic `tissue` values before annotation, so
these developmental labels SHOULD drop out before crosswalk lookup. The
draft block retains them so every drafted label is captured (NOT a
claim that the crosswalk is done — see "Status" at the top of this
file); the filter prevents their cells from polluting broad-bucket
counts.

**Granularity flags:**
- 6 LEC subtypes (LEC1-6) collapse to `lymphatic_endothelial` — fine
  bucket holds.
- 4 BEC subtypes collapse to `endothelial`.
- 8 enteric neuron Branch A/B labels + 4 glia labels all collapse to
  `mural/glia` + `enteroglial`. Significant developmental granularity
  loss but expected.
- 5 pericyte subtypes collapse into `fibroblast` per F8 design.
- 8 T-cell native labels collapse to `T cell` + `t_cell_broad`.
- Mesothelial cells (2 subtypes) are forced into `fibroblast` — no
  mesothelium bucket in 15-term vocab. Flag in any downstream paragraph
  that references HCA stromal counts.

**Structural extras (excluded from broad collapse):**
- `Erythroid` → `canonical_broad = (exclude)` — developmental, not in
  15-term vocab. Drop in QC; do NOT count as a structural zero.

**Uncertain buckets (marked `paper-supplement-uncertain`):**
`mast cell`, `monocyte/macrophage`, `dendritic cell`, `granulocyte` —
the myeloid compartment extraction was thin in the PMC excerpt; these
buckets are very likely populated by labels not in the draft block.
First HB obs enumeration should confirm and the draft block extended.

**Naming caveats (per row):**
- `Follicular dendritic cells (FDCs)` → `fibroblast` (stromal lineage,
  NOT myeloid DC despite the name). Easy confusion; called out in row
  notes.
- `Stromal 3` → `inflammatory_fibroblast` — Smillie's nomenclature
  treats S3 as inflammation-associated. HCA Gut S3 marker overlap is
  not confirmed; row note marks for verification.
- `Microfold (M) cells` → `colonocyte` broad as structural default; M
  cells are specialized antigen-sampling epithelium, not absorptive
  proper. No F8 bucket for them.

---

## Pan-GI Extended+ (Oliver et al., Nature 635:699, 2024)

**Sources:**
- Paper full-text: https://pmc.ncbi.nlm.nih.gov/articles/PMC11578898/
- Nature DOI: https://www.nature.com/articles/s41586-024-07571-1
- CELLxGENE deposit (h5ad on HB):
  `1dcf15ee-c103-4aaa-8b8c-0fc697fcccc8`
  (portal facets were unreachable from laptop fetch)

**Coverage:** ~45 of 136 paper-reported `level_3_annot` labels drafted
(~67% unseen at draft time). The fetched PMC text covered the
distinctive labels (INFLAREs, MGN, metaplastic Paneth, inflammatory
fibroblasts, the major immune subsets) but not the full enumeration —
many of the unseen ~91 are regional epithelial subtypes from non-colonic
GI tissue. Recommend dumping `adata.obs['level_3_annot'].unique()` on HB
and extending `_PANGI_DRAFT_FINE`.

**Region filter (CRITICAL):** Pan-GI is whole-GI (oral mucosa → stomach
→ SI → colon). ~30-40% of `level_3_annot` labels are non-colonic
(oral fibroblast, oesophagus fibroblast, gastric MGN, Brunner's gland,
surface foveolar, duodenal Paneth, etc.) and would pollute broad-bucket
counts. `load_pangi.py` already filters on `organ_unified` ∈ colonic
regions and `disease ∈ {normal, UC, IBD}` — these non-colonic labels
should drop before crosswalk lookup; the draft block retains them for
completeness.

**Granularity flags:**
- 17 T/NK fine subsets and 16 myeloid subsets (paper-reported) collapse
  to 4 broad buckets — same pattern as TAURUS/HCA.
- Pan-GI is the **only granulocyte anchor** across all five atlases.
  `Neutrophils` is the single source feeding `granulocyte` cross-atlas;
  Smillie / TAURUS / Garrido lose granulocytes to droplet-10x lysis.
  Pan-GI's granulocyte row is the cross-atlas pivot.
- Pan-GI's `Inflammatory fibroblasts` and `Oral mucosa fibroblasts`
  share the IBD-inflammatory signature per the paper — both map to
  `inflammatory_fibroblast`. Provides a clean cross-atlas anchor with
  Smillie's `Inflammatory Fibroblasts` and Garrido's
  `Inflammatory fibroblasts`.
- INFLAREs are NOVEL to this paper and the key UC-relevant finding. No
  direct counterpart in Smillie / Garrido / TAURUS / HCA. Tentative
  mapping to `paneth_like` is biologically reasonable (MUC6+, secretory,
  metaplastic) but the bucket will be Pan-GI-singleton at cross-atlas.

**Uncertain rows:**
- `INFLAREs` — paneth_like is the closest F8 bucket; biology is
  metaplastic foveolar / Brunner's-like, not true Paneth.
- `Surface foveolar cells` / `Surface foveolar-like cells` — gastric in
  healthy, appears in colon as metaplasia. `colonocyte` default but they
  are not absorptive; could argue `goblet` (secretory).
- `MGN cells` — same metaplastic-secretory ambiguity (paneth_like vs
  goblet).
- `Smooth muscle cells` — `mural/glia` broad is the only F8 home but
  it's a forced fit.
- `Deep crypt secretory (DCS)` — colonic functional equivalent of
  Paneth; `goblet` default with `paneth_like` as alternative.

**Integration caveat:** Pan-GI integrates Smillie + Kong donors per
DECISIONS (3/7). If the integrated h5ad re-labels original Smillie/Kong
donors with Pan-GI `level_3_annot`, treat as Pan-GI labels only — do
NOT double-count those donors as both atlases. `load_pangi.py` already
scans for the Smillie donor-ID regex (`SMILLIE_DONOR_REGEX`); a non-zero
overlap match should trigger the donor-overlap policy.

---

## Cross-atlas summary (likely structural-zero / single-source buckets)

| Bucket | Garrido | Smillie | TAURUS | HCA Gut | Pan-GI |
|---|---|---|---|---|---|
| `granulocyte` | yes (4 native) | structural zero | structural zero (biopsy + 10x) | uncertain (myeloid thin) | yes — **only anchor** |
| `mast cell` | yes (2 native) | yes (2 native) | uncertain | uncertain | yes |
| `mural/glia` (enteroglial) | yes (1 native) | yes (1 native) | uncertain (biopsy may not sample) | well-populated (developmental ENS focus) | yes |
| `lymphatic_endothelial` | yes (1 native) | (boundary case) | likely lumped into "Vascular" | yes (6 native: LEC1-6) | yes |
| `paneth_like` | study-specific (Paneth-like) | — | — | yes (Paneth cells) | yes (Paneth + Metaplastic Paneth + INFLAREs + MGN + DCS) |

**Implication for cross-atlas concordance:** granulocyte and lymphatic
endothelial buckets will be single-or-double-source rows in the 5×5
broad heatmap. The figure should mark these explicitly as
single-source rather than treat the missing rows as "low concordance."
This is a cross-atlas figure-construction note, not a flag — it
follows directly from the structural-zero rows already in the TSV.

---

## How to extend this on HB

When the first HB load enumerates obs for each atlas:

1. Dump the unique label set: `adata.obs['<fine_col>'].unique()`.
2. Diff against the drafted labels in `_TAURUS_DRAFT_FINE` /
   `_HCA_DRAFT_FINE` / `_PANGI_DRAFT_FINE`. Three buckets:
   - Drafted-and-confirmed: native label in obs matches the drafted
     spelling. Note row stays as-is; the `VERIFY ON LOAD:` note can
     be relaxed by removing the prefix.
   - Drafted-but-renamed: native label in obs is a spelling variant
     of a drafted label (e.g. `Plasma cell` vs `Plasma cells`). Update
     the drafted spelling; flag the rename in the row note.
   - Unseen: native label in obs that is NOT in the drafted block. Add
     a new row to the draft block with broad/fine mapping; remove from
     the `<pending-HB-unseen>` count.
3. Decrement the `_*_UNSEEN_FINE_EST` constant. When the count hits 0
   the residual skeleton row can be dropped (delete the
   `_build_paper_supplement_rows`'s residual-row emit for that atlas).
4. Convert `paper-supplement-uncertain` broad rows to
   `paper-supplement` or `paper-supplement-structural-zero` based on
   what obs shows.
5. Rebuild: `py.exe scripts/build_celltype_crosswalk.py`.
6. Log the update as a DECISIONS correction (one per atlas
   reconciliation pass) so the draft-to-confirmed transition is
   auditable.
