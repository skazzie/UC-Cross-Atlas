# Canonical broad-tier vocabulary — DRAFT (not locked)

> **Status: DRAFT — for red-line, not yet authoritative.**
> No `CANONICAL_BROAD` constant in `_shared/constants.py` yet. No
> DECISIONS.md lock. No loader retargets. Concordance scripts must not
> import from this document until it is promoted to LOCKED — see
> "Sequencing for the lock" below.

This document is the candidate single shared broad-tier vocabulary
across the five UC-Cross-Atlas atlases (Smillie, Garrido-Trigo,
TAURUS, HCA Gut, Pan-GI). It exists so the per-term Cell Ontology
assignments and the two outstanding biology calls (mural/glia,
enteroendocrine/tuft) can be reviewed and red-lined **before** any
concordance metric is computed — locking after seeing rho values would
make the regime-1-vs-2 comparison circular (fix-handoff item 4, P1).

The draft is built from the only two atlases whose fine tiers are
currently in hand (Garrido-Trigo and Smillie). HCA Gut, Pan-GI, and
TAURUS will contribute net-new fine→canonical maps once their fine
labels are enumerated on disk; per-term coverage verification (≥50
cells per atlas where biologically expected) happens at that point.

References: sanity-check fix handoff (2026-06-04, item 4);
`OPEN_FLAGS.md` F5 (broad-tier vocabulary size — partial resolution
proposed here); DECISIONS corrections (1/7), (2/7), (3/7), (5/7), (9),
(10); existing per-atlas maps in
`code/02_atlas_prep/load_garrido_trigo.py:FINE_TO_BROAD` and
`code/02_atlas_prep/load_smillie.py:FINE_TO_BROAD`.

---

## CL-subtree framing (departure from "single CL node")

Each canonical term maps to **one or more named Cell Ontology subtrees**
(a CL node plus its descendants), documented per row below. The
"each term = one CL node" idealization is dropped: biology won't honor
it for biologically coherent unions, and forcing it would either inflate
cardinality past the 10–15 target or compromise the auditability claim.
The unions used here (NK/ILC, monocyte/macrophage) span distinct CL
nodes but one coherent lineage; each is documented as such, not hidden.

**CL version pin (DONE 2026-06-06, DECISIONS correction 13).** Cell
Ontology release **2026-03-26** (versionIRI
`http://purl.obolibrary.org/obo/cl/releases/2026-03-26/cl.owl`); the
labels for the CL IDs referenced below are committed at
`data/reference/cl_terms_pinned.tsv` and the OWL itself is NOT
committed (re-downloadable from
`http://purl.obolibrary.org/obo/cl.owl`). The TSV is the
single source of truth for what each CL ID *meant* at the time we
recorded it; refresh procedure documented in the TSV header.

**CL pin caught six errors in the draft — all resolved as of
2026-06-07 (DECISIONS 22).** The 2026-03-26 lookup against the IDs
originally drafted produced three label-drift renames (same concept,
name updated upstream) and three flat-out wrong IDs. Label drifts:
doc updated. Wrong-ID rows: ontology-investigated against the pinned
release via OLS and resolved — two replacement IDs verified-existing
in the 2026-03-26 release and appended to `cl_terms_pinned.tsv` as
same-release post-hoc additions (DECISIONS 22), one wrong ID dropped
in favor of an already-pinned correct sibling.

| Drafted ID | Drafted label | Actual label in 2026-03-26 | Resolution |
|---|---|---|---|
| CL:1000347 | enterocyte of colon | colonocyte | label drift — doc updated |
| CL:0002204 | brush cell of intestine / tuft | tuft cell | label drift — doc updated |
| CL:0002138 | lymphatic endothelial cell | endothelial cell of lymphatic vessel | label drift — doc updated |
| ~~CL:1000280~~ | stem cell of intestine | (smooth muscle cell of colon) | **dropped** — CL:0002250 *intestinal crypt stem cell* (already pinned) covers the intended concept |
| ~~CL:0009039~~ | colon epithelial progenitor cell | (colon goblet cell) | **replaced by CL:0009010** *transit amplifying cell* (post-hoc pin add) |
| ~~CL:0002073~~ | enteric glial cell | (transitional myocyte) | **replaced by CL:4040002** *enteroglial cell* (exact synonym: "enteric glial cell"; post-hoc pin add) |

CL:0009010 last-modified 2024-04-03 and CL:4040002 last-modified
2023-04-03 (per OLS); both **directly verified present** in the local
65.88 MB 2026-03-26 OWL on 2026-06-07 (`owl:Class rdf:about=...CL_4040002`
and `...CL_0009010` declarations; labels match the tsv). The post-hoc
additions are same-release — the pin's release-date provenance is
unchanged.

v2 polish (deferred): CL:0009043 *intestinal crypt stem cell of colon*
is a colon-specific subtype of CL:0002250 that exists in the release;
not added to the v1 vocab because the broad-tier rows are anchored on
the parent CL:0002250 already and adding the colon-specific child
doesn't change broad-tier mapping.

---

## Resolved sub-decisions (in this DRAFT, pending Muskaan red-line)

1. **mural/glia → split.** The previous Garrido/Smillie shared vocab
   bundled pericyte and glia into one `mural/glia` term. Pericyte
   (mesenchymal, CL:0000669) and enteric glia (neural-crest,
   CL:0002073) are biologically unrelated, and the cross-atlas
   structure agrees: HCA Gut splits them across two top-level
   `category` values (Mesenchymal vs Neuronal). Resolution in this
   draft:
   - **`fibroblast` now includes pericyte.** Pericyte folds into the
     mesenchymal stromal term.
   - **`glia` becomes its own canonical term.** Small but real;
     populated in Garrido (`Glia`), Smillie (`Glia`), HCA (`Neuronal`
     category subset).
2. **Pericyte rationale.** Even though pericyte has its own CL node
   (CL:0000669), bundling it under `fibroblast` here is defensible
   because (a) the cross-atlas resolution is mesenchymal-tier across
   all candidates and (b) splitting pericyte into its own canonical
   term would push cardinality past the 10–15 target without serving
   any downstream comparison (pericyte cell counts are small in every
   atlas).

These resolve part of `OPEN_FLAGS.md` F5 ("Optional: fold Pericyte →
Fibroblast and reconsider Glia"). The remainder of F5 (final
cardinality after all five atlases are on disk) is still open.

---

## Open biology call (this DRAFT does not resolve)

**`enteroendocrine/tuft` — split or keep?** The two are distinct
lineages (hormone-secreting vs chemosensory), both rare. The current
union sits under one CL parent in some taxonomies (intestinal
epithelial cell, CL:0002563) but not in CL itself. Lower-stakes than
the mural/glia call because both populations are small in every atlas;
the split mostly matters if downstream concordance treats secretory vs
chemosensory differently.

- **Keep as union** (`enteroendocrine/tuft`, cardinality 15): preserves
  the 10–15 budget; treats rare epithelial sensors as one cross-atlas
  category.
- **Split** (`enteroendocrine` + `tuft`, cardinality 16): one term over
  budget, but biologically cleaner; both terms may fail the
  ≥50-cell-per-atlas gate and end up as structural zeros in some
  atlases anyway.

This DRAFT recommends **keep as union** pending coverage data from
HCA/Pan-GI/TAURUS; flip to split if both can clear ≥50 cells per term
per atlas, otherwise the split adds a term that no atlas can populate
densely.

---

## Candidate vocabulary (15 terms)

| # | Canonical term | CL subtree(s) | Composition / mapping notes |
|---|---|---|---|
| 1 | `colonocyte` | CL:1000347 (pinned label: **colonocyte**); CL:0011108 (colon epithelial cell) is parent | Absorptive epithelium. **Garrido (6):** Colonocyte 1/2, Inflammatory colonocyte, BEST4 OTOP2, Laminin colonocytes, PLCG2 colonocytes. **Smillie (5):** Enterocytes, Best4+ Enterocytes, Immature Enterocytes 1/2 (F7 pending), M cells (F7 pending, microfold). |
| 2 | `goblet` | CL:0000160 (goblet cell); CL:1000320 (large intestine goblet) | Mucus-secreting epithelium. **Garrido (3):** Goblet, Mature goblet, Paneth-like (currently here pending F4). **Smillie (2):** Goblet, Immature Goblet. |
| 3 | `enteroendocrine/tuft` *(OPEN: split or keep)* | CL:0000164 (enteroendocrine cell); CL:0002204 (pinned label: **tuft cell**) | Distinct lineages bundled as one rare-epithelial-sensor category — biology call open. **Garrido (2):** Enteroendocrine, Tuft cells. **Smillie (2):** Enteroendocrine, Tuft. |
| 4 | `epithelial progenitor` | CL:0002250 (intestinal crypt stem cell); CL:0009010 (transit amplifying cell) | Stem + transit-amplifying + cycling epithelial. **Garrido (2):** Cycling TA, Secretory progenitor. **Smillie (6):** Stem, TA 1, TA 2, Cycling TA, Secretory TA, Enterocyte Progenitors. CL:0009010 *transit amplifying cell* is the generic TA term; replaces the wrong-ID draft slots ~~CL:1000280~~ / ~~CL:0009039~~ per DECISIONS 22. Same-release children for v2 polish: CL:0009043 *intestinal crypt stem cell of colon* (colon-specific child of CL:0002250); CL:4047017 *transit amplifying cell of gut* (gut-specific child of CL:0009010, added 2024-09-24, verified in the pinned OWL). Both subsumed by their parents under subtree semantics — not added here. |
| 5 | `fibroblast` *(includes pericyte — mural/glia split)* | CL:0000057 (fibroblast) + descendants; **+ CL:0000669 (pericyte)** | Mesenchymal stromal. **Garrido (11):** S1, S1.2, S2a, S2b, S3, IER fibroblasts, Inflammatory fibroblasts, MT fibroblasts, Myofibroblasts, FRCs, Perycites. **Smillie (9):** WNT2B+ Fos-lo 1/2, WNT2B+ Fos-hi, WNT5B+ 1/2, RSPO3+, Inflammatory Fibroblasts, Myofibroblasts, Pericytes. |
| 6 | `endothelium` | CL:0000115 (endothelial cell) + descendants (incl. CL:0002138, pinned label: **endothelial cell of lymphatic vessel**) | Vascular + lymphatic. **Garrido (3):** Endothelium, Activated endothelium, Lymphatic endothelium. **Smillie (3):** Endothelial, Microvascular, Post-capillary Venules. |
| 7 | `glia` *(new — mural/glia split)* | CL:4040002 (enteroglial cell — exact synonym "enteric glial cell"); CL:0000125 (glial cell) is the verified general parent | Neural-crest, distinct from mesenchymal pericyte. **Garrido (1):** Glia. **Smillie (1):** Glia. HCA: maps from `Neuronal` category subset — confirms cross-atlas resolution. CL:4040002 replaces the wrong-ID draft slot ~~CL:0002073~~ per DECISIONS 22; OLS lookup against the live ontology confirms exact-synonym "enteric glial cell" and parent of CL:4047047 (type I enteric glial cell), with last-modification 2023-04-03 — predates the 2026-03-26 pin. |
| 8 | `T cell` | CL:0000084 (T cell) + descendants | **Garrido (15):** CD4 ANXA1, CD4 naive, CD8 CTL, CD8 CTL TRM, CD8 FGFBP2, Cycling T cells, DN EOMES, DN TNF, MT T cells (F2 pending), S1PR1 T cells, T cells CCL20, ThF, Tregs, gd IEL, MAIT. **Smillie (10):** CD4+ Memory, CD4+ Activated Fos-hi/lo, CD4+ PD1+, CD8+ LP, CD8+ IELs, CD8+ IL17+, Tregs, Cycling T, MT-hi (F2/F7 pending). |
| 9 | `NK/ILC` | CL:0000623 (natural killer cell); CL:0001065 (innate lymphoid cell) | **Coherent union** — NK is group-1 ILC under current taxonomy. Cross-atlas grouping preserves the biology. **Garrido (2):** NK, ILC3. **Smillie (2):** NKs, ILCs. |
| 10 | `B cell` | CL:0000236 (B cell) + descendants | **Garrido (7):** B cell, Memory B cell, Naive B cell, GC B cell, Cycling cells, Cycling cells 2, Cycling cells 3 (F3 pending). **Smillie (3):** Follicular, GC, Cycling B. |
| 11 | `plasma cell` | CL:0000786 (plasma cell) + descendants | **Garrido (15):** PC IER, PC immediate early response, PC IGLL5, PC IgA 1/2/3/4, PC IgA IgM, PC IgA Lambda 1, PC IgA heat shock 1/2, PC IgG 1/2, Plasmablast IgA Lambda 2, Plasmablast IgG, Plasmablast IgG Lambda. **Smillie (1):** Plasma. |
| 12 | `monocyte/macrophage` | CL:0000576 (monocyte); CL:0000235 (macrophage) | **Coherent union** — monocyte→macrophage continuum is one mononuclear-phagocyte lineage in tissue. **Garrido (7):** M0, M1 ACOD1, M1 CXCL5, M2, M2.2, IDA macrophage, Inflammatory monocytes, Cycling myeloid. **Smillie (3):** Macrophages, Inflammatory Monocytes, Cycling Monocytes. |
| 13 | `dendritic cell` | CL:0000451 (dendritic cell) + descendants | **Garrido (2):** DCs CCL22, DCs CD1c. **Smillie (2):** DC1, DC2. |
| 14 | `mast cell` | CL:0000097 (mast cell) + descendants | **Garrido (2):** Mast 1, Mast 2. **Smillie (2):** CD69+ Mast, CD69- Mast. |
| 15 | `granulocyte` | CL:0000094 (granulocyte) + descendants (neutrophil CL:0000775, eosinophil CL:0000771) | **Smillie structural zero** (taxonomy lacks granulocytes). **Garrido (4):** Neutrophil 1/2/3, Eosinophils. Document the Smillie zero explicitly in the coverage report; do not collapse the vocab to accommodate it. |

**Cardinality**: 15 with `enteroendocrine/tuft` kept as union; 16 if
split. Sits at the top of the 10–15 target either way; the split would
push one over but is defensible biology if downstream concordance
actually needs it.

---

## What this DRAFT does NOT contain

Deliberately deferred until HCA Gut, Pan-GI, and TAURUS are on disk
with fine tiers enumerated:

- The HCA `author_cell_type` (~120 fine labels) → canonical map.
- The Pan-GI `level_3_annot` (~70 fine labels) → canonical map.
- The TAURUS fine-tier → canonical map (loader skeleton in place per
  DECISIONS 16; LOW_TO_BROAD empty, populates on first run).
- Per-term, per-atlas coverage table (cell counts; ≥50-cell gate
  flagged per term per atlas; structural zeros documented explicitly).
- Per-atlas "canonical term this atlas cannot populate" log lines in
  loaders.
- Loader retargets:
  - Garrido `FINE_TO_BROAD`: rename `mural/glia` → split into
    `fibroblast` (pericyte) + new `glia`; otherwise no value changes.
  - Smillie `FINE_TO_BROAD`: same `mural/glia` → split.
  - HCA: net-new `author_cell_type` → canonical map.
  - Pan-GI: net-new `level_3_annot` → canonical map.

---

## Sequencing for the lock

1. Process the three remaining atlases (Smillie compute-node run;
   Garrido RAW.tar rewrite per correction 9 [DONE]; TAURUS implementation;
   HCA + Pan-GI already production loaders, no rewrite needed).
2. Enumerate fine labels from HCA `author_cell_type` and Pan-GI
   `level_3_annot` against the candidate vocabulary; build the two
   net-new maps.
3. Compute per-term, per-atlas coverage (cell counts ≥ 50). Document
   structural zeros (granulocyte in Smillie, plus anything else that
   surfaces).
4. Pin the CL release used for the subtree IDs above; commit the
   pinned `cl.owl` or term list at `data/reference/`.
5. Muskaan + reviewer red-line this document: resolve EE/tuft
   split-vs-keep; confirm or revise mural/glia split; confirm CL
   subtree assignments; confirm the 15 (or 16) terms.
6. **Then** promote: write `CANONICAL_BROAD` (and
   `CANONICAL_BROAD_CL` subtree map) into
   `code/_shared/constants.py`; retarget all four loaders to map their
   fine labels into `CANONICAL_BROAD`; add the coverage assertion
   ("emitted broad ⊆ CANONICAL_BROAD; structural-zero log for missing
   terms"); log a DECISIONS correction marking the lock; **only then**
   run any `06_concordance/` script.

The natural sequencing protects the one-way door: concordance can't run
until all atlases are processed, which can't happen until after the
lock — so the lock-before-concordance discipline holds without extra
gating.
