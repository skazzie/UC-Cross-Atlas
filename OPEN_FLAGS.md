# Open Flags — deferred decisions to commit to later

Running list of judgment calls and unresolved questions surfaced during the build but
not yet decided. The companion to DECISIONS.md: that file is *locked* decisions, this is
the *still-open* side. When a flag is resolved, log it in DECISIONS.md and strike it here.

_Last updated: June 4, 2026._

---

## F1 · UC tissue definition — inflamed vs non-inflamed
Smillie's `Health` field has **three** states: Healthy / Non-inflamed / Inflamed — UC
patients contribute both inflamed and non-inflamed biopsies. Decide whether the UC group
is inflamed-only, non-inflamed-only, or pooled — and make it **consistent across all
atlases** (confirm what Garrido's and Mennillo's "UC" tissue actually is). A mismatch
(e.g. Garrido = inflamed biopsies, Smillie = pooled) would confound concordance as
biology, not method.
- **Smillie structure (confirmed):** 12 healthy controls + 18 UC patients, paired —
  every UC patient contributed both an inflamed and a non-inflamed biopsy. Cells:
  110,110 Healthy / 125,119 Inflamed / 130,263 Non-inflamed (~366k total). UC-as-inflamed
  = 125k cells; UC-pooled = 255k from the same 18 patients (each present as 2 samples).
- **Bites at:** harmonization + concordance. The loader must *capture* `Health` now; the
  grouping decision is needed before the UC group is subset.

## F2 · Garrido QC-state labels (MT / heat-shock / IER)
Beyond Ribhi, Garrido has a second QC-state axis: `MT T cells`, `MT fibroblasts`
(mitochondrial-high), plus `PC IgA heat shock 1/2`, `PC immediate early response`,
`IER fibroblasts` (dissociation stress). 6 labels, cross-lineage, not true cell types.
Decide: extend the Ribhi collapse policy — marker-check, then collapse-to-parent or
exclude. **Currently unhandled by the Garrido loader.**
- **Smillie has the same axis:** `MT-hi` (442 cells) is the mitochondrial-high
  analogue. `load_smillie.py` maps it to `T cell` provisionally; the
  collapse-vs-exclude policy should cover both atlases uniformly.
- **Bites at:** crosswalk / loader (M2).

## F3 · Garrido lineage-ambiguous cycling labels
`Cycling cells` and `Cycling cells 3` have no resolvable lineage, so they can't harmonize
to any other atlas. Decide: marker-based reassignment or exclude before concordance.
- **Bites at:** crosswalk (M2).

## F4 · Garrido crosswalk REVIEW rows (12)
Biological calls for Muskaan, flagged REVIEW in `harmonization_crosswalk.csv`: Laminin /
PLCG2 colonocytes distinct or fold into Colonocyte; CD4 ANXA1 identity; S1PR1 T and
CCL20 T lineage; PC IGLL5 isotype grouping; merge Mature goblet into Goblet; etc.
- **Bites at:** crosswalk (M2).

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

## F7 · Smillie crosswalk REVIEW rows
Tentative `FINE_TO_BROAD` calls in `load_smillie.py` (marked `# REVIEW`), to
confirm with markers: `M cells` (microfold — mapped to colonocyte, no clean
broad home); `Immature Enterocytes 1/2` (→ colonocyte, progenitor/colonocyte
boundary); `Secretory TA` (→ epithelial progenitor, vs secretory/goblet
lineage); `MT-hi` (→ T cell, QC state — see F2). The loader runs with these in
place; they are not biology-locked.
- **Bites at:** crosswalk (M2).

> Note: gene-identifier harmonization is **not** a separate flag. Smillie ships
> HGNC symbols and Garrido/Pan-GI/HCA ship Ensembl, but every loader's final
> `ensembl_to_hgnc` step converges all atlases onto the NCBI-authoritative
> symbol set — that's the shared gene space. No symbol↔Ensembl map needed.

---

### Resolved (moved to DECISIONS.md)
_(none yet)_.
