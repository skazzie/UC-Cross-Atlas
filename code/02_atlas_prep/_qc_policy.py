"""Cross-atlas QC-state collapse + lineage-ambiguous exclusion rules.

Two rules, applied identically by every UC core loader (Garrido, Smillie,
TAURUS) so cross-atlas broad-tier concordance compares like with like and
doesn't manufacture discordance out of harmonization asymmetry. Both maps
live here (not in any one loader) precisely so the cross-atlas symmetry
is structural, not by-copy. The same drift risk that motivated
``_broad_vocab._BROAD_VOCAB`` (DECISIONS 20) applies here.

RULE 1 — ``QC_STATE_TO_PARENT``: collapse-to-parent for mitochondrial-high
(MT), heat-shock, and immediate-early-response (IER) labels. These are
real lineages stressed during dissociation; the original-author lineage
call is robust at the broad tier, so collapse-to-parent (not exclude)
keeps cells in the cohort and only removes QC noise from the fine tier.
Applied alongside ``RIBHI_TO_PARENT`` in each loader before
``FINE_TO_BROAD`` lookup. Atlases label stress states at different
granularity (Garrido breaks out MT-T / MT-fib / IER-fib / PC-heat-shock;
Smillie has only one ``MT-hi`` cluster in the Imm compartment); the
collapse harmonizes that asymmetry into the lineage tier everywhere.

RULE 2 — ``EXCLUDE_LINEAGE_AMBIGUOUS_FINE``: drop cycling clusters whose
fine label carries only proliferation markers (MKI67+ etc) and no stable
lineage call. Force-assigning these would create a cross-atlas discordance
that's pure harmonization artifact (e.g., Garrido→B cell vs Smillie→T
cell when each atlas's annotators happened to place them differently).
Symmetric exclusion is bounded — cycling-cell fractions are small in
every atlas — and biologically defensible: revisit in the marker-QC
step (MKI67 + lineage markers); reinstate per-cell only when markers and
cross-atlas placement agree.

References: DECISIONS (22) [this lock]; (20) [single-sourced vocab
pattern]; OPEN_FLAGS F2 (resolved here) / F3 (resolved here) / F4
(broad-tier portion resolved here).
"""

from __future__ import annotations


# RULE 1.
# Keys are raw fine labels as they appear in each atlas's native
# annotation (after _normalize_label). Values are the fine-tier parent
# the cell collapses into — that parent label MUST be a key in the
# loader's ``FINE_TO_BROAD`` (post-Ribhi-collapse universe), so the
# subsequent broad-tier lookup resolves.
QC_STATE_TO_PARENT: dict[str, str] = {
    # ---- Garrido-Trigo (Salas-lab 91-cluster annotation) ----
    # MT — mitochondrial-high; viability check is the deposit's own QC
    # filter (see DECISIONS 22 caveat).
    "MT T cells":                  "T",            # → "T cell" via FINE_TO_BROAD
    "MT fibroblasts":              "fibroblast",   # Ribhi-parent identity
    # IER — immediate early response (dissociation stress).
    "IER fibroblasts":             "fibroblast",
    # Heat shock — plasma cell stress states.
    "PC IgA heat shock 1":         "plasma cell",
    "PC IgA heat shock 2":         "plasma cell",
    "PC immediate early response": "plasma cell",
    # ---- Smillie (SCP259 51-cluster annotation) ----
    # ``MT-hi`` already maps directly to "T cell" at broad in the
    # Smillie loader's FINE_TO_BROAD; no fine-tier collapse needed
    # since it's a single-label, Imm-compartment-only cluster. Kept
    # out of this map by intent, not omission.
    # ---- TAURUS (Thomas/Dendrou 109 cell states) ----
    # Populated when the low-tier label set surfaces on first compute-node
    # run; any MT-/heat-shock-/IER-flavored low-tier label goes here.
}


# RULE 2.
# Fine labels (raw, post-_normalize_label) whose cells get dropped before
# the broad-tier mapping runs. Loaders log the drop count loudly so the
# exclusion is auditable.
EXCLUDE_LINEAGE_AMBIGUOUS_FINE: frozenset[str] = frozenset({
    # ---- Garrido-Trigo ----
    # Unprefixed cycling labels — no compartment tag, lineage call is
    # by-position in the Salas-lab clustering only. F3 in OPEN_FLAGS.
    "Cycling cells",
    "Cycling cells 2",
    "Cycling cells 3",
    # ---- Smillie ----
    # Smillie's cycling clusters are all compartment-prefixed
    # (``Cycling B``, ``Cycling T``, ``Cycling Monocytes``,
    # ``Cycling TA``) — lineage-tagged, kept by the rule. Empty by
    # intent.
    # ---- TAURUS ----
    # Populated when low-tier labels surface; any cycling cluster
    # without a lineage prefix goes here.
})
