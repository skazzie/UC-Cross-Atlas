"""Single source of the 15-term canonical broad vocabulary.

This module is imported by all UC core loaders (``load_smillie``,
``load_garrido_trigo``, ``load_taurus``) so the two-gate
canonical-vocab assertion (gate 1 at module load on ``FINE_TO_BROAD``
values, gate 2 at end-of-load on emitted ``obs['cell_type_broad']``)
measures *the same* vocab in every loader. Previously each loader
defined its own local frozenset; the three copies were identical but
duplicate-by-copy is exactly the drift risk this module exists to
prevent — if one loader's copy gets updated and the others don't,
``06_concordance``'s string-intersection over ``cell_type_broad``
would silently report vocab drift as biology.

Underscore prefix is intentional: this is **not** the locked public
``CANONICAL_BROAD`` of `canonical_broad_DRAFT.md`. It's the pre-lock
substrate used by atlas-prep only. When ``CANONICAL_BROAD`` locks
(with CL subtree IDs, deprecation pass, and PI v2 sign-off — see
DECISIONS 13 + 18(c)), the contents of this module promote to
``code/_shared/canonical_broad.py`` and become a public import for
downstream concordance code.

References: DECISIONS (13) [CL pin]; (16) [TAURUS swap]; (18)(c) [v2
escalation rolled back]; OPEN_FLAGS F5/F8.
"""

from __future__ import annotations

_BROAD_VOCAB: frozenset[str] = frozenset({
    "B cell",
    "NK/ILC",
    "T cell",
    "colonocyte",
    "dendritic cell",
    "endothelium",
    "enteroendocrine/tuft",
    "epithelial progenitor",
    "fibroblast",
    "goblet",
    "granulocyte",
    "mast cell",
    "monocyte/macrophage",
    "mural/glia",
    "plasma cell",
})
