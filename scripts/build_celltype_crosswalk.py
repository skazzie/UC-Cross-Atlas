"""Build data/atlases/celltype_crosswalk.tsv — the cross-atlas label
harmonization file.

Schema (tab-separated, one row per (atlas, native_label, tier)):

    atlas              one of {garrido, smillie, taurus, hca_gut, pangi}
    tier               broad | fine
    native_label       the label as it appears in the atlas's obs column
    native_col         the obs column name in the atlas
    post_qc_label      the label AFTER _qc_policy.py (Ribhi/QC collapse,
                       cycling exclusion); == native_label for atlases
                       that don't carry the QC-state axis
    canonical_broad    one of the 15 _BROAD_VOCAB terms; the row's
                       canonical-broad bucket (always non-empty for
                       broad-tier rows; populated for fine rows that
                       have a canonical broad ancestor)
    canonical_fine     one of the 14 candidate fine buckets from
                       OPEN_FLAGS F8 T2.5 (or empty if the native
                       label doesn't fit any of those 14)
    cl_anchor          CL term for the canonical_broad (per the pinned
                       2026-03-26 release; see canonical_broad_DRAFT.md)
    source             how this row was derived (loader-extracted,
                       paper-supplement, cellxgene-metadata, pending-HB)
    notes              free text (collapse target, exclusion reason,
                       study-specific tag, etc.)

This file IS the single source of truth for cross-atlas label
comparability. Anything that wants to compare cell-type-level
statistics across atlases must do so through the canonical_broad or
canonical_fine columns, not native_label.

Versioning: written with the script's commit SHA stamped in a header
line so a particular crosswalk version is reproducible.
"""

from __future__ import annotations

import argparse
import ast
import csv
import datetime as dt
import subprocess
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]


# ---- CL anchors per _BROAD_VOCAB term (from canonical_broad_DRAFT.md) -----
CL_ANCHORS: dict[str, str] = {
    "colonocyte":             "CL:1000347",
    "epithelial progenitor":  "CL:0002250,CL:0009010",
    "goblet":                 "CL:0000160",
    "enteroendocrine/tuft":   "CL:0000164,CL:0002204",
    "fibroblast":             "CL:0000057,CL:0000669",  # incl. pericyte
    "endothelium":            "CL:0000115,CL:0002138",
    "mural/glia":             "CL:4040002,CL:0000125",  # enteroglial cell (post DECISIONS 22)
    "T cell":                 "CL:0000084",
    "NK/ILC":                 "CL:0000623,CL:0001065",
    "B cell":                 "CL:0000236",
    "plasma cell":            "CL:0000786",
    "monocyte/macrophage":    "CL:0000576,CL:0000235",
    "dendritic cell":         "CL:0000451",
    "mast cell":              "CL:0000097",
    "granulocyte":            "CL:0000094,CL:0000775,CL:0000771",
}


# ---- Candidate canonical_fine buckets per OPEN_FLAGS F8 T2.5 --------------
# Maps a native (post-QC) Garrido/Smillie label to its canonical_fine
# bucket if it fits the 14-bucket F8 sketch. Native labels NOT in this
# map have canonical_fine = "" (study-specific subtype or unaligned).
CANONICAL_FINE_BUCKETS: dict[str, str] = {
    # absorptive enterocyte
    "Enterocytes":          "absorptive_enterocyte",
    "Best4+ Enterocytes":   "absorptive_enterocyte",
    "Colonocyte 1":         "absorptive_enterocyte",
    "Colonocyte 2":         "absorptive_enterocyte",
    "BEST4 OTOP2":          "absorptive_enterocyte",
    # immature enterocyte / TA-adjacent
    "Immature Enterocytes 1": "immature_enterocyte_ta",
    "Immature Enterocytes 2": "immature_enterocyte_ta",
    "Enterocyte Progenitors": "immature_enterocyte_ta",
    "Cycling TA":             "immature_enterocyte_ta",
    "Secretory progenitor":   "immature_enterocyte_ta",
    # crypt stem
    "Stem":                 "crypt_stem",
    # goblet
    "Goblet":               "goblet",
    "Immature Goblet":      "goblet",
    "Mature goblet":        "goblet",
    # Paneth-like (study-specific)
    "Paneth-like":          "paneth_like",
    # enteroendocrine
    "Enteroendocrine":      "enteroendocrine",
    # tuft
    "Tuft":                 "tuft",
    "Tuft cells":           "tuft",
    # fibroblast (stromal)
    "S1":                   "fibroblast_stromal",
    "S1.2":                 "fibroblast_stromal",
    "S2a":                  "fibroblast_stromal",
    "S2b":                  "fibroblast_stromal",
    "S3":                   "fibroblast_stromal",
    "WNT2B+ Fos-lo 1":      "fibroblast_stromal",
    "WNT2B+ Fos-lo 2":      "fibroblast_stromal",
    "WNT2B+ Fos-hi":        "fibroblast_stromal",
    "WNT5B+ 1":             "fibroblast_stromal",
    "WNT5B+ 2":             "fibroblast_stromal",
    "RSPO3+":               "fibroblast_stromal",
    # inflammatory fibroblast
    "Inflammatory fibroblasts": "inflammatory_fibroblast",
    "Inflammatory Fibroblasts": "inflammatory_fibroblast",
    # endothelial
    "Endothelium":          "endothelial",
    "Activated endothelium":"endothelial",
    "Endothelial":          "endothelial",
    "Microvascular":        "endothelial",
    # lymphatic endothelial
    "Lymphatic endothelium":"lymphatic_endothelial",
    # T cell (broad lymphoid placeholder)
    "CD4 ANXA1":            "t_cell_broad",
    "CD4 naive":            "t_cell_broad",
    "CD8 CTL":              "t_cell_broad",
    "CD8 CTL TRM":          "t_cell_broad",
    "CD8 FGFBP2":           "t_cell_broad",
    "CD4+ Memory":          "t_cell_broad",
    "CD4+ Activated Fos-hi":"t_cell_broad",
    "CD4+ Activated Fos-lo":"t_cell_broad",
    "CD4+ PD1+":            "t_cell_broad",
    "CD8+ LP":              "t_cell_broad",
    "CD8+ IELs":            "t_cell_broad",
    "CD8+ IL17+":           "t_cell_broad",
    "Tregs":                "t_cell_broad",
    "ThF":                  "t_cell_broad",
    "MAIT":                 "t_cell_broad",
    "gd IEL":               "t_cell_broad",
    "DN EOMES":             "t_cell_broad",
    "DN TNF":               "t_cell_broad",
    "S1PR1 T cells":        "t_cell_broad",
    "T cells CCL20":        "t_cell_broad",
    # plasma cell (isotype-collapsed)
    "Plasma":               "plasma_cell_collapsed",
    "PC IER":               "plasma_cell_collapsed",
    "PC IGLL5":             "plasma_cell_collapsed",
    "PC IgA 1":             "plasma_cell_collapsed",
    "PC IgA 2":             "plasma_cell_collapsed",
    "PC IgA 3":             "plasma_cell_collapsed",
    "PC IgA 4":             "plasma_cell_collapsed",
    "PC IgA IgM":           "plasma_cell_collapsed",
    "PC IgA Lambda 1":      "plasma_cell_collapsed",
    "PC IgG 1":             "plasma_cell_collapsed",
    "PC IgG 2":             "plasma_cell_collapsed",
    "Plasmablast IgA Lambda 2":"plasma_cell_collapsed",
    "Plasmablast IgG":      "plasma_cell_collapsed",
    "Plasmablast IgG Lambda":"plasma_cell_collapsed",
    # enteroglial
    "Glia":                 "enteroglial",
}


def _extract_fine_to_broad(loader_path: Path, dict_name: str) -> dict[str, str]:
    """Pull a FINE_TO_BROAD-style dict out of a loader file by AST walk."""
    src = loader_path.read_text(encoding="utf-8")
    mod = ast.parse(src)
    for node in ast.walk(mod):
        if isinstance(node, ast.Assign):
            for t in node.targets:
                if isinstance(t, ast.Name) and t.id == dict_name:
                    return ast.literal_eval(node.value)
    raise LookupError(f"{dict_name} not found in {loader_path}")


def _build_garrido_rows() -> list[dict]:
    """Garrido: 86 native fine -> 15 broad (post Ribhi-collapse + QC) +
    the 3 unprefixed Cycling labels excluded per DECISIONS 22.
    Source: code/02_atlas_prep/load_garrido_trigo.py FINE_TO_BROAD."""
    rows: list[dict] = []
    fine_to_broad = _extract_fine_to_broad(
        _REPO / "code/02_atlas_prep/load_garrido_trigo.py", "FINE_TO_BROAD"
    )
    # Native fine labels (the post-Ribhi-collapsed universe used by the loader).
    for native, broad in sorted(fine_to_broad.items()):
        # Skip the Ribhi-parent identity rows (epithelial, T, fibroblast,
        # mast, plasma cell) — these aren't true native labels, they're
        # collapse targets used internally. Mark them in notes.
        is_ribhi_parent = native in {"epithelial", "T", "fibroblast", "mast",
                                      "plasma cell"}
        rows.append({
            "atlas": "garrido",
            "tier": "fine",
            "native_label": native,
            "native_col": "annotation (GSE214695_cell_annotation.csv)",
            "post_qc_label": native,
            "canonical_broad": broad,
            "canonical_fine": CANONICAL_FINE_BUCKETS.get(native, ""),
            "cl_anchor": CL_ANCHORS.get(broad, ""),
            "source": "loader-extracted",
            "notes": ("Ribhi-parent identity row, target of RIBHI_TO_PARENT collapse"
                      if is_ribhi_parent else ""),
        })
    # Excluded unprefixed cycling labels (DECISIONS 22).
    for native in ("Cycling cells", "Cycling cells 2", "Cycling cells 3"):
        rows.append({
            "atlas": "garrido",
            "tier": "fine",
            "native_label": native,
            "native_col": "annotation (GSE214695_cell_annotation.csv)",
            "post_qc_label": "<EXCLUDED>",
            "canonical_broad": "<EXCLUDED>",
            "canonical_fine": "<EXCLUDED>",
            "cl_anchor": "",
            "source": "loader-extracted",
            "notes": "Lineage-ambiguous cycling cluster; dropped per DECISIONS 22(b) EXCLUDE_LINEAGE_AMBIGUOUS_FINE",
        })
    # Broad-tier rows: one per canonical broad term present in Garrido.
    for broad in sorted(set(fine_to_broad.values())):
        rows.append({
            "atlas": "garrido",
            "tier": "broad",
            "native_label": broad,
            "native_col": "cell_type_broad (loader-emitted)",
            "post_qc_label": broad,
            "canonical_broad": broad,
            "canonical_fine": "",
            "cl_anchor": CL_ANCHORS.get(broad, ""),
            "source": "loader-extracted",
            "notes": "",
        })
    return rows


def _build_smillie_rows() -> list[dict]:
    """Smillie: 51 native Cluster labels -> 14 broad (no granulocyte).
    Source: code/02_atlas_prep/load_smillie.py FINE_TO_BROAD."""
    rows: list[dict] = []
    fine_to_broad = _extract_fine_to_broad(
        _REPO / "code/02_atlas_prep/load_smillie.py", "FINE_TO_BROAD"
    )
    for native, broad in sorted(fine_to_broad.items()):
        notes = ""
        if native == "MT-hi":
            notes = "QC stress label; broad-mapped here, not collapsed via QC_STATE_TO_PARENT (single Imm-compartment cluster, see DECISIONS 22)"
        elif native in ("Cycling B", "Cycling T", "Cycling TA",
                        "Cycling Monocytes"):
            notes = "Lineage-prefixed cycling; KEPT per DECISIONS 22(b) (vs unprefixed-Cycling exclusion)"
        rows.append({
            "atlas": "smillie",
            "tier": "fine",
            "native_label": native,
            "native_col": "Cluster (all.meta2.txt)",
            "post_qc_label": native,
            "canonical_broad": broad,
            "canonical_fine": CANONICAL_FINE_BUCKETS.get(native, ""),
            "cl_anchor": CL_ANCHORS.get(broad, ""),
            "source": "loader-extracted",
            "notes": notes,
        })
    for broad in sorted(set(fine_to_broad.values())):
        rows.append({
            "atlas": "smillie",
            "tier": "broad",
            "native_label": broad,
            "native_col": "cell_type_broad (loader-emitted)",
            "post_qc_label": broad,
            "canonical_broad": broad,
            "canonical_fine": "",
            "cl_anchor": CL_ANCHORS.get(broad, ""),
            "source": "loader-extracted",
            "notes": "",
        })
    return rows


def _build_pending_hb_rows(
    atlas: str, fine_col: str, broad_col: str,
    expected_n_fine_est: str, paper_ref: str, notes: str,
) -> list[dict]:
    """Skeleton rows for atlases that live on HB only. Native labels
    are <pending-HB>; one row per canonical_broad documenting that
    the atlas IS expected to populate that bucket (or NOT, with
    structural-zero noted)."""
    rows: list[dict] = []
    # Skeleton fine row, listing the expected fine cardinality + source col.
    rows.append({
        "atlas": atlas,
        "tier": "fine",
        "native_label": f"<pending-HB: {expected_n_fine_est} labels>",
        "native_col": fine_col,
        "post_qc_label": "<pending-HB>",
        "canonical_broad": "<pending-HB>",
        "canonical_fine": "<pending-HB>",
        "cl_anchor": "",
        "source": "pending-HB",
        "notes": (f"{paper_ref}. {notes} Populate after first HB load by "
                  f"enumerating obs['{fine_col}'] and mapping each label "
                  "into _BROAD_VOCAB and canonical_fine buckets per "
                  "OPEN_FLAGS F8 T2.5."),
    })
    # Per-canonical-broad skeleton: 15 rows expected to be populated.
    for broad in sorted(CL_ANCHORS):
        rows.append({
            "atlas": atlas,
            "tier": "broad",
            "native_label": f"<pending-HB-mapping>",
            "native_col": broad_col,
            "post_qc_label": broad,
            "canonical_broad": broad,
            "canonical_fine": "",
            "cl_anchor": CL_ANCHORS[broad],
            "source": "pending-HB",
            "notes": notes,
        })
    return rows


# ---- Paper-supplement draft mappings (2026-06-27) ------------------------
# Each draft below is sourced from the published paper / supplementary
# tables / portal metadata for the corresponding atlas, NOT from the h5ad
# obs (which lives on HB only). Every drafted row is marked
# source='paper-supplement' and its notes column begins with
# 'VERIFY ON LOAD:' — the h5ad obs is authoritative if it disagrees with
# the paper, but a draft-to-verify is strictly safer than a blank-to-fill.
#
# Coverage is intentionally partial. Labels NOT in the draft remain as a
# single `<pending-HB-unseen: N labels>` skeleton row per atlas so the
# coverage gap is explicitly tracked, not implicitly omitted. The full
# narrative — including granularity flags, uncertain rows, and the
# research-agent source URLs — lives in
# `data/atlases/crosswalk_draft_notes.md` (companion human-readable file).

# TAURUS-IBD (Thomas et al., Nat Immunol 2024) — ~38 of 109 cell_state
# labels drafted from Extended Data Fig 1c-j dotplot legends (PMC11519010).
# Remaining ~71 cell_state labels stay `<pending-HB>` — Supplementary
# Table 1 (xlsx) and the Zenodo h5ad obs catalog were not reachable from
# the laptop fetch.
_TAURUS_DRAFT_FINE: tuple[tuple[str, str, str, str], ...] = (
    # (native_label, canonical_broad, canonical_fine, rationale)
    # T-cell compartment — all collapse to T cell (no F8 subset bucket)
    ("Th",                       "T cell", "t_cell_broad", "conventional CD4 helper"),
    ("Tfh",                      "T cell", "t_cell_broad", "follicular helper"),
    ("Tph",                      "T cell", "t_cell_broad", "peripheral helper, IBD-relevant"),
    ("Treg",                     "T cell", "t_cell_broad", "FOXP3+ regulatory"),
    ("TWIST1+ Treg",             "T cell", "t_cell_broad", "tissue Treg subset"),
    ("IFN-resp",                 "T cell", "t_cell_broad", "IFN-stimulated CD4"),
    ("CD8+ T",                   "T cell", "t_cell_broad", "CD8 conventional"),
    ("GZMK int",                 "T cell", "t_cell_broad", "CD8 effector intermediate"),
    ("GZMK hi",                  "T cell", "t_cell_broad", "CD8 effector high"),
    ("FGFBP2+",                  "T cell", "t_cell_broad", "cytotoxic CD8 / NK-like"),
    ("CTLA4 hi TIGIT hi",        "T cell", "t_cell_broad", "exhausted CD8"),
    ("gamma-delta T",            "T cell", "t_cell_broad", "unconventional T"),
    ("MAIT",                     "T cell", "t_cell_broad", "invariant MR1-restricted"),
    # NK / ILC
    ("NK",                       "NK/ILC", "", "classical NK"),
    ("ILC",                      "NK/ILC", "", "innate lymphoid"),
    # B / plasma
    ("Naive B",                  "B cell", "", ""),
    ("Memory B",                 "B cell", "", ""),
    ("IFN-resp B",               "B cell", "", "IFN-stim B"),
    ("GC B",                     "B cell", "", "germinal centre"),
    ("Plasmablast",              "plasma cell", "plasma_cell_collapsed", "proliferative"),
    ("Plasma cell",              "plasma cell", "plasma_cell_collapsed", ""),
    ("IgG+ Plasma cell",         "plasma cell", "plasma_cell_collapsed", "IBD class-switched"),
    # Myeloid
    ("Monocyte",                 "monocyte/macrophage", "", "classical mono"),
    ("S100A8/9 hi TNF hi IL6+",  "monocyte/macrophage", "", "inflammatory mono — IBD signature"),
    ("C1Q hi IL1B lo",           "monocyte/macrophage", "", "resident macrophage"),
    ("C1Q hi IL1B hi",           "monocyte/macrophage", "", "inflammatory macrophage"),
    ("Macrophage",               "monocyte/macrophage", "", ""),
    ("DC",                       "dendritic cell", "", ""),
    ("LAMP3+ IL1B+ DC",          "dendritic cell", "", "mature/activated DC"),
    ("pDC",                      "dendritic cell", "", "plasmacytoid"),
    # Stromal
    ("Fibroblast",               "fibroblast", "fibroblast_stromal", ""),
    ("THY1+ FAP+ PDPN+",         "fibroblast", "inflammatory_fibroblast", "IBD inflammatory fibroblast signature"),
    ("Pericyte",                 "fibroblast", "", "F8 collapses pericyte into fibroblast bucket"),
    ("Vascular",                 "endothelium", "endothelial", "endothelial cells; LEC/BEC likely lumped"),
    # Epithelial — ileal
    ("Enterocyte",               "colonocyte", "absorptive_enterocyte",
     "SI absorptive — canonical_broad reuses 'colonocyte' as the absorptive bucket"),
    ("TA (ileal)",               "epithelial progenitor", "immature_enterocyte_ta", "transit-amplifying"),
    ("Undiff (ileal)",           "epithelial progenitor", "immature_enterocyte_ta", "undifferentiated TA-adjacent"),
    ("Goblet (ileal)",           "goblet", "goblet", ""),
    ("Tuft",                     "enteroendocrine/tuft", "tuft", ""),
    ("EEC (ileal)",              "enteroendocrine/tuft", "enteroendocrine", ""),
    # Epithelial — colonic
    ("Colonocyte",               "colonocyte", "absorptive_enterocyte", ""),
    ("LGR5+ Stem",               "epithelial progenitor", "crypt_stem", "LGR5+ crypt-base stem"),
    ("TA (colonic)",             "epithelial progenitor", "immature_enterocyte_ta", ""),
    ("Undiff (colonic)",         "epithelial progenitor", "immature_enterocyte_ta", ""),
    ("Goblet (colonic)",         "goblet", "goblet", ""),
    ("EEC (colonic)",            "enteroendocrine/tuft", "enteroendocrine", ""),
    ("CT Goblet",                "goblet", "goblet", "crypt-top goblet, UC-relevant"),
)
_TAURUS_UNSEEN_FINE_EST: int = 64  # 109 paper-reported - 45 drafted above
_TAURUS_STRUCTURAL_ZEROS: frozenset[str] = frozenset({
    "granulocyte",   # mucosal-biopsy droplet 10x lyses granulocytes (cf. Smillie)
})
_TAURUS_UNCERTAIN_BROAD: frozenset[str] = frozenset({
    "mast cell",          # no mast label surfaced in ED Fig 1; may sit inside myeloid
    "mural/glia",         # no enteric glia label surfaced; biopsy may not sample submucosa
})

# HCA Gut (Elmentaite et al., Nature 597:250) — ~110 of ~120 author_cell_type
# labels drafted from PMC8426186 full-text. Region-filter to adult colon
# BEFORE scoring (per F1 + per the loader's tissue filter) so developmental
# / non-colon labels do not pollute broad-bucket counts.
_HCA_DRAFT_FINE: tuple[tuple[str, str, str, str], ...] = (
    # Epithelial
    ("Stem cells",                        "epithelial progenitor", "crypt_stem", "LGR5+ stem"),
    ("Proximal progenitors",              "epithelial progenitor", "crypt_stem", "regional stem; developmental"),
    ("Distal progenitors",                "epithelial progenitor", "crypt_stem", "regional stem; developmental"),
    ("Transit-amplifying (TA)",           "epithelial progenitor", "immature_enterocyte_ta", ""),
    ("CLDN10 cells",                      "epithelial progenitor", "", "pancreatic-progenitor-like; ambiguous"),
    ("Enterocytes",                       "colonocyte", "absorptive_enterocyte", "SI absorptive"),
    ("Colonocytes",                       "colonocyte", "absorptive_enterocyte", ""),
    ("BEST4 enterocytes",                 "colonocyte", "absorptive_enterocyte", "BEST4+ subtype"),
    ("BEST2+ goblet cells",               "goblet", "goblet", ""),
    ("Goblet cells",                      "goblet", "goblet", ""),
    ("Paneth cells",                      "goblet", "paneth_like", "secretory; F8 paneth_like bucket"),
    ("Tuft cells",                        "enteroendocrine/tuft", "tuft", ""),
    ("Microfold (M) cells",               "colonocyte", "", "specialized antigen-sampling; no F8 bucket"),
    ("Enteroendocrine cells",             "enteroendocrine/tuft", "enteroendocrine", ""),
    ("Enterochromaffin (EC) cells",       "enteroendocrine/tuft", "enteroendocrine", ""),
    ("M/X cells",                         "enteroendocrine/tuft", "enteroendocrine", ""),
    ("D cells",                           "enteroendocrine/tuft", "enteroendocrine", "SST+"),
    ("beta cells",                        "enteroendocrine/tuft", "enteroendocrine", "INS+; FLAG developmental"),
    ("L cells",                           "enteroendocrine/tuft", "enteroendocrine", "GCG+"),
    ("N cells",                           "enteroendocrine/tuft", "enteroendocrine", "NTS+"),
    ("K cells",                           "enteroendocrine/tuft", "enteroendocrine", "GIP+"),
    ("I cells",                           "enteroendocrine/tuft", "enteroendocrine", "CCK+"),
    ("NPW-EC cells",                      "enteroendocrine/tuft", "enteroendocrine", ""),
    ("TAC1-EC cells",                     "enteroendocrine/tuft", "enteroendocrine", ""),
    ("NEUROG3+ progenitors",              "enteroendocrine/tuft", "", "EEC progenitor"),
    # Endothelial
    ("Arterial endothelial cells",        "endothelium", "endothelial", ""),
    ("Venous endothelial cells",          "endothelium", "endothelial", ""),
    ("Capillary endothelial cells",       "endothelium", "endothelial", ""),
    ("Arterial capillaries",              "endothelium", "endothelial", ""),
    ("LEC1",                              "endothelium", "lymphatic_endothelial", "lymphatic subset"),
    ("LEC2",                              "endothelium", "lymphatic_endothelial", "lymphatic subset"),
    ("LEC3",                              "endothelium", "lymphatic_endothelial", "lymphatic subset"),
    ("LEC4",                              "endothelium", "lymphatic_endothelial", "lymphatic subset"),
    ("LEC5",                              "endothelium", "lymphatic_endothelial", "lymphatic subset"),
    ("LEC6",                              "endothelium", "lymphatic_endothelial", "lymphatic subset"),
    # Neural
    ("Enteric neural crest cells (ENCCs)","mural/glia", "enteroglial", "neural crest progenitor; developmental"),
    ("Neuroblasts",                       "mural/glia", "enteroglial", "developing neuron"),
    ("Branch A1 (iMN)",                   "mural/glia", "enteroglial", "enteric neuron lineage"),
    ("Branch A2 (IPAN/IN)",               "mural/glia", "enteroglial", "enteric neuron lineage"),
    ("Branch A3 (IPAN/IN)",               "mural/glia", "enteroglial", "enteric neuron lineage"),
    ("Branch A4 (IN)",                    "mural/glia", "enteroglial", "enteric neuron lineage"),
    ("Branch B1 (immature eMN)",          "mural/glia", "enteroglial", "enteric neuron lineage"),
    ("Branch B2 (eMN)",                   "mural/glia", "enteroglial", "enteric neuron lineage"),
    ("Branch B3 (IPAN)",                  "mural/glia", "enteroglial", "enteric neuron lineage"),
    ("Glia 1 (DHH+)",                     "mural/glia", "enteroglial", ""),
    ("Glia 2 (ELN+)",                     "mural/glia", "enteroglial", ""),
    ("Glia 3 (BCAN+)",                    "mural/glia", "enteroglial", ""),
    ("Differentiating glia (COL20A1+)",   "mural/glia", "enteroglial", ""),
    # Mesenchymal
    ("Mesoderm 1",                        "fibroblast", "", "developmental progenitor"),
    ("Mesoderm 2",                        "fibroblast", "", "developmental progenitor"),
    ("Stromal 1",                         "fibroblast", "fibroblast_stromal", "S1 crypt/villus tip"),
    ("Stromal 2",                         "fibroblast", "fibroblast_stromal", "S2 crypt base"),
    ("Stromal 3",                         "fibroblast", "inflammatory_fibroblast",
     "S3 inflammation-associated per Smillie nomenclature; VERIFY marker overlap"),
    ("Stromal 4",                         "fibroblast", "fibroblast_stromal", "S4 submucosal"),
    ("T reticular cells",                 "fibroblast", "fibroblast_stromal", "reticular network"),
    ("Follicular dendritic cells (FDCs)", "fibroblast", "fibroblast_stromal",
     "stromal lineage NOT myeloid DC; confusable naming"),
    ("FMO2 stromal cells",                "fibroblast", "fibroblast_stromal", ""),
    ("Myofibroblasts",                    "fibroblast", "", "no F8 myofib bucket"),
    ("Cycling myofibroblasts",            "fibroblast", "", ""),
    ("Smooth muscle cells",               "mural/glia", "", "no F8 SM bucket; broad mural"),
    ("Interstitial cells of Cajal (ICC)", "mural/glia", "", "interstitial pacemaker"),
    ("Immature pericytes",                "fibroblast", "", "broad collapses pericyte"),
    ("Contractile pericytes",             "fibroblast", "", ""),
    ("Angiogenic pericytes",              "fibroblast", "", ""),
    ("CD36+ pericytes",                   "fibroblast", "", ""),
    ("Mature pericytes",                  "fibroblast", "", ""),
    ("Mesothelial cells",                 "fibroblast", "", "serosal; no mesothelium bucket — forced fit"),
    ("RGS5+ mesothelial cells",           "fibroblast", "", "serosal; forced fit"),
    ("Mesenchymal lymphoid tissue organizers (mLTo)",
                                          "fibroblast", "fibroblast_stromal", "developmental organizer"),
    # T lymphoid
    ("CD4 T cells",                       "T cell", "t_cell_broad", ""),
    ("SELL+ CD4 T cells",                 "T cell", "t_cell_broad", "naive-like CD4"),
    ("T regulatory cells (Treg)",         "T cell", "t_cell_broad", "FOXP3+"),
    ("CD8 T cells",                       "T cell", "t_cell_broad", ""),
    ("Tissue-resident CD8 T cells",       "T cell", "t_cell_broad", "Trm CD8"),
    ("TCRalpha-beta+ T cells",            "T cell", "t_cell_broad", "conventional alpha-beta"),
    ("TCRgamma-delta+ T cells",           "T cell", "t_cell_broad", "unconventional gamma-delta"),
    ("Cycling T cells",                   "T cell", "t_cell_broad", "cycling-T-lineage; NOT unprefixed-cycling"),
    # B lymphoid
    ("Common lymphoid progenitor (CLP)",  "B cell", "", "progenitor; developmental"),
    ("Pro-B cells",                       "B cell", "", "developmental"),
    ("Pre-B cells",                       "B cell", "", "developmental"),
    ("Immature B cells",                  "B cell", "", ""),
    ("Naive B cells",                     "B cell", "", ""),
    ("Memory B cells",                    "B cell", "", ""),
    ("FCRL4+ memory B cells",             "B cell", "", "tissue-resident memory"),
    ("Cycling B cells",                   "B cell", "", "cycling-B-lineage"),
    ("IgM plasma cells",                  "plasma cell", "plasma_cell_collapsed", ""),
    ("IgA plasma cells",                  "plasma cell", "plasma_cell_collapsed", ""),
    ("IgG plasma cells",                  "plasma cell", "plasma_cell_collapsed", ""),
    # ILC / NK
    ("ILCPs",                             "NK/ILC", "", "ILC progenitor"),
    ("NCR+ ILC3",                         "NK/ILC", "", ""),
    ("NCR- ILC3",                         "NK/ILC", "", ""),
    ("LTi-like ILC3",                     "NK/ILC", "", "lymphoid tissue inducer"),
    ("NK cells",                          "NK/ILC", "", "classical NK"),
    ("Adult ILC3",                        "NK/ILC", "", ""),
    # Erythroid — STRUCTURAL EXTRA, not in 15-term vocab; drop
    ("Erythroid",                         "(exclude)", "", "developmental; drop in QC (not in 15-term vocab)"),
)
_HCA_UNSEEN_FINE_EST: int = 30  # ~120 paper-reported - ~90 drafted; PMC excerpt thin on myeloid
_HCA_STRUCTURAL_ZEROS: frozenset[str] = frozenset()  # all 15 likely populated by full obs
_HCA_UNCERTAIN_BROAD: frozenset[str] = frozenset({
    "mast cell",          # myeloid extraction was thin in PMC excerpt
    "monocyte/macrophage",
    "dendritic cell",
    "granulocyte",
})

# Pan-GI Extended+ (Oliver et al., Nature 635:699, 2024) — ~45 of 136
# level_3_annot labels drafted from PMC11578898 (Fig 1-4 captions). Spans
# whole-GI: region-filter to colon BEFORE scoring; ~30-40% of labels are
# non-colon (oral mucosa, oesophagus, gastric, duodenal) and pollute
# broad-bucket counts if not filtered.
_PANGI_DRAFT_FINE: tuple[tuple[str, str, str, str], ...] = (
    # Epithelial — note INFLAREs + MGN are pyloric/Brunner's-like, appear
    # in UC/Crohn's diseased colon as metaplasia.
    ("INFLAREs",                          "goblet", "paneth_like",
     "metaplastic secretory MUC6+; NOVEL to Pan-GI — no direct counterpart in Smillie/Garrido/TAURUS"),
    ("MGN cells",                         "goblet", "paneth_like", "mucous gland neck; metaplastic"),
    ("Surface foveolar cells",            "colonocyte", "", "gastric surface; non-colon, region-filter"),
    ("Surface foveolar-like cells",       "colonocyte", "", "metaplastic; UNCERTAIN broad"),
    ("Paneth cells",                      "goblet", "paneth_like", "native SI Paneth"),
    ("Metaplastic Paneth cells",          "goblet", "paneth_like", "colonic metaplastic Paneth, UC-relevant"),
    ("Goblet cells",                      "goblet", "goblet", ""),
    ("BEST4 enterocytes",                 "colonocyte", "absorptive_enterocyte", "BEST4+ subtype"),
    ("Colonocytes",                       "colonocyte", "absorptive_enterocyte", ""),
    ("LGR5+ stem cells",                  "epithelial progenitor", "crypt_stem", ""),
    ("Transit amplifying (TA)",           "epithelial progenitor", "immature_enterocyte_ta", ""),
    ("Deep crypt secretory (DCS)",        "goblet", "", "colonic functional equivalent of Paneth"),
    ("Tuft cells",                        "enteroendocrine/tuft", "tuft", ""),
    ("Enteroendocrine cells",             "enteroendocrine/tuft", "enteroendocrine", "subtypes not enumerated in fetched text"),
    # Mesenchymal
    ("Crypt fibroblasts (PI16+)",         "fibroblast", "fibroblast_stromal", "crypt-base PI16+"),
    ("Lamina propria fibroblasts (ADAMDEC1+)",
                                          "fibroblast", "fibroblast_stromal", "LP ADAMDEC1+"),
    ("Villus fibroblasts (F3+)",          "fibroblast", "fibroblast_stromal", ""),
    ("Inflammatory fibroblasts",          "fibroblast", "inflammatory_fibroblast", "IBD-associated; clean cross-atlas anchor"),
    ("Oral mucosa fibroblasts",           "fibroblast", "inflammatory_fibroblast",
     "appear metaplastically in UC/Crohn's colon per paper"),
    ("Oesophagus fibroblasts",            "fibroblast", "fibroblast_stromal", "non-colon, region-filter"),
    ("Rectum fibroblasts",                "fibroblast", "fibroblast_stromal", ""),
    ("Smooth muscle cells",               "mural/glia", "", "no F8 SM bucket"),
    # Endothelial
    ("Venous endothelial (ACKR1+)",       "endothelium", "endothelial", ""),
    ("Arterial endothelial",              "endothelium", "endothelial", ""),
    ("Capillary endothelial",             "endothelium", "endothelial", ""),
    ("Lymphatic endothelial",             "endothelium", "lymphatic_endothelial", ""),
    # Neural
    ("Glia (subtypes)",                   "mural/glia", "enteroglial", "placeholder; subtypes not enumerated in fetched text"),
    ("Enteric neurons (subtypes)",        "mural/glia", "enteroglial", "placeholder; subtypes not enumerated"),
    # T/NK
    ("CD4 Th17",                          "T cell", "t_cell_broad", ""),
    ("CD4 Treg",                          "T cell", "t_cell_broad", ""),
    ("CD4 TEM",                           "T cell", "t_cell_broad", "effector memory"),
    ("CD4 TRM",                           "T cell", "t_cell_broad", "tissue-resident memory"),
    ("CD8 T",                             "T cell", "t_cell_broad", ""),
    ("gamma-delta T",                     "T cell", "t_cell_broad", ""),
    ("MAIT",                              "T cell", "t_cell_broad", ""),
    ("CD56bright NK",                     "NK/ILC", "", ""),
    ("ILC3",                              "NK/ILC", "", ""),
    # Myeloid
    ("Macrophages",                       "monocyte/macrophage", "", ""),
    ("LYVE1+ macrophages",                "monocyte/macrophage", "", "tissue-resident"),
    ("Dendritic cell subtypes",           "dendritic cell", "", "placeholder; expect cDC1/cDC2/migratory/pDC"),
    ("Mast cells",                        "mast cell", "", ""),
    ("Neutrophils",                       "granulocyte", "", "only granulocyte anchor across all five atlases"),
    # B / plasma
    ("Progenitor B cells",                "B cell", "", ""),
    ("Mature B cells",                    "B cell", "", ""),
    ("IgA plasma cells",                  "plasma cell", "plasma_cell_collapsed", ""),
    ("IgA2 plasma cells",                 "plasma cell", "plasma_cell_collapsed", ""),
    ("IgM plasma cells",                  "plasma cell", "plasma_cell_collapsed", ""),
)
_PANGI_UNSEEN_FINE_EST: int = 91  # ~136 paper-reported - 45 drafted
_PANGI_STRUCTURAL_ZEROS: frozenset[str] = frozenset()  # whole-GI atlas; all 15 likely populated
_PANGI_UNCERTAIN_BROAD: frozenset[str] = frozenset()


def _build_paper_supplement_rows(
    atlas: str, fine_col: str, broad_col: str,
    paper_ref: str, draft_fine: tuple[tuple[str, str, str, str], ...],
    unseen_fine_est: int, structural_zeros: frozenset[str],
    uncertain_broad: frozenset[str],
) -> list[dict]:
    """Emit paper-supplement-sourced draft rows for an HB-bound atlas.

    Each drafted native fine label becomes a 'fine' tier row with
    source='paper-supplement' and notes prefixed 'VERIFY ON LOAD:'. A
    single residual `<pending-HB-unseen: N labels>` row tracks the
    coverage gap so unseen labels are explicit, not implicit.

    For broad tier, one row per _BROAD_VOCAB term: 'paper-supplement'
    if the draft populates it (with the per-bucket native_label count),
    'paper-supplement-structural-zero' if the atlas is known not to
    produce that bucket (e.g. granulocyte for biopsy droplet-10x), or
    'paper-supplement-uncertain' if the fetched text did not surface
    enough evidence to call populated-vs-zero.
    """
    rows: list[dict] = []

    # Fine rows — one per drafted native_label
    for native, broad, fine, rationale in draft_fine:
        cl = CL_ANCHORS.get(broad, "")
        rows.append({
            "atlas": atlas,
            "tier": "fine",
            "native_label": native,
            "native_col": fine_col,
            "post_qc_label": native,
            "canonical_broad": broad,
            "canonical_fine": fine,
            "cl_anchor": cl,
            "source": "paper-supplement",
            "notes": "VERIFY ON LOAD: " + (rationale or "draft from paper/portal"),
        })

    # Residual unseen-labels skeleton row
    rows.append({
        "atlas": atlas,
        "tier": "fine",
        "native_label": f"<pending-HB-unseen: ~{unseen_fine_est} labels>",
        "native_col": fine_col,
        "post_qc_label": "<pending-HB>",
        "canonical_broad": "<pending-HB>",
        "canonical_fine": "<pending-HB>",
        "cl_anchor": "",
        "source": "pending-HB",
        "notes": (f"{paper_ref}. Drafted rows above cover the subset reachable "
                  "from paper supplement / portal at draft time; remaining labels "
                  "must be enumerated from obs on first HB load and added by "
                  "extending the draft block in scripts/build_celltype_crosswalk.py. "
                  "See data/atlases/crosswalk_draft_notes.md for the per-atlas "
                  "granularity flags and uncertain rows."),
    })

    # Per-canonical-broad rows: which buckets the draft populates, with
    # source provenance + a count of contributing fine labels.
    per_broad_count: dict[str, int] = {}
    for _native, broad, _fine, _ratio in draft_fine:
        if broad in CL_ANCHORS:
            per_broad_count[broad] = per_broad_count.get(broad, 0) + 1
    for broad in sorted(CL_ANCHORS):
        cl = CL_ANCHORS[broad]
        if broad in structural_zeros:
            source = "paper-supplement-structural-zero"
            note = (f"VERIFY ON LOAD: {paper_ref}. Atlas not expected to populate "
                    f"this bucket (technical / cohort reason — see "
                    "crosswalk_draft_notes.md).")
        elif per_broad_count.get(broad, 0) > 0:
            source = "paper-supplement"
            note = (f"VERIFY ON LOAD: {paper_ref}. Drafted from "
                    f"{per_broad_count[broad]} native fine label(s) in the draft block.")
        elif broad in uncertain_broad:
            source = "paper-supplement-uncertain"
            note = (f"VERIFY ON LOAD: {paper_ref}. Fetched paper/portal text did not "
                    "surface a label for this bucket; could not call "
                    "populated-vs-zero. Confirm from obs on first HB load.")
        else:
            source = "paper-supplement-uncertain"
            note = (f"VERIFY ON LOAD: {paper_ref}. No drafted fine label maps to "
                    "this bucket from the fetched paper/portal subset; status "
                    "(populated by an unseen label vs structural zero) is unknown "
                    "until obs enumeration on first HB load.")
        rows.append({
            "atlas": atlas,
            "tier": "broad",
            "native_label": "<draft-from-paper>",
            "native_col": broad_col,
            "post_qc_label": broad,
            "canonical_broad": broad,
            "canonical_fine": "",
            "cl_anchor": cl,
            "source": source,
            "notes": note,
        })
    return rows


def _build_taurus_rows() -> list[dict]:
    return _build_paper_supplement_rows(
        atlas="taurus",
        fine_col="cell_state",
        broad_col="low (mapped via LOW_TO_BROAD)",
        paper_ref=("Thomas et al., Nat Immunol 25:2152-2165 (2024); "
                   "Zenodo v3 10.5281/zenodo.14007626; draft from PMC11519010 "
                   "Extended Data Fig 1c-j (2026-06-27)"),
        draft_fine=_TAURUS_DRAFT_FINE,
        unseen_fine_est=_TAURUS_UNSEEN_FINE_EST,
        structural_zeros=_TAURUS_STRUCTURAL_ZEROS,
        uncertain_broad=_TAURUS_UNCERTAIN_BROAD,
    )


def _build_hca_gut_rows() -> list[dict]:
    return _build_paper_supplement_rows(
        atlas="hca_gut",
        fine_col="author_cell_type",
        broad_col="category (mapped via category -> _BROAD_VOCAB)",
        paper_ref=("Elmentaite et al., Nature 597:250 (2021); CELLxGENE deposit "
                   "f34d2b82-9265-4a73-bda4-852933bf2a8d; draft from PMC8426186 "
                   "full-text (2026-06-27)"),
        draft_fine=_HCA_DRAFT_FINE,
        unseen_fine_est=_HCA_UNSEEN_FINE_EST,
        structural_zeros=_HCA_STRUCTURAL_ZEROS,
        uncertain_broad=_HCA_UNCERTAIN_BROAD,
    )


def _build_pangi_rows() -> list[dict]:
    return _build_paper_supplement_rows(
        atlas="pangi",
        fine_col="level_3_annot",
        broad_col="level_2_annot (mapped via level_2_annot -> _BROAD_VOCAB)",
        paper_ref=("Oliver et al., Nature 635:699 (2024); CELLxGENE deposit "
                   "1dcf15ee-c103-4aaa-8b8c-0fc697fcccc8; draft from PMC11578898 "
                   "Fig 1-4 captions (2026-06-27). Pan-GI integrates "
                   "Smillie+Kong donors — donor-overlap scan flagged in DECISIONS."),
        draft_fine=_PANGI_DRAFT_FINE,
        unseen_fine_est=_PANGI_UNSEEN_FINE_EST,
        structural_zeros=_PANGI_STRUCTURAL_ZEROS,
        uncertain_broad=_PANGI_UNCERTAIN_BROAD,
    )


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=_REPO
        ).decode().strip()
    except Exception:
        return "UNKNOWN"


HEADER_COMMENT = """# UC-Cross-Atlas cell-type label crosswalk
#
# Generated by scripts/build_celltype_crosswalk.py — rerun on any change
# to _broad_vocab.py / loader FINE_TO_BROAD maps / canonical_broad_DRAFT.md
# CL pin. Do not hand-edit; this file is regenerable.
#
# Tier 'fine' rows = the post-QC fine labels emitted into obs['cell_type_fine']
# by each loader. Tier 'broad' rows = the 15 _BROAD_VOCAB terms; one row per
# (atlas, broad) documenting CL anchor + whether the atlas populates that
# bucket.
#
# Atlases that live on HB only (TAURUS, HCA Gut, Pan-GI) carry
# paper-supplement-sourced DRAFT rows from 2026-06-27 (per DECISIONS 31),
# marked source='paper-supplement' and notes prefixed 'VERIFY ON LOAD:'.
# A single residual `<pending-HB-unseen: ~N labels>` row per atlas tracks
# the coverage gap for native labels that were NOT reachable from the
# paper / portal at draft time. Companion narrative (granularity flags,
# uncertain rows, source URLs): data/atlases/crosswalk_draft_notes.md.
#
# build_sha: %s
# build_date: %s
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--out", type=Path,
        default=_REPO / "data/atlases/celltype_crosswalk.tsv",
    )
    args = parser.parse_args()

    rows = []
    rows.extend(_build_garrido_rows())
    rows.extend(_build_smillie_rows())
    rows.extend(_build_taurus_rows())
    rows.extend(_build_hca_gut_rows())
    rows.extend(_build_pangi_rows())

    fields = ["atlas", "tier", "native_label", "native_col",
              "post_qc_label", "canonical_broad", "canonical_fine",
              "cl_anchor", "source", "notes"]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8", newline="") as fh:
        fh.write(HEADER_COMMENT % (_git_sha(), dt.date.today().isoformat()))
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                            lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    n_per_atlas: dict[str, int] = {}
    n_pending: dict[str, int] = {}
    n_paper: dict[str, int] = {}
    for r in rows:
        n_per_atlas[r["atlas"]] = n_per_atlas.get(r["atlas"], 0) + 1
        if r["source"] == "pending-HB":
            n_pending[r["atlas"]] = n_pending.get(r["atlas"], 0) + 1
        if r["source"].startswith("paper-supplement"):
            n_paper[r["atlas"]] = n_paper.get(r["atlas"], 0) + 1

    print(f"[build_celltype_crosswalk] wrote {args.out} ({len(rows)} rows)")
    for a in sorted(n_per_atlas):
        pend = n_pending.get(a, 0)
        paper = n_paper.get(a, 0)
        if paper > 0:
            status = (f"DRAFT (paper-supplement; verify against obs) — "
                      f"{paper} drafted rows, {pend} pending-HB-unseen")
        elif pend > 0:
            status = f"{pend} pending-HB"
        else:
            status = "complete (loader-extracted)"
        print(f"  {a:10s}  {n_per_atlas[a]:4d} rows  ({status})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
