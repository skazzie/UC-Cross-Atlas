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


def _build_taurus_rows() -> list[dict]:
    return _build_pending_hb_rows(
        atlas="taurus",
        fine_col="cell_state",
        broad_col="low (mapped via LOW_TO_BROAD)",
        expected_n_fine_est="~109",
        paper_ref="Thomas et al., Nat Immunol 25:2152-2165 (2024); Zenodo v3 10.5281/zenodo.14007626",
        notes=("LOW_TO_BROAD ships EMPTY in load_taurus.py; gate (2) "
               "fails loud on first HB run with the actual low-tier "
               "label set, then this crosswalk row should be replaced "
               "with the populated mapping."),
    )


def _build_hca_gut_rows() -> list[dict]:
    return _build_pending_hb_rows(
        atlas="hca_gut",
        fine_col="author_cell_type",
        broad_col="category",
        expected_n_fine_est="~120",
        paper_ref="Elmentaite et al., Nature 597:250 (2021); CELLxGENE deposit f34d2b82",
        notes=("Native 'category' values (Mesenchymal / Neuronal / "
               "Epithelial / Lymphoid / Myeloid / etc.) do NOT match "
               "_BROAD_VOCAB; load_hca_gut.py passes them through "
               "unmapped. Net-new category -> _BROAD_VOCAB mapping "
               "needs to be authored when first HB load enumerates "
               "category labels."),
    )


def _build_pangi_rows() -> list[dict]:
    return _build_pending_hb_rows(
        atlas="pangi",
        fine_col="level_3_annot",
        broad_col="level_2_annot",
        expected_n_fine_est="~70",
        paper_ref="Oliver et al. 2024 (Pan-GI Extended+); CELLxGENE deposit 1dcf15ee",
        notes=("Native 'level_2_annot' values do NOT match _BROAD_VOCAB; "
               "load_pangi.py passes them through unmapped. Net-new "
               "level_2_annot -> _BROAD_VOCAB mapping needs to be "
               "authored when first HB load enumerates the labels. "
               "Pan-GI integrates Smillie+Kong donors — donor-overlap "
               "scan flagged in DECISIONS."),
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
# Atlases that live on HB only (TAURUS, HCA Gut, Pan-GI) get skeleton rows
# with native_label='<pending-HB>' — these MUST be replaced with the actual
# native label set after the first HB load.
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
    for r in rows:
        n_per_atlas[r["atlas"]] = n_per_atlas.get(r["atlas"], 0) + 1
        if r["source"] == "pending-HB":
            n_pending[r["atlas"]] = n_pending.get(r["atlas"], 0) + 1

    print(f"[build_celltype_crosswalk] wrote {args.out} ({len(rows)} rows)")
    for a in sorted(n_per_atlas):
        pend = n_pending.get(a, 0)
        status = f"{pend} pending-HB" if pend else "complete"
        print(f"  {a:10s}  {n_per_atlas[a]:4d} rows  ({status})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
