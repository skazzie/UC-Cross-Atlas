"""Driver: load Pan-GI Extended+, remap cell_type_broad to canonical
vocab, precompute neighbors, and write h5ad + partial covariate file.

Source pin: CELLxGENE deposit 1dcf15ee-c103-4aaa-8b8c-0fc697fcccc8.h5ad;
see ``code/02_atlas_prep/atlas_schemas.md``.

Applies every fix worked out during the TAURUS debug loop:

- ``load_pangi.load()`` already reads backed + materializes only the
  DISEASE_KEEP x ORGAN_KEEP subset (fits 128 GB VM; unbacked read of
  the full 1.6M-cell h5ad OOMs a 64 GB box).
- Two-gate canonical-vocab assertion around ``PANGI_TO_BROAD``: gate 1
  at module import (every mapped value must be in ``_BROAD_VOCAB``);
  gate 2 at run-time (every ``level_2_annot`` label must have a
  mapping). ``PANGI_TO_BROAD`` ships EMPTY in v0; gate 2 raises on
  first run with the label list so it can be populated iteratively
  (one commit per Muskaan biology pass), same pattern as
  ``MAJOR_TO_BROAD`` in ``load_taurus.py``.
- scDRS-replica pre-filter (``min_genes=250`` / ``min_cells=50``)
  applied here so ``obsp['connectivities']`` on the shipped h5ad
  survives ``scdrs perform-downstream --flag-filter-data False``.
  ``uns['scdrs_prefilter']`` stamped for provenance.
- Cov file matches ``smillie_covariates.tsv`` / the fixed
  ``taurus_covariates.tsv``: ``const``, ``log_n_genes``,
  ``log_n_counts``, then explicit ``sample_<id>`` float 0.0/1.0
  columns built with ``pd.get_dummies(...).astype(float)``. Categorical
  columns are NOT written — scDRS's ``category2dummy`` produces bool
  dummies which upcast ``df_cov.values`` to object and crash
  ``np.linalg.solve``. An all-numeric assertion runs before write.
- ``donor`` dropped from cov: samples nest within donors, so donor
  dummies are a linear combination of sample dummies and joint dummies
  + const are rank-deficient. Nesting gate runs before write and
  fails loud if the assumption is violated. If ``--sample-col`` is
  omitted and no finer column is found, sample defaults to
  ``donorID_unified`` (sample == donor); the drop is then vacuous.
- NO disease / health covariate. Pan-GI's ``disease`` col MUST NOT
  enter the covariate — that would scrub the very signal scDRS tests
  (no-disease-covariate rule).

Reference: ``run_taurus_load.py`` — mirror driver structure.
"""

from __future__ import annotations

import argparse
import logging

import numpy as np
import pandas as pd
import scanpy as sc

from load_pangi import load
from _broad_vocab import _BROAD_VOCAB

logger = logging.getLogger(__name__)

# scDRS 1.0.2 load-time filter defaults (scdrs/util.py:82-88). Match
# here so the on-disk cell set == what --flag-filter-data True would
# produce.
SCDRS_MIN_GENES_PER_CELL = 250
SCDRS_MIN_CELLS_PER_GENE = 50

# Pan-GI's donor identifier per load_pangi.load(). The nesting gate
# below groups by this. If Pan-GI ever migrates to a different donor
# column, update in one place.
DONOR_COL = "donorID_unified"

# Sample-level column candidates. Pan-GI CELLxGENE obs does not
# guarantee a finer-than-donor sample identifier; auto-detect falls
# back to DONOR_COL (sample == donor, drop-donor is vacuous). Override
# with --sample-col if the deposit ships a library / biopsy field.
_SAMPLE_COL_CANDIDATES = ("sample_id", "sample", "library_id", "library_uuid")

# Map from Pan-GI level_2_annot labels into the canonical broad vocab
# (single-sourced from _broad_vocab._BROAD_VOCAB — same 15 terms as
# Garrido / Smillie / TAURUS so 06_concordance's string-intersection
# on cell_type_broad aligns across atlases).
#
# 43 level_2_annot labels enumerated on the 128 GB VM 2026-08-09.
#
# Scope + collapse notes:
# - No "mast cell" bucket at Pan-GI level_2 (would live at level_3
#   under myeloid); mast cell broad tier is therefore not estimable
#   from Pan-GI (fine tier still works if level_3 has it).
# - All conventional T (Conventional_CD4 / CD8 / Treg) + cycling T/NK
#   + Unconventional_T/ILC collapse to "T cell"; only pure NK goes to
#   "NK/ILC". Rationale: mirrors TAURUS's split where innate lymphoid
#   cells sit in NK/ILC and everything else T-lineage goes to T cell.
# - Erythrocyte and Megakaryocyte/platelet have NO natural bucket in
#   _BROAD_VOCAB (all 15 vocab terms are mucosal-resident lineages).
#   Mapping to "granulocyte" as the closest blood-cell bucket; these
#   are typically droplet-contamination in colonic biopsies and could
#   alternatively be excluded via QC on the next pass.
# - Basal / Keratinocyte / Salivary/Submucosal_glands / Mesothelium /
#   Melanocyte / Myoblast/myocyte are non-colonic lineages that should
#   be ~empty after load_pangi's ORGAN_KEEP filter; mappings are kept
#   because the categorical dtype's `.categories` still lists the label
#   even at zero cells and gate (2) checks the full label set.
# - "Unknown" has no biological broad bucket. Mapping to "fibroblast"
#   is a PLACEHOLDER; preferred handling is to exclude these cells
#   before broad-tier concordance (candidate: add a
#   PANGI_EXCLUDE_LEVEL_2 set + filter in main() on the next pass).
PANGI_TO_BROAD: dict[str, str] = {
    # B / plasma lineage
    "B_plasma":                    "plasma cell",
    "Immature_B":                  "B cell",
    "Mature_B":                    "B cell",
    # T / NK / ILC
    "Conventional_CD4":            "T cell",
    "Conventional_CD8":            "T cell",
    "Treg":                        "T cell",
    "Cycling_T/NK":                "T cell",
    "Unconventional_T/ILC":        "T cell",
    "NK":                          "NK/ILC",
    # Myeloid
    "Macrophage":                  "monocyte/macrophage",
    "Monocyte":                    "monocyte/macrophage",
    "DC":                          "dendritic cell",
    "Granulocyte":                 "granulocyte",
    # Non-mucosal blood-lineage — see scope note (contamination-prone)
    "Erythrocyte":                 "granulocyte",
    "Megakaryocyte/platelet":      "granulocyte",
    # Colonic epithelium
    "Absorptive":                  "colonocyte",
    "Microfold":                   "colonocyte",  # M cells; FAE-specialized
    "Secretory":                   "goblet",
    "Enteroendocrine":             "enteroendocrine/tuft",
    "Epithelial_progenitor":       "epithelial progenitor",
    "Epithelial_stem":             "epithelial progenitor",
    "Transit_amplifying":          "epithelial progenitor",
    "Cycling_epithelia":           "epithelial progenitor",
    # Non-colonic epithelium (should be ~empty post-ORGAN_KEEP)
    "Basal":                       "epithelial progenitor",
    "Keratinocyte":                "epithelial progenitor",
    "Salivary/Submucosal_glands":  "goblet",       # mucus-secretory glands
    "Mesothelium":                 "fibroblast",   # mesodermal, imperfect fit
    # Endothelia
    "Vascular_endothelia":         "endothelium",
    "Lymphatic_endothelia":        "endothelium",
    "Cycling_endothelia":          "endothelium",
    # Stroma / fibroblast
    "Fibroblast":                  "fibroblast",
    "Myofibroblast":               "fibroblast",
    "Lymphoid_stromal_cell":       "fibroblast",   # FRC-like
    "Mesoderm":                    "fibroblast",   # mesenchymal default
    # Mural / glia / neural
    "Pericyte":                    "mural/glia",
    "Smooth_muscle":               "mural/glia",
    "Myoblast/myocyte":            "mural/glia",
    "Intestinal_Cajal_cell":       "mural/glia",   # gut pacemaker
    "Glia":                        "mural/glia",
    "Neuron":                      "mural/glia",
    "Neuron_progenitor":           "mural/glia",
    "Melanocyte":                  "mural/glia",   # neural-crest lineage
    # Ambiguous — see 'Unknown' note above
    "Unknown":                     "fibroblast",
}

# Gate (1): every value the map ships must be in the canonical vocab.
# Vacuously true while the map is empty; protects the day we start
# filling it in.
_unmapped_broad_out = set(PANGI_TO_BROAD.values()) - _BROAD_VOCAB
if _unmapped_broad_out:
    raise ValueError(
        f"run_pangi_load.PANGI_TO_BROAD ships broad values outside "
        f"_BROAD_VOCAB: {sorted(_unmapped_broad_out)}. Typo on the "
        f"value side of the map; see canonical_broad_DRAFT.md."
    )
del _unmapped_broad_out


def _raw_counts_matrix(adata):
    """Depth covariates (log_n_counts, log_n_genes) MUST be computed
    from raw UMIs. Prefer ``layers['counts']``, then ``adata.raw.X``,
    then ``adata.X`` (with a loud warning if X is likely normalized).

    load_pangi.load() as of this commit does not populate a counts
    layer; CELLxGENE deposits typically ship raw counts in ``.raw.X``
    and log1p(CP10k) in ``.X``. We surface the fallback so a future
    schema change is visible in the driver log rather than silently
    miscounting depth.
    """
    if "counts" in adata.layers:
        logger.info("[driver] depth from layers['counts']")
        return adata.layers["counts"]
    if getattr(adata, "raw", None) is not None:
        logger.info("[driver] depth from adata.raw.X (no layers['counts'])")
        return adata.raw.X
    logger.warning(
        "[driver] no layers['counts'] and no adata.raw; using adata.X "
        "for depth. If X is log-normalized this MISCOUNTS depth — fix "
        "upstream in load_pangi.load()."
    )
    return adata.X


def _resolve_sample_col(obs: pd.DataFrame, explicit: str | None) -> str:
    if explicit:
        if explicit not in obs.columns:
            raise SystemExit(
                f"[driver] --sample-col={explicit!r} not in obs. "
                f"Available: {sorted(obs.columns)}"
            )
        return explicit
    for c in _SAMPLE_COL_CANDIDATES:
        if c in obs.columns:
            return c
    return DONOR_COL


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--h5ad", required=True,
                   help="Pan-GI CELLxGENE h5ad "
                        "(1dcf15ee-c103-4aaa-8b8c-0fc697fcccc8.h5ad)")
    # Defaults are cwd-relative (no ../../ prefix): if the driver is
    # invoked from the repo root (VM convention), scratch/data/atlases/
    # resolves correctly. Override if you invoke from elsewhere.
    p.add_argument("--out-h5ad",
                   default="scratch/data/atlases/pangi.h5ad")
    p.add_argument("--out-cov",
                   default="scratch/data/atlases/pangi.cov.tsv")
    p.add_argument(
        "--sample-col", default=None,
        help=(f"Sample-level obs column for cov dummies. Defaults to "
              f"auto-detect from {_SAMPLE_COL_CANDIDATES}, else "
              f"{DONOR_COL!r} (sample==donor)."),
    )
    a = p.parse_args()

    adata = load(a.h5ad)  # apply_v1_filter=True; DISEASE_KEEP x ORGAN_KEEP
    n_broad_raw = int(adata.obs["cell_type_broad"].nunique())
    logger.info(
        "[driver] loaded: %d cells x %d genes; %d unique level_2_annot labels",
        adata.n_obs, adata.n_vars, n_broad_raw,
    )

    # ---- Cell-type remap (gate 2) ----
    raw_broad = adata.obs["cell_type_broad"].astype(str)
    unmapped = sorted(set(raw_broad.unique()) - set(PANGI_TO_BROAD))
    if unmapped:
        raise KeyError(
            f"Pan-GI driver: {len(unmapped)} level_2_annot labels have no "
            f"PANGI_TO_BROAD entry. Extend the map "
            f"(run_pangi_load.PANGI_TO_BROAD). Unmapped labels (full list): "
            f"{unmapped}"
        )
    broad = raw_broad.map(PANGI_TO_BROAD)
    outside = set(broad.dropna().unique()) - _BROAD_VOCAB
    if outside:
        raise ValueError(
            f"Pan-GI driver: emitted cell_type_broad values "
            f"{sorted(outside)} not in _BROAD_VOCAB. Fix PANGI_TO_BROAD "
            f"value side; see canonical_broad_DRAFT.md."
        )
    adata.obs["cell_type_broad"] = broad.astype("category")
    logger.info(
        "[driver] cell_type_broad remapped: %d unique canonical labels",
        int(adata.obs["cell_type_broad"].nunique()),
    )

    # ---- scDRS-replica pre-filter (matches add_neighbors.py so
    #      perform-downstream can be invoked with --flag-filter-data False
    #      without a cell-set mismatch). filter_cells counts nonzero-X
    #      genes per cell; log1p(0)==0 so the cell set is identical to
    #      what scDRS would produce on raw counts.
    n0_cells, n0_genes = adata.n_obs, adata.n_vars
    sc.pp.filter_cells(adata, min_genes=SCDRS_MIN_GENES_PER_CELL)
    sc.pp.filter_genes(adata, min_cells=SCDRS_MIN_CELLS_PER_GENE)
    print(
        f"[driver] scDRS-replica filter: "
        f"{n0_cells}->{adata.n_obs} cells "
        f"(min_genes={SCDRS_MIN_GENES_PER_CELL}), "
        f"{n0_genes}->{adata.n_vars} genes "
        f"(min_cells={SCDRS_MIN_CELLS_PER_GENE})"
    )

    # ---- PCA + kNN so obsp['connectivities'] ships with the h5ad and
    #      scdrs perform-downstream skips its ~45-min rebuild. Params
    #      match scDRS group-analysis defaults (knn_n_pcs=20,
    #      knn_n_neighbors=15).
    sc.pp.pca(adata, n_comps=20)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20)
    assert "connectivities" in adata.obsp, (
        "sc.pp.neighbors did not populate obsp['connectivities']"
    )

    adata.uns["scdrs_prefilter"] = {
        "min_genes_per_cell": SCDRS_MIN_GENES_PER_CELL,
        "min_cells_per_gene": SCDRS_MIN_CELLS_PER_GENE,
        "n_pcs": 20,
        "n_neighbors": 15,
        "invoke_perform_downstream_with":
            "--flag-filter-data False --flag-raw-count False",
    }

    adata.write_h5ad(a.out_h5ad)
    print(
        f"[driver] wrote {a.out_h5ad}: {adata.n_obs} cells x "
        f"{adata.n_vars} genes; obsp['connectivities'] populated "
        f"(n_pcs=20, n_neighbors=15)."
    )

    # ---- Cov file build ----
    sample_col_name = _resolve_sample_col(adata.obs, a.sample_col)
    logger.info(
        "[driver] sample_col=%r  donor_col=%r", sample_col_name, DONOR_COL,
    )

    # Sample-in-donor nesting gate. Required so that dropping donor
    # from the cov is safe — sample dummies must fully absorb donor-
    # level variation. If sample_col == donor_col the check is vacuous
    # but harmless (every sample maps to exactly one donor by def).
    nest_df = pd.DataFrame({
        "sample": adata.obs[sample_col_name].astype(str).values,
        "donor":  adata.obs[DONOR_COL].astype(str).values,
    })
    donors_per_sample = nest_df.groupby("sample")["donor"].nunique()
    multi_donor_samples = donors_per_sample[donors_per_sample > 1]
    if len(multi_donor_samples) > 0:
        raise SystemExit(
            f"[driver] Sample-in-donor nesting violated: "
            f"{len(multi_donor_samples)} sample(s) span multiple donors "
            f"(sample_col={sample_col_name!r}, donor_col={DONOR_COL!r}). "
            f"Top 10 offenders (sample -> n_donors): "
            f"{multi_donor_samples.head(10).to_dict()}. "
            f"Fix: pass --sample-col=<truly-biopsy-level-column>, or "
            f"fall back to {DONOR_COL!r} (sample==donor)."
        )
    print(
        f"[driver] sample-nesting gate passed: "
        f"{donors_per_sample.size} samples across "
        f"{nest_df['donor'].nunique()} donors, each sample -> one donor."
    )

    X = _raw_counts_matrix(adata)
    n_counts = np.asarray(X.sum(axis=1)).ravel()
    n_genes  = np.asarray((X > 0).sum(axis=1)).ravel()

    # Pre-expand sample into explicit float 0.0/1.0 dummy columns
    # (sample_<id>). Matches smillie_covariates.tsv structure and
    # bypasses scDRS's category2dummy bool-dummy path that upcasts
    # df_cov.values to object dtype.
    sample_values = adata.obs[sample_col_name].astype(str).values
    sample_dummies = pd.get_dummies(
        pd.Series(sample_values, index=adata.obs_names, name="sample"),
        prefix="sample",
    ).astype(float)

    depth = pd.DataFrame({
        "const":        1.0,
        "log_n_genes":  np.log1p(n_genes),
        "log_n_counts": np.log1p(n_counts),
    }, index=adata.obs_names)
    cov = pd.concat([depth, sample_dummies], axis=1)
    cov.index.name = "cell"

    # Belt-and-suspenders: every column must be a numeric dtype so
    # scDRS's np.asarray(df_cov.values, dtype=float) doesn't hit
    # dtype('O'). Regression here would be caught in the driver, not
    # deep inside scDRS.
    non_numeric = [c for c in cov.columns
                   if not pd.api.types.is_numeric_dtype(cov[c])]
    if non_numeric:
        raise SystemExit(
            f"[driver] cov has non-numeric columns after pre-expansion: "
            f"{non_numeric}. scDRS will crash on df_cov.values upcast."
        )

    cov.to_csv(a.out_cov, sep="\t")
    print(
        f"[driver] wrote {a.out_cov}: {len(cov.columns)} numeric cols "
        f"(const + 2 depth + {sample_dummies.shape[1]} sample_<id> "
        f"float dummies from obs['{sample_col_name}']) — donor DROPPED "
        f"(absorbed by sample dummies); no disease/health per "
        f"no-disease-covariate rule."
    )


if __name__ == "__main__":
    main()
