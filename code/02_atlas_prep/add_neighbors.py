"""Retrofit script: pre-filter + precompute PCA + kNN graph on an atlas h5ad.

The naive version (just PCA + neighbors) did not skip the ~45-min rebuild
inside ``scdrs perform-downstream`` because scDRS's load chain
(``sc.pp.filter_cells`` → ``filter_genes`` → ``normalize_per_cell`` →
``log1p``, all via ``scdrs.util.load_h5ad`` at ``scdrs/util.py:82-88``)
was silently dropping ``obsp['connectivities']`` between the disk read
and the ``if "connectivities" not in adata.obsp`` guard at
``bin/scdrs:669``. Cell-count mismatch wasn't the cause — the guard only
checks presence, and ``downstream_group_analysis`` intersects cells via
``set(adata.obs_names) & set(df_full_score.index)`` at
``scdrs/method.py:748`` so no strict size-match is required. The actual
culprit was the ``_inplace_subset_obs`` / view→copy chain across four
sequential scanpy calls on a 365k×365k sparse obsp.

The clean fix, baked in here:

1. Replicate scDRS's exact load-time filter (``min_genes=250``,
   ``min_cells=50``) in this script BEFORE building the graph, so the
   on-disk cell set matches what scDRS's ``flag_filter_data=True`` would
   produce.
2. Run PCA + neighbors on the filtered adata and save.
3. Invoke perform-downstream with ``--flag-filter-data=False`` AND
   ``--flag-raw-count=False``. ``scdrs.util.load_h5ad`` then degrades to
   a bare ``read_h5ad`` — no scanpy chain, obsp survives verbatim, guard
   at ``bin/scdrs:669`` does not trigger. ``downstream_group_analysis``
   does not touch adata.X, so X's normalization state is irrelevant to
   the group-analysis output (``method.py:705-741``).

Parameters match scDRS group-analysis defaults (``knn_n_pcs=20``,
``knn_n_neighbors=15``).

Preserves in place: ``layers['counts']``, ``obs`` (filtered), ``var``
(filtered). Adds: ``obsm['X_pca']``, ``varm['PCs']``, ``uns['pca']``,
``obsp['connectivities']``, ``obsp['distances']``, ``uns['neighbors']``.
Also stamps ``uns['scdrs_prefilter']`` so future changes to scDRS
defaults are catchable.

Usage
-----
    python add_neighbors.py --h5ad path/to/atlas.h5ad

Then invoke scDRS as:

    scdrs perform-downstream \\
        --h5ad-file path/to/atlas.h5ad \\
        --score-file <...>/UC.full_score.gz \\
        --out-folder <...> \\
        --group-analysis cell_type_broad,cell_type_fine \\
        --flag-filter-data False \\
        --flag-raw-count False
"""
import argparse
import scanpy as sc

# scDRS 1.0.2 load defaults — hardcoded so a version bump that changes
# them is caught by comparing this stamp against scdrs.util.load_h5ad.
SCDRS_MIN_GENES_PER_CELL = 250
SCDRS_MIN_CELLS_PER_GENE = 50


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--h5ad", required=True, help="atlas h5ad; re-saved in place")
    p.add_argument("--n-pcs", type=int, default=20, help="PCs (scDRS default)")
    p.add_argument("--n-neighbors", type=int, default=15, help="kNN k (scDRS default)")
    a = p.parse_args()

    adata = sc.read_h5ad(a.h5ad)
    print(f"[add_neighbors] loaded {a.h5ad}: {adata.n_obs} cells x {adata.n_vars} genes")

    # Pre-flight checks (fail before the ~10-min PCA on Smillie-scale input).
    if "counts" not in adata.layers:
        raise KeyError(
            "layers['counts'] missing — loader invariant violated; refusing to "
            "run so the raw-count layer isn't silently lost on re-save."
        )

    # Replicate scDRS's flag_filter_data=True behavior so the saved cell
    # set matches what scDRS would have filtered to. filter_cells counts
    # genes with nonzero X per cell; identical result whether X is raw
    # or log1p (log1p(0)=0), so we don't need to touch X first.
    n0_cells, n0_genes = adata.n_obs, adata.n_vars
    sc.pp.filter_cells(adata, min_genes=SCDRS_MIN_GENES_PER_CELL)
    sc.pp.filter_genes(adata, min_cells=SCDRS_MIN_CELLS_PER_GENE)
    print(
        f"[add_neighbors] scDRS-replica filter: "
        f"{n0_cells}→{adata.n_obs} cells (min_genes={SCDRS_MIN_GENES_PER_CELL}), "
        f"{n0_genes}→{adata.n_vars} genes (min_cells={SCDRS_MIN_CELLS_PER_GENE})"
    )

    if a.n_pcs > min(adata.shape) - 1:
        raise ValueError(
            f"--n-pcs {a.n_pcs} > min(adata.shape)-1 = {min(adata.shape) - 1}"
        )

    sc.pp.pca(adata, n_comps=a.n_pcs)
    sc.pp.neighbors(adata, n_neighbors=a.n_neighbors, n_pcs=a.n_pcs)
    assert "connectivities" in adata.obsp, "sc.pp.neighbors did not populate obsp"

    # Provenance stamp — lets future audits catch a scDRS default drift.
    adata.uns["scdrs_prefilter"] = {
        "min_genes_per_cell": SCDRS_MIN_GENES_PER_CELL,
        "min_cells_per_gene": SCDRS_MIN_CELLS_PER_GENE,
        "n_pcs": a.n_pcs,
        "n_neighbors": a.n_neighbors,
        "invoke_perform_downstream_with": "--flag-filter-data False --flag-raw-count False",
    }

    adata.write_h5ad(a.h5ad)
    print(
        f"[add_neighbors] re-saved {a.h5ad}: "
        f"{adata.n_obs} cells x {adata.n_vars} genes, "
        f"obsp['connectivities'] populated "
        f"(n_pcs={a.n_pcs}, n_neighbors={a.n_neighbors})."
    )
    print(
        "[add_neighbors] IMPORTANT: invoke scdrs perform-downstream with "
        "`--flag-filter-data False --flag-raw-count False` so scDRS's "
        "load chain doesn't strip obsp."
    )


if __name__ == "__main__":
    main()
