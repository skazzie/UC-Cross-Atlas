"""Retrofit script: precompute PCA + kNN graph on an existing atlas h5ad.

Without this, ``scdrs perform-downstream --group-analysis`` rebuilds the
graph on every invocation (~45 min on Smillie-scale atlases; prints
``connectivities not found in adata.obsp; run sc.pp.neighbors first``).

Parameters match scDRS group-analysis defaults (``knn_n_pcs=20``,
``knn_n_neighbors=15``) so ``adata.obsp['connectivities']`` is populated
identically to what scDRS would build itself, and scDRS then skips the
rebuild via the ``if "connectivities" not in adata.obsp`` guard
(scdrs_cli.py:669).

Preserves in place: ``layers['counts']``, ``obs``, ``var``. Adds:
``obsm['X_pca']``, ``varm['PCs']``, ``uns['pca']``, ``obsp['connectivities']``,
``obsp['distances']``, ``uns['neighbors']``.

Usage
-----
    python add_neighbors.py --h5ad path/to/atlas.h5ad
"""
import argparse
import scanpy as sc


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
    if a.n_pcs > min(adata.shape) - 1:
        raise ValueError(
            f"--n-pcs {a.n_pcs} > min(adata.shape)-1 = {min(adata.shape) - 1}"
        )

    sc.pp.pca(adata, n_comps=a.n_pcs)
    sc.pp.neighbors(adata, n_neighbors=a.n_neighbors, n_pcs=a.n_pcs)
    assert "connectivities" in adata.obsp, "sc.pp.neighbors did not populate obsp"

    adata.write_h5ad(a.h5ad)
    print(
        f"[add_neighbors] re-saved {a.h5ad} with obsp['connectivities'] "
        f"(n_pcs={a.n_pcs}, n_neighbors={a.n_neighbors})"
    )


if __name__ == "__main__":
    main()
