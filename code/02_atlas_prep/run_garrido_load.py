"""Driver: load Garrido-Trigo, write h5ad + covariate file.

Loader output (verified 2026-10-07 against VM run in CONTINUITY):
  29,675 cells x 22,414 genes x 12 donors

sex is NOT available in GSE214695 (not in RAW.tar or annotation CSV);
cov file omits it. Source per-donor sex externally for DECISIONS-locked
runs. See module docstring in ``load_garrido_trigo.py``.

Covariates (per 2026-10-07 directive): float 0/1 dummies pre-expanded,
sample kept, donor dropped, disease status never included. For Garrido
donor == sample 1:1 per CONTINUITY, so sample dummies carry the
donor-level batch structure.
"""
import argparse, numpy as np, pandas as pd, scanpy as sc
from load_garrido_trigo import load

# scDRS filter constants — stamped for provenance and so a future scDRS
# default drift is catchable by comparing against scdrs.util.load_h5ad's
# constants. Not applied here; see note in the stamp.
SCDRS_MIN_GENES_PER_CELL = 250
SCDRS_MIN_CELLS_PER_GENE = 50


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tar", required=True)
    p.add_argument("--csv", required=True)
    p.add_argument("--out-h5ad", required=True)
    p.add_argument("--out-cov", required=True)
    a = p.parse_args()

    adata = load(a.tar, a.csv)                      # apply_v1_filter=True default
    print(f"[driver] loader output: {adata.n_obs} cells x {adata.n_vars} genes")

    # Precompute PCA + kNN so obsp['connectivities'] ships with the h5ad
    # and scdrs perform-downstream skips its ~45-min rebuild. Params
    # match scDRS group-analysis defaults (knn_n_pcs=20, knn_n_neighbors=15).
    sc.pp.pca(adata, n_comps=20)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20)
    assert "connectivities" in adata.obsp, "sc.pp.neighbors did not populate obsp"

    # Provenance stamp. NOTE: unlike TAURUS/Pan-GI/HCA loaders, Garrido's
    # loader does NOT apply sc.pp.filter_cells / sc.pp.filter_genes here —
    # the published Garrido matrix is pre-QC'd and the gate is on the
    # loader-raw shape (29,675 x 22,414). The scDRS stage asserts the
    # post-filter shape (28,706 x 14,578) separately. scdrs is still
    # invoked with --flag_filter_data=False (loader kNN is authoritative).
    adata.uns["scdrs_prefilter"] = {
        "min_genes_per_cell": SCDRS_MIN_GENES_PER_CELL,
        "min_cells_per_gene": SCDRS_MIN_CELLS_PER_GENE,
        "n_pcs": 20,
        "n_neighbors": 15,
        "filter_applied_in_loader": False,
        "invoke_perform_downstream_with":
            "--flag-filter-data False --flag-raw-count False",
    }

    adata.write_h5ad(a.out_h5ad)
    print(f"[driver] wrote {a.out_h5ad}: {adata.n_obs} cells x {adata.n_vars} genes")

    X = adata.layers["counts"]                      # raw counts for depth covars
    n_counts = np.asarray(X.sum(axis=1)).ravel()
    n_genes  = np.asarray((X > 0).sum(axis=1)).ravel()

    # Pre-expand sample as float 0/1 dummies; drop donor (redundant with
    # sample on Garrido 1:1). drop_first=True avoids collinearity with
    # the const column.
    sample_dummies = pd.get_dummies(
        adata.obs["sample"].astype(str),
        prefix="sample",
        drop_first=True,
    ).astype(float)
    cov = pd.concat(
        [
            pd.DataFrame(
                {
                    "const": 1.0,
                    "log_n_genes": np.log1p(n_genes),
                    "log_n_counts": np.log1p(n_counts),
                },
                index=adata.obs_names,
            ),
            sample_dummies,
        ],
        axis=1,
    )
    cov.index.name = "cell"
    cov.to_csv(a.out_cov, sep="\t")
    print(
        f"[driver] wrote {a.out_cov}: const + log_n_genes + log_n_counts + "
        f"{sample_dummies.shape[1]} sample dummies (sex UNAVAILABLE; donor dropped)"
    )


if __name__ == "__main__":
    main()
