"""Driver: load Garrido-Trigo, write h5ad + partial covariate file.
sex is NOT available in GSE214695 (not in RAW.tar or annotation CSV);
cov file omits it. Source per-donor sex externally for DECISIONS-locked runs.
"""
import argparse, numpy as np, pandas as pd, scanpy as sc
from load_garrido_trigo import load

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tar", required=True)
    p.add_argument("--csv", required=True)
    p.add_argument("--out-h5ad", required=True)
    p.add_argument("--out-cov", required=True)
    a = p.parse_args()

    adata = load(a.tar, a.csv)                      # apply_v1_filter=True default
    adata.write_h5ad(a.out_h5ad)
    print(f"[driver] wrote {a.out_h5ad}: {adata.n_obs} cells x {adata.n_vars} genes")

    X = adata.layers["counts"]                      # raw counts for depth covars
    n_counts = np.asarray(X.sum(axis=1)).ravel()
    n_genes  = np.asarray((X > 0).sum(axis=1)).ravel()
    cov = pd.DataFrame({
        "const": 1,
        "log_n_genes":  np.log1p(n_genes),
        "log_n_counts": np.log1p(n_counts),
        "donor":  adata.obs["donor"].astype(str).values,
        "sample": adata.obs["sample"].astype(str).values,
    }, index=adata.obs_names)
    cov.index.name = "cell"
    cov.to_csv(a.out_cov, sep="\t")                 # sex omitted — see module docstring
    print(f"[driver] wrote {a.out_cov}: {list(cov.columns)} (sex UNAVAILABLE)")

if __name__ == "__main__":
    main()
