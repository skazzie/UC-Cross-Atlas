"""Driver: load Smillie (SCP259), write h5ad + covariate file.
sex NOT in SCP259 metadata (header: Subject/Health/Location/Sample) -> omitted.
Covariates: log_n_genes, log_n_counts, donor, sample (sample != donor here).
Dummies pre-numericized to float64 (scdrs 1.0.2 + modern pandas bool-dummy bug).
"""
import argparse, numpy as np, pandas as pd, scanpy as sc
from load_smillie import load

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--scp259-dir", required=True)
    p.add_argument("--out-h5ad", required=True)
    p.add_argument("--out-cov", required=True)
    a = p.parse_args()

    adata = load(a.scp259_dir)                       # apply_v1_filter=True default
    adata.write_h5ad(a.out_h5ad)
    print(f"[driver] wrote {a.out_h5ad}: {adata.n_obs} cells x {adata.n_vars} genes")
    print(f"[driver] donors={adata.obs['donor'].nunique()} "
          f"samples={adata.obs['sample'].nunique() if 'sample' in adata.obs else 'NA'}")

    X = adata.layers["counts"] if "counts" in adata.layers else adata.X
    n_counts = np.asarray(X.sum(axis=1)).ravel()
    n_genes  = np.asarray((X > 0).sum(axis=1)).ravel()
    cov = pd.DataFrame({
        "const": 1.0,
        "log_n_genes":  np.log1p(n_genes),
        "log_n_counts": np.log1p(n_counts),
        "donor":  adata.obs["donor"].astype(str).values,
        "sample": adata.obs["sample"].astype(str).values,
    }, index=adata.obs_names)
    # pre-numericize categorical dummies to float64 (scdrs 1.0.2 needs numeric cov matrix)
    cov = pd.get_dummies(cov, columns=["donor", "sample"], drop_first=True).astype(np.float64)
    cov.index.name = "cell"
    cov.to_csv(a.out_cov, sep="\t")
    print(f"[driver] wrote {a.out_cov}: {len(cov.columns)} cols (sex UNAVAILABLE)")

if __name__ == "__main__":
    main()
