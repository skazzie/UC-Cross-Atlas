"""Dump an h5ad atlas to flat files that R can read without zellkonverter.

Windows-laptop workaround for the failed `zellkonverter` install: that R
package wants Rtools 4.5 to source-build the basilisk + reticulate
machinery, and Rtools is a separate Windows installer we don't auto-pull.
This script bypasses the whole problem by dumping the h5ad into a
directory of standard files:

  <out_dir>/
    X.mtx              log1p(CP10k) sparse matrix in MatrixMarket format
                       (cells x genes — R-side script transposes to
                       genes x cells for SCE convention)
    counts.mtx         raw integer counts sparse matrix, same shape
    cells.tsv          obs metadata (cell barcode index + all obs cols)
    genes.tsv          var metadata (gene symbol index + var cols)
    layout.json        shapes + which-axis-is-cells, for the R loader

Then `code/04_seismic/sce_from_flat.R` reads these and builds a
SingleCellExperiment with assays `counts` + `logcounts` and colData /
rowData matching the AnnData obs / var. seismicGWAS consumes that
SCE directly via calc_specificity().

Usage:
  py.exe scripts/h5ad_to_flat.py data/atlases/garrido.h5ad data/atlases/garrido_flat/
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("h5ad", type=Path)
    p.add_argument("out_dir", type=Path)
    args = p.parse_args()

    print(f"[h5ad_to_flat] reading {args.h5ad}", flush=True)
    adata = ad.read_h5ad(args.h5ad)
    print(f"[h5ad_to_flat] shape: {adata.n_obs} cells x {adata.n_vars} genes", flush=True)

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # X — log-normalized expression (cells x genes).
    X = adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    print(f"[h5ad_to_flat] writing X.mtx ({X.nnz:,} nnz)", flush=True)
    scipy.io.mmwrite(str(args.out_dir / "X.mtx"), X, field="real")

    # counts layer (raw integer counts) if present.
    if "counts" in adata.layers:
        counts = adata.layers["counts"]
        if not sp.issparse(counts):
            counts = sp.csr_matrix(counts)
        print(f"[h5ad_to_flat] writing counts.mtx ({counts.nnz:,} nnz)", flush=True)
        scipy.io.mmwrite(str(args.out_dir / "counts.mtx"), counts, field="integer")
    else:
        print("[h5ad_to_flat] WARNING: no 'counts' layer — raw counts not exported")

    # obs (cell metadata).
    obs = adata.obs.copy()
    obs.index.name = "cell_id"
    print(f"[h5ad_to_flat] writing cells.tsv ({obs.shape})", flush=True)
    obs.to_csv(args.out_dir / "cells.tsv", sep="\t", index=True)

    # var (gene metadata) — var_names go into 'gene' column; var index used as primary key.
    var = adata.var.copy()
    var.index.name = "gene"
    print(f"[h5ad_to_flat] writing genes.tsv ({var.shape})", flush=True)
    var.to_csv(args.out_dir / "genes.tsv", sep="\t", index=True)

    # Layout sidecar — tells the R loader the orientation + assay names available.
    layout = {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "X_orientation": "cells_x_genes",   # both .mtx files are cells x genes
        "X_assay_name": "logcounts",         # what to call X when building SCE
        "has_counts_layer": "counts" in adata.layers,
        "obs_columns": list(obs.columns),
        "var_columns": list(var.columns),
    }
    (args.out_dir / "layout.json").write_text(json.dumps(layout, indent=2))
    print(f"[h5ad_to_flat] wrote {args.out_dir / 'layout.json'}", flush=True)
    print("[h5ad_to_flat] done", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
