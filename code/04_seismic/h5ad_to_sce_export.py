"""Export an h5ad atlas to flat files so R can build a SingleCellExperiment
WITHOUT zellkonverter/basilisk.

zellkonverter::readH5AD routes through basilisk, which ignores
RETICULATE_PYTHON and tries to source-build Python 3.14 — that fails on
our sudo-less GCP VM. This script sidesteps the whole chain: read the
h5ad in the conda env with anndata, drop flat files, and hand off to
sce_from_export.R (pure R, no Python called).

Our atlas h5ads carry raw counts in layers['counts'] and log1p(CP10k)
in X. We export BOTH:
  - counts.mtx from layers['counts'] (raw), so a `counts` assay is
    available for anything that expects raw.
  - logcounts.mtx from X (log1p(CP10k)), so seismic's calc_specificity
    (which defaults to assay_name="logcounts") consumes the SAME
    normalization our scDRS runs use. Keeping the two methods on
    identical input means cross-method concordance reflects method
    differences, not normalization differences.
We refuse to fall back to X for counts (would silently feed
log-normalized values into anything expecting raw).

Outputs (in --out-dir):
  counts.mtx     scipy MatrixMarket, cells x genes, raw integer counts
                 (from layers['counts'])
  logcounts.mtx  scipy MatrixMarket, cells x genes, log1p(CP10k) values
                 (from adata.X)
  genes.tsv      one gene symbol per line, no header (matches barcodes.tsv)
  barcodes.tsv   one cell id per line, no header, row-aligned with counts.mtx
  obs.tsv        cell metadata with header; first column 'cell_id' is the
                 barcode. cell_type_broad, cell_type_fine, donor required
                 (seismic keys off these).

Usage:
  python code/04_seismic/h5ad_to_sce_export.py \\
      --h5ad data/atlases/garrido_trigo.h5ad \\
      --out-dir data/atlases/garrido_trigo_export/
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import scipy.io
import scipy.sparse as sp


REQUIRED_OBS = ("cell_type_broad", "cell_type_fine", "donor")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--h5ad", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    args = p.parse_args()

    print(f"[h5ad_to_sce_export] reading {args.h5ad}", flush=True)
    adata = ad.read_h5ad(args.h5ad)
    print(
        f"[h5ad_to_sce_export] shape: {adata.n_obs} cells x {adata.n_vars} genes",
        flush=True,
    )

    if "counts" not in adata.layers:
        raise SystemExit(
            "ERROR: layers['counts'] missing — refusing to export X as counts "
            "(seismic requires RAW counts; X is log1p-normalized in our atlases)."
        )

    missing_obs = [c for c in REQUIRED_OBS if c not in adata.obs.columns]
    if missing_obs:
        raise SystemExit(
            f"ERROR: obs missing required columns: {missing_obs}. "
            f"Available: {list(adata.obs.columns)}"
        )

    args.out_dir.mkdir(parents=True, exist_ok=True)

    counts = adata.layers["counts"]
    if not sp.issparse(counts):
        counts = sp.csr_matrix(counts)
    print(
        f"[h5ad_to_sce_export] writing counts.mtx ({counts.nnz:,} nnz, "
        f"cells x genes)",
        flush=True,
    )
    scipy.io.mmwrite(str(args.out_dir / "counts.mtx"), counts, field="integer")

    # X is log1p(CP10k) in our atlases — the same normalization scDRS
    # consumes. Exporting it as logcounts.mtx lets seismic's
    # calc_specificity(assay_name="logcounts") run on identical input,
    # so cross-method concordance isolates method effects from
    # normalization effects.
    logcounts = adata.X
    if not sp.issparse(logcounts):
        logcounts = sp.csr_matrix(logcounts)
    print(
        f"[h5ad_to_sce_export] writing logcounts.mtx ({logcounts.nnz:,} nnz, "
        f"cells x genes)",
        flush=True,
    )
    scipy.io.mmwrite(str(args.out_dir / "logcounts.mtx"), logcounts, field="real")

    barcodes_path = args.out_dir / "barcodes.tsv"
    barcodes_path.write_text("\n".join(adata.obs_names.astype(str)) + "\n")
    print(f"[h5ad_to_sce_export] wrote {barcodes_path} ({adata.n_obs} lines)", flush=True)

    genes_path = args.out_dir / "genes.tsv"
    genes_path.write_text("\n".join(adata.var_names.astype(str)) + "\n")
    print(f"[h5ad_to_sce_export] wrote {genes_path} ({adata.n_vars} lines)", flush=True)

    obs = adata.obs.copy()
    obs.index.name = "cell_id"
    obs_path = args.out_dir / "obs.tsv"
    obs.to_csv(obs_path, sep="\t", index=True)
    print(f"[h5ad_to_sce_export] wrote {obs_path} ({obs.shape})", flush=True)

    print("[h5ad_to_sce_export] done", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
