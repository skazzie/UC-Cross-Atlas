#!/usr/bin/env python3
"""MAGMA gene-property regression on Smillie at broad tier under de Lange.

Spec: see ./README.md and DECISIONS.md.

Per broad cell type:
  1. mean log-normalized expression per gene
  2. OLS of MAGMA gene-Z ~ mean expression + confounders
     (gene_length_log, gene_gene_ld, transcript_count from the MAGMA TSV)
  3. BH-FDR across cell types

Output TSV columns:
  cell_type, coefficient, se, pvalue, fdr, n_genes
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

LOGGER = logging.getLogger(Path(__file__).stem)

_REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO / "code"))

from _shared.constants import SEED  # noqa: E402
from _shared.git import git_sha  # noqa: E402


def require_path(p: Path, descr: str) -> None:
    if not p.exists():
        LOGGER.error("Missing input: %s (%s)", p, descr)
        raise SystemExit(2)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--atlas", default="smillie",
                        help="Atlas slug (default: smillie)")
    parser.add_argument("--h5ad-path", type=Path, required=True)
    parser.add_argument("--magma-z", type=Path, required=True,
                        help="MAGMA gene-Z TSV with columns: gene, z, "
                             "gene_length_log, gene_gene_ld, transcript_count")
    parser.add_argument("--tier", default="broad", choices=["broad", "fine"])
    parser.add_argument("--gwas", default="delange")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s: %(message)s",
    )

    require_path(args.h5ad_path, "atlas h5ad")
    require_path(args.magma_z, "MAGMA gene-Z TSV")

    t0 = time.time()

    import numpy as np
    import pandas as pd
    import scanpy as sc
    import statsmodels.api as sm
    from statsmodels.stats.multitest import multipletests

    LOGGER.info("Loading atlas %s from %s", args.atlas, args.h5ad_path)
    adata = sc.read_h5ad(args.h5ad_path)
    ct_col = f"cell_type_{args.tier}"
    if ct_col not in adata.obs.columns:
        LOGGER.error("obs has no '%s'; available: %s",
                     ct_col, list(adata.obs.columns))
        raise SystemExit(2)
    LOGGER.info("Atlas: %d cells x %d genes", adata.n_obs, adata.n_vars)

    LOGGER.info("Loading MAGMA gene-Z from %s", args.magma_z)
    magma = pd.read_csv(args.magma_z, sep="\t")
    required = {"gene", "z", "gene_length_log", "gene_gene_ld", "transcript_count"}
    missing = required - set(magma.columns)
    if missing:
        LOGGER.error("MAGMA TSV missing columns: %s", missing)
        raise SystemExit(2)
    magma = magma.dropna(subset=["z"]).copy()
    LOGGER.info("MAGMA: %d genes", len(magma))

    # Mean log-normalized expression per gene per cell type (broad).
    LOGGER.info("Computing per-cell-type mean expression on %s", ct_col)
    labels = adata.obs[ct_col].astype("category")
    X = adata.X
    if hasattr(X, "toarray"):
        # Sparse: compute group means via a one-hot accumulator to avoid
        # materializing the dense matrix.
        cats = labels.cat.categories
        means = np.zeros((len(cats), adata.n_vars), dtype=np.float32)
        for i, cat in enumerate(cats):
            mask = (labels == cat).to_numpy()
            sub = X[mask]
            sums = np.asarray(sub.sum(axis=0)).ravel()
            means[i] = sums / max(mask.sum(), 1)
    else:
        cats = labels.cat.categories
        means = np.zeros((len(cats), adata.n_vars), dtype=np.float32)
        for i, cat in enumerate(cats):
            mask = (labels == cat).to_numpy()
            means[i] = np.asarray(X[mask]).mean(axis=0)

    gene_names = pd.Index(adata.var_names)
    common = gene_names.intersection(magma["gene"])
    LOGGER.info("Genes in both atlas and MAGMA: %d", len(common))
    if len(common) < 100:
        LOGGER.error("Too few overlapping genes (%d); check HGNC remap.",
                     len(common))
        raise SystemExit(2)

    gene_pos = {g: i for i, g in enumerate(gene_names)}
    magma_common = magma[magma["gene"].isin(common)].copy()
    magma_common = magma_common.set_index("gene").loc[list(common)]
    cols = [gene_pos[g] for g in common]

    rows = []
    for i, cat in enumerate(cats):
        expr = means[i, cols]
        X_design = pd.DataFrame({
            "expr": expr,
            "gene_length_log": magma_common["gene_length_log"].to_numpy(),
            "gene_gene_ld":    magma_common["gene_gene_ld"].to_numpy(),
            "transcript_count": magma_common["transcript_count"].to_numpy(),
        })
        X_design = sm.add_constant(X_design)
        y = magma_common["z"].to_numpy()
        model = sm.OLS(y, X_design, missing="drop").fit()
        rows.append({
            "cell_type": str(cat),
            "coefficient": float(model.params["expr"]),
            "se":          float(model.bse["expr"]),
            "pvalue":      float(model.pvalues["expr"]),
            "n_genes":     int(model.nobs),
        })

    result = pd.DataFrame(rows)
    _, fdr, _, _ = multipletests(result["pvalue"], method="fdr_bh")
    result["fdr"] = fdr
    result["tool_version"] = "statsmodels-ols"
    result["git_sha"] = git_sha()
    result["atlas"] = args.atlas
    result["tier"] = args.tier
    result["gwas"] = args.gwas

    # Sanity track per README: UC-relevant lineages should appear in top 5.
    KNOWN_UC_RELEVANT_BROAD = {"T cell", "Plasma", "Myeloid", "Enterocyte"}
    top_5 = (
        result.nlargest(5, "coefficient")["cell_type"].astype(str).tolist()
    )
    overlap = {ct for ct in top_5
               if any(known.lower() in ct.lower()
                      for known in KNOWN_UC_RELEVANT_BROAD)}
    LOGGER.info("UC-relevant cell types in top 5: %s", sorted(overlap))
    if not overlap:
        LOGGER.warning(
            "No UC-relevant cell types in top 5; pipeline may have a bug "
            "(informational per README sanity-track caveat)."
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.out, sep="\t", index=False)
    LOGGER.info("Wrote %s (%d rows) in %.1fs",
                args.out, len(result), time.time() - t0)
    return 0


if __name__ == "__main__":
    sys.exit(main())
