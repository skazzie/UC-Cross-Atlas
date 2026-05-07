"""
Build the headline cross-atlas concordance table from cell-type test
statistic files.

Per the locked v1 plan (PLAN.md, DECISIONS.md):

- Headline metric is Spearman rho on cell-type-level Z-scores (scDRS) /
  regression coefficients (seismicGWAS), NOT p-values.
- Bootstrap 95% CIs (1000 iters, percentile, seed=42) on every reported
  rho.
- Concordance computed on the shared cell-type intersection per pair,
  with min cell-count threshold (>=50 cells in BOTH atlases).
- Kappa with marginal-saturation contingency: report kappa @ FDR<0.01
  if >=80% saturation at FDR<0.05.
- Top-k Jaccard at k=5,10 (broad) or k=5,10,20 (fine).

Input format: one TSV per (atlas, method, GWAS, tier) with columns
    cell_type   score   pval   qval   n_cells
where `score` is the larger-is-stronger statistic for that method
(scDRS: mean per-cell Z within cell type; seismicGWAS: regression
coefficient). `n_cells` is the cell count for that cell type in that
atlas (used by the min-cell-count filter).

Output: a long-format CSV where each row is one (atlas_a, atlas_b,
method, gwas, tier) comparison with Spearman rho + bootstrap 95% CI,
top-k Jaccard, Cohen's kappa with saturation flag, marginals, and
exclusion counts.

Usage:
    python compute_concordance.py \\
        --input \\
            results/scdrs/smillie_delange_broad/UC.cell_type.tsv:smillie:scdrs:delange:broad \\
            results/scdrs/kong_delange_broad/UC.cell_type.tsv:kong:scdrs:delange:broad \\
            results/scdrs/mennillo_delange_broad/UC.cell_type.tsv:mennillo:scdrs:delange:broad \\
            ... \\
        --out results/concordance/cross_atlas_table.csv
"""

import argparse
import itertools
from collections import defaultdict
from pathlib import Path

import pandas as pd

from metrics import DEFAULT_MIN_CELLS, DEFAULT_SEED, concordance


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--input", nargs="+", required=True,
        help="One or more entries of the form PATH:ATLAS:METHOD:GWAS:TIER",
    )
    p.add_argument("--out", required=True, help="Output CSV path")
    p.add_argument("--fdr-threshold", type=float, default=0.05)
    p.add_argument("--fdr-threshold-strict", type=float, default=0.01,
                   help="Saturation contingency threshold (PLAN.md, DECISIONS.md).")
    p.add_argument("--min-cells", type=int, default=DEFAULT_MIN_CELLS)
    p.add_argument("--bootstrap-iter", type=int, default=1000)
    p.add_argument("--seed", type=int, default=DEFAULT_SEED)
    p.add_argument("--cell-type-col", default="cell_type")
    p.add_argument("--score-col", default="score",
                   help="Larger-is-stronger statistic column. scDRS: mean per-cell Z; "
                        "seismicGWAS: regression coefficient.")
    p.add_argument("--pval-col", default="pval")
    p.add_argument("--qval-col", default="qval")
    p.add_argument("--n-cells-col", default="n_cells",
                   help="Column with the cell count per cell type. Used for min-cell filter.")
    return p.parse_args()


def load(path, cell_type_col, score_col, pval_col, qval_col, n_cells_col):
    df = pd.read_csv(path, sep="\t")
    for col in (cell_type_col, score_col, pval_col, qval_col):
        if col not in df.columns:
            raise SystemExit(f"{path}: missing column '{col}'. Have: {list(df.columns)}")
    scores = dict(zip(df[cell_type_col], df[score_col]))
    pvals = dict(zip(df[cell_type_col], df[pval_col]))
    qvals = dict(zip(df[cell_type_col], df[qval_col]))
    if n_cells_col in df.columns:
        n_cells = dict(zip(df[cell_type_col], df[n_cells_col]))
    else:
        n_cells = None
    return scores, pvals, qvals, n_cells


def main():
    args = parse_args()

    # data[(method, gwas, tier)][atlas] = (scores, pvals, qvals, n_cells)
    data = defaultdict(dict)
    for entry in args.input:
        parts = entry.split(":")
        if len(parts) != 5:
            raise SystemExit(f"Bad --input entry (need PATH:ATLAS:METHOD:GWAS:TIER): {entry}")
        path, atlas, method, gwas, tier = parts
        data[(method, gwas, tier)][atlas] = load(
            path, args.cell_type_col, args.score_col, args.pval_col,
            args.qval_col, args.n_cells_col,
        )
        print(f"[concordance] loaded {atlas}/{method}/{gwas}/{tier}: "
              f"{len(data[(method, gwas, tier)][atlas][0])} cell types", flush=True)

    rows = []
    for (method, gwas, tier), atlases in data.items():
        is_fine = tier == "fine"
        for a, b in itertools.combinations(sorted(atlases), 2):
            scores_a, pvals_a, qvals_a, counts_a = atlases[a]
            scores_b, pvals_b, qvals_b, counts_b = atlases[b]
            res = concordance(
                scores_a, scores_b, qvals_a, qvals_b,
                fdr_threshold=args.fdr_threshold,
                fdr_threshold_strict=args.fdr_threshold_strict,
                cell_counts_a=counts_a, cell_counts_b=counts_b,
                min_cells=args.min_cells,
                n_bootstrap=args.bootstrap_iter,
                seed=args.seed,
                is_fine_tier=is_fine,
                larger_is_stronger=True,
            )
            rows.append({
                "method": method,
                "gwas": gwas,
                "tier": tier,
                "atlas_a": a,
                "atlas_b": b,
                "n_common": res.n_common,
                "n_excluded_low_count": res.excluded_low_count,
                "spearman_rho": res.spearman_rho,
                "spearman_p": res.spearman_p,
                "spearman_ci_lo": res.spearman_ci_lo,
                "spearman_ci_hi": res.spearman_ci_hi,
                "jaccard_top5": res.jaccard_top5,
                "jaccard_top10": res.jaccard_top10,
                "jaccard_top20": res.jaccard_top20,
                "kappa": res.kappa,
                "kappa_threshold": res.kappa_threshold,
                "kappa_saturation_contingency": res.kappa_threshold_used_due_to_saturation,
                "n_sig_a": res.n_sig_a,
                "n_sig_b": res.n_sig_b,
            })

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)
    print(f"[concordance] wrote {out} ({len(df)} rows)", flush=True)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
