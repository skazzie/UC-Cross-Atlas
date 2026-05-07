"""
Build the headline cross-atlas concordance table from cell-type p-value files.

Input: one TSV per (atlas, method) with columns
    cell_type<tab>pval<tab>qval
(scDRS .scdrs_group output can be massaged into this with a one-liner; same for
seismicGWAS and scPagwas.)

Output: a long-format CSV where each row is one (atlas_a, atlas_b, method)
comparison with Spearman ρ + 95% CI, top-5/top-10 Jaccard, Cohen's κ, and
significance marginals — i.e., the headline table from spec §2.5.

Usage:
    python compute_concordance.py \
        --input results/scdrs/smillie/UC.cell_type.tsv:smillie:scdrs \
                results/scdrs/kong/UC.cell_type.tsv:kong:scdrs \
                results/scdrs/mennillo/UC.cell_type.tsv:mennillo:scdrs \
                results/seismic/smillie/UC.cell_type.tsv:smillie:seismic \
                ... \
        --out results/concordance/headline_table.csv
"""

import argparse
import itertools
from collections import defaultdict
from pathlib import Path

import pandas as pd

from metrics import bootstrap_spearman_ci, concordance


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", nargs="+", required=True,
                   help="One or more entries of the form PATH:ATLAS:METHOD")
    p.add_argument("--out", required=True, help="Output CSV path")
    p.add_argument("--fdr-threshold", type=float, default=0.05)
    p.add_argument("--bootstrap-iter", type=int, default=1000)
    p.add_argument("--cell-type-col", default="cell_type")
    p.add_argument("--pval-col", default="pval")
    p.add_argument("--qval-col", default="qval")
    return p.parse_args()


def load(path, cell_type_col, pval_col, qval_col):
    df = pd.read_csv(path, sep="\t")
    for col in (cell_type_col, pval_col, qval_col):
        if col not in df.columns:
            raise SystemExit(f"{path}: missing column '{col}'. Have: {list(df.columns)}")
    pvals = dict(zip(df[cell_type_col], df[pval_col]))
    qvals = dict(zip(df[cell_type_col], df[qval_col]))
    return pvals, qvals


def main():
    args = parse_args()

    # data[method][atlas] = (pvals, qvals)
    data: dict[str, dict[str, tuple[dict, dict]]] = defaultdict(dict)
    for entry in args.input:
        parts = entry.split(":")
        if len(parts) != 3:
            raise SystemExit(f"Bad --input entry (need PATH:ATLAS:METHOD): {entry}")
        path, atlas, method = parts
        data[method][atlas] = load(path, args.cell_type_col, args.pval_col, args.qval_col)
        print(f"[concordance] loaded {atlas}/{method}: {len(data[method][atlas][0])} cell types", flush=True)

    rows = []
    for method, atlases in data.items():
        for a, b in itertools.combinations(sorted(atlases), 2):
            pa, qa = atlases[a]
            pb, qb = atlases[b]
            res = concordance(pa, pb, qa, qb, fdr_threshold=args.fdr_threshold)
            ci_lo, ci_hi = bootstrap_spearman_ci(pa, pb, n_iter=args.bootstrap_iter)
            rows.append({
                "method": method,
                "atlas_a": a,
                "atlas_b": b,
                "n_common": res.n_common,
                "spearman_rho": res.spearman_rho,
                "spearman_p": res.spearman_p,
                "spearman_ci_lo": ci_lo,
                "spearman_ci_hi": ci_hi,
                "jaccard_top5": res.jaccard_top5,
                "jaccard_top10": res.jaccard_top10,
                "kappa": res.kappa,
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
