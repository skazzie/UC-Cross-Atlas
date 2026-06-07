"""
Convert a GWAS summary statistics file into the two MAGMA inputs:
  - <prefix>.snp.loc : rsid, chr, bp  (whitespace-separated, no header)
  - <prefix>.pval    : SNP, P, N      (whitespace-separated, with header)

Locked v1 GWAS in this pipeline:
  - de Lange 2017 UC (primary)              GWAS Catalog GCST004131
  - Liu 2023 multi-ancestry IBD (UC arm)    cross-GWAS sensitivity
  - Trubetskoy 2022 schizophrenia            negative control on Smillie

Pre-committed filters applied (DECISIONS.md):
  - autosomes only (chr 1-22), X chromosome excluded
  - drop rows with missing rsid/chr/bp/p
  - drop p outside (0, 1]
  - drop MAF < 0.01 (if FRQ column is provided)
  - drop INFO < 0.6 (if INFO column is provided)

The Liu 2023 download may include a per-SNP N column; if so, pass
--col-n. If absent, pass --n-fixed and document the fixed-N approximation
in DECISIONS.md and Methods.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

AUTOSOMES = {str(c) for c in range(1, 23)}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="GWAS summary statistics (.tsv/.tsv.gz)")
    p.add_argument("--out-prefix", required=True, help="Output prefix; writes <prefix>.snp.loc and <prefix>.pval")
    p.add_argument("--col-snp", default="hm_rsid")
    p.add_argument("--col-chr", default="hm_chrom")
    p.add_argument("--col-bp", default="hm_pos")
    p.add_argument("--col-p", default="p_value")
    p.add_argument("--col-n", default=None, help="Column with per-SNP N. If omitted, --n-fixed is used.")
    p.add_argument("--n-fixed", type=int, default=None, help="Single N value applied to every SNP if --col-n is absent.")
    p.add_argument("--col-frq", default=None, help="Allele frequency column (optional, used for MAF filter)")
    p.add_argument("--col-info", default=None, help="Imputation INFO column (optional, used for INFO filter)")
    p.add_argument("--maf-min", type=float, default=0.01)
    p.add_argument("--info-min", type=float, default=0.60)
    p.add_argument("--sep", default="\t")
    p.add_argument(
        "--keep-non-autosomes",
        action="store_true",
        help="Disable the autosome-only filter (default: drop X/Y/MT/non-1-22). "
             "Locked v1 keeps autosomes only; only set this for exploratory analyses.",
    )
    p.add_argument(
        "--n-min-frac",
        type=float,
        default=0.67,
        help="Drop SNPs whose per-SNP N is below this fraction of max(N). "
             "Matches LDSC munge_sumstats default (0.67); discards the noisy "
             "low-N tail from meta-analyses where some variants are tested in "
             "small sub-cohorts. No-op when --n-fixed is used (every SNP shares "
             "the same N). Set to 0 to disable. See DECISIONS correction 26.",
    )
    p.add_argument(
        "--lambda-gc-out",
        default=None,
        help="If set, compute genomic inflation factor lambda_GC and write to this path "
             "as a single-line TSV: gwas\\tlambda_gc\\tn_snps. Per DECISIONS.md, "
             "lambda_GC > 1.1 should be flagged for revision response.",
    )
    return p.parse_args()


def compute_lambda_gc(pvals):
    """Genomic inflation factor: median of chi-sq statistics divided by chi-sq median (0.4549)."""
    pvals = np.asarray(pvals, dtype=float)
    pvals = pvals[(pvals > 0) & (pvals <= 1)]
    if len(pvals) == 0:
        return float("nan")
    chisq = -2 * np.log(pvals)  # rough chi-sq under null; works for any df
    # The convention is to use 1-df chi-sq from p-values directly via the inverse survival fn.
    from scipy.stats import chi2
    chisq = chi2.isf(pvals, df=1)
    return float(np.median(chisq) / chi2.ppf(0.5, df=1))


def main():
    args = parse_args()
    if args.col_n is None and args.n_fixed is None:
        sys.exit("Either --col-n or --n-fixed must be provided.")

    print(f"[prepare_gwas] reading {args.input}", flush=True)
    df = pd.read_csv(args.input, sep=args.sep, low_memory=False)
    n0 = len(df)
    print(f"[prepare_gwas] {n0:,} rows on input", flush=True)

    required = [args.col_snp, args.col_chr, args.col_bp, args.col_p]
    missing = [c for c in required if c not in df.columns]
    if missing:
        sys.exit(f"Missing required columns: {missing}\nAvailable: {list(df.columns)}")

    df = df.dropna(subset=required)
    df = df[(df[args.col_p] > 0) & (df[args.col_p] <= 1)]

    if args.col_frq and args.col_frq in df.columns:
        frq = df[args.col_frq].astype(float)
        maf = frq.where(frq <= 0.5, 1 - frq)
        df = df[maf >= args.maf_min]
        print(f"[prepare_gwas] MAF >= {args.maf_min}: {len(df):,} rows", flush=True)

    if args.col_info and args.col_info in df.columns:
        df = df[df[args.col_info].astype(float) >= args.info_min]
        print(f"[prepare_gwas] INFO >= {args.info_min}: {len(df):,} rows", flush=True)

    df[args.col_chr] = df[args.col_chr].astype(str).str.replace("^chr", "", regex=True)

    if not args.keep_non_autosomes:
        n_before = len(df)
        df = df[df[args.col_chr].isin(AUTOSOMES)]
        print(f"[prepare_gwas] autosomes only (chr 1-22): {len(df):,} rows "
              f"(dropped {n_before - len(df):,} non-autosome SNPs)", flush=True)

    snp_loc = df[[args.col_snp, args.col_chr, args.col_bp]].copy()
    snp_loc.columns = ["SNP", "CHR", "BP"]
    snp_loc = snp_loc.drop_duplicates(subset=["SNP"])

    pval = df[[args.col_snp, args.col_p]].copy()
    pval.columns = ["SNP", "P"]
    if args.col_n and args.col_n in df.columns:
        pval["N"] = df[args.col_n].astype(int).values
        if args.n_min_frac > 0:
            n_max = int(pval["N"].max())
            n_threshold = int(n_max * args.n_min_frac)
            n_before = len(pval)
            pval = pval[pval["N"] >= n_threshold]
            n_dropped = n_before - len(pval)
            print(f"[prepare_gwas] N filter: drop SNPs with N < {n_threshold:,} "
                  f"({args.n_min_frac:.2f} x max={n_max:,}); kept {len(pval):,} of "
                  f"{n_before:,} ({n_dropped:,} dropped, {100*n_dropped/n_before:.1f}%)",
                  flush=True)
            # Keep snp_loc in sync with pval after the N filter.
            snp_loc = snp_loc[snp_loc["SNP"].isin(pval["SNP"])]
    else:
        pval["N"] = args.n_fixed
        # --n-fixed path: every SNP shares one N, so the LDSC-style low-N tail
        # filter is a no-op and is skipped silently.
    pval = pval.drop_duplicates(subset=["SNP"])

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    snp_loc_path = Path(str(out_prefix) + ".snp.loc")
    pval_path = Path(str(out_prefix) + ".pval")

    snp_loc.to_csv(snp_loc_path, sep="\t", index=False, header=False)
    pval.to_csv(pval_path, sep="\t", index=False)

    print(f"[prepare_gwas] wrote {snp_loc_path} ({len(snp_loc):,} SNPs)", flush=True)
    print(f"[prepare_gwas] wrote {pval_path} ({len(pval):,} SNPs)", flush=True)
    print(f"[prepare_gwas] dropped {n0 - len(snp_loc):,} rows during QC", flush=True)

    if args.lambda_gc_out:
        lam = compute_lambda_gc(pval["P"].values)
        out = Path(args.lambda_gc_out)
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w") as fh:
            fh.write("gwas\tlambda_gc\tn_snps\n")
            fh.write(f"{out_prefix.name}\t{lam:.4f}\t{len(pval)}\n")
        flag = " *** > 1.1 — flag for revision response ***" if lam > 1.1 else ""
        print(f"[prepare_gwas] lambda_GC = {lam:.4f}{flag}", flush=True)
        print(f"[prepare_gwas] wrote {out}", flush=True)


if __name__ == "__main__":
    main()
