"""
Convert a GWAS summary statistics file into the two MAGMA inputs:
  - <prefix>.snp.loc : rsid, chr, bp  (whitespace-separated, no header)
  - <prefix>.pval    : SNP, P, N      (whitespace-separated, with header)

The de Lange 2017 UC summary stats from the GWAS Catalog and IIBDGC use slightly
different column names. Pass the column mapping explicitly with --col-* flags;
defaults match the GWAS Catalog harmonized format.

Filters applied: drop rows with missing rsid/chr/bp/p, MAF < 0.01 (if FRQ column
is provided), INFO < 0.6 (if INFO column is provided), p outside (0, 1].
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


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
    return p.parse_args()


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

    snp_loc = df[[args.col_snp, args.col_chr, args.col_bp]].copy()
    snp_loc.columns = ["SNP", "CHR", "BP"]
    snp_loc = snp_loc.drop_duplicates(subset=["SNP"])

    pval = df[[args.col_snp, args.col_p]].copy()
    pval.columns = ["SNP", "P"]
    if args.col_n and args.col_n in df.columns:
        pval["N"] = df[args.col_n].astype(int).values
    else:
        pval["N"] = args.n_fixed
    pval = pval.drop_duplicates(subset=["SNP"])

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    snp_loc_path = out_prefix.with_suffix(out_prefix.suffix + ".snp.loc") if out_prefix.suffix else Path(str(out_prefix) + ".snp.loc")
    pval_path = out_prefix.with_suffix(out_prefix.suffix + ".pval") if out_prefix.suffix else Path(str(out_prefix) + ".pval")

    snp_loc.to_csv(snp_loc_path, sep="\t", index=False, header=False)
    pval.to_csv(pval_path, sep="\t", index=False)

    print(f"[prepare_gwas] wrote {snp_loc_path} ({len(snp_loc):,} SNPs)", flush=True)
    print(f"[prepare_gwas] wrote {pval_path} ({len(pval):,} SNPs)", flush=True)
    print(f"[prepare_gwas] dropped {n0 - len(snp_loc):,} rows during QC", flush=True)


if __name__ == "__main__":
    main()
