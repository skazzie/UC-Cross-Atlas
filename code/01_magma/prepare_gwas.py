"""
Convert a GWAS summary statistics file into the two MAGMA inputs:
  - <prefix>.snp.loc : rsid, chr, bp  (whitespace-separated, no header)
  - <prefix>.pval    : SNP, P, N      (whitespace-separated, with header)

Locked v1 GWAS in this pipeline:
  - de Lange 2017 UC (primary)              GWAS Catalog GCST004133
  - Liu 2023 multi-ancestry IBD (UC arm)    cross-GWAS sensitivity
  - Trubetskoy 2022 schizophrenia            negative control on Smillie

Pre-committed filters applied (DECISIONS.md):
  - autosomes only (chr 1-22), X chromosome excluded
  - drop rows with missing rsid/chr/bp/p
  - drop p outside (0, 1]
  - drop MAF < 0.01 (if FRQ column is provided)
  - drop INFO < 0.6 (if INFO column is provided)

Build-safe mode (DECISIONS correction 2026-10-07):
  Pass --bim <plink.bim> to re-anchor every SNP to the LD-panel's
  rsID+chr+bp. This drops SNPs not in the panel and overwrites any
  build-mismatched chr/bp from the sumstats (de Lange's harmonised and
  Liu's raw ship GRCh38 positions, which do not match the GRCh37
  g1000_eur panel — SNPs get mis-assigned to neighbouring genes and
  canonical UC loci like IL23R / PTPN22 end up with noise Z-scores).
  --bim-match=rsid        join sumstats[col-snp] to bim[rsID]
  --bim-match=chrpos      parse chr:pos from sumstats[col-snp] (handles
                          'chr:pos', 'chr:pos_A1_A2', 'chr_pos_A1_A2'),
                          join sumstats[(chr,bp)] to bim[(chr,bp)]
  In both cases, the output rsID/chr/bp are bim's. For raw de Lange the
  SNP id is MarkerName='chr:pos_A1_A2' build37-native (bim-match=chrpos
  is a true-position join). For Liu use harmonised rsid + bim-match=rsid
  (Liu's raw variant_id is build38 chr:pos and does not safely join on
  position).

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
        "--bim",
        default=None,
        help="PLINK .bim path. If set, re-anchor every SNP to the panel's "
             "rsID+chr+bp. Overrides --col-chr/--col-bp from sumstats.",
    )
    p.add_argument(
        "--bim-match",
        default="rsid",
        choices=["rsid", "chrpos"],
        help="Join key when --bim is set. 'rsid' = sumstats --col-snp -> "
             "bim rsID. 'chrpos' = parse chr:pos from sumstats --col-snp "
             "(handles ':' and '_' separators) and join on (chr,bp).",
    )
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
        "--skip-comments",
        action="store_true",
        help="PGC sumstats VCF v1.0 convention: skip leading lines starting "
             "with '##' (metadata), then read the column-header line. If the "
             "header itself starts with a single '#', that prefix is stripped "
             "from the first column name. Use for Trubetskoy SCZ figshare "
             "deposit and any other PGC-style file. See DECISIONS correction 28.",
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


def _drop_ambiguous_chrpos(df, chr_col, bp_col, side_label):
    """Drop every row whose (chr,bp) appears more than once in df.
    This is the cheap guard rather than resolving multi-allelics: at an
    ambiguous site we don't know which record MAGMA should use, so drop
    all records at that site. Returns (df_filtered, n_sites_dropped).
    """
    key = pd.MultiIndex.from_arrays([df[chr_col].values, df[bp_col].values])
    counts = key.value_counts()
    amb_sites = counts[counts > 1].index
    if len(amb_sites) == 0:
        return df, 0
    mask = key.isin(amb_sites)
    n_dropped_rows = int(mask.sum())
    df_out = df[~mask].reset_index(drop=True)
    print(
        f"[prepare_gwas] dedup {side_label}: dropped {len(amb_sites):,} "
        f"ambiguous chr:pos sites ({n_dropped_rows:,} rows)",
        flush=True,
    )
    return df_out, int(len(amb_sites))


def _bim_reanchor(df, args, n_raw):
    """Join sumstats to bim and overwrite args.col_snp/col_chr/col_bp
    with bim's values. Returns the inner-joined df. Mutates args.

    Ambiguous (multi-allelic) chr:pos sites are dropped on both sides
    before the merge, per 2026-10-07 directive — simpler than resolving
    on alleles (MAGMA gene-analysis consumes SNP/P/N only; alleles can't
    change a gene Z) and avoids using de Lange's non-strand-resolved
    METAL alleles. Also reports (i) overall merge rate vs raw input,
    (ii) ambiguous-site counts; hard-fails if merge rate < 90%.
    """
    bim = pd.read_csv(
        args.bim,
        sep=r"\s+",
        header=None,
        names=["CHR_BIM", "RSID_BIM", "CM_BIM", "BP_BIM", "A1_BIM", "A2_BIM"],
        dtype={"CHR_BIM": str},
    )
    bim["CHR_BIM"] = bim["CHR_BIM"].str.replace("^chr", "", regex=True)
    print(f"[prepare_gwas] bim: {len(bim):,} variants from {args.bim}", flush=True)

    # Dedup bim upfront — multi-allelic positions appear as multiple rows
    # in the bim and would create multi-match ambiguity on any join.
    bim, n_bim_amb = _drop_ambiguous_chrpos(bim, "CHR_BIM", "BP_BIM", "bim")

    if args.bim_match == "chrpos":
        # Parse chr:pos from sumstats --col-snp. Handles 'chr:pos',
        # 'chr:pos_A1_A2' (de Lange MarkerName) and 'chr_pos_A1_A2'
        # (Liu variant_id). The last format uses '_' for all separators;
        # the first two use ':' between chr and pos. Normalize both.
        snp = df[args.col_snp].astype(str).str.replace("^chr", "", regex=True)
        # If a ':' is present, split on it first; else fall back to '_'.
        has_colon = snp.str.contains(":")
        chr_tok = np.where(has_colon, snp.str.split(":").str[0],
                           snp.str.split("_").str[0])
        pos_rest = np.where(has_colon, snp.str.split(":").str[1],
                            snp.str.split("_", n=1).str[1])
        pos_tok = pd.Series(pos_rest).astype(str).str.split("_").str[0]
        df = df.assign(_SS_CHR=pd.Series(chr_tok).values,
                       _SS_BP=pd.to_numeric(pos_tok, errors="coerce"))
        df = df.dropna(subset=["_SS_BP"])
        df["_SS_BP"] = df["_SS_BP"].astype(int)
        print(f"[prepare_gwas] bim-match=chrpos: parsed {len(df):,} chr:pos "
              f"from {args.col_snp!r} (sample: {df.iloc[0][args.col_snp]!r})",
              flush=True)
        # Dedup sumstats on parsed chr:pos
        df, n_ss_amb = _drop_ambiguous_chrpos(df, "_SS_CHR", "_SS_BP", "sumstats")
        df = df.merge(
            bim[["RSID_BIM", "CHR_BIM", "BP_BIM"]],
            left_on=["_SS_CHR", "_SS_BP"],
            right_on=["CHR_BIM", "BP_BIM"],
            how="inner",
        )
    else:
        # rsID join. Dedup sumstats on (col-chr, col-bp) if those columns
        # exist (Liu harmonised ships chromosome + base_pair_location);
        # otherwise just the implicit uniqueness of the rsID key.
        print(f"[prepare_gwas] bim-match=rsid: joining sumstats {args.col_snp!r} "
              f"to bim rsID", flush=True)
        if args.col_chr in df.columns and args.col_bp in df.columns:
            df, n_ss_amb = _drop_ambiguous_chrpos(
                df, args.col_chr, args.col_bp, "sumstats")
        else:
            n_ss_amb = 0
            print(f"[prepare_gwas] dedup sumstats: skipped (no chr/bp cols)",
                  flush=True)
        df = df.merge(
            bim[["RSID_BIM", "CHR_BIM", "BP_BIM"]],
            left_on=args.col_snp,
            right_on="RSID_BIM",
            how="inner",
        )

    # Overwrite SNP/CHR/BP with bim's values. Downstream code reads from
    # args.col_snp/args.col_chr/args.col_bp — update those pointers.
    df[args.col_snp] = df["RSID_BIM"].values
    _ss_chr_col = "_SS_FINAL_CHR"
    _ss_bp_col = "_SS_FINAL_BP"
    df[_ss_chr_col] = df["CHR_BIM"].values
    df[_ss_bp_col] = df["BP_BIM"].values
    args.col_chr = _ss_chr_col
    args.col_bp = _ss_bp_col
    n_after = len(df)
    rate = 100.0 * n_after / max(n_raw, 1)
    print(f"[prepare_gwas] after bim-join: {n_after:,} rows ({rate:.1f}% of "
          f"{n_raw:,} raw sumstats rows); ambiguous sites dropped — "
          f"bim: {n_bim_amb:,}, sumstats: {n_ss_amb:,}", flush=True)
    if rate < 90.0:
        sys.exit(
            f"[prepare_gwas] FAIL: merge rate {rate:.1f}% below 90% gate — "
            f"something is wrong with the sumstats <-> bim correspondence "
            f"(wrong build, wrong rsID column, truncated bim, ...). "
            f"Investigate before proceeding."
        )
    return df


def main():
    args = parse_args()
    if args.col_n is None and args.n_fixed is None:
        sys.exit("Either --col-n or --n-fixed must be provided.")

    print(f"[prepare_gwas] reading {args.input}", flush=True)
    # PGC sumstats VCF v1.0 convention: ## metadata lines, then a single-'#'
    # column header line. pd.read_csv's `comment='#'` would eat the header too,
    # so we count the ## prefix and pass skiprows instead, then strip a leading
    # '#' from the first column name if it survived. See DECISIONS 28.
    skiprows = 0
    if args.skip_comments:
        import gzip as _gzip
        opener = _gzip.open if str(args.input).endswith(".gz") else open
        # Explicit utf-8 + errors='replace': PGC sumstats VCF metadata lines
        # contain accented author names (e.g. acknowledgments listing), and
        # Windows default cp1252 chokes on byte 0x9d. errors='replace' is
        # safe because we only USE the ## prefix detection on these lines,
        # not their content.
        with opener(args.input, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("##"):
                    skiprows += 1
                else:
                    break
        print(f"[prepare_gwas] PGC mode: skipping {skiprows} leading '##' lines",
              flush=True)
    df = pd.read_csv(args.input, sep=args.sep, low_memory=False, skiprows=skiprows)
    if args.skip_comments and len(df.columns) > 0 and str(df.columns[0]).startswith("#"):
        new_first = str(df.columns[0]).lstrip("#")
        print(f"[prepare_gwas] PGC mode: stripped '#' from header column 0: "
              f"{df.columns[0]!r} -> {new_first!r}", flush=True)
        df.columns = [new_first] + list(df.columns[1:])
    n0 = len(df)
    print(f"[prepare_gwas] {n0:,} rows on input", flush=True)

    # Column requirement check — when --bim is set, chr/bp come from bim
    # so sumstats need only --col-snp + --col-p.
    if args.bim:
        required = [args.col_snp, args.col_p]
    else:
        required = [args.col_snp, args.col_chr, args.col_bp, args.col_p]
    missing = [c for c in required if c not in df.columns]
    if missing:
        sys.exit(f"Missing required columns: {missing}\nAvailable: {list(df.columns)}")

    df = df.dropna(subset=required)
    df = df[(df[args.col_p] > 0) & (df[args.col_p] <= 1)]

    # Build-safe re-anchor to the LD panel's rsID+chr+bp. Must happen
    # before autosome filter + snp.loc output so args.col_chr/col_bp
    # point at the bim-derived columns.
    if args.bim:
        df = _bim_reanchor(df, args, n_raw=n0)

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
