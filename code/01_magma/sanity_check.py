"""
Sanity-check a MAGMA .genes.out file against expected top-gene patterns
for the trait under analysis.

Locked v1 expected patterns (PLAN.md §2.1, DECISIONS.md):

- de Lange 2017 UC / Liu 2023 UC — hard gate 2026-10-07 (replaces the
  earlier 15-gene / top-20 / min-2 panel, whose expected values and
  PLCL1-Z=9.81 target were derived from a build-mismatch artifact that
  mis-assigned SNPs to large neighbouring genes):
    panel = {IL23R, IL10, TNFSF15, NKX2-3, STAT3, JAK2, FCGR2A, SMAD3}
    at least 6 of 8 must appear in the top 150 non-MHC genes.
  MHC region (chr6:28-34Mb) is excluded from the ranking because UC
  MAGMA output is dominated by MHC and the canonical non-MHC panel
  would otherwise be crowded out of the top of the list.
  PLCL1 and PTPN22 are intentionally NOT in the panel:
    PLCL1 is not an established UC locus; its old rank-1 was a
      large-gene build-mismatch artifact.
    PTPN22 is a Crohn's-specific locus (rs2476601 UC OR 0.98 p=0.88,
      CD OR 0.81 p=7.4e-6); bottom-of-ranking in a UC-only MAGMA
      is correct behaviour.

- Trubetskoy 2022 schizophrenia (negative control): top genes should be
  brain-related (e.g. CACNA1C, GRIN2A, DRD2, FURIN, SP4). Legacy top-20
  / min-2 gate retained — SCZ's MHC burden is cleaner.

Exit code 0 = looks OK, 1 = something looks wrong.
"""

import argparse
import sys

import pandas as pd

# UC canonical panel (2026-10-07 revision). Match on uppercase SYMBOL.
UC_PANEL = {
    "IL23R", "IL10", "TNFSF15", "NKX2-3",
    "STAT3", "JAK2", "FCGR2A", "SMAD3",
}

SCZ_PANEL = {
    "CACNA1C", "GRIN2A", "DRD2", "FURIN", "SP4",
    "GRIA3", "TCF4", "ZNF804A", "CACNB2", "GABBR2",
    "RBFOX1", "NRGN", "SETD1A",
}

# MHC region (GRCh37). Any gene whose window overlaps this is excluded
# from the UC top-N ranking. Standard scDRS/MAGMA convention.
MHC_CHR = "6"
MHC_START = 28_000_000
MHC_END = 34_000_000


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genes-out", required=True)
    p.add_argument("--gene-loc", required=True, help="MAGMA gene location file with Entrez->symbol mapping")
    p.add_argument(
        "--trait-class",
        choices=["uc", "schizophrenia"],
        default="uc",
    )
    p.add_argument(
        "--top-n", type=int, default=None,
        help="Top-N by ZSTAT to search (default: 150 for UC, 20 for SCZ)",
    )
    p.add_argument(
        "--min-hits", type=int, default=None,
        help="Minimum panel genes required in top-N (default: 6 for UC, 2 for SCZ)",
    )
    p.add_argument(
        "--mhc-exclude", type=lambda s: s.lower() in ("true", "1", "yes"),
        default=None,
        help="Exclude MHC region (chr6:28-34Mb) from the ranking before "
             "taking top-N (default: True for UC, False for SCZ)",
    )
    return p.parse_args()


def main():
    args = parse_args()
    if args.trait_class == "uc":
        panel = UC_PANEL
        top_n = args.top_n if args.top_n is not None else 150
        min_hits = args.min_hits if args.min_hits is not None else 6
        mhc_exclude = args.mhc_exclude if args.mhc_exclude is not None else True
    else:
        panel = SCZ_PANEL
        top_n = args.top_n if args.top_n is not None else 20
        min_hits = args.min_hits if args.min_hits is not None else 2
        mhc_exclude = args.mhc_exclude if args.mhc_exclude is not None else False

    loc = pd.read_csv(
        args.gene_loc, sep=r"\s+", header=None,
        names=["GENE_ID", "GLCHR", "GSTART", "GEND", "STRAND", "SYMBOL"],
    )
    loc["GENE_ID"] = loc["GENE_ID"].astype(str)
    loc["GLCHR"] = loc["GLCHR"].astype(str)
    sym_map = loc.set_index("GENE_ID")["SYMBOL"].str.upper()
    chr_map = loc.set_index("GENE_ID")["GLCHR"]
    start_map = loc.set_index("GENE_ID")["GSTART"]
    end_map = loc.set_index("GENE_ID")["GEND"]

    genes = pd.read_csv(args.genes_out, sep=r"\s+")
    genes["GENE"] = genes["GENE"].astype(str)
    genes["SYMBOL"] = genes["GENE"].map(sym_map)
    genes["GLCHR"] = genes["GENE"].map(chr_map)
    genes["GSTART"] = genes["GENE"].map(start_map)
    genes["GEND"] = genes["GENE"].map(end_map)

    if mhc_exclude:
        is_mhc = (
            (genes["GLCHR"] == MHC_CHR)
            & (genes["GSTART"] < MHC_END)
            & (genes["GEND"] > MHC_START)
        )
        n_mhc = int(is_mhc.sum())
        print(f"[sanity_check] excluding {n_mhc} MHC genes "
              f"(chr{MHC_CHR}:{MHC_START // 10 ** 6}-{MHC_END // 10 ** 6}Mb)")
        pool = genes[~is_mhc]
    else:
        pool = genes

    top = pool.sort_values("ZSTAT", ascending=False).head(top_n)

    print(
        f"Top {top_n} genes by MAGMA Z-score "
        f"(trait_class={args.trait_class}, mhc_exclude={mhc_exclude}):"
    )
    print(top[["SYMBOL", "GENE", "NSNPS", "ZSTAT", "P"]].head(20).to_string(index=False))
    if top_n > 20:
        print(f"  ... ({top_n - 20} more rows omitted)")
    print()

    top_symbols = set(top["SYMBOL"].dropna())
    hits = top_symbols & panel
    missing = panel - hits
    print(
        f"Panel genes found in top {top_n} non-MHC: "
        f"{sorted(hits) or 'NONE'}"
    )
    print(f"({len(hits)} of {len(panel)} panel; threshold = {min_hits})")
    if missing:
        print(f"Panel genes NOT in top {top_n}: {sorted(missing)}")

    if len(hits) < min_hits:
        print("\nFAIL: too few panel genes in top hits.", file=sys.stderr)
        if args.trait_class == "uc":
            print(
                "Likely causes: genome-build mismatch, wrong --col-* mapping in "
                "prepare_gwas.py, LD reference ancestry mismatch, or wrong arm "
                "of Liu 2023 (must be UC arm only, not CD or IBD-combined). "
                "Verify prepare_gwas.py merge-rate was >=90% and ambiguous-site "
                "dedup counts are not pathological.",
                file=sys.stderr,
            )
        else:
            print(
                "Likely causes: genome-build mismatch, wrong --col-* mapping, "
                "or schizophrenia GWAS download is corrupted/wrong file.",
                file=sys.stderr,
            )
        sys.exit(1)

    print("\nPASS: panel hits in top meet threshold.")


if __name__ == "__main__":
    main()
