"""
Sanity-check a MAGMA .genes.out file for UC.

Spec section 2.2: the top-20 genes by MAGMA Z-score for de Lange UC should
include known IBD/UC risk genes. If they don't, something upstream is broken
(usually a genome-build mismatch or column-misread).

Exit code 0 = looks OK, 1 = something looks wrong.
"""

import argparse
import sys

import pandas as pd

UC_EXPECTED = {
    "IL23R", "JAK2", "TYK2", "STAT3", "SMAD3",
    "IL12B", "NKX2-3", "ATG16L1", "IRGM", "CARD9",
    "FCGR2A", "PTPN22", "IL10", "TNFSF15", "IL2RA",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genes-out", required=True)
    p.add_argument("--gene-loc", required=True, help="MAGMA gene location file with Entrez->symbol mapping")
    p.add_argument("--top-n", type=int, default=20)
    p.add_argument("--min-hits", type=int, default=2,
                   help="Minimum number of expected UC genes that must appear in top-N for the check to pass")
    return p.parse_args()


def main():
    args = parse_args()

    sym_map = pd.read_csv(args.gene_loc, sep=r"\s+", header=None,
                          names=["GENE_ID", "CHR", "START", "END", "STRAND", "SYMBOL"])
    sym_map["GENE_ID"] = sym_map["GENE_ID"].astype(str)
    sym_map = sym_map.set_index("GENE_ID")["SYMBOL"].str.upper()

    genes = pd.read_csv(args.genes_out, sep=r"\s+")
    genes["GENE"] = genes["GENE"].astype(str)
    genes["SYMBOL"] = genes["GENE"].map(sym_map)
    top = genes.sort_values("ZSTAT", ascending=False).head(args.top_n)

    print(f"Top {args.top_n} genes by MAGMA Z-score:")
    print(top[["SYMBOL", "GENE", "NSNPS", "ZSTAT", "P"]].to_string(index=False))
    print()

    top_symbols = set(top["SYMBOL"].dropna())
    hits = top_symbols & UC_EXPECTED
    print(f"Expected UC risk genes found in top {args.top_n}: {sorted(hits) or 'NONE'}")
    print(f"({len(hits)} of {len(UC_EXPECTED)} expected; threshold = {args.min_hits})")

    if len(hits) < args.min_hits:
        print("\nFAIL: too few known UC risk genes in top hits.", file=sys.stderr)
        print("Likely causes: genome-build mismatch, wrong --col-* mapping in prepare_gwas.py,", file=sys.stderr)
        print("or LD reference ancestry mismatch.", file=sys.stderr)
        sys.exit(1)

    print("\nPASS: top hits look UC-like.")


if __name__ == "__main__":
    main()
