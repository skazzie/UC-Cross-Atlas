"""
Sanity-check a MAGMA .genes.out file against expected top-gene patterns
for the trait under analysis.

Locked v1 expected patterns (PLAN.md §2.1, DECISIONS.md):

- de Lange 2017 UC: top 20 should include several of IL23R, JAK2, TYK2,
  STAT3, SMAD3, IL12B, NKX2-3, IRGM, CARD9. NOD2 is weaker than for CD.
- Liu 2023 multi-ancestry UC: broadly similar top genes; if top genes
  are unrelated to immunity, the wrong ancestry arm or CD/IBD-combined
  arm has been grabbed by mistake.
- Trubetskoy 2022 schizophrenia (negative control): top genes should be
  brain-related (e.g. CACNA1C, GRIN2A, DRD2, FURIN, SP4). MHC will not
  appear because MHC genes are excluded from the .genes.out by upstream
  autosome filter + the downstream make_scdrs_gs.py exclusion. (This
  script reads .genes.out which still contains MHC; MHC presence in this
  file is fine — it gets filtered later.)

Exit code 0 = looks OK, 1 = something looks wrong.
"""

import argparse
import sys

import pandas as pd

EXPECTED = {
    "uc": {
        "IL23R", "JAK2", "TYK2", "STAT3", "SMAD3",
        "IL12B", "NKX2-3", "ATG16L1", "IRGM", "CARD9",
        "FCGR2A", "PTPN22", "IL10", "TNFSF15", "IL2RA",
    },
    "schizophrenia": {
        "CACNA1C", "GRIN2A", "DRD2", "FURIN", "SP4",
        "GRIA3", "TCF4", "ZNF804A", "CACNB2", "GABBR2",
        "RBFOX1", "NRGN", "SETD1A",
    },
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genes-out", required=True)
    p.add_argument("--gene-loc", required=True, help="MAGMA gene location file with Entrez->symbol mapping")
    p.add_argument(
        "--trait-class",
        choices=sorted(EXPECTED),
        default="uc",
        help="Which expected-gene panel to compare against. 'uc' for de Lange 2017 + Liu 2023; "
             "'schizophrenia' for Trubetskoy 2022 negative control.",
    )
    p.add_argument("--top-n", type=int, default=20)
    p.add_argument(
        "--min-hits", type=int, default=2,
        help="Minimum number of expected genes that must appear in top-N for the check to pass",
    )
    return p.parse_args()


def main():
    args = parse_args()
    expected = EXPECTED[args.trait_class]

    sym_map = pd.read_csv(
        args.gene_loc, sep=r"\s+", header=None,
        names=["GENE_ID", "CHR", "START", "END", "STRAND", "SYMBOL"],
    )
    sym_map["GENE_ID"] = sym_map["GENE_ID"].astype(str)
    sym_map = sym_map.set_index("GENE_ID")["SYMBOL"].str.upper()

    genes = pd.read_csv(args.genes_out, sep=r"\s+")
    genes["GENE"] = genes["GENE"].astype(str)
    genes["SYMBOL"] = genes["GENE"].map(sym_map)
    top = genes.sort_values("ZSTAT", ascending=False).head(args.top_n)

    print(f"Top {args.top_n} genes by MAGMA Z-score (trait_class={args.trait_class}):")
    print(top[["SYMBOL", "GENE", "NSNPS", "ZSTAT", "P"]].to_string(index=False))
    print()

    top_symbols = set(top["SYMBOL"].dropna())
    hits = top_symbols & expected
    print(
        f"Expected {args.trait_class} genes found in top {args.top_n}: "
        f"{sorted(hits) or 'NONE'}"
    )
    print(f"({len(hits)} of {len(expected)} expected; threshold = {args.min_hits})")

    if len(hits) < args.min_hits:
        print("\nFAIL: too few expected genes in top hits.", file=sys.stderr)
        if args.trait_class == "uc":
            print(
                "Likely causes: genome-build mismatch, wrong --col-* mapping in "
                "prepare_gwas.py, LD reference ancestry mismatch, or wrong arm "
                "of Liu 2023 (must be UC arm only, not CD or IBD-combined).",
                file=sys.stderr,
            )
        else:
            print(
                "Likely causes: genome-build mismatch, wrong --col-* mapping, "
                "or schizophrenia GWAS download is corrupted/wrong file.",
                file=sys.stderr,
            )
        sys.exit(1)

    print("\nPASS: top hits match expected pattern.")


if __name__ == "__main__":
    main()
