"""
Convert MAGMA .genes.out to an scDRS .gs file.

scDRS expects, per row: TRAIT <tab> GENESET, where GENESET is a comma-separated
list of GENE:WEIGHT entries. The convention from the scDRS paper is to take the
top 1000 genes by MAGMA Z-score and weight them by Z.

Genes are emitted as gene SYMBOLS (uppercase). MAGMA's .genes.out uses Entrez
IDs by default; pass --gene-loc (the same NCBI37.3.gene.loc used for annotation)
to map back to symbols.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genes-out", required=True, help="MAGMA .genes.out file")
    p.add_argument("--gene-loc", required=True, help="MAGMA gene location file (e.g. NCBI37.3.gene.loc) used to map Entrez->symbol")
    p.add_argument("--trait", default="UC", help="Trait label written into the .gs file")
    p.add_argument("--top-n", type=int, default=1000)
    p.add_argument("--out", required=True, help="Output .gs path")
    return p.parse_args()


def load_gene_loc(path):
    # NCBI37.3.gene.loc columns: GENE_ID CHR START END STRAND SYMBOL
    df = pd.read_csv(path, sep=r"\s+", header=None,
                     names=["GENE_ID", "CHR", "START", "END", "STRAND", "SYMBOL"])
    df["GENE_ID"] = df["GENE_ID"].astype(str)
    return df.set_index("GENE_ID")["SYMBOL"].str.upper()


def main():
    args = parse_args()

    genes = pd.read_csv(args.genes_out, sep=r"\s+")
    if "ZSTAT" not in genes.columns or "GENE" not in genes.columns:
        sys.exit(f"Expected GENE and ZSTAT columns in {args.genes_out}; got {list(genes.columns)}")
    genes["GENE"] = genes["GENE"].astype(str)

    sym_map = load_gene_loc(args.gene_loc)
    genes["SYMBOL"] = genes["GENE"].map(sym_map)

    n_unmapped = genes["SYMBOL"].isna().sum()
    if n_unmapped:
        print(f"[make_scdrs_gs] dropping {n_unmapped:,} genes with no symbol mapping", flush=True)
        genes = genes.dropna(subset=["SYMBOL"])

    genes = genes.dropna(subset=["ZSTAT"])
    genes = genes.sort_values("ZSTAT", ascending=False).drop_duplicates("SYMBOL")
    top = genes.head(args.top_n).copy()
    top = top[top["ZSTAT"] > 0]

    if len(top) < args.top_n:
        print(f"[make_scdrs_gs] WARNING: only {len(top)} positive-Z genes available (asked for top {args.top_n})", flush=True)

    geneset = ",".join(f"{sym}:{z:.4f}" for sym, z in zip(top["SYMBOL"], top["ZSTAT"]))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("TRAIT\tGENESET\n")
        fh.write(f"{args.trait}\t{geneset}\n")

    print(f"[make_scdrs_gs] wrote {out} ({len(top)} genes, trait={args.trait})", flush=True)


if __name__ == "__main__":
    main()
