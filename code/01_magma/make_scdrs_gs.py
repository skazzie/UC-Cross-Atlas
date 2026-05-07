"""
Convert MAGMA .genes.out to an scDRS .gs file (and a long-format gene-Z
table for seismicGWAS), with MHC-region exclusion locked in.

scDRS expects, per row: TRAIT <tab> GENESET, where GENESET is a
comma-separated list of GENE:WEIGHT entries. The convention from the scDRS
paper is to take the top 1000 genes by MAGMA Z-score and weight them by Z.

Genes are emitted as gene SYMBOLS (uppercase). MAGMA's .genes.out uses
Entrez IDs by default; pass --gene-loc (the same NCBI37.3.gene.loc used
for annotation) to map back to symbols.

MHC region exclusion (DECISIONS.md):
  Per the locked v1 plan, MHC genes (chr 6: 28,477,797-33,448,354 in
  GRCh37) are excluded from the scDRS top-1000 gene set and from the
  seismicGWAS gene-Z-score table. UC has the strongest GWAS signals in
  the MHC region; without exclusion, scDRS scores can be dominated by
  MHC-driven cells (mostly antigen-presenting cells), confounding
  cross-atlas comparison. Use --keep-mhc only for the supplementary
  MHC-included sensitivity analysis on Smillie x de Lange.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

# IMGT/HLA region in GRCh37 coordinates
MHC_CHR = "6"
MHC_START = 28_477_797
MHC_END = 33_448_354


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genes-out", required=True, help="MAGMA .genes.out file")
    p.add_argument(
        "--gene-loc",
        required=True,
        help="MAGMA gene location file (e.g. NCBI37.3.gene.loc) used to map "
             "Entrez->symbol AND to identify MHC-region genes by coordinates.",
    )
    p.add_argument("--trait", default="UC", help="Trait label written into the .gs file")
    p.add_argument("--top-n", type=int, default=1000)
    p.add_argument("--out", required=True, help="Output .gs path")
    p.add_argument(
        "--out-zscore-table",
        default=None,
        help="Optional output path for a long-format gene-Z-score table "
             "(SYMBOL, ENTREZ, CHR, START, END, ZSTAT, P) used by "
             "seismicGWAS. MHC genes are also excluded from this table "
             "unless --keep-mhc is set.",
    )
    p.add_argument(
        "--keep-mhc",
        action="store_true",
        help="Disable MHC-region gene exclusion. Locked v1 default is to "
             "exclude MHC; only use this for the MHC-included sensitivity "
             "run on Smillie x de Lange.",
    )
    return p.parse_args()


def load_gene_loc(path):
    """NCBI37.3.gene.loc columns: GENE_ID CHR START END STRAND SYMBOL"""
    df = pd.read_csv(
        path, sep=r"\s+", header=None,
        names=["GENE_ID", "CHR", "START", "END", "STRAND", "SYMBOL"],
    )
    df["GENE_ID"] = df["GENE_ID"].astype(str)
    df["CHR"] = df["CHR"].astype(str)
    df["SYMBOL"] = df["SYMBOL"].str.upper()
    return df


def is_mhc_gene(chr_, start, end):
    """A gene overlaps the MHC region if its body intersects MHC_START..MHC_END on chr 6."""
    if str(chr_) != MHC_CHR:
        return False
    return not (end < MHC_START or start > MHC_END)


def main():
    args = parse_args()

    genes = pd.read_csv(args.genes_out, sep=r"\s+")
    if "ZSTAT" not in genes.columns or "GENE" not in genes.columns:
        sys.exit(f"Expected GENE and ZSTAT columns in {args.genes_out}; got {list(genes.columns)}")
    genes["GENE"] = genes["GENE"].astype(str)

    gene_loc = load_gene_loc(args.gene_loc)
    sym_map = gene_loc.set_index("GENE_ID")["SYMBOL"]
    chr_map = gene_loc.set_index("GENE_ID")["CHR"]
    start_map = gene_loc.set_index("GENE_ID")["START"]
    end_map = gene_loc.set_index("GENE_ID")["END"]

    genes["SYMBOL"] = genes["GENE"].map(sym_map)
    genes["CHR"] = genes["GENE"].map(chr_map)
    genes["START"] = genes["GENE"].map(start_map)
    genes["END"] = genes["GENE"].map(end_map)

    n_unmapped = genes["SYMBOL"].isna().sum()
    if n_unmapped:
        print(f"[make_scdrs_gs] dropping {n_unmapped:,} genes with no symbol mapping", flush=True)
        genes = genes.dropna(subset=["SYMBOL", "CHR", "START", "END"])

    genes = genes.dropna(subset=["ZSTAT"])

    if not args.keep_mhc:
        is_mhc = [
            is_mhc_gene(c, s, e)
            for c, s, e in zip(genes["CHR"], genes["START"], genes["END"])
        ]
        n_mhc = sum(is_mhc)
        genes = genes[[not flag for flag in is_mhc]]
        print(
            f"[make_scdrs_gs] MHC exclusion (chr {MHC_CHR}: "
            f"{MHC_START:,}-{MHC_END:,}): dropped {n_mhc:,} MHC genes",
            flush=True,
        )
    else:
        print(
            "[make_scdrs_gs] WARNING: --keep-mhc set; MHC genes retained. "
            "Locked v1 default is to exclude MHC.",
            flush=True,
        )

    genes = genes.sort_values("ZSTAT", ascending=False).drop_duplicates("SYMBOL")
    top = genes.head(args.top_n).copy()
    top = top[top["ZSTAT"] > 0]

    if len(top) < args.top_n:
        print(
            f"[make_scdrs_gs] WARNING: only {len(top)} positive-Z genes available "
            f"(asked for top {args.top_n})",
            flush=True,
        )

    geneset = ",".join(f"{sym}:{z:.4f}" for sym, z in zip(top["SYMBOL"], top["ZSTAT"]))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("TRAIT\tGENESET\n")
        fh.write(f"{args.trait}\t{geneset}\n")
    print(f"[make_scdrs_gs] wrote {out} ({len(top)} genes, trait={args.trait})", flush=True)

    if args.out_zscore_table:
        ztable = genes[["SYMBOL", "GENE", "CHR", "START", "END", "ZSTAT", "P"]].rename(
            columns={"GENE": "ENTREZ"}
        )
        ztable_path = Path(args.out_zscore_table)
        ztable_path.parent.mkdir(parents=True, exist_ok=True)
        ztable.to_csv(ztable_path, sep="\t", index=False)
        print(
            f"[make_scdrs_gs] wrote {ztable_path} ({len(ztable):,} genes, "
            f"for seismicGWAS)",
            flush=True,
        )


if __name__ == "__main__":
    main()
