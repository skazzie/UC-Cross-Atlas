#!/usr/bin/env bash
# Run MAGMA annotation + gene-based test.
#
# Locked v1 (DECISIONS.md):
#   - 10 kb upstream + 10 kb downstream window — passed as the 5th
#     positional argument. MAGMA's --annotate default is "0,0"
#     (intragenic only), which produces a gene-length-biased ranking
#     (large genes capture more signal than short ones). scDRS's
#     published settings use 10,10. The window arg here is MANDATORY;
#     the script refuses to run with an empty value.
#   - 1000G EUR LD reference for both UC GWAS (acknowledged approximate
#     for multi-ancestry Liu 2023; documented in Methods)
#   - autosomes only — enforced upstream by prepare_gwas.py (which drops
#     non-1-22 SNPs by default)
#   - MHC region exclusion is applied LATER by make_scdrs_gs.py to the
#     downstream gene set; MAGMA itself sees MHC SNPs to maintain the
#     LD null model
#
# Run this once per GWAS:
#   - de Lange 2017 UC (primary)
#   - Liu 2023 multi-ancestry UC arm (cross-GWAS sensitivity)
#   - Trubetskoy 2022 schizophrenia (negative control on Smillie)
#
# Usage:
#   ./run_magma.sh <snp_loc> <pval_file> <gene_loc> <bfile_prefix> <window_kb> <out_prefix>
#
# Example (de Lange UC):
#   ./run_magma.sh data/gwas/uc_delange.snp.loc data/gwas/uc_delange.pval \
#                  data/reference/NCBI37.3.gene.loc \
#                  data/reference/g1000_eur \
#                  10 results/magma/uc_delange_10kb

set -euo pipefail

if [[ $# -ne 6 ]]; then
    echo "Usage: $0 <snp_loc> <pval_file> <gene_loc> <bfile_prefix> <window_kb> <out_prefix>" >&2
    exit 1
fi

SNP_LOC="$1"
PVAL_FILE="$2"
GENE_LOC="$3"
BFILE="$4"
WINDOW_KB="$5"
OUT_PREFIX="$6"

if [[ -z "$WINDOW_KB" || ! "$WINDOW_KB" =~ ^[0-9]+$ ]]; then
    echo "ERROR: window_kb must be a positive integer (got: '${WINDOW_KB}'). " \
         "MAGMA default 0,0 produces gene-length-biased rankings." >&2
    exit 1
fi

mkdir -p "$(dirname "$OUT_PREFIX")"

# NOTE: window=${WINDOW_KB},${WINDOW_KB} is passed EXPLICITLY. MAGMA default
# is 0,0 (intragenic SNPs only), which is NOT what scDRS expects. See
# header comment and DECISIONS.md.
echo "[magma] step 1: annotate (window=${WINDOW_KB},${WINDOW_KB} — EXPLICIT, not default)"
magma --annotate "window=${WINDOW_KB},${WINDOW_KB}" \
      --snp-loc "$SNP_LOC" \
      --gene-loc "$GENE_LOC" \
      --out "${OUT_PREFIX}_annot"

echo "[magma] step 2: gene-based test"
magma --bfile "$BFILE" \
      --pval "$PVAL_FILE" use=SNP,P ncol=N \
      --gene-annot "${OUT_PREFIX}_annot.genes.annot" \
      --out "${OUT_PREFIX}"

echo "[magma] done. Output: ${OUT_PREFIX}.genes.out"
