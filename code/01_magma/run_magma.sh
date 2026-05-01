#!/usr/bin/env bash
# Run MAGMA annotation + gene-based test for the UC GWAS.
#
# Usage:
#   ./run_magma.sh <snp_loc> <pval_file> <gene_loc> <bfile_prefix> <window_kb> <out_prefix>
#
# Example (10 kb default):
#   ./run_magma.sh data/gwas/uc.snp.loc data/gwas/uc.pval \
#                  data/reference/NCBI37.3.gene.loc \
#                  data/reference/g1000_eur \
#                  10 results/magma/uc_10kb
#
# Example (50 kb sensitivity, week 7.5):
#   ./run_magma.sh data/gwas/uc.snp.loc data/gwas/uc.pval \
#                  data/reference/NCBI37.3.gene.loc \
#                  data/reference/g1000_eur \
#                  50 results/magma/uc_50kb

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

mkdir -p "$(dirname "$OUT_PREFIX")"

echo "[magma] step 1: annotate (window=${WINDOW_KB},${WINDOW_KB})"
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
