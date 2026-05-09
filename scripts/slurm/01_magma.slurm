#!/usr/bin/env bash
#SBATCH --job-name=ucc_magma
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=128x24
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=PUT_YOUR_EMAIL_HERE
# ^^ Edit --partition / --mail-user, OR override at submit time:
#    sbatch -p $UCC_PARTITION --mail-user=$UCC_EMAIL ...
#
# Runs the full MAGMA pipeline for one GWAS:
#   prepare_gwas -> annotate -> gene-test -> sanity_check -> .gs + gene-Z table
#   (+ MHC-included sensitivity .gs, only for de Lange)
#
# Submit one of:
#   sbatch --export=ALL,GWAS=delange scripts/slurm/01_magma.slurm
#   sbatch --export=ALL,GWAS=liu     scripts/slurm/01_magma.slurm
#   sbatch --export=ALL,GWAS=scz     scripts/slurm/01_magma.slurm

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "$SCRIPT_DIR/../config.sh"

# Activate conda env
if command -v module >/dev/null 2>&1; then
    module load miniconda3 || module load anaconda3 || true
fi
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$UCC_CONDA_ENV"

# Make MAGMA binary discoverable
export PATH="$UCC_DATA/reference:$PATH"

cd "$UCC_REPO"

# ---- Per-GWAS parameters -------------------------------------------------

case "${GWAS:-}" in
  delange)
    INPUT="$UCC_DATA/gwas/uc_delange_harmonised.tsv.gz"
    OUT_PREFIX="$UCC_DATA/gwas/uc_delange"
    N_FIXED=45975                # 12,366 cases + 33,609 controls
    TRAIT_CLASS=uc
    TRAIT=UC
    ;;
  liu)
    INPUT="$UCC_DATA/gwas/uc_liu.tsv.gz"
    OUT_PREFIX="$UCC_DATA/gwas/uc_liu"
    # FILL IN combined N from Liu 2023 paper at M1, OR pass --col-n if
    # the file includes a per-SNP N column (preferred — DECISIONS.md).
    N_FIXED=999999
    TRAIT_CLASS=uc
    TRAIT=UC
    ;;
  scz)
    INPUT="$UCC_DATA/gwas/scz_trubetskoy.tsv.gz"
    OUT_PREFIX="$UCC_DATA/gwas/scz_trubetskoy"
    # 76,755 cases + 53,880 controls; verify exact N at M1
    N_FIXED=130635
    TRAIT_CLASS=schizophrenia
    TRAIT=SCZ
    ;;
  *)
    echo "ERROR: set GWAS=delange | liu | scz" >&2
    exit 1
    ;;
esac

echo "[01_magma] GWAS=$GWAS  input=$INPUT"

# ---- 1. Reformat (autosomes, MAF/INFO QC, lambda_GC) ---------------------

python code/01_magma/prepare_gwas.py \
    --input "$INPUT" \
    --out-prefix "$OUT_PREFIX" \
    --col-snp hm_rsid --col-chr hm_chrom --col-bp hm_pos \
    --col-p p_value --n-fixed "$N_FIXED" \
    --lambda-gc-out "$UCC_RESULTS/magma/${GWAS}_lambda_gc.tsv"

# ---- 2. MAGMA annotate + gene-test ---------------------------------------

bash code/01_magma/run_magma.sh \
    "${OUT_PREFIX}.snp.loc" "${OUT_PREFIX}.pval" \
    "$UCC_DATA/reference/NCBI37.3.gene.loc" \
    "$UCC_DATA/reference/g1000_eur" \
    10 "$UCC_RESULTS/magma/${GWAS}_10kb"

# ---- 3. Sanity check (top genes match expected pattern) ------------------

python code/01_magma/sanity_check.py \
    --genes-out "$UCC_RESULTS/magma/${GWAS}_10kb.genes.out" \
    --gene-loc "$UCC_DATA/reference/NCBI37.3.gene.loc" \
    --trait-class "$TRAIT_CLASS"

# ---- 4. .gs file + gene-Z table (MHC excluded) ---------------------------

python code/01_magma/make_scdrs_gs.py \
    --genes-out "$UCC_RESULTS/magma/${GWAS}_10kb.genes.out" \
    --gene-loc "$UCC_DATA/reference/NCBI37.3.gene.loc" \
    --trait "$TRAIT" --top-n 1000 \
    --out "$UCC_DATA/gwas/${GWAS}_top1000.gs" \
    --out-zscore-table "$UCC_DATA/gwas/${GWAS}_gene_z.tsv"

# ---- 5. MHC-included sensitivity gene set (de Lange only) ----------------

if [ "$GWAS" = "delange" ]; then
    python code/01_magma/make_scdrs_gs.py \
        --genes-out "$UCC_RESULTS/magma/${GWAS}_10kb.genes.out" \
        --gene-loc "$UCC_DATA/reference/NCBI37.3.gene.loc" \
        --trait "$TRAIT" --top-n 1000 \
        --out "$UCC_DATA/gwas/${GWAS}_top1000_with_mhc.gs" \
        --keep-mhc
fi

echo "[01_magma] DONE for GWAS=$GWAS"
