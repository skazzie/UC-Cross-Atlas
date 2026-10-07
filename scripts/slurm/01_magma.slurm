#!/usr/bin/env bash
#SBATCH --job-name=ucc_magma
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=128x24
#SBATCH --account=128x24
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mukhinda@ucsc.edu
#
# Runs the full MAGMA pipeline for one GWAS:
#   prepare_gwas -> annotate -> gene-test -> sanity_check -> .gs + gene-Z table
#   (+ MHC-included sensitivity .gs, only for de Lange)
#
# Submit one of:
#   sbatch --export=ALL,GWAS=delange scripts/slurm/01_magma.slurm
#   sbatch --export=ALL,GWAS=liu     scripts/slurm/01_magma.slurm
#   sbatch --export=ALL,GWAS=scz     scripts/slurm/01_magma.slurm
#
# All UC inputs go through --bim re-anchoring to the g1000_eur GRCh37 panel
# (DECISIONS 2026-10-07) so build-mismatched SNPs can't be assigned to
# the wrong gene. See prepare_gwas.py --bim/--bim-match docstring.

set -euo pipefail

# Under SLURM, BASH_SOURCE points at the script copy in /var/spool/slurmd,
# so prefer SLURM_SUBMIT_DIR (= dir sbatch was invoked from = repo root).
if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    # shellcheck disable=SC1091
    source "$SLURM_SUBMIT_DIR/scripts/config.sh"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    # shellcheck disable=SC1091
    source "$SCRIPT_DIR/../config.sh"
fi

# Activate conda env
if command -v module >/dev/null 2>&1; then
    module load miniconda3/3.13
    module load git/2.52.0 || true
fi
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$UCC_CONDA_ENV"

# Make MAGMA binary discoverable
export PATH="$UCC_DATA/reference:$PATH"

cd "$UCC_REPO"

BIM="$UCC_DATA/reference/g1000_eur.bim"
BIM_MATCH=""          # set per-GWAS below

# ---- Per-GWAS parameters -------------------------------------------------

case "${GWAS:-}" in
  delange)
    # Raw de Lange 2017 deposit: build37-native, no rsID column. The SNP
    # identifier is MarkerName = 'chr:pos_A1_A2'. Parse chr:pos and join
    # to g1000_eur.bim on (chr,bp) — this is a build37-build37 join, so
    # positions align exactly with the LD panel. Harmonised file is
    # GRCh38 and misassigns SNPs under the GRCh37 gene-loc (see 797529
    # post-mortem: IL23R/PTPN22 got noise-Z, PLCL1 Z=9.59 vs expected 9.81).
    INPUT="$UCC_DATA/gwas/uc_delange_GCST004133.txt.gz"
    OUT_PREFIX="$UCC_DATA/gwas/uc_delange"
    N_FIXED=45975                # 12,366 cases + 33,609 controls
    COL_SNP=MarkerName
    COL_P=P.value
    BIM_MATCH=chrpos
    TRAIT_CLASS=uc
    TRAIT=UC
    ;;
  liu)
    # Liu 2023 harmonised deposit: ships a real rsID column
    # (lowercase `rsid`, e.g. rs1274919517). The raw file's variant_id
    # is chr_pos_A1_A2 in GRCh38 and does not safely join by position;
    # harmonised rsID -> bim rsID gives a build-independent anchor.
    INPUT="$UCC_DATA/gwas/uc_liu_harmonised.tsv.gz"
    OUT_PREFIX="$UCC_DATA/gwas/uc_liu"
    N_FIXED=375508
    COL_SNP=rsid
    COL_P=p_value
    BIM_MATCH=rsid
    TRAIT_CLASS=uc
    TRAIT=UC
    ;;
  scz)
    INPUT="$UCC_DATA/gwas/scz_trubetskoy.tsv.gz"
    OUT_PREFIX="$UCC_DATA/gwas/scz_trubetskoy"
    # 76,755 cases + 53,880 controls; verify exact N at M1
    N_FIXED=130635
    COL_SNP=FIXME; COL_P=FIXME
    BIM_MATCH=rsid
    TRAIT_CLASS=schizophrenia
    TRAIT=SCZ
    ;;
  *)
    echo "ERROR: set GWAS=delange | liu | scz" >&2
    exit 1
    ;;
esac

echo "[01_magma] GWAS=$GWAS  input=$INPUT  bim-match=$BIM_MATCH"

# ---- 1. Reformat (autosomes, MAF/INFO QC, lambda_GC) ---------------------
#
# --bim re-anchors every SNP to g1000_eur's rsID+chr+bp (build37); SNPs
# not in the panel are dropped. See prepare_gwas.py --bim/--bim-match.

python code/01_magma/prepare_gwas.py \
    --input "$INPUT" \
    --out-prefix "$OUT_PREFIX" \
    --col-snp "$COL_SNP" --col-p "$COL_P" \
    --n-fixed "$N_FIXED" \
    --bim "$BIM" --bim-match "$BIM_MATCH" \
    --lambda-gc-out "$UCC_MAGMA_OUT/${GWAS}_lambda_gc.tsv"

# ---- 2. MAGMA annotate + gene-test ---------------------------------------

bash code/01_magma/run_magma.sh \
    "${OUT_PREFIX}.snp.loc" "${OUT_PREFIX}.pval" \
    "$UCC_DATA/reference/NCBI37.3.gene.loc" \
    "$UCC_DATA/reference/g1000_eur" \
    10 "$UCC_MAGMA_OUT/${GWAS}_10kb"

# ---- 3. Sanity check (top genes match expected pattern) ------------------

python code/01_magma/sanity_check.py \
    --genes-out "$UCC_MAGMA_OUT/${GWAS}_10kb.genes.out" \
    --gene-loc "$UCC_DATA/reference/NCBI37.3.gene.loc" \
    --trait-class "$TRAIT_CLASS"

# ---- 4. .gs file + gene-Z table (MHC excluded) ---------------------------

python code/01_magma/make_scdrs_gs.py \
    --genes-out "$UCC_MAGMA_OUT/${GWAS}_10kb.genes.out" \
    --gene-loc "$UCC_DATA/reference/NCBI37.3.gene.loc" \
    --trait "$TRAIT" --top-n 1000 \
    --out "$UCC_DATA/gwas/${GWAS}_top1000.gs" \
    --out-zscore-table "$UCC_DATA/gwas/${GWAS}_gene_z.tsv"

# ---- 5. MHC-included sensitivity gene set (de Lange only) ----------------

if [ "$GWAS" = "delange" ]; then
    python code/01_magma/make_scdrs_gs.py \
        --genes-out "$UCC_MAGMA_OUT/${GWAS}_10kb.genes.out" \
        --gene-loc "$UCC_DATA/reference/NCBI37.3.gene.loc" \
        --trait "$TRAIT" --top-n 1000 \
        --out "$UCC_DATA/gwas/${GWAS}_top1000_with_mhc.gs" \
        --keep-mhc
fi

echo "[01_magma] DONE for GWAS=$GWAS"
