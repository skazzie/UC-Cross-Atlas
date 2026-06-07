#!/usr/bin/env bash
#SBATCH --job-name=ucc_scdrs
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=128x24
#SBATCH --account=128x24
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amoli@ucsc.edu
#
# For Pan-GI (~1.1M cells) override --mem and --partition at submit:
#    sbatch -p 256x44 --mem=48G --time=06:00:00 \
#           --export=ALL,ATLAS=pangi,GWAS=delange,TIER_LIST=broad \
#           scripts/slurm/03_scdrs_compute.slurm
#
# Locked configuration (DECISIONS.md):
#   - all-cells policy (no disease-status subsetting)
#   - covariates: log_n_genes, log_n_counts, donor, sample, sex
#   - top-1000 MAGMA gene set (MHC-excluded by default)
#   - random seed = 42 for headline runs
#   - n_ctrl = 1000 Monte Carlo null draws
#
# Submit examples:
#   sbatch --export=ALL,ATLAS=smillie,GWAS=delange       scripts/slurm/03_scdrs_compute.slurm
#   sbatch --export=ALL,ATLAS=garrido_trigo,GWAS=delange scripts/slurm/03_scdrs_compute.slurm
#   sbatch --export=ALL,ATLAS=taurus,GWAS=delange        scripts/slurm/03_scdrs_compute.slurm
#   sbatch --export=ALL,ATLAS=smillie,GWAS=liu      scripts/slurm/03_scdrs_compute.slurm
#   ...
#
# Optional env vars:
#   GS_SUFFIX=_with_mhc   # for the MHC-included sensitivity (de Lange + Smillie only)
#   SEED=42               # override for test-retest seeds 1, 2, 3
#   TIER_LIST="broad fine"

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "$SCRIPT_DIR/../config.sh"

if command -v module >/dev/null 2>&1; then
    module load miniconda3 || module load anaconda3 || true
fi
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$UCC_CONDA_ENV"

cd "$UCC_REPO"

# ---- Required parameters -------------------------------------------------

: "${ATLAS:?set ATLAS=smillie|garrido_trigo|taurus|hca_gut|pangi}"
: "${GWAS:?set GWAS=delange|liu|scz}"

GS_SUFFIX="${GS_SUFFIX:-}"
SEED="${SEED:-42}"
TIER_LIST="${TIER_LIST:-broad fine}"

H5AD="$UCC_DATA/atlases/${ATLAS}.h5ad"
COV="$UCC_DATA/atlases/${ATLAS}_covariates.tsv"
GS="$UCC_DATA/gwas/${GWAS}_top1000${GS_SUFFIX}.gs"
OUT="$UCC_RESULTS/scdrs/${ATLAS}_${GWAS}${GS_SUFFIX}_seed${SEED}"

mkdir -p "$OUT"

for f in "$H5AD" "$COV" "$GS"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: missing input $f" >&2
        exit 1
    fi
done

echo "[03_scdrs] atlas=$ATLAS  gwas=$GWAS  seed=$SEED  out=$OUT"

# ---- 1. compute-score (the heavy step) -----------------------------------

# NOTE (DECISIONS.md correction 2026-05-20 (5/7)): CELLxGENE atlases ship
# log-normalized X with no raw layer. v1 policy locks --flag-raw-count False
# uniformly across all atlases. Smillie SCP259, Garrido RAW.tar, and TAURUS
# Zenodo all ship raw counts; their loaders re-normalize via log1p(CP10k)
# to keep input distributions comparable. (Per DECISIONS 16, TAURUS-IBD
# replaced the previously planned Mennillo as the third UC core atlas.)
scdrs compute-score \
    --h5ad-file "$H5AD" \
    --gs-file "$GS" \
    --cov-file "$COV" \
    --flag-filter-data True \
    --flag-raw-count False \
    --n-ctrl 1000 \
    --random-seed "$SEED" \
    --out-folder "$OUT/"

# ---- 2. Group analysis at each requested tier ----------------------------

# scDRS writes <TRAIT>.full_score.gz; figure out the trait label.
SCORE_FILE=$(ls "$OUT"/*.full_score.gz | head -1)

for tier in $TIER_LIST; do
    GROUP_OUT="${OUT}_${tier}"
    mkdir -p "$GROUP_OUT"
    scdrs perform-downstream \
        --h5ad-file "$H5AD" \
        --score-file "$SCORE_FILE" \
        --group-analysis "cell_type_${tier}" \
        --out-folder "$GROUP_OUT/"
done

# ---- 3. Aggregate per-null-draw cell-type Z (for Brown's method, M5) -----

# Saves a feather indexed by (null_draw_idx, cell_type) per tier.
# This is a thin wrapper script; if it doesn't exist yet at M3, this
# block is a no-op and you can backfill.
if [ -f code/03_scdrs/aggregate_null_draws.py ]; then
    for tier in $TIER_LIST; do
        python code/03_scdrs/aggregate_null_draws.py \
            --scdrs-out-folder "$OUT" \
            --h5ad-file "$H5AD" \
            --tier "$tier" \
            --out "$OUT/null_aggregations_${tier}.feather"
    done
else
    echo "[03_scdrs] note: code/03_scdrs/aggregate_null_draws.py not present yet"
    echo "[03_scdrs]       (M3 work — Brown's empirical correlation matrix needs this)"
fi

echo "[03_scdrs] DONE: $ATLAS x $GWAS  seed=$SEED"
