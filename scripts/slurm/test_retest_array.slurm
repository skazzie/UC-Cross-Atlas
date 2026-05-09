#!/usr/bin/env bash
#SBATCH --job-name=ucc_test_retest
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --partition=128x24
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=PUT_YOUR_EMAIL_HERE
#SBATCH --array=1-3
# ^^ 3 seeds (1, 2, 3). Test-retest gates per DECISIONS.md:
#    scDRS pass >= 0.9; seismicGWAS pass >= 0.999.
#
# Run only on de Lange (DECISIONS.md). Submit one array per (atlas, method):
#   sbatch --export=ALL,ATLAS=smillie,METHOD=scdrs    scripts/slurm/test_retest_array.slurm
#   sbatch --export=ALL,ATLAS=smillie,METHOD=seismic  scripts/slurm/test_retest_array.slurm
#   ... and similarly for kong, mennillo.
#
# Total: 3 atlases x 2 methods x 3 seeds = 18 runs, all under de Lange.

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

: "${ATLAS:?set ATLAS=smillie|kong|mennillo}"
: "${METHOD:?set METHOD=scdrs|seismic}"
GWAS=delange

SEED="${SLURM_ARRAY_TASK_ID:-1}"

echo "[test_retest] atlas=$ATLAS method=$METHOD seed=$SEED"

case "$METHOD" in
  scdrs)
    H5AD="$UCC_DATA/atlases/${ATLAS}.h5ad"
    COV="$UCC_DATA/atlases/${ATLAS}_covariates.tsv"
    GS="$UCC_DATA/gwas/${GWAS}_top1000.gs"
    OUT="$UCC_RESULTS/scdrs/test_retest/${ATLAS}_${GWAS}_seed${SEED}"
    mkdir -p "$OUT"

    scdrs compute-score \
        --h5ad-file "$H5AD" \
        --gs-file "$GS" \
        --cov-file "$COV" \
        --flag-filter-data True \
        --flag-raw-count True \
        --n-ctrl 1000 \
        --random-seed "$SEED" \
        --out-folder "$OUT/"

    SCORE_FILE=$(ls "$OUT"/*.full_score.gz | head -1)
    for tier in broad fine; do
        scdrs perform-downstream \
            --h5ad-file "$H5AD" \
            --score-file "$SCORE_FILE" \
            --group-analysis "cell_type_${tier}" \
            --out-folder "${OUT}_${tier}/"
    done
    ;;

  seismic)
    # seismicGWAS is deterministic given fixed labels + expression matrix.
    # Test-retest "seed" really only varies whatever stochastic seeding
    # exists in calc_specificity / get_ct_trait_associations. Pass through.
    SCE="$UCC_DATA/atlases/${ATLAS}_sce.rds"
    GENE_Z="$UCC_DATA/gwas/${GWAS}_gene_z.tsv"
    OUT="$UCC_RESULTS/seismic/test_retest/${ATLAS}_${GWAS}_seed${SEED}"
    mkdir -p "$OUT"

    DRIVER="code/04_seismic/run_seismic.R"
    if [ ! -f "$DRIVER" ]; then
        echo "[test_retest] note: $DRIVER not present yet — write it in M3." >&2
        exit 0
    fi

    Rscript "$DRIVER" \
        --sce "$SCE" \
        --gene-z "$GENE_Z" \
        --tiers broad,fine \
        --out-dir "$OUT" \
        --seed "$SEED" \
        --permute 0
    ;;

  *)
    echo "ERROR: METHOD must be scdrs|seismic" >&2
    exit 1
    ;;
esac

echo "[test_retest] DONE atlas=$ATLAS method=$METHOD seed=$SEED"
