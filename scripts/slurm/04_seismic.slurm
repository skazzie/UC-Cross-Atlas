#!/usr/bin/env bash
#SBATCH --job-name=ucc_seismic
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --partition=128x24
#SBATCH --account=128x24
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amoli@ucsc.edu
#
# seismicGWAS regime-1 pipeline — one (atlas, GWAS, tier) per task.
#
# Default array covers all 12 UC-trio jobs:
#   ATLAS in {smillie, garrido_trigo, mennillo}
#   GWAS  in {delange, liu}
#   TIER  in {broad, fine}
#
# Submit the full sweep:
#   sbatch --array=0-11 scripts/slurm/04_seismic.slurm
#
# Or a single job (broad atlas comparator, no array):
#   sbatch --export=ALL,ATLAS=hca_gut,GWAS=delange,TIER=broad scripts/slurm/04_seismic.slurm
#
# Permutations are on by default (M=1000 per task, ~30 min). Disable with
#   PERMUTE=0 for a quick smoke check.

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

# ---- Resolve (ATLAS, GWAS, TIER) from array index OR env -----------------

ATLASES=(smillie smillie smillie smillie garrido_trigo garrido_trigo garrido_trigo garrido_trigo mennillo mennillo mennillo mennillo)
GWASES=(delange  delange  liu      liu      delange       delange       liu           liu           delange  delange  liu      liu)
TIERS=(broad     fine     broad    fine     broad         fine          broad         fine          broad    fine     broad    fine)

if [ -n "${SLURM_ARRAY_TASK_ID:-}" ] && [ -z "${ATLAS:-}" ]; then
    IDX="$SLURM_ARRAY_TASK_ID"
    ATLAS="${ATLASES[$IDX]}"
    GWAS="${GWASES[$IDX]}"
    TIER="${TIERS[$IDX]}"
fi

: "${ATLAS:?set ATLAS=smillie|garrido_trigo|mennillo|hca_gut|pangi (or submit with --array=0-11)}"
: "${GWAS:?set GWAS=delange|liu}"
: "${TIER:=broad}"
PERMUTE="${PERMUTE:-1}"
SEED="${SEED:-42}"

H5AD="$UCC_DATA/atlases/${ATLAS}.h5ad"
GENE_Z="$UCC_DATA/gwas/${GWAS}_gene_z.tsv"
OUT="$UCC_RESULTS/seismic"

mkdir -p "$OUT" "$OUT/permutations"

for f in "$H5AD" "$GENE_Z"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: missing input $f" >&2
        exit 1
    fi
done

echo "[04_seismic] atlas=$ATLAS gwas=$GWAS tier=$TIER permute=$PERMUTE"

DRIVER="code/04_seismic/run_seismic.R"
if [ ! -f "$DRIVER" ]; then
    echo "ERROR: $DRIVER missing" >&2
    exit 1
fi

PERM_ARG=()
if [ "$PERMUTE" != "1" ]; then
    PERM_ARG=(--skip-permutations)
fi

Rscript "$DRIVER" \
    --atlas "$ATLAS" \
    --gwas "$GWAS" \
    --tier "$TIER" \
    --h5ad-path "$H5AD" \
    --magma-z "$GENE_Z" \
    --out-dir "$OUT" \
    --permutations 1000 \
    --seed "$SEED" \
    "${PERM_ARG[@]}"

echo "[04_seismic] DONE ($ATLAS $GWAS $TIER)"
