#!/usr/bin/env bash
#SBATCH --job-name=ucc_donor_loo
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --partition=128x24
#SBATCH --account=128x24
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=amoli@ucsc.edu
# ^^ Set the array bounds at submit time to match donor count:
#    Smillie ~30, garrido_trigo 12, Mennillo ~10-15. Per DECISIONS.md:
#       --array=0-29      for Smillie
#       --array=0-11      for garrido_trigo  (small-cohort guard applies; see below)
#       --array=0-N       for Mennillo (verify N at M1 after subsetting)
#
# Locked scope (DECISIONS.md):
#   - broad-tier headline metrics under de Lange only
#   - both methods (scDRS + seismicGWAS)
#   - Liu = point estimates only (no LOO)
#   - garrido_trigo: enforce >=5 cases per group (6 UC + 6 healthy total);
#     LOO iterations that drop below this are skipped.
#
# Submit examples:
#   sbatch --array=0-29 \
#       --export=ALL,ATLAS=smillie,METHOD=scdrs \
#       scripts/slurm/donor_loo_array.slurm
#
#   sbatch --array=0-11 \
#       --export=ALL,ATLAS=garrido_trigo,METHOD=seismic \
#       scripts/slurm/donor_loo_array.slurm

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

: "${ATLAS:?set ATLAS=smillie|garrido_trigo|mennillo}"
: "${METHOD:?set METHOD=scdrs|seismic}"
GWAS=delange
TIER=broad

DONOR_IDX="${SLURM_ARRAY_TASK_ID:?donor_loo expects to run as a SLURM array}"

DONOR_META="$UCC_DATA/atlases/donor_metadata/${ATLAS}_donor_metadata.csv"
if [ ! -f "$DONOR_META" ]; then
    echo "ERROR: missing $DONOR_META — generate during M1 atlas prep." >&2
    exit 1
fi

# Pull the donor_id at row $DONOR_IDX (0-based, skipping header).
DROPPED_DONOR=$(awk -F, -v i="$DONOR_IDX" 'NR > 1 { if (NR - 2 == i) print $1; }' "$DONOR_META")
if [ -z "$DROPPED_DONOR" ]; then
    echo "ERROR: could not resolve donor at index $DONOR_IDX in $DONOR_META" >&2
    exit 1
fi

echo "[donor_loo] atlas=$ATLAS method=$METHOD dropping_donor=$DROPPED_DONOR"

# garrido_trigo-specific guard: skip iterations that violate >=5 cases per group.
if [ "$ATLAS" = "garrido_trigo" ]; then
    GROUP=$(awk -F, -v d="$DROPPED_DONOR" 'NR > 1 && $1 == d { print $4 }' "$DONOR_META")
    if [ -z "$GROUP" ]; then
        echo "[donor_loo] no disease_status for $DROPPED_DONOR; skipping" >&2
        exit 0
    fi
    REMAINING_IN_GROUP=$(awk -F, -v g="$GROUP" -v d="$DROPPED_DONOR" \
        'NR > 1 && $4 == g && $1 != d { c++ } END { print c+0 }' "$DONOR_META")
    if [ "$REMAINING_IN_GROUP" -lt 5 ]; then
        echo "[donor_loo] garrido_trigo: dropping $DROPPED_DONOR would leave $REMAINING_IN_GROUP" \
             "in group=$GROUP; skipping per DECISIONS.md."
        exit 0
    fi
fi

# Driver script that re-runs the chosen method on the atlas with one
# donor masked out. Lives next to the method runner.
case "$METHOD" in
  scdrs)
    DRIVER="code/03_scdrs/run_donor_loo.py"
    ;;
  seismic)
    DRIVER="code/04_seismic/run_donor_loo.R"
    ;;
  *)
    echo "ERROR: METHOD must be scdrs|seismic" >&2
    exit 1
    ;;
esac

if [ ! -f "$DRIVER" ]; then
    echo "[donor_loo] note: $DRIVER not present yet (M3 work). Skipping." >&2
    exit 0
fi

OUT="$UCC_RESULTS/${METHOD}/donor_loo/${ATLAS}_${GWAS}_${TIER}/dropped_${DROPPED_DONOR}"
mkdir -p "$OUT"

if [ "$METHOD" = "scdrs" ]; then
    python "$DRIVER" \
        --atlas "$ATLAS" \
        --gwas "$GWAS" \
        --tier "$TIER" \
        --drop-donor "$DROPPED_DONOR" \
        --out-folder "$OUT"
else
    Rscript "$DRIVER" \
        --atlas "$ATLAS" \
        --gwas "$GWAS" \
        --tier "$TIER" \
        --drop-donor "$DROPPED_DONOR" \
        --out-dir "$OUT"
fi

echo "[donor_loo] DONE atlas=$ATLAS method=$METHOD dropped=$DROPPED_DONOR"
