#!/usr/bin/env bash
#SBATCH --job-name=ucc_scdrs
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=128x24
#SBATCH --account=128x24
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mukhinda@ucsc.edu
#
# scDRS compute-score + perform-downstream for one (ATLAS, GWAS, MHC) run.
#
# Hummingbird facts baked in (confirmed 2026-10-03):
#   - account and partition are both "128x24" (they genuinely share the name).
#   - 256x44 is the 247 GB partition — override --partition at submit for the
#     largest atlases (Pan-GI, HCA Gut); 128x24 at --mem=64G covers everything
#     else (Pan-GI compute-score peaked ~44 GB RSS).
#   - Time limits are infinite on every partition. 12 h is a backfill-friendly
#     cap (TAURUS 2x2 took ~4.5 h on GCP).
#   - Only miniconda3/3.13 and git/2.52.0 are available as modules; no system
#     R, no system MAGMA.
#
# Locked configuration (DECISIONS.md):
#   - all-cells policy (no disease-status subsetting)
#   - covariates: const + log_n_genes + log_n_counts + pre-expanded sample
#     dummies (donor DROPPED as collinear with sample; sex ABSENT from SCP259
#     and Garrido metadata). HCA Gut swaps sample for {assay, batch, Fraction,
#     sex} dummies — see run_hca_load.py.
#   - top-1000 MAGMA gene set (MHC-excluded by default; _with_mhc for the MHC
#     sensitivity)
#   - random seed = 42 for headline runs; override SEED=1|2|3 for test-retest
#   - n_ctrl = 1000 Monte Carlo null draws
#
# Submit examples:
#   sbatch --export=ALL,ATLAS=smillie,GWAS=delange       scripts/slurm/03_scdrs_compute.slurm
#   sbatch --export=ALL,ATLAS=garrido_trigo,GWAS=delange scripts/slurm/03_scdrs_compute.slurm
#   sbatch --export=ALL,ATLAS=taurus,GWAS=delange        scripts/slurm/03_scdrs_compute.slurm
#   sbatch -p 256x44 --export=ALL,ATLAS=pangi,GWAS=delange   scripts/slurm/03_scdrs_compute.slurm
#   sbatch -p 256x44 --export=ALL,ATLAS=hca_gut,GWAS=delange scripts/slurm/03_scdrs_compute.slurm
#   sbatch --export=ALL,ATLAS=smillie,GWAS=delange,GS_SUFFIX=_with_mhc \
#          scripts/slurm/03_scdrs_compute.slurm

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "$SCRIPT_DIR/../config.sh"

if command -v module >/dev/null 2>&1; then
    module load miniconda3/3.13
    module load git/2.52.0
fi
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$UCC_CONDA_ENV"

cd "$UCC_REPO"

# ---- Required parameters -------------------------------------------------

: "${ATLAS:?set ATLAS=smillie|garrido_trigo|taurus|hca_gut|pangi}"
: "${GWAS:?set GWAS=delange|liu|scz}"

GS_SUFFIX="${GS_SUFFIX:-}"
SEED="${SEED:-42}"

# MHC tag used in the output dir name. GS_SUFFIX is empty for MHC-excluded
# (the default; what seismic is paired against) and "_with_mhc" for the
# MHC-included sensitivity. Translate to a short, unambiguous dir tag.
if [ -z "$GS_SUFFIX" ]; then
    MHC="excl"
else
    MHC="incl"
fi

H5AD="$UCC_DATA/atlases/${ATLAS}.h5ad"
COV="$UCC_DATA/atlases/${ATLAS}_covariates.tsv"
GS="$UCC_DATA/gwas/${GWAS}_top1000${GS_SUFFIX}.gs"

# Standardised output naming (replaces per-atlas historical accidents like
# "garrido_delange_smoketest_group" and "smillie_delange_excl_fullcov_group").
# SCDRS_GROUP_DIRS in code/08_cross_method/run_scdrs_seismic_concordance.py
# must match this pattern.
OUT="$UCC_RESULTS/scdrs/${ATLAS}_${GWAS}_${MHC}"
GROUP_OUT="$UCC_RESULTS/scdrs/${ATLAS}_${GWAS}_${MHC}_group"

# scDRS will not create its own output dirs — it dies after computing if
# the folder is missing. Create both up front.
mkdir -p "$OUT" "$GROUP_OUT"

for f in "$H5AD" "$COV" "$GS"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: missing input $f" >&2
        exit 1
    fi
done

echo "[03_scdrs] atlas=$ATLAS  gwas=$GWAS  mhc=$MHC  seed=$SEED"
echo "[03_scdrs]   score dir: $OUT"
echo "[03_scdrs]   group dir: $GROUP_OUT"

# ---- 1. compute-score ----------------------------------------------------
#
# scdrs 1.0.2 CLI is positional with underscore flags here. Loaders bake in
# the scDRS-replica prefilter + sc.pp.pca(n_comps=20) + sc.pp.neighbors
# (n_neighbors=15, n_pcs=20) and stamp uns['scdrs_prefilter']; passing
# --flag_filter_data=True drops obsp and forces a ~45 min kNN rebuild per
# run (and OOMs). --flag_raw_count=False is uniform across atlases per the
# v1 policy (CELLxGENE deposits ship log-normalised; the raw-count atlases'
# loaders re-normalise so inputs are comparable).

scdrs compute-score \
    "$H5AD" human "$GS" human "$OUT/" \
    --cov_file="$COV" \
    --flag_filter_data=False \
    --flag_raw_count=False \
    --n_ctrl=1000

# ---- 2. perform-downstream (one call, both tiers) ------------------------
#
# perform-downstream takes DASHED flags on scdrs 1.0.2 — this asymmetry
# with compute-score (underscores) is how the CLI actually behaves, not a
# typo. Pass both tiers comma-separated to --group-analysis; one scDRS
# invocation is faster and avoids re-reading the h5ad per tier.

SCORE_FILE=$(ls "$OUT"/*.full_score.gz | head -1)

scdrs perform-downstream \
    "$H5AD" "$SCORE_FILE" "$GROUP_OUT/" \
    --group-analysis=cell_type_broad,cell_type_fine \
    --flag-filter-data=False \
    --flag-raw-count=False

# ---- 3. Aggregate per-null-draw cell-type Z (for Brown's method, M5) -----
#
# Saves a feather indexed by (null_draw_idx, cell_type) per tier. Thin
# wrapper; no-op if the script hasn't landed yet.
if [ -f code/03_scdrs/aggregate_null_draws.py ]; then
    for tier in broad fine; do
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

echo "[03_scdrs] DONE: $ATLAS x $GWAS  mhc=$MHC  seed=$SEED"
