#!/usr/bin/env bash
#SBATCH --job-name=ucc_atlas_prep
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=128x24
#SBATCH --account=128x24
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mukhinda@ucsc.edu
#
# Build per-atlas <ATLAS>.h5ad + <ATLAS>_covariates.tsv using the
# loaders in code/02_atlas_prep/. The output filenames match the
# convention 03_scdrs_compute.slurm / 04_seismic.slurm consume
# ($UCC_DATA/atlases/${ATLAS}.h5ad and ${ATLAS}_covariates.tsv),
# so downstream jobs find them without further coordination.
#
# Submit:
#   sbatch --export=ALL,ATLAS=garrido_trigo \
#          [--export=ALL,ATLAS=...,EXPECT_CELLS=N,EXPECT_GENES=M,EXPECT_DONORS=K] \
#          scripts/slurm/02_atlas_prep.slurm
#
# ATLAS values:
#   garrido_trigo   — GSE214695 RAW.tar + annotation CSV (both on disk)
#   smillie         — SCP259 directory tree under scratch/data/atlases/smillie/
#   taurus          — Zenodo pooled h5ad TAURUS_raw_counts_annotated_final.h5ad
#   pangi           — CELLxGENE 757945c8-a916-431d-aceb-1afbc80a7c55.h5ad
#   hca_gut         — CELLxGENE f34d2b82-9265-4a73-bda4-852933bf2a8d.h5ad
#
# Optional hard gates (set to enforce — exit non-zero on mismatch,
# so --dependency=afterok:<this> downstream chains refuse to fire on
# a wrong-shape load):
#   EXPECT_CELLS, EXPECT_GENES, EXPECT_DONORS

set -euo pipefail
set -x

module load miniconda3/3.13

# Prefer SLURM_SUBMIT_DIR; under SLURM, BASH_SOURCE points at the batch
# script copy in /var/spool/slurmd (see 797217 post-mortem).
if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    # shellcheck disable=SC1091
    source "$SLURM_SUBMIT_DIR/scripts/config.sh"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    # shellcheck disable=SC1091
    source "$SCRIPT_DIR/../config.sh"
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$UCC_CONDA_ENV"

cd "$UCC_REPO"

: "${ATLAS:?set ATLAS=garrido_trigo|smillie|taurus|pangi|hca_gut}"

H5AD="$UCC_DATA/atlases/${ATLAS}.h5ad"
COV="$UCC_DATA/atlases/${ATLAS}_covariates.tsv"

echo "[02_atlas_prep] atlas=$ATLAS  -> $H5AD  +  $COV"

case "$ATLAS" in
  garrido_trigo)
    # Loader gate — matches CONTINUITY VM-run shape (verified 2026-10-07).
    # Garrido's loader does NOT apply sc.pp.filter_* (see run_garrido_load.py);
    # the published matrix is pre-QC'd and the gate is on loader-raw shape.
    EXPECT_CELLS="${EXPECT_CELLS:-29675}"
    EXPECT_GENES="${EXPECT_GENES:-22414}"
    EXPECT_DONORS="${EXPECT_DONORS:-12}"
    python code/02_atlas_prep/run_garrido_load.py \
        --tar "$UCC_DATA/atlases/GSE214695_RAW.tar" \
        --csv "$UCC_DATA/atlases/GSE214695_cell_annotation.csv.gz" \
        --out-h5ad "$H5AD" \
        --out-cov "$COV"
    ;;
  smillie)
    python code/02_atlas_prep/run_smillie_load.py \
        --scp259-dir "$UCC_DATA/atlases/smillie/SCP259" \
        --out-h5ad "$H5AD" \
        --out-cov "$COV"
    ;;
  taurus)
    python code/02_atlas_prep/run_taurus_load.py \
        --h5ad "$UCC_DATA/atlases/TAURUS_raw_counts_annotated_final.h5ad" \
        --out-h5ad "$H5AD" \
        --out-cov "$COV"
    ;;
  pangi)
    python code/02_atlas_prep/run_pangi_load.py \
        --h5ad "$UCC_DATA/atlases/757945c8-a916-431d-aceb-1afbc80a7c55.h5ad" \
        --out-h5ad "$H5AD" \
        --out-cov "$COV"
    ;;
  hca_gut)
    python code/02_atlas_prep/run_hca_load.py \
        --h5ad "$UCC_DATA/atlases/f34d2b82-9265-4a73-bda4-852933bf2a8d.h5ad" \
        --out-h5ad "$H5AD" \
        --out-cov "$COV"
    ;;
  *)
    echo "ERROR: unknown ATLAS=$ATLAS" >&2
    exit 1
    ;;
esac

echo "[02_atlas_prep] load DONE for $ATLAS"

# Optional hard-gate on cells/genes/donors. Any set -> all checked.
# Downstream chains should use --dependency=afterok:<this>; a hard-gate
# miss returns non-zero and prevents fan-out to a wrong-shape atlas.
if [ -n "${EXPECT_CELLS:-}${EXPECT_GENES:-}${EXPECT_DONORS:-}" ]; then
    python - "$H5AD" "${EXPECT_CELLS:-}" "${EXPECT_GENES:-}" "${EXPECT_DONORS:-}" <<'PYEOF'
import sys, anndata as ad
h5ad, ec, eg, ed = sys.argv[1:5]
a = ad.read_h5ad(h5ad)
nc, ng = a.n_obs, a.n_vars
if "donor" in a.obs:
    nd = int(a.obs["donor"].nunique())
elif "donor_id" in a.obs:
    nd = int(a.obs["donor_id"].nunique())
else:
    nd = -1
print(f"[02_atlas_prep] actual: cells={nc}  genes={ng}  donors={nd}")
fails = []
if ec and int(ec) != nc:
    fails.append(f"cells {nc} != expected {ec}")
if eg and int(eg) != ng:
    fails.append(f"genes {ng} != expected {eg}")
if ed and int(ed) != nd:
    fails.append(f"donors {nd} != expected {ed}")
if fails:
    print("[02_atlas_prep] HARD-GATE FAIL: " + "; ".join(fails), file=sys.stderr)
    sys.exit(2)
print("[02_atlas_prep] HARD-GATE OK")
PYEOF
fi

echo "[02_atlas_prep] $ATLAS DONE"
