#!/usr/bin/env bash
# One-time environment setup on Hummingbird.
#
# Run from a login node (NOT a compute node), takes ~10-20 minutes the
# first time. Idempotent — re-running updates the conda env.
#
# Usage:
#   bash scripts/setup_env.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "$SCRIPT_DIR/config.sh"

echo "[setup_env] repo:    $UCC_REPO"
echo "[setup_env] scratch: $UCC_SCRATCH"
echo "[setup_env] env:     $UCC_CONDA_ENV"

# ---- 1. Modules ----------------------------------------------------------
# Adjust to whatever `module avail` shows on Hummingbird. miniconda3 is
# the most reliable conda provider on academic clusters; if Hummingbird
# uses anaconda3 or a different module name, change accordingly.

echo "[setup_env] loading modules"
if command -v module >/dev/null 2>&1; then
    module load miniconda3 || module load anaconda3 || true
fi

if ! command -v conda >/dev/null 2>&1; then
    echo "[setup_env] ERROR: conda not on PATH after module load." >&2
    echo "[setup_env] Run 'module avail' on Hummingbird and edit this" >&2
    echo "[setup_env] script to load the right conda module." >&2
    exit 1
fi

# Make 'conda activate' work in this non-interactive shell.
source "$(conda info --base)/etc/profile.d/conda.sh"

# ---- 2. Conda env --------------------------------------------------------

echo "[setup_env] creating/updating conda env: $UCC_CONDA_ENV"
if conda env list | awk '{print $1}' | grep -qx "$UCC_CONDA_ENV"; then
    echo "[setup_env] env exists; updating"
    conda env update -n "$UCC_CONDA_ENV" -f "$SCRIPT_DIR/environment.yml" --prune
else
    conda env create -n "$UCC_CONDA_ENV" -f "$SCRIPT_DIR/environment.yml"
fi

conda activate "$UCC_CONDA_ENV"
echo "[setup_env] active env: $(python -c 'import sys; print(sys.prefix)')"

# ---- 3. R packages from GitHub (not on conda) ----------------------------

echo "[setup_env] installing R packages from GitHub (seismicGWAS, scOntoMatch)"
Rscript "$SCRIPT_DIR/install_r_deps.R"

# ---- 4. Repo as editable Python package ----------------------------------

echo "[setup_env] pip install -e ."
cd "$UCC_REPO"
pip install -e ".[atlas,dev]"

# ---- 5. Smoke test -------------------------------------------------------

echo "[setup_env] smoke test: scDRS import + concordance pytest"
python -c "import scdrs; print('scdrs OK,', scdrs.__version__)"
pytest code/06_concordance/test_metrics.py -q

echo
echo "[setup_env] DONE."
echo
echo "Next steps:"
echo "  1. bash scripts/download_refs.sh         # auto-fetch 1000G + MAGMA + de Lange"
echo "  2. (manually fetch atlases + Liu + SCZ; see download_refs.sh notes)"
echo "  3. sbatch scripts/slurm/01_magma.slurm   # once GWAS files are in place"
