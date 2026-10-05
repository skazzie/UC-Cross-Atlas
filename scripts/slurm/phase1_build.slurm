#!/usr/bin/env bash
#SBATCH --job-name=ucc_phase1_env_v2
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=128x24
#SBATCH --account=128x24
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mukhinda@ucsc.edu

set -euo pipefail
set -x

module load miniconda3/3.13
module load git/2.52.0

CONDA_BASE="$(conda info --base)"
# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"

ENV_NAME="uc-cross-atlas"

# Nuke any partial env from the failed 791207 build so the solver can't
# silently preserve bad versions.
if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    echo "[phase1] removing stale env $ENV_NAME before clean rebuild"
    conda env remove -y -n "$ENV_NAME"
fi

echo "[phase1] creating env $ENV_NAME (python=3.10)"
conda create -y -n "$ENV_NAME" python=3.10

conda activate "$ENV_NAME"
echo "[phase1] active python: $(which python)  -> $(python --version)"

if ! command -v mamba >/dev/null 2>&1; then
    echo "[phase1] installing mamba in the env"
    conda install -y -n "$ENV_NAME" -c conda-forge mamba
fi

# ============================================================================
# Python side. numpy pinned <2.0 (ABI break on 2.0 has burned this project).
# ============================================================================
echo "[phase1] python-side packages via mamba (conda-forge)"
mamba install -y -n "$ENV_NAME" -c conda-forge \
    'numpy<2.0' \
    scanpy \
    anndata \
    pandas \
    scipy \
    h5py \
    pyarrow \
    scikit-learn \
    numba \
    igraph \
    leidenalg \
    pip

# scdrs 1.0.2 is required — 1.0.4 is GitHub-only and has a different CLI.
# EmpiricalBrownsMethod is NOT on PyPI (confirmed in code/07_regime2_meta/
# run_brown.py:11-14 — Brown's method is implemented inline). R side has
# bioconductor-empiricalbrownsmethod.
echo "[phase1] scdrs 1.0.2 via python -m pip"
python -m pip install --no-cache-dir 'scdrs==1.0.2'

# ============================================================================
# R side. r-remotes is the missing-piece fix for job 791207 — devtools 2.5+
# delegates install_github to remotes and errors out if remotes isn't present.
# Do NOT install zellkonverter — basilisk tries to source-build Python
# without root and ignores RETICULATE_PYTHON. The flat-export chain in
# code/04_seismic/ exists specifically to bypass it.
# ============================================================================
echo "[phase1] R + Bioconductor via mamba (prebuilt binaries from conda-forge + bioconda)"
mamba install -y -n "$ENV_NAME" -c conda-forge -c bioconda \
    r-base \
    r-devtools \
    r-remotes \
    r-optparse \
    r-matrix \
    r-data.table \
    r-dplyr \
    r-biocmanager \
    bioconductor-singlecellexperiment \
    bioconductor-summarizedexperiment \
    bioconductor-s4vectors \
    bioconductor-empiricalbrownsmethod

# seismicGWAS via remotes directly, skipping the devtools delegation layer
# that failed 791207. r-pak is a fallback (linked in by r-devtools).
echo "[phase1] seismicGWAS from GitHub via remotes::install_github"
if ! R --vanilla -e 'remotes::install_github("ylaboratory/seismic", upgrade="never")'; then
    echo "[phase1] remotes install failed; falling back to pak"
    R --vanilla -e 'pak::pak("ylaboratory/seismic")'
fi

# ============================================================================
# Verification
# ============================================================================
echo "[phase1] ================= VERIFICATION ================="
which scdrs
scdrs --help 2>&1 | head -5
python -c "import scanpy, anndata, numpy, sklearn, numba, igraph, leidenalg; print('numpy:', numpy.__version__, '| scanpy:', scanpy.__version__, '| anndata:', anndata.__version__)"
R --vanilla -e 'suppressPackageStartupMessages({
  library(seismicGWAS); library(remotes); library(optparse);
  library(Matrix); library(SingleCellExperiment);
  library(data.table); library(dplyr)
})
cat("seismicGWAS version:", as.character(packageVersion("seismicGWAS")), "\n")
print(sessionInfo())'

echo "[phase1] DONE"
