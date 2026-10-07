#!/usr/bin/env bash
#SBATCH --job-name=ucc_phase2_download
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=128x24
#SBATCH --account=128x24
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mukhinda@ucsc.edu
#
# Wraps scripts/download_refs.sh in a batch job so it survives SSH drops
# and so MAGMA jobs can chain via --dependency=afterok:<this_jobid>.
# Downloads 1000G EUR LD panel, MAGMA binary, NCBI37.3 gene.loc,
# Cell Ontology OWL, and the de Lange / Liu / Yengo / SCZ GWAS files.

set -euo pipefail
set -x

module load miniconda3/3.13
module load git/2.52.0

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

# Activate env so `unzip` (installed via mamba 2026-10-07) is on PATH —
# compute nodes don't have /usr/bin/unzip; job 797221 died on this.
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$UCC_CONDA_ENV"

cd "$UCC_REPO"
bash scripts/download_refs.sh

echo "[phase2_download] DONE"
