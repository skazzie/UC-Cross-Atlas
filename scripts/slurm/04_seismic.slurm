#!/usr/bin/env bash
#SBATCH --job-name=ucc_seismic
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=128x24
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=PUT_YOUR_EMAIL_HERE
#
# seismicGWAS regime-1 pipeline for one (atlas, GWAS) at both tiers.
# Computes specificity (GWAS-independent, cached per atlas), runs
# get_ct_trait_associations() for the requested GWAS, and saves the
# specificity scores in long format for v2 reuse.
#
# For Brown's method support, also runs M=1000 gene-Z permutations
# (~30 min per (atlas, GWAS, tier)) when PERMUTE=1.
#
# Submit examples:
#   sbatch --export=ALL,ATLAS=smillie,GWAS=delange scripts/slurm/04_seismic.slurm
#   sbatch --export=ALL,ATLAS=kong,GWAS=delange,PERMUTE=1 scripts/slurm/04_seismic.slurm

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

: "${ATLAS:?set ATLAS=smillie|kong|mennillo|hca_gut|pangi}"
: "${GWAS:?set GWAS=delange|liu}"
PERMUTE="${PERMUTE:-0}"

SCE="$UCC_DATA/atlases/${ATLAS}_sce.rds"
GENE_Z="$UCC_DATA/gwas/${GWAS}_gene_z.tsv"
OUT="$UCC_RESULTS/seismic/${ATLAS}_${GWAS}"
SPEC_LONG="$UCC_REPO/code/04_seismic/specificity_long/${ATLAS}_specificity.feather"

mkdir -p "$OUT" "$(dirname "$SPEC_LONG")"

for f in "$SCE" "$GENE_Z"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: missing input $f" >&2
        exit 1
    fi
done

echo "[04_seismic] atlas=$ATLAS  gwas=$GWAS  permute=$PERMUTE"

# Driver script TBD in M3 — once written, this block invokes it.
DRIVER="code/04_seismic/run_seismic.R"
if [ ! -f "$DRIVER" ]; then
    cat <<EOF
[04_seismic] note: $DRIVER does not exist yet.
This SLURM template assumes a small R driver of the form:

    Rscript code/04_seismic/run_seismic.R \\
        --sce $SCE \\
        --gene-z $GENE_Z \\
        --tiers broad,fine \\
        --out-dir $OUT \\
        --specificity-long $SPEC_LONG \\
        --permute $PERMUTE \\
        --m-permutations 1000

Write that driver during M3 (see code/04_seismic/README.md). For now
this job exits successfully but does no work.
EOF
    exit 0
fi

Rscript "$DRIVER" \
    --sce "$SCE" \
    --gene-z "$GENE_Z" \
    --tiers broad,fine \
    --out-dir "$OUT" \
    --specificity-long "$SPEC_LONG" \
    --permute "$PERMUTE" \
    --m-permutations 1000

echo "[04_seismic] DONE"
