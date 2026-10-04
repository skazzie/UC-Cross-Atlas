#!/usr/bin/env bash
#SBATCH --job-name=ucc_seismic
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --partition=128x24
#SBATCH --account=128x24
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mukhinda@ucsc.edu
#
# seismicGWAS regime-1 pipeline — one atlas per task via the FLAT-EXPORT
# chain. Three steps:
#   (1) python h5ad_to_sce_export.py   — anndata reads h5ad in-env and
#       drops counts.mtx / logcounts.mtx / genes.tsv / barcodes.tsv /
#       obs.tsv (expect ~340M nonzero for a 230k-cell atlas → the --mem=64G
#       request lives here).
#   (2) Rscript sce_from_export.R      — pure Matrix::readMM, builds a
#       SingleCellExperiment and saves ${ATLAS}_sce.rds.
#   (3) Rscript run_seismic.R          — loops over {delange, liu} x
#       {broad, fine} = 4 inner runs, consuming the .rds from (2).
#
# WHY the flat-export detour: run_seismic.R's --h5ad-path path loads via
# zellkonverter::readH5AD, which drags in basilisk, which ignores
# RETICULATE_PYTHON and tries to source-build Python 3.14 without root.
# That failure blocked seismic from June until August on GCP. We NEVER
# pass --h5ad-path from this slurm — always --sce-rds pointing at (2).
#
# --magma-z MUST be the symbol-keyed ${GWAS}_gene_z.tsv (SYMBOL + ZSTAT
# columns) written by code/01_magma/make_scdrs_gs.py. Do NOT pass
# .genes.out directly — MAGMA keys those on Entrez IDs, which the
# specificity matrix doesn't carry, producing "1 overlapping gene" and
# an all-NA output.
#
# Memory peaks in step (1) (anndata + the mtx write). seismic itself and
# the SCE build are modest. Per-atlas run time: ~1 h for the four
# analytic inner runs; permutations extend it to roughly 2–3 h depending
# on atlas size — the --time=04:00:00 cap has headroom.
#
# Permutations default ON (PERMUTE=1), with --permutations 1000 to match
# the locked SEISMIC_N_PERMUTATIONS in code/_shared/constants.py. The
# GCP runs never set this flag, so the single FDR-significant seismic
# result in the project (TAURUS monocyte/macrophage, FDR 0.012) has only
# analytic p-values and can't carry a headline claim yet. Set PERMUTE=0
# for a fast smoke check.
#
# Submit examples:
#   sbatch --array=0-4                           scripts/slurm/04_seismic.slurm
#   sbatch --export=ALL,ATLAS=taurus             scripts/slurm/04_seismic.slurm
#   sbatch -p 256x44 --export=ALL,ATLAS=pangi    scripts/slurm/04_seismic.slurm
#   sbatch --export=ALL,ATLAS=smillie,PERMUTE=0  scripts/slurm/04_seismic.slurm   # smoke

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

# seismicGWAS pulls in reticulate indirectly through some deps; pin
# the Python it would see at our env's interpreter so no code path
# accidentally wakes basilisk looking for an outside Python.
export RETICULATE_PYTHON="$CONDA_PREFIX/bin/python"

cd "$UCC_REPO"

# ---- Resolve ATLAS from array index or env -------------------------------

ATLASES=(smillie garrido_trigo taurus pangi hca_gut)

if [ -n "${SLURM_ARRAY_TASK_ID:-}" ] && [ -z "${ATLAS:-}" ]; then
    ATLAS="${ATLASES[$SLURM_ARRAY_TASK_ID]}"
fi

: "${ATLAS:?set ATLAS=smillie|garrido_trigo|taurus|pangi|hca_gut (or submit with --array=0-4)}"
PERMUTE="${PERMUTE:-1}"
SEED="${SEED:-42}"

H5AD="$UCC_DATA/atlases/${ATLAS}.h5ad"
EXPORT_DIR="$UCC_DATA/seismic_export/${ATLAS}"
SCE_RDS="$UCC_DATA/seismic_export/${ATLAS}_sce.rds"
OUT="$UCC_RESULTS/seismic"

mkdir -p "$EXPORT_DIR" "$OUT"
if [ "$PERMUTE" = "1" ]; then
    mkdir -p "$OUT/permutations"
fi

if [ ! -f "$H5AD" ]; then
    echo "ERROR: missing input $H5AD" >&2
    exit 1
fi

echo "[04_seismic] atlas=$ATLAS permute=$PERMUTE seed=$SEED"
echo "[04_seismic]   h5ad:       $H5AD"
echo "[04_seismic]   export dir: $EXPORT_DIR"
echo "[04_seismic]   sce rds:    $SCE_RDS"
echo "[04_seismic]   out dir:    $OUT"

# ---- 1. h5ad -> flat export (python / anndata) ---------------------------
#
# ~340M-nonzero mtx for TAURUS (~230k cells); Pan-GI will be larger.
# This is the memory peak of the pipeline — downstream R work on the
# .rds is modest.

python code/04_seismic/h5ad_to_sce_export.py \
    --h5ad "$H5AD" \
    --out-dir "$EXPORT_DIR"

# ---- 2. flat export -> SingleCellExperiment .rds (pure R) ----------------
#
# No zellkonverter, no basilisk — Matrix::readMM reads counts.mtx and
# logcounts.mtx from (1). Saves one .rds we reuse across every (gwas,
# tier) inner run below.

Rscript code/04_seismic/sce_from_export.R \
    --export-dir "$EXPORT_DIR" \
    --out-rds "$SCE_RDS"

# ---- 3. seismic regression: 2 GWAS x 2 tiers from the shared .rds --------
#
# Pass --sce-rds (not --h5ad-path — the latter wakes zellkonverter and
# dies on basilisk). --magma-z must be the SYMBOL-keyed _gene_z.tsv;
# .genes.out would key on Entrez and overlap 1 gene with the specificity
# matrix. --run-permutations enables the M=1000 shuffle path per
# SEISMIC_N_PERMUTATIONS; omit for a quick smoke (set PERMUTE=0).

PERM_ARGS=()
if [ "$PERMUTE" = "1" ]; then
    PERM_ARGS=(--run-permutations --permutations 1000)
else
    echo "[04_seismic] PERMUTE=0 — skipping --run-permutations (analytic p-values only)."
fi

for GWAS in delange liu; do
    GENE_Z="$UCC_DATA/gwas/${GWAS}_gene_z.tsv"
    if [ ! -f "$GENE_Z" ]; then
        echo "ERROR: missing MAGMA gene-Z $GENE_Z — run 01_magma.slurm for GWAS=$GWAS first" >&2
        exit 1
    fi
    for TIER in broad fine; do
        echo "[04_seismic] run: atlas=$ATLAS gwas=$GWAS tier=$TIER"
        Rscript code/04_seismic/run_seismic.R \
            --sce-rds "$SCE_RDS" \
            --magma-z "$GENE_Z" \
            --atlas "$ATLAS" \
            --gwas "$GWAS" \
            --tier "$TIER" \
            --out-dir "$OUT" \
            --seed "$SEED" \
            "${PERM_ARGS[@]}"
    done
done

echo "[04_seismic] DONE: $ATLAS (4 inner runs)"
