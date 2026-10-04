# Source this file to set common paths and SLURM parameters.
#
# Edit the values below for your Hummingbird account, OR override at the
# command line by exporting UCC_* env vars before sourcing.
#
# Get your partition / account values from:
#     sacctmgr show assoc user=$USER format=Cluster,Partition,Account,QOS

# ---- Project locations ---------------------------------------------------

# The cloned repo on Hummingbird. Lives on scratch (not $HOME) because
# $HOME has no space for atlas data and Hummingbird's convention is to
# work out of /hb/scratch/$USER/. Verified 2026-10-03 (rebuild handoff).
export UCC_REPO="${UCC_REPO:-/hb/scratch/mukhinda/UC-Cross-Atlas}"

# Big input files (atlases, references, GWAS summary stats) live under
# scratch/ inside the repo tree — matches the "scratch/data/..." relative
# paths baked into drivers and the handoff. Override UCC_SCRATCH if you
# want to put the data tree somewhere else.
export UCC_SCRATCH="${UCC_SCRATCH:-$UCC_REPO/scratch}"
export UCC_DATA="$UCC_SCRATCH/data"
export UCC_LOGS="$UCC_SCRATCH/logs"

# Analysis outputs (scDRS group tables, seismic TSVs, concordance CSV)
# live at REPO-ROOT results/ — that's what
# code/08_cross_method/run_scdrs_seismic_concordance.py resolves relative
# to the repo root (`results/scdrs/{group_dir}/...`,
# `results/seismic/{atlas}_{gwas}_{tier}.tsv`,
# `results/concordance/scdrs_vs_seismic.csv`), and .gitignore's
# `results/*` rule assumes the folder is at repo root (keeping
# `results/.gitkeep` + `results/MANIFEST.md` tracked). Don't move this
# into scratch — concordance will stop finding the files.
export UCC_RESULTS="${UCC_RESULTS:-$UCC_REPO/results}"

# MAGMA intermediates (`.genes.out`, `.lambda_gc.tsv`, `.snp.loc` etc.)
# are big and never read back from results/ — only `make_scdrs_gs.py`
# consumes them, and the downstream artifacts (`${gwas}_top1000.gs`,
# `${gwas}_gene_z.tsv`) land in `$UCC_DATA/gwas/`. Keep MAGMA output on
# scratch so it doesn't bloat the committed-folder tree.
export UCC_MAGMA_OUT="${UCC_MAGMA_OUT:-$UCC_SCRATCH/results/magma}"

# ---- SLURM parameters ----------------------------------------------------

# General-purpose partition. Hummingbird verified via `sinfo`:
#   128x24*    up infinite 18 nodes 122880 MB
#   96x24gpu4  up infinite  1 node   94200 MB
#   256x44     up infinite  1 node  252952 MB
export UCC_PARTITION="${UCC_PARTITION:-128x24}"

# High-memory partition for Pan-GI (~30 GB RAM, ~1.1M cells).
export UCC_PARTITION_HIGHMEM="${UCC_PARTITION_HIGHMEM:-256x44}"

# SLURM account (verified via `sacctmgr show assoc user=$USER`).
export UCC_ACCOUNT="${UCC_ACCOUNT:-128x24}"

# Email for SLURM notifications.
export UCC_EMAIL="${UCC_EMAIL:-mukhinda@ucsc.edu}"

# Conda env name.
export UCC_CONDA_ENV="${UCC_CONDA_ENV:-uc-cross-atlas}"

# ---- Bootstrap directories on first source -------------------------------

mkdir -p \
    "$UCC_SCRATCH" \
    "$UCC_DATA" \
    "$UCC_DATA/gwas" \
    "$UCC_DATA/atlases" \
    "$UCC_DATA/atlases/donor_metadata" \
    "$UCC_DATA/reference" \
    "$UCC_MAGMA_OUT" \
    "$UCC_RESULTS" \
    "$UCC_RESULTS/scdrs" \
    "$UCC_RESULTS/seismic" \
    "$UCC_RESULTS/concordance" \
    "$UCC_RESULTS/regime2" \
    "$UCC_RESULTS/figures" \
    "$UCC_LOGS"

# ---- Helpers -------------------------------------------------------------

# Optional --account flag for sbatch. Empty if not set.
ucc_account_flag() {
    if [ -n "${UCC_ACCOUNT:-}" ]; then
        echo "--account=$UCC_ACCOUNT"
    fi
}
