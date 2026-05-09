# Source this file to set common paths and SLURM parameters.
#
# Edit the values below for your Hummingbird account, OR override at the
# command line by exporting UCC_* env vars before sourcing.
#
# Get your partition / account values from:
#     sacctmgr show assoc user=$USER format=Cluster,Partition,Account,QOS

# ---- Project locations ---------------------------------------------------

# The cloned repo on Hummingbird (small files, $HOME quota is fine).
export UCC_REPO="${UCC_REPO:-$HOME/uc-cross-atlas}"

# Big files (atlases, references, results) live on scratch, not $HOME.
# Hummingbird scratch path; verify the convention for your group:
#   /hb/scratch/$USER/...   per-user
#   /hb/groups/<group>/...  shared
export UCC_SCRATCH="${UCC_SCRATCH:-/hb/scratch/$USER/uc-cross-atlas}"
export UCC_DATA="$UCC_SCRATCH/data"
export UCC_RESULTS="$UCC_SCRATCH/results"
export UCC_LOGS="$UCC_SCRATCH/logs"

# ---- SLURM parameters ----------------------------------------------------

# General-purpose partition. Verify with `sinfo`. Common Hummingbird
# names include "128x24", "Instruction", "240x44"; yours may differ.
export UCC_PARTITION="${UCC_PARTITION:-128x24}"

# High-memory partition for Pan-GI (~30 GB RAM, ~1.1M cells).
export UCC_PARTITION_HIGHMEM="${UCC_PARTITION_HIGHMEM:-256x44}"

# Account / project code if your cluster requires one. Leave blank to omit.
export UCC_ACCOUNT="${UCC_ACCOUNT:-}"

# Email for SLURM notifications.
export UCC_EMAIL="${UCC_EMAIL:-$USER@ucsc.edu}"

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
    "$UCC_RESULTS" \
    "$UCC_RESULTS/magma" \
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
