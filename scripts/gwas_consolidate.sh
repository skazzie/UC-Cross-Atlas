#!/usr/bin/env bash
# One-shot consolidation: gzip Liu, download harmonised de Lange.
# Prep step before resubmitting MAGMA under the correct filenames.
set -euo pipefail

GWAS_DIR="/mnt/beegfs/cluster/scratch/mukhinda/UC-Cross-Atlas/scratch/data/gwas"
cd "$GWAS_DIR"

# --- Liu: produce uc_liu.tsv.gz from the already-downloaded raw .tsv ---
if [ -f uc_liu.tsv.gz ]; then
    echo "[consolidate] uc_liu.tsv.gz already present, skipping gzip"
elif [ -f uc_liu_GCST90446794.tsv ]; then
    echo "[consolidate] gzip -> uc_liu.tsv.gz (from uc_liu_GCST90446794.tsv)"
    gzip -c uc_liu_GCST90446794.tsv > uc_liu.tsv.gz.partial
    mv uc_liu.tsv.gz.partial uc_liu.tsv.gz
else
    echo "[consolidate] ERROR: no Liu source file to compress" >&2
    exit 1
fi

# --- de Lange: fetch harmonised build (hm_* columns, build37-consistent) ---
# NOTE 2026-10-07: 01_magma.slurm now reads the RAW build37 deposit
# uc_delange_GCST004133.txt.gz (already on disk from phase2) and bim-joins
# via chr:pos — see prepare_gwas.py --bim/--bim-match. The harmonised file
# below is kept for provenance / debugging but is no longer the MAGMA input.
DELANGE_URL="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004133/harmonised/28067908-GCST004133-EFO_0000729.h.tsv.gz"
if [ -f uc_delange_harmonised.tsv.gz ]; then
    echo "[consolidate] uc_delange_harmonised.tsv.gz already present, skipping"
else
    echo "[consolidate] fetching harmonised de Lange -> uc_delange_harmonised.tsv.gz"
    curl --retry 5 --continue-at - --fail -L "$DELANGE_URL" -o uc_delange_harmonised.tsv.gz.partial
    mv uc_delange_harmonised.tsv.gz.partial uc_delange_harmonised.tsv.gz
fi

# --- Liu: fetch harmonised deposit (has `rsid` column for bim-join) ---
# Raw Liu's variant_id is chr_pos_A1_A2 in GRCh38 — no rsID, so there's
# no safe build-37 join path from the raw file. Harmonised ships both
# `variant_id` and a real `rsid` column (e.g. rs1274919517). 01_magma.slurm
# bim-joins on rsid -> g1000_eur. See DECISIONS 2026-10-07.
LIU_HARM_URL="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90446001-GCST90447000/GCST90446794/harmonised/GCST90446794.h.tsv.gz"
if [ -f uc_liu_harmonised.tsv.gz ]; then
    echo "[consolidate] uc_liu_harmonised.tsv.gz already present, skipping"
else
    echo "[consolidate] fetching harmonised Liu -> uc_liu_harmonised.tsv.gz"
    curl --retry 5 --continue-at - --fail -L "$LIU_HARM_URL" -o uc_liu_harmonised.tsv.gz.partial
    mv uc_liu_harmonised.tsv.gz.partial uc_liu_harmonised.tsv.gz
fi

echo "[consolidate] DONE  --  gwas dir:"
ls -la "$GWAS_DIR"

# Quick header sanity on harmonised file. SIGPIPE guard: zcat | head
# returns 141 under set -euo pipefail (head closes the pipe before zcat
# finishes), misleading SLURM into marking a good step as failed. Mask
# zcat's "Broken pipe" write with 2>/dev/null and swallow the pipe exit
# with || true.
echo "[consolidate] harmonised de Lange header:"
zcat uc_delange_harmonised.tsv.gz 2>/dev/null | head -1 || true
