#!/usr/bin/env bash
# One-shot TAURUS pooled h5ad download from Zenodo v3 (10.5281/zenodo.14007626).
# md5 of pooled file (per DECISIONS 16 pin, 2026-06-06):
#   c1bd13b92cacb164a401c6c4a4e7912c
set -euo pipefail

OUT_DIR="/mnt/beegfs/cluster/scratch/mukhinda/UC-Cross-Atlas/scratch/data/atlases"
OUT="$OUT_DIR/TAURUS_raw_counts_annotated_final.h5ad"
URL="https://zenodo.org/records/14007626/files/TAURUS_raw_counts_annotated_final.h5ad?download=1"
EXPECTED_MD5="c1bd13b92cacb164a401c6c4a4e7912c"

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

echo "[taurus_download] URL: $URL"
echo "[taurus_download] -> $OUT"

if [ -f "$OUT" ]; then
    echo "[taurus_download] already present; verifying md5"
else
    curl --retry 5 --continue-at - --fail -L "$URL" -o "$OUT.partial"
    mv "$OUT.partial" "$OUT"
fi

echo "[taurus_download] computing md5"
GOT_MD5="$(md5sum "$OUT" | awk '{print $1}')"
echo "[taurus_download] got $GOT_MD5"
echo "[taurus_download] exp $EXPECTED_MD5"
if [ "$GOT_MD5" != "$EXPECTED_MD5" ]; then
    echo "[taurus_download] HARD-GATE FAIL: md5 mismatch" >&2
    exit 2
fi
echo "[taurus_download] md5 OK  --  size:"
ls -la "$OUT"
echo "[taurus_download] DONE"
