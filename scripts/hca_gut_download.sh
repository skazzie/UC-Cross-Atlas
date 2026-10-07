#!/usr/bin/env bash
# One-shot HCA Gut Cell Atlas download from CELLxGENE (asset UUID
# f34d2b82-9265-4a73-bda4-852933bf2a8d).
set -euo pipefail

OUT_DIR="/mnt/beegfs/cluster/scratch/mukhinda/UC-Cross-Atlas/scratch/data/atlases"
OUT="$OUT_DIR/f34d2b82-9265-4a73-bda4-852933bf2a8d.h5ad"
URL="https://datasets.cellxgene.cziscience.com/f34d2b82-9265-4a73-bda4-852933bf2a8d.h5ad"

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

echo "[hca_gut_download] URL: $URL"
echo "[hca_gut_download] -> $OUT"

if [ -f "$OUT" ]; then
    echo "[hca_gut_download] already present, verifying size"
else
    curl --retry 5 --continue-at - --fail -L "$URL" -o "$OUT.partial"
    mv "$OUT.partial" "$OUT"
fi

ls -la "$OUT"
echo "[hca_gut_download] DONE"
