#!/usr/bin/env bash
# One-shot Pan-GI download from CELLxGENE (asset UUID
# 757945c8-a916-431d-aceb-1afbc80a7c55). The deposit UUID 404s on
# datasets.cellxgene.cziscience.com; the asset UUID is the correct
# download target (per DECISIONS 2026-10-03).
set -euo pipefail

OUT_DIR="/mnt/beegfs/cluster/scratch/mukhinda/UC-Cross-Atlas/scratch/data/atlases"
OUT="$OUT_DIR/757945c8-a916-431d-aceb-1afbc80a7c55.h5ad"
URL="https://datasets.cellxgene.cziscience.com/757945c8-a916-431d-aceb-1afbc80a7c55.h5ad"

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

echo "[pangi_download] URL: $URL"
echo "[pangi_download] -> $OUT"

if [ -f "$OUT" ]; then
    echo "[pangi_download] already present, verifying size"
else
    curl --retry 5 --continue-at - --fail -L "$URL" -o "$OUT.partial"
    mv "$OUT.partial" "$OUT"
fi

ls -la "$OUT"
echo "[pangi_download] DONE"
