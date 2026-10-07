#!/usr/bin/env bash
# One-shot SCP259 Smillie download. Token in URL expires in minutes, so
# this script is submitted once at paste-time and not reused.
set -euo pipefail

cd /mnt/beegfs/cluster/scratch/mukhinda/UC-Cross-Atlas/scratch/data/atlases/smillie

echo "[smillie_fetch] fetching curl config from SCP"
curl --retry 3 --fail "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP259&auth_code=3k2sHfmx&directory=all&context=study" -o cfg.txt

echo "[smillie_fetch] cfg.txt size: $(wc -c < cfg.txt) bytes, head:"
head -20 cfg.txt

echo "[smillie_fetch] downloading via curl -K cfg.txt"
curl --retry 3 -K cfg.txt

echo "[smillie_fetch] cleaning up cfg.txt"
rm cfg.txt

echo "[smillie_fetch] DONE  --  tree:"
ls -laR | head -80
du -sh .
