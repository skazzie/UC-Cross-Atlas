#!/usr/bin/env bash
# Auto-fetch reference data and the GWAS files that are publicly
# downloadable without an account. Manual steps (atlases, registered
# GWAS) are listed at the end.
#
# Usage (from a Hummingbird login node, after setup_env.sh):
#   bash scripts/download_refs.sh
#
# Idempotent: existing files are skipped.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "$SCRIPT_DIR/config.sh"

cd "$UCC_DATA"

# Helper: download $1 to $2 if $2 doesn't already exist.
_fetch() {
    local url="$1"
    local out="$2"
    if [ -f "$out" ]; then
        echo "[download] skipping (exists): $out"
        return
    fi
    echo "[download] $out"
    mkdir -p "$(dirname "$out")"
    if ! wget -q --show-progress -O "$out.partial" "$url"; then
        echo "[download] WARNING: $url failed. Skipping; you'll need to fetch manually." >&2
        rm -f "$out.partial"
        return
    fi
    mv "$out.partial" "$out"
}

# ---- 1. 1000 Genomes Phase 3 EUR (MAGMA bfile) ---------------------------

if [ ! -f reference/g1000_eur.bed ]; then
    _fetch \
        "https://vu.data.surf.nl/s/VZNByNwpD8qqINe/download" \
        reference/g1000_eur.zip
    if [ -f reference/g1000_eur.zip ]; then
        unzip -o reference/g1000_eur.zip -d reference/
        rm reference/g1000_eur.zip
    fi
fi

# ---- 2. MAGMA binary -----------------------------------------------------

if [ ! -x reference/magma ]; then
    # Version may need updating; check https://ctg.cncr.nl/software/magma
    _fetch \
        "https://ctg.cncr.nl/software/MAGMA/prog/magma_v1.10.zip" \
        reference/magma.zip
    if [ -f reference/magma.zip ]; then
        mkdir -p reference/_magma_unpacked
        unzip -o reference/magma.zip -d reference/_magma_unpacked/
        cp reference/_magma_unpacked/magma reference/magma
        chmod +x reference/magma
        rm -rf reference/_magma_unpacked reference/magma.zip
    fi
fi

# ---- 3. NCBI37.3 gene location file --------------------------------------

if [ ! -f reference/NCBI37.3.gene.loc ]; then
    _fetch \
        "https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI37.3.zip" \
        reference/ncbi37_aux.zip
    if [ -f reference/ncbi37_aux.zip ]; then
        unzip -o -j reference/ncbi37_aux.zip -d reference/
        rm reference/ncbi37_aux.zip
    fi
fi

# ---- 4. Cell Ontology OWL ------------------------------------------------

_fetch \
    "http://purl.obolibrary.org/obo/cl.owl" \
    reference/cl.owl

# Pin the release date in DECISIONS.md per plan.
if [ -f reference/cl.owl ]; then
    cl_date=$(stat -c %y reference/cl.owl 2>/dev/null | cut -d' ' -f1 || date +%Y-%m-%d)
    echo "[download] Cell Ontology OWL downloaded $cl_date — record this in DECISIONS.md"
fi

# ---- 5. de Lange 2017 UC summary statistics (GWAS Catalog GCST004131) ----
# The harmonised file is hosted under the FTP tree. The exact filename
# can change with re-harmonisation; verify at:
#   https://www.ebi.ac.uk/gwas/studies/GCST004131
# A common pattern is:
#   ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/.../GCST004131/

cat <<'NOTE'

[download] de Lange 2017 UC (GCST004131): the harmonised .tsv.gz lives
under the GWAS Catalog FTP tree. The URL changes with re-harmonisation,
so this script does not auto-fetch. Visit:

    https://www.ebi.ac.uk/gwas/studies/GCST004131

and download the harmonised summary statistics to:

    $UCC_DATA/gwas/uc_delange_harmonised.tsv.gz

NOTE

# ---- 6. Manual-only items ------------------------------------------------

cat <<EOF

[download] Auto-fetchable refs done. Remaining manual steps for M1:

  Atlases (need cellxgene Discover or Single Cell Portal accounts)
  ----------------------------------------------------------------

  1. Smillie 2019 (SCP259) — https://singlecell.broadinstitute.org/
     Download to: \$UCC_DATA/atlases/smillie.h5ad
     (Free Single Cell Portal account required.)

  2. Garrido-Trigo 2023 UC subset (GSE214695 / CELLxGENE) —
     https://datasets.cellxgene.cziscience.com/b1a62801-f509-45f8-b55f-533fbb7e7800.h5ad
     UC subset filter applied by code/02_atlas_prep/load_garrido_trigo.py.
     Save to:
       \$UCC_DATA/atlases/garrido_trigo.h5ad

  3. Mennillo 2024 — verify GEO accession at M1; subset to
     pre-treatment baseline samples per DECISIONS.md before saving:
       \$UCC_DATA/atlases/mennillo.h5ad

  4. HCA Gut Cell Atlas (Elmentaite 2021) —
     https://cellxgene.cziscience.com/  (search "HCA gut")
     Filter to large-intestine cells, download as:
       \$UCC_DATA/atlases/hca_gut.h5ad

  5. Pan-GI (Oliver 2024) — same site, large-intestine subset:
       \$UCC_DATA/atlases/pangi.h5ad
     (~30 GB; use a tmux session, this takes a while.)

  GWAS (registration may be required)
  -----------------------------------

  6. Liu 2023 multi-ancestry IBD (UC arm) — Nat Genet 55:796-806.
     Check journal supplement and the GWAS Catalog. Save as:
       \$UCC_DATA/gwas/uc_liu.tsv.gz
     Verify per-SNP N column at M1.

  7. Trubetskoy 2022 schizophrenia (negative control) —
     https://www.med.unc.edu/pgc/download-results/
     Save as:
       \$UCC_DATA/gwas/scz_trubetskoy.tsv.gz

After every download, sanity-check column names with:
    zcat <file> | head -2

Then proceed to:
    sbatch --export=GWAS=delange scripts/slurm/01_magma.slurm
    sbatch --export=GWAS=liu     scripts/slurm/01_magma.slurm
    sbatch --export=GWAS=scz     scripts/slurm/01_magma.slurm

EOF
