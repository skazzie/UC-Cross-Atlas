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
        "https://vu.data.surf.nl/s/lxDgt2dNdNr6DYt/download" \
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
        "https://vu.data.surf.nl/s/Pj2orwuF2JYyKxq/download" \
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

# ---- 5. de Lange 2017 UC summary statistics (GWAS Catalog GCST004133) ----
# NOTE 2026-06-06: this script previously listed GCST004131 here, but that
# is the IBD-combined arm (25,042 cases). The UC-only arm is GCST004133
# (12,366 UC cases / 33,609 controls). See DECISIONS correction (14).
# Schema (captured 2026-06-06): MarkerName, Allele1, Allele2, Effect,
# StdErr, P.value, Direction, HetISq, HetChiSq, HetDf, HetPVal,
# Pval_IBDseq, Pval_IIBDGC, Pval_GWAS3, Min_single_cohort_pval. NO
# per-SNP N — fixed N = 45,975 for all variants.

_fetch \
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004133/uc_build37_45975_20161107.txt.gz" \
    gwas/uc_delange_GCST004133.txt.gz

# ---- 5b. Liu 2023 multi-ancestry UC (GWAS Catalog GCST90446794) ----------
# Liu 2023 Nat Genet 55:796 — UC arm of the multi-ancestry IBD analysis.
# Schema (captured 2026-06-06): chromosome, base_pair_location,
# effect_allele, other_allele, beta, standard_error,
# effect_allele_frequency, p_value, variant_id, FreqSE, MinFreq, MaxFreq,
# Direction, HetISq, HetChiSq, HetDf, HetPVal. NO per-SNP N — fixed N =
# 375,508 (22,318 EAS + 353,190 EUR). Build GRCh38, 1-based.
# File is 2.49 GB uncompressed; allow time. See DECISIONS correction (14).

_fetch \
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90446001-GCST90447000/GCST90446794/GCST90446794.tsv" \
    gwas/uc_liu_GCST90446794.tsv

# ---- 5c. Yengo 2022 height EUR (GWAS Catalog GCST90245992) ---------------
# Positive-control GWAS. EUR-ancestry subset with full p-values.
# Schema (captured 2026-06-06): chromosome, base_pair_location,
# effect_allele, other_allele, beta, standard_error,
# effect_allele_frequency, p_value, variant_id, n. HAS per-SNP N column.
# Build GRCh37. 1,597,374 European samples. ~95 MB uncompressed.

_fetch \
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90245001-GCST90246000/GCST90245992/GCST90245992_buildGRCh37.tsv" \
    gwas/height_yengo_GCST90245992.tsv

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

  GWAS (Liu, de Lange, Yengo are now auto-fetched above; only SCZ is manual)
  --------------------------------------------------------------------------

  6. Trubetskoy 2022 schizophrenia (negative control) — Nature 604:502.
     GWAS Catalog accession GCST90128471, but full summary statistics are
     NOT hosted there (fullPvalueSet = False). Download from PGC instead:
       https://www.med.unc.edu/pgc/download-results/
     Registration required; the PGC3 SCZ file is the multi-ancestry meta
     (53,386 EUR + 14,004 EAS + 6,152 AA + 1,234 Latino cases). Save as:
       \$UCC_DATA/gwas/scz_trubetskoy.tsv.gz
     See DECISIONS correction (14).

After every download, sanity-check column names with:
    zcat <file> | head -2

Per-SNP N column status (DECISIONS 14):
  - Yengo height (GCST90245992):           HAS per-SNP 'n' column.
  - de Lange UC (GCST004133):              NO per-SNP N; fixed = 45,975.
  - Liu UC (GCST90446794):                 NO per-SNP N; fixed = 375,508.
  - Trubetskoy SCZ (PGC):                  check on download.

Then proceed to:
    sbatch --export=GWAS=delange scripts/slurm/01_magma.slurm
    sbatch --export=GWAS=liu     scripts/slurm/01_magma.slurm
    sbatch --export=GWAS=yengo   scripts/slurm/01_magma.slurm
    sbatch --export=GWAS=scz     scripts/slurm/01_magma.slurm

EOF
