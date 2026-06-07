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

# ---- 5d. Trubetskoy 2022 SCZ EUR (PGC3 wave3) ---------------------------
# Negative-control GWAS. EUR-ancestry subset, PGC3 wave 3.
# Hosted on **public figshare** (DOI 10.6084/m9.figshare.19426775.v7,
# CC BY 4.0, is_public=true, is_embargoed=false). No registration
# required — the "PGC Terms and Conditions" mentioned on
# https://pgc.unc.edu/for-researchers/download-results/ are a
# Fort Lauderdale honor-code researcher pledge, not an access gate.
# The publication embargo cleared when Trubetskoy 2022 (Nature 604:502)
# landed; free for genome-wide analyses now. See DECISIONS (19) for the
# (18)(d) retraction.
#
# File: PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz, 240 MB,
# md5 6ebe2376f5cda972d37efa0f214c4df0.
#
# Format: PGC sumstats VCF v1.0 (PGC-internal extension of VCF —
# non-standard; munge step needs to extract the 14 columns from rows
# rather than treat as standard GWAS-Catalog SSF). 14 columns captured
# 2026-06-06: chr, rsid, pos, A1, A2, frq_A, frq_U, info, beta, se, p,
# n_case_cohort (52017+1369), n_ctrl_cohort (75889+1369), n_eff_per_snp.
# Build GRCh37. **HAS per-SNP N_eff** in the last column.

_fetch \
    "https://ndownloader.figshare.com/files/34517828" \
    gwas/scz_trubetskoy_eur_PGC3_v3.vcf.tsv.gz

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

  3. TAURUS-IBD (Thomas 2024) — Zenodo deposit; replaces the previously
     planned Mennillo 2024 per DECISIONS correction (16). Pinned to Zenodo
     v3: 10.5281/zenodo.14007626 (TAURUS_raw_counts_annotated_final.h5ad
     md5 c1bd13b92cacb164a401c6c4a4e7912c, 12.7 GB). CC-BY-4.0, public,
     27.9 GB total across the pooled file + per-lineage h5ads. Subset to UC donors
     only + a single time-point per donor on load. Save as:
       \$UCC_DATA/atlases/taurus.h5ad

  4. HCA Gut Cell Atlas (Elmentaite 2021) —
     https://cellxgene.cziscience.com/  (search "HCA gut")
     Filter to large-intestine cells, download as:
       \$UCC_DATA/atlases/hca_gut.h5ad

  5. Pan-GI (Oliver 2024) — same site, large-intestine subset:
       \$UCC_DATA/atlases/pangi.h5ad
     (~30 GB; use a tmux session, this takes a while.)

  GWAS (Liu, de Lange, Yengo are now auto-fetched above; only SCZ is manual)
  --------------------------------------------------------------------------

  6. Trubetskoy 2022 SCZ is now auto-fetched above (5d). The previous
     belief that PGC required multi-day registration was wrong — the
     summary stats are public on figshare, no registration. See
     DECISIONS (19).

After every download, sanity-check column names with:
    zcat <file> | head -2

Per-SNP N column status (DECISIONS 14, 19):
  - Yengo height (GCST90245992):           HAS per-SNP 'n' column.
  - de Lange UC (GCST004133):              NO per-SNP N; fixed = 45,975
                                           (N_eff ≈ 36,168 — DECISIONS 16).
  - Liu UC (GCST90446794):                 NO per-SNP N; fixed = 375,508
                                           (N_eff ≈ 87,242 — DECISIONS 16).
  - Trubetskoy SCZ EUR (PGC3 figshare):    HAS per-SNP n_eff (last column).
                                           Cohort: 53,386 cases (incl. trios)
                                           + 77,258 controls; N_eff ≈ 123,000.

Then proceed to:
    sbatch --export=GWAS=delange scripts/slurm/01_magma.slurm
    sbatch --export=GWAS=liu     scripts/slurm/01_magma.slurm
    sbatch --export=GWAS=yengo   scripts/slurm/01_magma.slurm
    sbatch --export=GWAS=scz     scripts/slurm/01_magma.slurm

EOF
