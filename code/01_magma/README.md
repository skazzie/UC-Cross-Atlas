# 01_magma — GWAS preprocessing → gene scores → scDRS / seismicGWAS inputs

Pipeline for converting GWAS summary statistics into:

- a top-1000 MAGMA gene set (`.gs` file) for scDRS
- a long-format gene-Z-score table for seismicGWAS

Locked v1 commitments (see DECISIONS.md):

- **Four GWAS**: de Lange 2017 UC (primary), Liu 2023 multi-ancestry UC
  arm (cross-GWAS sensitivity), Yengo 2022 height EUR (positive control),
  Trubetskoy 2022 schizophrenia (negative control).
- **Per-GWAS N is N_eff, not total** (DECISIONS 23(a) + 27(a)) —
  case-control studies pass `--n-fixed <N_eff>` where
  `N_eff = 4·cases·ctrls/(cases+ctrls)`; using total N over-states
  MAGMA's per-SNP precision weight (~27% inflation for de Lange).
- **MHC region excluded** (chr 6: 28,477,797–33,448,354 GRCh37) from the
  scDRS top-1000 gene set and the seismicGWAS Z-score table. Sensitivity
  re-run with MHC retained: `--keep-mhc` flag, applied to Smillie ×
  de Lange only.
- **Autosomes only** (chr 1–22). X chromosome excluded.
- **LD reference**: 1000G EUR for de Lange and Yengo; Liu's ancestry-LD
  choice has three options (1000G EUR vs trans-ancestry vs per-ancestry
  split + meta) — see DECISIONS 27(c). Default-conservative recommendation
  is 1000G EUR with mis-specification documented in Methods.
- **λ_GC** reported per GWAS (`results/magma/<trait>_lambda_gc.tsv`).
  The rule-of-thumb `≤ 1.1` is the WRONG check at large N (polygenic
  signal inflates λ_GC at high power; Bulik-Sullivan 2015) — the
  stratification gate is the **LDSC intercept**, tracked as
  OPEN_FLAGS F10 (pre-narrative, not pre-figure).
- **`--n-min-frac` low-N tail filter** (DECISIONS 26(a)): the munge
  drops SNPs below 0.67 · max(N) when per-SNP N is used (LDSC
  convention). No-op on the `--n-fixed` path.

## One-time setup

Download into `data/reference/`:

- **MAGMA binary** — https://ctg.cncr.nl/software/magma (Linux x86_64).
  Put on `$PATH`.
- **NCBI37.3.gene.loc** — ships with MAGMA. Protein-coding genes with
  Entrez IDs and symbols.
- **g1000_eur** — 1000 Genomes Phase 3 EUR reference panel (`.bed/.bim/.fam`
  from the MAGMA site, ~1 GB).

Download into `data/gwas/` (most are auto-fetched by
`scripts/download_refs.sh`; see DECISIONS 14 + 19):

- **de Lange 2017 UC** — GWAS Catalog **`GCST004133`** (UC-only;
  per DECISIONS 14, **not** `GCST004131` which is IBD-combined).
  Use the harmonized variant
  `harmonised/28067908-GCST004133-EFO_0000729.h.tsv.gz` (`hm_*`
  columns; ~299 MB), NOT the raw deposit
  `uc_build37_45975_20161107.txt.gz` (which ships `MarkerName /
  P.value` with no `chr`/`bp` columns and would need a separate
  rsid → position lookup).
- **Liu 2023 multi-ancestry UC** — GWAS Catalog `GCST90446794` (UC
  arm of *Nat Genet* 55:796). NO per-SNP N column (DECISIONS 14);
  use `--n-fixed 87249` (N_eff per (27)(a), not total 375,508). The
  2.49 GB deposit is gated on the ancestry-LD decision (DECISIONS
  27(c)) — don't download until the call lands.
- **Yengo 2022 height EUR** — GWAS Catalog `GCST90245992` (positive
  control). HAS per-SNP `n` column (range 344 → 1.6M per
  DECISIONS 25(b); `--n-min-frac` filters the low-N tail).
- **Trubetskoy 2022 schizophrenia** — figshare DOI
  `10.6084/m9.figshare.19426775.v7` (DECISIONS 19; no PGC
  registration). Negative control on Smillie at broad tier. **Format
  caveat**: PGC sumstats VCF v1.0 (`##` metadata header + single-`#`
  column header), which the current `prepare_gwas.py` cannot parse
  as-is — a `--comment-char` extension is required (flagged in
  DECISIONS 24(d)).

## Run the pipeline (per GWAS)

### de Lange 2017 UC (primary)

```bash
# 1. Reformat into MAGMA inputs (autosome filter + MAF/INFO QC + lambda_GC).
#    N=36160 is the case-control N_eff (DECISIONS 27(a)), NOT total 45,975.
#    Harmonized file ships hm_* columns matching prepare_gwas defaults.
python code/01_magma/prepare_gwas.py \
    --input data/gwas/uc_delange_GCST004133.h.tsv.gz \
    --out-prefix data/gwas/uc_delange \
    --col-snp hm_rsid --col-chr hm_chrom --col-bp hm_pos \
    --col-p p_value --n-fixed 36160 \
    --lambda-gc-out results/magma/uc_delange_lambda_gc.tsv

# 2. MAGMA annotate + gene-test (10 kb window)
bash code/01_magma/run_magma.sh \
    data/gwas/uc_delange.snp.loc data/gwas/uc_delange.pval \
    data/reference/NCBI37.3.gene.loc data/reference/g1000_eur \
    10 results/magma/uc_delange_10kb

# 3. Sanity-check top genes for UC pattern
python code/01_magma/sanity_check.py \
    --genes-out results/magma/uc_delange_10kb.genes.out \
    --gene-loc data/reference/NCBI37.3.gene.loc \
    --trait-class uc

# 4. Build the scDRS .gs file (top 1000 by Z-score, MHC-excluded)
#    AND the long-format gene-Z table for seismicGWAS
python code/01_magma/make_scdrs_gs.py \
    --genes-out results/magma/uc_delange_10kb.genes.out \
    --gene-loc data/reference/NCBI37.3.gene.loc \
    --trait UC --top-n 1000 \
    --out data/gwas/uc_delange_top1000.gs \
    --out-zscore-table data/gwas/uc_delange_gene_z.tsv

# 5. MHC-included sensitivity gene set (for Smillie x de Lange supplementary run)
python code/01_magma/make_scdrs_gs.py \
    --genes-out results/magma/uc_delange_10kb.genes.out \
    --gene-loc data/reference/NCBI37.3.gene.loc \
    --trait UC --top-n 1000 \
    --out data/gwas/uc_delange_top1000_with_mhc.gs \
    --keep-mhc
```

### Yengo 2022 height EUR (positive control)

```bash
# Yengo ships a per-SNP n column; pass --col-n n and let the LDSC-style
# --n-min-frac filter drop the low-N tail (DECISIONS 26(a)). Range
# 344 to 1.6M reflects rare variants tested in small sub-cohorts —
# real data shape, not a munge bug.
python code/01_magma/prepare_gwas.py \
    --input data/gwas/yengo_height_GCST90245992.tsv \
    --out-prefix data/gwas/yengo_height \
    --col-snp variant_id --col-chr chromosome --col-bp base_pair_location \
    --col-p p_value --col-n n --col-frq effect_allele_frequency \
    --lambda-gc-out results/magma/yengo_height_lambda_gc.tsv

bash code/01_magma/run_magma.sh \
    data/gwas/yengo_height.snp.loc data/gwas/yengo_height.pval \
    data/reference/NCBI37.3.gene.loc data/reference/g1000_eur \
    10 results/magma/yengo_height_10kb

# Positive control: top genes should be height-associated (HMGA1,
# HMGA2, GDF5, ACAN, ADAMTSL3, etc.) — sanity check uses the height
# trait class (add to sanity_check.py if not present).

python code/01_magma/make_scdrs_gs.py \
    --genes-out results/magma/yengo_height_10kb.genes.out \
    --gene-loc data/reference/NCBI37.3.gene.loc \
    --trait HEIGHT --top-n 1000 \
    --out data/gwas/yengo_height_top1000.gs
```

### Liu 2023 UC arm (cross-GWAS sensitivity)

```bash
# Gated on the ancestry-LD decision (DECISIONS 27(c)). Defaults below
# assume Option 1 (1000G EUR LD) — the default-conservative
# recommendation. If Option 2 (combined trans-ancestry panel) or
# Option 3 (per-ancestry split + meta) lands, the run_magma.sh bfile
# argument changes, not the prepare_gwas invocation.
#
# Schema per scripts/download_refs.sh (DECISIONS 14):
#   chromosome, base_pair_location, effect_allele, other_allele,
#   beta, standard_error, effect_allele_frequency, p_value, variant_id
# NO per-SNP N column; --n-fixed = N_eff = 87,249 (precise:
#   4 * 23252 * 352256 / 375508 = 87,248.81), per DECISIONS 27(a).
python code/01_magma/prepare_gwas.py \
    --input data/gwas/uc_liu_GCST90446794.tsv \
    --out-prefix data/gwas/uc_liu \
    --col-snp variant_id --col-chr chromosome --col-bp base_pair_location \
    --col-p p_value --n-fixed 87249 \
    --col-frq effect_allele_frequency \
    --lambda-gc-out results/magma/uc_liu_lambda_gc.tsv

bash code/01_magma/run_magma.sh \
    data/gwas/uc_liu.snp.loc data/gwas/uc_liu.pval \
    data/reference/NCBI37.3.gene.loc data/reference/g1000_eur \
    10 results/magma/uc_liu_10kb

python code/01_magma/sanity_check.py \
    --genes-out results/magma/uc_liu_10kb.genes.out \
    --gene-loc data/reference/NCBI37.3.gene.loc \
    --trait-class uc

python code/01_magma/make_scdrs_gs.py \
    --genes-out results/magma/uc_liu_10kb.genes.out \
    --gene-loc data/reference/NCBI37.3.gene.loc \
    --trait UC --top-n 1000 \
    --out data/gwas/uc_liu_top1000.gs \
    --out-zscore-table data/gwas/uc_liu_gene_z.tsv
```

### Trubetskoy 2022 schizophrenia (negative control)

**Blocked on `prepare_gwas.py` extension** (DECISIONS 24(d)). The
figshare deposit is PGC sumstats VCF v1.0:

- `##`-prefixed metadata header lines at the top of the file
- Single-`#`-prefixed column header line
- Tab-separated data rows with 14 columns including `n_eff_per_snp`

`prepare_gwas.py` calls `pd.read_csv` without `comment=` or
`skiprows=` handling, so the file does not parse as-is. The patch is
a small CLI flag addition (e.g., `--comment-char "##"` mapping
appropriately) plus single-`#` header handling. Not implemented yet.

Once the script supports the format, the invocation will be:

```bash
# Use n_eff_per_snp (last column) — already trial-effective N, per
# DECISIONS 19 / 21(b). PGC3 EUR cohort: 53,386 cases + 77,258 ctrls.
python code/01_magma/prepare_gwas.py \
    --input data/gwas/scz_trubetskoy_eur_PGC3_v3.vcf.tsv.gz \
    --out-prefix data/gwas/scz_trubetskoy \
    --comment-char "##" \
    --col-snp rsid --col-chr chr --col-bp pos \
    --col-p p --col-n n_eff_per_snp \
    --lambda-gc-out results/magma/scz_trubetskoy_lambda_gc.tsv

bash code/01_magma/run_magma.sh \
    data/gwas/scz_trubetskoy.snp.loc data/gwas/scz_trubetskoy.pval \
    data/reference/NCBI37.3.gene.loc data/reference/g1000_eur \
    10 results/magma/scz_trubetskoy_10kb

# Sanity check uses the schizophrenia trait class
python code/01_magma/sanity_check.py \
    --genes-out results/magma/scz_trubetskoy_10kb.genes.out \
    --gene-loc data/reference/NCBI37.3.gene.loc \
    --trait-class schizophrenia

# scDRS .gs for the negative control (MHC-excluded per locked policy)
python code/01_magma/make_scdrs_gs.py \
    --genes-out results/magma/scz_trubetskoy_10kb.genes.out \
    --gene-loc data/reference/NCBI37.3.gene.loc \
    --trait SCZ --top-n 1000 \
    --out data/gwas/scz_trubetskoy_top1000.gs
```

## Notes on column names

The `--col-*` defaults match the GWAS Catalog **harmonised** format
(`hm_rsid`, `hm_chrom`, `hm_pos`, `p_value`) — verified to drop in
against the de Lange harmonised deposit. Non-harmonised deposits use
different names (`chromosome` / `base_pair_location` / `variant_id`
for Liu and Yengo) — open the file first and pass the right column
names. If the file has per-SNP N, use `--col-n` (the LDSC-style
`--n-min-frac 0.67` filter applies); otherwise pass
`--n-fixed <N_eff>`.

## What "looks right"

After step 3, the top 20 should include several of: IL23R, JAK2, TYK2,
STAT3, SMAD3, IL12B, NKX2-3, IRGM, CARD9. NOD2 will be lower for UC than
for Crohn's — expected. If the top list is unrecognizable, suspect a build
mismatch (summary stats must be GRCh37 to match `NCBI37.3.gene.loc` and
`g1000_eur`).

For Liu 2023, the top genes should be broadly similar to de Lange. If they
are not, you may have grabbed the CD or IBD-combined arm by mistake (the
correct UC arm is `GCST90446794`, NOT IBD-combined).

For Yengo height (positive control), the top genes should be growth /
chondrogenesis-associated (HMGA1, HMGA2, GDF5, ACAN, ADAMTSL3, LCORL,
LIN28B). λ_GC will be very high (~5) — that's polygenic signal at
N=1.6M, not stratification; see DECISIONS 25(b)–(c).

For schizophrenia, the top genes should be brain-related (CACNA1C, GRIN2A,
DRD2, etc.).

## What's been run (state as of DECISIONS 27)

| GWAS | Download | `prepare_gwas.py` | `run_magma.sh` | `.gs` |
|---|---|---|---|---|
| de Lange UC (primary) | ✓ laptop | ✓ N_eff=36,160 | — HB | — HB |
| Yengo height (pos ctrl) | ✓ laptop | ✓ filtered tail | — HB | — HB |
| Liu UC (cross-GWAS) | gated on (27)(c) | gated | gated | gated |
| Trubetskoy SCZ (neg ctrl) | auto-fetchable | blocked on script ext (24)(d) | blocked | blocked |

Intermediates at `data/gwas/{uc_delange,yengo_height}.{snp.loc,pval}`
are ready to ship to Hummingbird for the `run_magma.sh` step.
