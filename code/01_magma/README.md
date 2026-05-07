# 01_magma — GWAS preprocessing → gene scores → scDRS / seismicGWAS inputs

Pipeline for converting GWAS summary statistics into:

- a top-1000 MAGMA gene set (`.gs` file) for scDRS
- a long-format gene-Z-score table for seismicGWAS

Locked v1 commitments (see DECISIONS.md):

- **Three GWAS**: de Lange 2017 UC (primary), Liu 2023 multi-ancestry UC arm
  (cross-GWAS sensitivity), Trubetskoy 2022 schizophrenia (negative control).
- **MHC region excluded** (chr 6: 28,477,797–33,448,354 GRCh37) from the
  scDRS top-1000 gene set and the seismicGWAS Z-score table. Sensitivity
  re-run with MHC retained: `--keep-mhc` flag, applied to Smillie ×
  de Lange only.
- **Autosomes only** (chr 1–22). X chromosome excluded.
- **LD reference**: 1000G EUR for both UC GWAS (approximate for multi-
  ancestry Liu 2023 — documented in Methods).
- **λ_GC ≤ 1.1** verified in M1 for both GWAS, written to
  `results/magma/lambda_gc.tsv`.

## One-time setup

Download into `data/reference/`:

- **MAGMA binary** — https://ctg.cncr.nl/software/magma (Linux x86_64).
  Put on `$PATH`.
- **NCBI37.3.gene.loc** — ships with MAGMA. Protein-coding genes with
  Entrez IDs and symbols.
- **g1000_eur** — 1000 Genomes Phase 3 EUR reference panel (`.bed/.bim/.fam`
  from the MAGMA site, ~1 GB).

Download into `data/gwas/`:

- **de Lange 2017 UC** — GWAS Catalog `GCST004131` (UC-only) harmonized
  `.tsv.gz`.
- **Liu 2023 multi-ancestry IBD** — *Nat Genet* 55:796–806. UC arm only.
  Verify per-SNP N column at M1; if absent, document fixed-N approximation.
- **Trubetskoy 2022 schizophrenia** — *Nature* 604:502–508. Used as
  negative control on Smillie at broad tier.

## Run the pipeline (per GWAS)

### de Lange 2017 UC (primary)

```bash
# 1. Reformat into MAGMA inputs (autosome filter + MAF/INFO QC + lambda_GC)
python code/01_magma/prepare_gwas.py \
    --input data/gwas/GCST004131_harmonised.tsv.gz \
    --out-prefix data/gwas/uc_delange \
    --col-snp hm_rsid --col-chr hm_chrom --col-bp hm_pos \
    --col-p p_value --n-fixed 45975 \
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

### Liu 2023 UC arm (cross-GWAS sensitivity)

```bash
python code/01_magma/prepare_gwas.py \
    --input data/gwas/<liu_uc_file>.tsv.gz \
    --out-prefix data/gwas/uc_liu \
    --col-snp <rsid_col> --col-chr <chr_col> --col-bp <bp_col> \
    --col-p <p_col> --col-n <n_col> \
    --lambda-gc-out results/magma/uc_liu_lambda_gc.tsv

# (If Liu lacks per-SNP N, swap --col-n for --n-fixed <combined_N> and
#  document fixed-N approximation in DECISIONS.md.)

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

```bash
python code/01_magma/prepare_gwas.py \
    --input data/gwas/<scz_file>.tsv.gz \
    --out-prefix data/gwas/scz_trubetskoy \
    --col-snp <rsid_col> --col-chr <chr_col> --col-bp <bp_col> \
    --col-p <p_col> --n-fixed <combined_N>

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

The `--col-*` defaults match the GWAS Catalog harmonized format
(`hm_rsid`, `hm_chrom`, `hm_pos`, `p_value`). The IIBDGC portal version
uses different names — open the file first and pass the right column
names. If the file has per-SNP N, use `--col-n`; otherwise pass
`--n-fixed`.

## What "looks right"

After step 3, the top 20 should include several of: IL23R, JAK2, TYK2,
STAT3, SMAD3, IL12B, NKX2-3, IRGM, CARD9. NOD2 will be lower for UC than
for Crohn's — expected. If the top list is unrecognizable, suspect a build
mismatch (summary stats must be GRCh37 to match `NCBI37.3.gene.loc` and
`g1000_eur`).

For Liu 2023, the top genes should be broadly similar to de Lange. If they
are not, you may have grabbed the CD or IBD-combined arm by mistake.

For schizophrenia, the top genes should be brain-related (CACNA1C, GRIN2A,
DRD2, etc.).
