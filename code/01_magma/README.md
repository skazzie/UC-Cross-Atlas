# 01_magma — GWAS preprocessing → gene scores → scDRS gene set

Pipeline for converting the de Lange 2017 UC summary statistics into the
top-1000 MAGMA gene set that scDRS, scPagwas, and seismicGWAS all consume.

## One-time setup

Download into `data/reference/`:

- **MAGMA binary** — https://ctg.cncr.nl/software/magma (Linux x86_64). Put on `$PATH`.
- **NCBI37.3.gene.loc** — ships with MAGMA. Protein-coding genes with Entrez IDs and symbols.
- **g1000_eur** — 1000 Genomes Phase 3 EUR reference panel (`.bed/.bim/.fam` from the MAGMA site, ~1 GB).

Download into `data/gwas/`:

- **de Lange 2017 UC summary stats** — GWAS Catalog `GCST004131` (UC-only). Harmonized `.tsv.gz`.

## Run the pipeline

```bash
# 1. Reformat GWAS into MAGMA inputs (drops MAF<0.01, INFO<0.6 if columns present)
python code/01_magma/prepare_gwas.py \
    --input data/gwas/GCST004131_harmonised.tsv.gz \
    --out-prefix data/gwas/uc \
    --col-snp hm_rsid --col-chr hm_chrom --col-bp hm_pos \
    --col-p p_value --n-fixed 45975

# 2. MAGMA annotate + gene-test (10 kb window, primary)
bash code/01_magma/run_magma.sh \
    data/gwas/uc.snp.loc data/gwas/uc.pval \
    data/reference/NCBI37.3.gene.loc data/reference/g1000_eur \
    10 results/magma/uc_10kb

# 3. Sanity-check top genes (must include some of IL23R/JAK2/TYK2/STAT3/...)
python code/01_magma/sanity_check.py \
    --genes-out results/magma/uc_10kb.genes.out \
    --gene-loc data/reference/NCBI37.3.gene.loc

# 4. Build the scDRS .gs file (top 1000 by Z-score)
python code/01_magma/make_scdrs_gs.py \
    --genes-out results/magma/uc_10kb.genes.out \
    --gene-loc data/reference/NCBI37.3.gene.loc \
    --trait UC --top-n 1000 \
    --out data/gwas/uc_top1000.gs

# 5. Week 7.5 sensitivity: same pipeline with 50 kb window
bash code/01_magma/run_magma.sh \
    data/gwas/uc.snp.loc data/gwas/uc.pval \
    data/reference/NCBI37.3.gene.loc data/reference/g1000_eur \
    50 results/magma/uc_50kb
python code/01_magma/make_scdrs_gs.py \
    --genes-out results/magma/uc_50kb.genes.out \
    --gene-loc data/reference/NCBI37.3.gene.loc \
    --trait UC --top-n 1000 \
    --out data/gwas/uc_top1000_50kb.gs
```

## Notes on column names

The `--col-*` defaults match the GWAS Catalog harmonized format (`hm_rsid`,
`hm_chrom`, `hm_pos`, `p_value`). The IIBDGC portal version uses different
names — open the file first and pass the right column names. If the file has
per-SNP N, use `--col-n`; if not, `--n-fixed 45975` (≈ 12,366 cases + 33,609
controls for de Lange UC).

## What "looks right"

After step 3, the top 20 should include several of: IL23R, JAK2, TYK2, STAT3,
SMAD3, IL12B, NKX2-3, IRGM, CARD9. NOD2 will be lower for UC than for Crohn's;
that is expected. If the top list is unrecognizable, suspect a build mismatch
(summary stats must be GRCh37 to match `NCBI37.3.gene.loc` and `g1000_eur`).
