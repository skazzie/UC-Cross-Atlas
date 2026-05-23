# 12_positive_control — sanity scaffolding GWAS × atlas controls

Locked-core sanity controls from PLAN.md §"Sanity scaffolding":

## scDRS positive control: Tabula Muris × Yengo 2022 height GWAS

Confirms scDRS *can* detect signal (modulo mouse-human ortholog mapping
caveat). 1 compute-score run + group analysis.

Expected: cell types involved in skeletal/musculoskeletal development
(chondrocytes, osteoblasts) light up with FDR < 0.05.

## scDRS negative control: Trubetskoy 2022 schizophrenia × Smillie at broad tier

Confirms scDRS *doesn't manufacture* signal. 1 compute-score run + group
analysis.

Expected: no colon cell type achieves FDR < 0.05 for schizophrenia. If colon
cell types do enrich, the pipeline is producing spurious signal — debug
before proceeding.

**Caveat on negative-control purity:** schizophrenia hits include
MHC-region signal and complement-pathway genes, which can theoretically
enrich for antigen-presenting / immune cell types in non-brain tissues.
Mitigated by (a) excluding MHC genes from all scDRS gene sets including
schizophrenia (per the MHC policy in DECISIONS.md) and (b) framing
"absence of enrichment" rather than "exact null" as the operational
expectation. A perfectly clean negative control (e.g., LDL on
non-hepatocyte cells) is deferred to revision if reviewers ask.

## scDRS Smillie × de Lange MHC-included sensitivity

One supplementary scDRS run with MHC genes retained, to confirm MHC
inclusion shifts antigen-presenting cell rankings as expected. Compares
against the MHC-excluded headline.

## Stretch #3: LDL × Tabula Sapiens hepatocytes

Conditional, only if locked core is clean and stretches are green-lit at
M6. Closes the human-only-code-path gap. ~1 day.

## Output

- `results/sanity/tabula_muris_height/` — positive control
- `results/sanity/smillie_schizophrenia/` — negative control
- `results/sanity/smillie_delange_mhc_included/` — MHC sensitivity
- `results/sanity/ts_hepatocytes_ldl/` — stretch #3 (if activated)

## Driver script

`code/12_positive_control/run_controls.py`. One control type per
invocation. CLI examples:

```bash
# Positive: Tabula Muris x height
python code/12_positive_control/run_controls.py \
    --control positive \
    --atlas tabula_muris \
    --h5ad-path data/atlases/tabula_muris.h5ad \
    --gwas-gene-set results/magma/height_yengo_top1k.gs \
    --out-dir results/sanity/tabula_muris_height/

# Negative: Smillie x schizophrenia
python code/12_positive_control/run_controls.py \
    --control negative \
    --atlas smillie \
    --h5ad-path data/atlases/smillie.h5ad \
    --gwas-gene-set results/magma/scz_trubetskoy_top1k.gs \
    --out-dir results/sanity/smillie_schizophrenia/

# MHC sensitivity: Smillie x de Lange (MHC retained)
python code/12_positive_control/run_controls.py \
    --control mhc-sensitivity \
    --atlas smillie \
    --h5ad-path data/atlases/smillie.h5ad \
    --gwas-gene-set results/magma/delange_top1k_mhc_included.gs \
    --baseline-group-file results/scdrs/smillie_delange/cell_type_broad/smillie_delange.scdrs_group \
    --out-dir results/sanity/smillie_delange_mhc_included/
```

Warn-don't-fail semantics for the assertions, matching the README's
sanity-track scope (the user decides whether to investigate). Runs on
login node when the atlas is small; submit to SLURM otherwise.
