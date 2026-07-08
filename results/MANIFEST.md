# results/MANIFEST.md — analysis output ledger

Ground-truth inventory of which analyses have SAVED outputs. `results/` is
gitignored, so this manifest is the tracked record of what exists and where.
Keep it current after every run. Last updated: 2026-07-07.

Legend: [x] file confirmed | [~] ran but output not saved / must redo |
[ ] not run

Machines: LAPTOP = `C:\Users\muska\UC-Cross-Atlas`;
VM = `ucca-compute` GCP (CANONICAL going forward); HB = Hummingbird (NO results).

---

## MAGMA gene-level (`results/magma/`)

| GWAS | .genes.out | .gs built | where |
|------|-----------|-----------|-------|
| de Lange (GCST004133 harmonised) | [x] | [x] std + [x] with_mhc | LAPTOP + VM |
| Liu (GCST90446794 harmonised)    | [x] | [x] std + [x] with_mhc | VM (2026-07-07) |
| SCZ (neg control)                | [x] | [x] | (June) |
| Yengo height (pos control)       | [x] | [x] | (June) |

Note: de Lange used hm_rsid/hm_chrom/hm_pos; Liu used rsid/chromosome/
base_pair_location (different harmonised schema — per-trait COL_* vars added
to 01_magma.slurm). Liu sanity passed only at top-n=100 (HLA-saturated).

---

## scDRS (`results/scdrs/<atlas>_<gwas>[...]/`)

### Garrido-Trigo — DONE (scDRS)
| run | broad | fine | cov set | machine |
|-----|-------|------|---------|---------|
| garrido x delange (headline seed42) | [x] | [x] | (June; cov set unverified this session) | LAPTOP |
| garrido x delange MHC-included      | [x] | [x] | (June) | LAPTOP |
| garrido x delange no-cov (diag)     | [x] | [x] | none | LAPTOP |
| garrido x height (pos ctrl)         | [x] | [ ] | (June) | LAPTOP |
| garrido x SCZ (neg ctrl)            | [x] | [ ] | (June) | LAPTOP |
| garrido x delange MHC-excl (smoketest) | [x] | [ ] | log_n_genes,log_n_counts,donor | VM 07-07 |
| garrido x delange MHC-incl          | [x] | [ ] | same | VM 07-07 |
| garrido x Liu MHC-excl              | [x] | [ ] | same | VM 07-07 |
| garrido x Liu MHC-incl              | [x] | [ ] | same | VM 07-07 |

CAVEAT: VM Garrido runs use depth+donor cov only (sex unavailable). NOT the
locked headline cov config. Laptop seed42 cov set not re-verified this session.

### Other atlases — NOT RUN
| atlas | scDRS |
|-------|-------|
| Smillie (PRIMARY, account-gated) | [ ] |
| TAURUS   | [ ] |
| HCA Gut  | [ ] |
| Pan-GI   | [ ] |

---

## seismicGWAS (`results/seismic/`)

| atlas | seismic | status |
|-------|---------|--------|
| Garrido | [~] | RAN June (cross-method rho 0.60-0.79 reported in chat) but NO OUTPUT FILE SAVED. MUST RE-RUN with outputs saved. |
| Smillie / TAURUS / HCA / Pan-GI | [ ] | not run |

Prereq: confirm `seismicGWAS` R package installed on VM; need `<atlas>_sce.rds`
inputs.

---

## Downstream analyses — NONE RUN

| step | dir | status |
|------|-----|--------|
| Concordance (cross-atlas)   | 06_concordance      | [ ] |
| Cross-method (scDRS vs seismic) | 08_cross_method | [ ] |
| Cross-GWAS (de Lange vs Liu)| 09_cross_gwas       | [ ] |
| Broad-atlas HCA             | 10_broad_atlas_hca  | [ ] |
| Broad-atlas Pan-GI          | 11_broad_atlas_pangi| [ ] |
| Test-retest (seeds 1,2,3)   | —                   | [ ] |
| Regime-2 / Brown's          | 07_regime2_meta     | [ ] null-draw feather not generated |
| MAGMA gene-property         | 05_magma_geneprop   | [ ] not confirmed |

---

## Headline confirmed values (Garrido, cell_type_broad, VM 2026-07-07)

de Lange and Liu, {MHC-excluded, MHC-included}, assoc_mcz per cell type:

| cell type            | dL excl | dL incl | Liu excl | Liu incl |
|----------------------|---------|---------|----------|----------|
| B cell               | -0.30   | +7.29   | +0.59    | +7.44    |
| dendritic cell       | -0.31   | +6.30   | +0.75    | +6.32    |
| monocyte/macrophage  | +0.76   | +4.98   | +0.63    | +4.28    |
| epithelial progenitor| +3.26   | +1.55   | +3.92    | +2.17    |
| colonocyte           | +2.62   | +2.62   | +3.21    | +1.89    |
| T cell               | -2.64   | -2.64   | -2.18    | -2.39    |

(Positive mcz = enriched; these are depth+donor-corrected, sex-unavailable —
NOT the locked headline cov config. Treat as validation, not final numbers.)
