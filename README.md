# Cross-Atlas Reproducibility of GWAS-Driven Cell-Type Prioritization in Ulcerative Colitis

A 12-month, two-author dry-lab project at UC Santa Cruz, 2026. v1 of a
two-paper research line (v2: UC + CD + RA cross-tissue generalization,
2027–2028).

## Overview

We compare two operational regimes for GWAS-driven cell-type prioritization
on UC colon mucosa:

- **Regime 1 (single-atlas):** run scDRS / seismicGWAS on whichever atlas
  you have.
- **Regime 2 (per-atlas-then-meta-analyze):** run each method per atlas,
  then combine cell-type p-values across atlases via Brown's method
  (correlation-aware Fisher's combination).

We compare these regimes side-by-side on three independent UC colon
atlases (Smillie 2019, Kong 2023, Mennillo 2024) against two UC GWAS
(de Lange 2017, Liu 2023 multi-ancestry), using two methods (scDRS,
seismicGWAS) at two cell-type granularity tiers (broad ~10–15 cell types;
fine ~30–50 cell states). We also compare to two broad multi-tissue gut
atlases — HCA Gut Cell Atlas (external single-atlas reference; nested
within Pan-GI as ~25% of Extended+, per DECISIONS.md correction 3/7) and
Pan-GI (integration-pipeline-robustness comparator with known donor
overlap).

The full plan is in [`docs/uc-cross-atlas-v1-plan.pdf`](docs/uc-cross-atlas-v1-plan.pdf);
[`PLAN.md`](PLAN.md) is a tight operational summary.
[`DECISIONS.md`](DECISIONS.md) is the append-only log of locked
methodological commitments.

## Authors

This is a two-author project. Division of labor is documented below; the
manuscript drafter is first author (or co-first if contributions are
genuinely indistinguishable). Revisit at month 4.

- Author A: `[name, email, role]`
- Author B: `[name, email, role]`

PI mentor: `[name]`. PI sign-off on locked-core scope, stretch ladder,
and v2 trajectory at M0.

### Division of labor (sketch — finalize in M0)

- **Joint:** cell-type harmonization, regime 2 implementation, manuscript.
- **Method-axis split (natural with two methods in core):**
  - scDRS pipeline (compute-score, group analysis, null aggregation,
    test-retest, donor-LOO): `[author]`
  - seismicGWAS pipeline (specificity, regressions, gene-Z permutations,
    test-retest, donor-LOO): `[author]`
- **Other:** MAGMA pipeline (joint), broad-atlas comparators (joint),
  figures (joint).

## Status

`[M0 / M1 / ...]`. See [`PLAN.md`](PLAN.md) §"Three explicit fork points"
for descope decisions.

## Repository structure

```
.
├── DECISIONS.md                # locked methodological commitments (append-only)
├── PLAN.md                     # operational summary
├── README.md                   # this file
├── docs/
│   ├── uc-cross-atlas-v1-plan.pdf
│   └── uc-cross-atlas-pre-project-curriculum.pdf
├── data/
│   ├── gwas/                   # de Lange + Liu 2023 + Trubetskoy
│   ├── atlases/                # Smillie, Kong-UC, Mennillo, HCA Gut, Pan-GI
│   │   └── donor_metadata/     # per-atlas donor_metadata.csv (v2-setup pattern)
│   └── reference/              # 1000G EUR, CL ontology, HGNC pin
├── code/
│   ├── 01_magma/               # ×2 GWAS + negative control (MHC-excluded, autosomes only)
│   ├── 02_atlas_prep/          # download, harmonize, HGNC remap, donor-attribution
│   ├── 03_scdrs/               # regime 1 (×2 GWAS, ×2 tiers)
│   ├── 04_seismic/             # regime 1 (×2 GWAS, ×2 tiers)
│   │   └── specificity_long/   # serialized specificity for v2 reuse
│   ├── 05_magma_geneprop/      # Watanabe-style sanity track
│   ├── 06_concordance/         # 5 axes; metrics library + cross-atlas driver
│   ├── 07_regime2_meta/        # Brown's method (analytical, de Lange only)
│   ├── 08_cross_method/        # scDRS vs seismicGWAS within-atlas
│   ├── 09_cross_gwas/          # de Lange vs Liu within-atlas
│   ├── 10_broad_atlas_hca/     # HCA Gut comparator (independent)
│   ├── 11_broad_atlas_pangi/   # Pan-GI comparator (with/without donor-overlap)
│   ├── 12_positive_control/    # Tabula Muris × height + schizophrenia negative
│   └── 13_figures/
├── results/                    # generated outputs (not tracked)
├── notebooks/
└── manuscript/
```

Each `code/<step>/` has its own README with run instructions.

## Methods locked in v1

- **scDRS** (Zhang 2022) and **seismicGWAS** (Lai 2025).
- **scPagwas is not in v1** (compute infeasible).
- **Regime 3** (scANVI integration) is deferred to v2.

## Reproducing the analysis

Prerequisites:

- Python 3.10 (pinned for scDRS 1.0.4 compatibility).
- R 4.3+ for seismicGWAS (`devtools::install_github("ylaboratory/seismicGWAS")`)
  and Brown's method (`EmpiricalBrownsMethod::kostsMethod()`).
- MAGMA precompiled binary from CTGLab.
- plink 1.9.

```bash
# Python deps
pip install -e .[atlas,dev]

# R deps — see code/04_seismic/README.md
```

Per-step instructions live in each `code/<step>/README.md`.

## Reading list

See PLAN.md §"Reading list" or PDF Part 6. Must-read in M1:

- Watanabe 2019, Skene 2018, Bryois 2020, Zhang 2022 (scDRS),
  Lai 2025 (seismicGWAS), Lakkis 2024, Smillie 2019, de Lange 2017,
  Liu 2023.

## Contact

`[Author A email]` / `[Author B email]`
