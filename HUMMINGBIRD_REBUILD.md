# UC-Cross-Atlas: Hummingbird Rebuild Handoff

_2026-10-03 · @So_

The GCP VM was deleted for non-payment and all results were lost; the
full repo survives on GitHub at `db3148c`. This is the complete
instruction set for rebuilding on UCSC Hummingbird, including every bug
fixed during the GCP sessions so none has to be re-discovered.

---

## What survived, what was lost

All code survived. The repo is intact on GitHub through commit `db3148c`,
including every loader, driver, analysis script, `DECISIONS.md` with the
canonical atlas-roster header, and `CONTINUITY.md`. The hard debugging
work is preserved.

| Asset | State | Cost to rebuild |
| --- | --- | --- |
| Repo code (all loaders, drivers, SLURM templates) | Safe on GitHub at `db3148c` | None — `git clone` |
| `DECISIONS.md`, `CONTINUITY.md` | Safe on GitHub | None |
| conda env (scdrs, scanpy, R, seismicGWAS) | Lost | 1–2 h, mostly unattended |
| MAGMA outputs + `*_gene_z.tsv` + `.gs` gene sets | Lost | ~1 h, gates everything |
| Raw atlas downloads (~24 GB across 4 atlases) | Lost | ~1 h of network |
| Processed `.h5ad` + covariate files | Lost | ~1–2 h of loader runs |
| All scDRS group analyses, seismic TSVs, concordance CSV | Lost | ~25 h of compute, unattended |

The 175 MB results archive was built on the VM but never downloaded
before the project was deleted.

---

## Hummingbird environment — verified facts

All of the following was checked directly on the cluster on 2026-10-03.
Do not assume anything beyond it.

| Fact | Value |
| --- | --- |
| Login | `ssh mukhinda@hb.ucsc.edu` — requires UCSC VPN (`vpn.ucsc.edu`, Cisco AnyConnect + Duo) from off campus |
| Repo location | `/hb/scratch/mukhinda/UC-Cross-Atlas` (already cloned) |
| Data location | `/hb/scratch/mukhinda/...` — scratch, not home |
| SLURM account | `128x24` (confirmed via `sacctmgr`; account and partition genuinely share the name) |
| Partitions | `128x24` (120 GB, 24 cores), `256x44` (247 GB, 44 cores), `96x24gpu4` |
| Time limits | infinite on all three partitions |
| Core cap | 72 cores simultaneous (open access) |
| Modules present | `miniconda3/3.13`, `git/2.52.0`, `plink`, `gatk`, `cellranger` |
| Modules absent | **No R. No MAGMA.** Both must be installed into the conda env or fetched as binaries |
| Login-node network | Working (`curl` to Zenodo returns 200), so downloads happen on the head node |
| Data transfer host | `hbfeeder.ucsc.edu` for `sftp`/`scp`, to keep load off the head node |

- `git` is not on the default PATH — `module load git/2.52.0` first, in
  both interactive sessions and job scripts.
- Hummingbird does not back up user data. Build an `rsync` of `results/`
  back to the laptop into the workflow from day one; the whole results
  tree was only 175 MB.

---

## Non-negotiable lessons

These each cost hours or days on GCP and are already fixed in the
committed code. Preserve every one; do not "simplify" them away.

1. **Covariate files must contain pre-expanded float dummies, never a
   raw categorical column.** Handing scDRS a categorical `sample` column
   makes its `category2dummy` emit **bool** dummies; `df_cov.values`
   then upcasts to `object` and `np.linalg.solve` dies with
   `Cannot cast dtype('O') to float64`. Cov format is `const`,
   `log_n_genes`, `log_n_counts`, `sample_<id>...` all float64. This
   crashed TAURUS three times.
2. **Drop collinear donor columns.** Where samples nest within donors,
   donor dummies plus sample dummies plus the intercept give a
   rank-deficient design. TAURUS: rank 54 of 77. Keep `sample`, drop
   `donor`, and assert nesting before writing.
3. **Never use disease, health or inflammation status as a covariate.**
   It regresses out the signal being tested.
4. **`--flag_filter_data=False` everywhere.** Loaders apply the
   scDRS-replica pre-filter (`min_genes=250`, `min_cells=50`)
   themselves, then bake in `sc.pp.pca(n_comps=20)` +
   `sc.pp.neighbors(n_neighbors=15, n_pcs=20)` and stamp
   `uns['scdrs_prefilter']`. With `True`, scDRS's load chain does four
   in-place subset copies that **drop `obsp`**, forcing a ~45-minute
   kNN rebuild per run and causing three OOM crashes.
5. **The scdrs 1.0.2 CLI is positional with underscore flags**, not the
   dashed keyword form in the old SLURM templates. See Phase 4 for the
   exact invocations.
6. **seismicGWAS runs through the flat-export bypass, not
   `zellkonverter`.** `readH5AD` pulls in basilisk, which tries to
   compile Python 3.14 from source and fails without root; basilisk
   ignores `RETICULATE_PYTHON`. The working chain is
   `h5ad_to_sce_export.py` → `sce_from_export.R` (pure
   `Matrix::readMM`) → `run_seismic.R`. This was blocked from June
   until August.
7. **seismic API specifics**: `get_ct_trait_associations()` takes **no**
   `confounders` argument (it corrects gene length and LD internally);
   pass `magma_gene_col="SYMBOL"`, `magma_z_col="ZSTAT"`; it returns
   only `cell_type, pvalue, FDR`; it needs a `logcounts` assay, not
   just `counts`; and no feather cache — the round-trip loses gene
   rownames and yields "1 overlapping gene".
8. **Coerce object-dtype obs columns to str before `write_h5ad`.** h5py
   cannot serialise mixed/NaN object columns; this killed the first
   Pan-GI load on `obs['batch']` after all the real work had succeeded.
9. **Always pass explicit `--out-h5ad` / `--out-cov`.** The drivers'
   `../../` defaults resolve outside the repo and fail at the write
   step, after the full load has run.
10. **`git stash && git pull --rebase && git stash drop` at the start
    of every session.** Unstaged changes silently block pulls, so the
    VM ran stale code at least four times.

The MAGMA `.genes.out` files are keyed on **Entrez IDs**, not symbols.
seismic consumes the symbol-keyed `scratch/data/gwas/<gwas>_gene_z.tsv`
written by `make_scdrs_gs.py`. Passing `.genes.out` directly produces
the "1 overlapping gene" failure.

---

## Phase 1 — environment bootstrap

This is the only phase that needs real attention; everything after it
is mechanical. Run it on the login node.

Use the cluster's `miniconda3/3.13` module rather than installing
miniforge, unless it proves unusable. The GCP env was Python 3.10 with
`scdrs 1.0.2`, and that version pairing matters — `scdrs 1.0.4` is
GitHub-only and the CLI differs.

```
module load miniconda3/3.13 git/2.52.0
cd /hb/scratch/mukhinda/UC-Cross-Atlas
conda create -y -n uc-cross-atlas python=3.10
conda activate uc-cross-atlas
mamba install -y -c conda-forge scanpy anndata numpy'<2.0' pandas scipy h5py pyarrow
python -m pip install scdrs==1.0.2
```

Keep `numpy<2.0`. A 1.26→2.0 upgrade caused an ABI break on a previous
machine.

The R side is the hard part:

```
mamba install -y -c conda-forge r-base r-devtools r-optparse r-matrix \
  bioconductor-singlecellexperiment bioconductor-summarizedexperiment
R -e 'devtools::install_github("ylaboratory/seismic")'
```

Note the repo is `seismic` but the installed package is `seismicGWAS`.
Do **not** install `zellkonverter` — see lesson 6.

Verify before moving on:

```
which scdrs && scdrs --help | head -5
R -e 'library(seismicGWAS); packageVersion("seismicGWAS")'
R -e 'library(optparse); library(Matrix); library(SingleCellExperiment)'
```

If any R install fails against cluster toolchains, `hb-team@ucsc.edu`
handles system-wide install requests.

---

## Phase 2 — GWAS and MAGMA

Everything downstream consumes MAGMA output, so this runs before any
atlas work. All of it was lost with the VM.

Needed artifacts, all under `scratch/`:
- `scratch/data/reference/` — MAGMA binary, `g1000_eur` LD panel,
  `NCBI37.3.gene.loc`
- `scratch/results/magma/delange_10kb.genes.out` and `liu_10kb.genes.out`
- `scratch/data/gwas/delange_gene_z.tsv`, `liu_gene_z.tsv` — symbol-keyed,
  `SYMBOL + ZSTAT`, consumed by seismic
- `scratch/data/gwas/delange_top1000.gs` and
  `delange_top1000_with_mhc.gs` — consumed by scDRS

The CTG-NL hosting moved; `scripts/download_refs.sh` already carries the
working replacement URLs on `vu.data.surf.nl`. Run it rather than
hunting links.

GWAS sources: de Lange 2017 is **GCST004133** (not GCST004131 —
corrected in DECISIONS); Liu 2023 needs its UC arm separated from the
CD and IBD-combined arms.

```
bash scripts/download_refs.sh
python code/01_magma/prepare_gwas.py                              # per GWAS
sbatch --export=ALL,GWAS=delange scripts/slurm/01_magma.slurm
sbatch --export=ALL,GWAS=liu     scripts/slurm/01_magma.slurm
python code/01_magma/make_scdrs_gs.py                             # writes .gs + *_gene_z.tsv
```

Sanity-check before proceeding — this passed on GCP and should again:

```
python code/01_magma/sanity_check.py \
  --genes-out scratch/results/magma/liu_10kb.genes.out \
  --gene-loc scratch/data/reference/NCBI37.3.gene.loc \
  --trait-class uc --top-n 100
```

Expected: `PASS`, with `PLCL1` top (Z≈9.81) and `FCGR2A` + `JAK2` among
the expected UC genes in the top 100.

---

## Phase 3 — atlases and loaders

Four atlases were working on GCP; all four loaders now run first-try, so
this is download plus execution. Download on the login node into
`scratch/data/atlases/`.

| Atlas | Source | Size | Loader status |
| --- | --- | --- | --- |
| Smillie 2019 | Single Cell Portal **SCP259** (bulk download, short-lived auth token) | ~5 GB | `run_smillie_load.py` — committed at `db3148c` |
| Garrido-Trigo 2023 | GEO **GSE214695** — `GSE214695_RAW.tar` + `GSE214695_cell_annotation.csv.gz` | ~2 GB | `run_garrido_load.py` |
| TAURUS-IBD | Zenodo **10.5281/zenodo.14007626**, file `TAURUS_raw_counts_annotated_final.h5ad`, **v3 only** | 12.7 GB | `run_taurus_load.py` |
| Pan-GI Extended+ | CELLxGENE `https://datasets.cellxgene.cziscience.com/757945c8-a916-431d-aceb-1afbc80a7c55.h5ad` | 11 GB | `run_pangi_load.py` |

**TAURUS**: use **v3** (posted 2024-10-30). Earlier revisions ship
incorrect donor/sample metadata and fail the loader's 22-donor gate.
md5 is `c1bd13b92cacb164a401c6c4a4e7912c`; verify it.

**Pan-GI**: the download UUID (`757945c8…`) differs from the deposit ID
originally pinned in the loader docstring (`1dcf15ee…`). The download
UUID is the correct one — reconciled 2026-10-03 in `load_pangi.py` and
`DECISIONS.md`.

Loader invocations, all with explicit output paths:

```
conda activate uc-cross-atlas
export RETICULATE_PYTHON="$CONDA_PREFIX/bin/python"

python code/02_atlas_prep/run_taurus_load.py \
  --h5ad scratch/data/atlases/TAURUS_raw_counts_annotated_final.h5ad \
  --out-h5ad scratch/data/atlases/taurus.h5ad \
  --out-cov scratch/data/atlases/taurus_covariates.tsv

python code/02_atlas_prep/run_pangi_load.py \
  --h5ad scratch/data/atlases/pangi_extended_plus.h5ad \
  --out-h5ad scratch/data/atlases/pangi.h5ad \
  --out-cov scratch/data/atlases/pangi_covariates.tsv
```

Run the big loads on `256x44`, or interactively on a compute node —
Pan-GI peaked around 44 GB RSS during scDRS and the backed read of 1.6 M
cells is the memory-heavy step.

**HCA Gut** now has a driver: `run_hca_load.py` (added 2026-10-03,
mirrors `run_pangi_load.py`; cov = const + depth + float dummies for
`assay`, `batch`, `Fraction`, `sex`; donor dropped as collinear with
batch under the nesting gate).

---

## Phase 4 — SLURM script rewrites

The five templates in `scripts/slurm/` were written in May, before
anything had actually run. They encode assumptions the GCP sessions
disproved and would fail on submission. Fix them before submitting
anything.

### `03_scdrs_compute.slurm` — six defects (fixed 2026-10-03)

| Defect | Previous | Correct |
| --- | --- | --- |
| CLI form | `--h5ad-file X --gs-file Y --cov-file Z --out-folder` | positional, underscore flags — see below |
| Filter flag | `--flag-filter-data True` | `False` — `True` destroys the baked kNN graph |
| Tier handling | loops `perform-downstream` once per tier | one call, `--group-analysis=cell_type_broad,cell_type_fine` |
| Output dirs | not created before the run | `mkdir -p` both score and group dirs; scDRS will not create them and dies after computing |
| Memory | `--mem=24G` | `--mem=64G` — Pan-GI peaked ~44 GB RSS |
| Walltime | `--time=04:00:00` | `--time=12:00:00` — TAURUS's 2×2 took ~4.5 h; partitions allow infinite, but shorter requests backfill better |

The invocation that actually worked:

```
scdrs compute-score \
  "$H5AD" human "$GS" human "$OUT/" \
  --cov_file="$COV" \
  --flag_filter_data=False --flag_raw_count=False --n_ctrl=1000

scdrs perform-downstream \
  "$H5AD" "$OUT/UC.full_score.gz" "${OUT}_group/" \
  --group-analysis=cell_type_broad,cell_type_fine \
  --flag-filter-data=False --flag-raw-count=False
```

Note the asymmetry: `compute-score` takes underscores,
`perform-downstream` takes dashes on its flags. That is how scdrs 1.0.2
actually behaves.

### Output directory naming

The GCP directory names were historical accidents
(`garrido_delange_smoketest_group`, `smillie_delange_excl_fullcov_group`,
`taurus_delange_excl_group`). Standardised 2026-10-03 on
`${ATLAS}_${GWAS}_${MHC}_group`; `SCDRS_GROUP_DIRS` in
`code/08_cross_method/run_scdrs_seismic_concordance.py` matches.

### Other fixes

- `--mail-user=amoli@ucsc.edu` → `mukhinda@ucsc.edu` in all five
  templates.
- The header comment claiming covariates include `donor`, `sample`,
  `sex` is stale — donor is dropped as collinear, and sex is absent from
  both SCP259 and Garrido metadata.
- `module load miniconda3 || module load anaconda3` should be
  `module load miniconda3/3.13`, and `git/2.52.0` added where the
  scripts touch git.
- Check `scripts/config.sh` points `UCC_REPO` at
  `/hb/scratch/mukhinda/UC-Cross-Atlas` and `UCC_CONDA_ENV` at
  `uc-cross-atlas`.

---

## Phase 5 — submission plan

The 72-core cap allows several atlases to run at once, which the serial
GCP VM could not. Submit the scDRS jobs together and walk away.

```
# scDRS: 4 atlases x de Lange x {MHC-excluded, MHC-included}
for atlas in smillie garrido_trigo taurus pangi; do
  for suffix in "" "_with_mhc"; do
    sbatch --export=ALL,ATLAS=$atlas,GWAS=delange,GS_SUFFIX=$suffix \
           scripts/slurm/03_scdrs_compute.slurm
  done
done
```

Budget 4–6 h per job. HCA at ~400 k cells is the largest and will want
`256x44`.

Seismic is a three-step chain per atlas and is much cheaper — roughly
1 h for all four runs of one atlas:

```
python code/04_seismic/h5ad_to_sce_export.py \
  --h5ad scratch/data/atlases/${ATLAS}.h5ad \
  --out-dir scratch/data/seismic_export/${ATLAS}

Rscript code/04_seismic/sce_from_export.R \
  --export-dir scratch/data/seismic_export/${ATLAS} \
  --out-rds scratch/data/seismic_export/${ATLAS}_sce.rds

for gwas in delange liu; do for tier in broad fine; do
  Rscript code/04_seismic/run_seismic.R \
    --sce-rds scratch/data/seismic_export/${ATLAS}_sce.rds \
    --magma-z scratch/data/gwas/${gwas}_gene_z.tsv \
    --atlas ${ATLAS} --gwas ${gwas} --tier ${tier} \
    --out-dir results/seismic
done; done
```

Then concordance, which is seconds once the inputs exist:

```
python code/08_cross_method/run_scdrs_seismic_concordance.py \
  --out results/concordance/scdrs_vs_seismic.csv
```

Monitoring is `squeue -u $USER` and, after the fact,
`sacct -j <id> --format=JobID,State,Elapsed,MaxRSS,ExitCode` or
`seff <id>`. None of the marker-file and `pgrep` machinery from GCP is
needed — SLURM records exit codes and peak memory itself.

---

## Verification targets

These are the exact figures the GCP runs produced. A rebuilt pipeline
that reproduces them is faithful; a divergence means something in the
chain changed.

### Loader outputs

| Atlas | Cells after pre-filter | Genes | Donors | Cov columns |
| --- | --- | --- | --- | --- |
| Smillie | 348,044 (from 365,492) | 15,207 | 30 | 164 |
| Garrido-Trigo | 28,706 | 14,578 | 12 | — |
| TAURUS | 230,130 | 17,698 | 22 UC, 52 samples | 55 |
| Pan-GI | 203,007 (from 1,596,200) | 15,985 | 67 | 70 |

Pan-GI also logs `Elmentaite2021 cells (HCA Gut overlap): 0` and
`STUDY_EXCLUDE=['Huang2019'] dropped 22626 cells`, and reports depth
from `adata.raw.X`. HGNC canonical-hit survival should be 5/5 for every
atlas.

### Cross-method concordance — 8 rows, all positive and significant

| Atlas | GWAS | Tier | n | Spearman ρ | p |
| --- | --- | --- | --- | --- | --- |
| smillie | delange | broad | 14 | 0.675 | 0.008 |
| smillie | delange | fine | 51 | 0.405 | 0.003 |
| garrido_trigo | delange | broad | 15 | 0.525 | 0.044 |
| garrido_trigo | delange | fine | 75 | 0.480 | 1.3e-5 |
| garrido_trigo | liu | broad | 15 | 0.696 | 0.004 |
| garrido_trigo | liu | fine | 75 | 0.370 | 0.001 |
| taurus | delange | broad | 11 | 0.709 | 0.015 |
| taurus | delange | fine | 38 | 0.576 | 1.6e-4 |

The pattern: broad ρ is higher (fewer cell types), fine ρ is lower but
with far smaller p-values across more types, so fine tier is the more
robust evidence.

### Headline biology to reproduce

**TAURUS de Lange broad, MHC-excluded**: only colonocyte is significant
(`assoc_mcz` 4.18, p=0.002); everything else null. **MHC-included**: the
ranking inverts to immune — dendritic cell 8.19, B cell 6.97,
monocyte/macrophage 5.92, colonocyte drops to 2.52, T cell goes
negative. This reproduces the MHC fragility seen in Garrido on an
independent cohort.

**TAURUS seismic de Lange broad** ranks monocyte/macrophage first at
FDR 0.012 — the only FDR-significant seismic result in the project —
while scDRS ranks colonocyte first. High ρ with a different top hit:
the methods agree on overall ordering, not on the head.

---

## Open items and decisions

### Atlas roster — settled, do not re-litigate

Five atlases: **Smillie, Garrido-Trigo, TAURUS** (the three core UC
atlases in `UC_ATLASES`) plus **HCA Gut** and **Pan-GI** as broad
comparators. Mennillo was dropped on 2026-07-01 because its entire UC
cohort is on treatment, leaving an empty treatment-naive subset;
TAURUS-IBD replaced it. Kong 2023 is Crohn's and was corrected to
Garrido-Trigo. `DECISIONS.md` now carries a canonical roster header at
the top; stale Mennillo and Kong references below it are historical
record.

### Code tasks outstanding

- [x] Write `run_hca_load.py` — the driver for `load_hca_gut.py`,
      applying all ten lessons above. DECISIONS locks HCA covariates as
      `assay, batch, Fraction, sex`, and the v1 filter chain lives in
      the loader. _(landed 2026-10-03)_
- [x] Rewrite the `03_scdrs_compute.slurm` template per Phase 4.
      _(landed 2026-10-03)_
- [ ] Rewrite the remaining four SLURM templates per Phase 4
      (`01_magma.slurm`, `04_seismic.slurm`,
      `donor_loo_array.slurm`, `test_retest_array.slurm`) — same
      `module load miniconda3/3.13` + `mail-user=mukhinda@ucsc.edu`
      fixes.
- [x] Reconcile the Pan-GI deposit UUID in `load_pangi.py` and
      DECISIONS. _(landed 2026-10-03)_
- [ ] Update `scripts/config.sh`: `UCC_REPO=/hb/scratch/mukhinda/UC-Cross-Atlas`,
      `UCC_EMAIL=mukhinda@ucsc.edu`.
- [ ] Update `results/MANIFEST.md` once results exist again.

### Scientific decisions pending

**Seismic permutations.** `constants.py` locks
`SEISMIC_N_PERMUTATIONS=1000`, but every seismic run so far used
`--run-permutations FALSE`, so the p-values are analytic. The single
FDR-significant result in the project — TAURUS monocyte/macrophage at
0.012 — cannot carry a headline claim until a permutation run backs it.

**The Liu runs.** scDRS on Liu exists only for Garrido. Smillie×Liu,
TAURUS×Liu and Pan-GI×Liu are each a 2×2 job, roughly 14 h in total,
and they are what the cross-GWAS concordance axis needs. Worth deciding
deliberately whether that axis is headline evidence or a robustness
check that Garrido alone can carry with a documented limitation.

**Pan-GI batch correction is coarse.** Its `sample_col` and `donor_col`
both resolve to `donorID_unified` — 67 samples across 67 donors, one to
one — because the atlas carries no biopsy-level identifier. Within-donor
batch effects are uncorrected. A methods-section caveat, not a blocker.

### Remaining compute

Roughly 35–45 hours after the rebuild, nearly all of it unattended: HCA
load and both methods, the three Liu scDRS runs, seismic permutations,
test-retest across seeds 1–3, Brown's null-draw feather (~5 h), and the
concordance axes (fast). Actual attention time is maybe 4–6 hours.

The binding constraint remains the manuscript, not the compute. It
should start running in parallel with the last of these jobs rather
than waiting for them.
