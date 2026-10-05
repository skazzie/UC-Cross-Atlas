# CONTINUITY — UC-Cross-Atlas

Session handoff / project-state-of-truth. Purpose: stop re-deriving what's
already done and stop assuming things are done that aren't. Last updated:
2026-10-05.

This file is the human-readable status narrative. For the exact
(atlas x method x GWAS x tier) file inventory, see `results/MANIFEST.md`.

---

## Canonical machine (DECISION 2026-10-03, superseding 2026-07-07)

**UCSC Hummingbird is the canonical compute + results machine going
forward.** The GCP VM `ucca-compute` was deleted for non-payment; all VM
results (scDRS group analyses, seismic TSVs, concordance CSV, MAGMA
intermediates, processed `.h5ad` files, raw atlas downloads) were lost
with it. The 175 MB results archive was built on the VM but never
downloaded before deletion. Repo code survived intact on GitHub through
`db3148c`. The GCP trial clock is gone.

Hummingbird layout (verified 2026-10-03, repo at `db3148c`):
- Repo: `/hb/scratch/mukhinda/UC-Cross-Atlas` (cloned, pulled)
- Data/results root: `/hb/scratch/mukhinda/UC-Cross-Atlas/scratch`
- Account/partition: `128x24` (both; verified via `sacctmgr`)
- Modules: `miniconda3/3.13` + `git/2.52.0`; **no system R, no system MAGMA** — both installed into the conda env or fetched as binaries
- 72-core simultaneous cap, infinite partition walltime

See `HUMMINGBIRD_REBUILD.md` for the complete rebuild plan and the
non-negotiable lessons preserved from GCP sessions.

### Phase 1 — conda env rebuilt on Hummingbird (2026-10-05)

Job `794999` COMPLETED in 44:38 (`128x24`, 16 GB, 4 cpus). Env name
`uc-cross-atlas`. Verified: `scdrs 1.0.2`, `numpy 1.26.4`, `scanpy
1.11.5`, `anndata 0.11.4`, `seismicGWAS 1.0.0`, R 4.5.3 with `remotes`,
`optparse`, `Matrix`, `SingleCellExperiment`, `SummarizedExperiment`,
`S4Vectors`, `data.table`, `dplyr`. `scripts/environment.yml` and
`scripts/requirements.txt` re-exported from the live env.

What broke in the first attempt (job `791207`) and the fix:
`devtools::install_github()` on devtools >= 2.5 delegates to the
`remotes` package and errors if it's not installed. The fix was to add
`r-remotes` to the mamba install list and switch to
`remotes::install_github("ylaboratory/seismic", upgrade="never")`
directly, bypassing the devtools delegation layer entirely.

The pip-only Brown's method dependency (`EmpiricalBrownsMethod`) was
NOT installed — it does not exist on PyPI (see
`code/07_regime2_meta/run_brown.py:11-14`); the algorithm is implemented
inline in `run_brown.py`. The R binary from `bioconductor-empiricalbrownsmethod`
IS installed for the few places that still shell out to R for Brown's.

---

## What is GENUINELY done (file-confirmed)

### MAGMA gene-level (`01_magma`) — 4 GWAS
- de Lange UC (GCST004133, harmonised, N_fixed=45975) — genes.out exists
- Liu UC (GCST90446794, harmonised, N_fixed=375508) — run 2026-07-07 on VM
- SCZ (negative control) — .gs built => MAGMA ran (June)
- Yengo height (positive control) — .gs built => MAGMA ran (June)

### scDRS on Garrido-Trigo (`03_scdrs`) — ONE atlas only
Laptop-confirmed outputs (June):
- garrido x de Lange, headline (seed42): broad + fine
- garrido x de Lange, MHC-included: broad + fine
- garrido x de Lange, no-cov (diagnostic): broad + fine
- garrido x height (pos control): broad
- garrido x SCZ (neg control): broad

VM outputs (2026-07-07, depth+donor covariates only — see cov note):
- garrido x de Lange, MHC-excluded + MHC-included (2x2 baseline)
- garrido x Liu, MHC-excluded + MHC-included

### Key scientific finding (file-backed, 2026-07-07)
MHC exclusion INVERTS the top cell-type compartment on Garrido, and this
REPLICATES across both UC GWAS (de Lange + Liu):
- MHC-excluded: epithelial (colonocyte, epithelial progenitor) at top;
  APCs null/negative; T cell strongly negative.
- MHC-included: B cell / dendritic / monocyte-macrophage at top
  (mcz ~4.3-7.4); epithelium marginal; T cell still negative.
- de Lange and Liu APC enrichment near-identical when MHC included
  (B cell +7.29 vs +7.44) despite divergent MAGMA top-hits => cell-type
  attribution is robust to GWAS choice.
Implication: the locked v1 MHC-EXCLUDED default would report an EPITHELIAL
UC signal when the dominant biology is antigen-presenting immune cells.
ESCALATE to M-level: report MHC-included, or both, for cross-atlas cell-type
conclusions. (Logged in DECISIONS 2026-07-07.)

---

## What is UNCERTAIN / must be REDONE

### Garrido x seismicGWAS — RAN but OUTPUTS NOT SAVED. MUST RE-RUN.
The June session reported cross-method Spearman correlations (0.60-0.79)
for Garrido scDRS-vs-seismic across MHC variants. Those NUMBERS exist in
chat history, but NO seismic output files exist on the laptop or HB. A
result with a number but no artifact is not reproducible and cannot go in
the paper. This run must be redone and its outputs saved to
`results/seismic/`.

Prerequisite to verify before re-running: is the `seismicGWAS` R package
actually installed in the VM env? The VM env has R 4.5.3 +
SingleCellExperiment + zellkonverter, but seismicGWAS itself (GitHub
`ylaboratory/seismicGWAS`) was NOT confirmed installed. Check:
```
conda activate uc-cross-atlas && R -e 'library(seismicGWAS); packageVersion("seismicGWAS")'
```
If missing, install before "re-run seismic" is actionable. seismicGWAS also
needs SingleCellExperiment `.rds` inputs per atlas (`<atlas>_sce.rds`) — not
confirmed present on the VM.

---

## What is NOT started (no files anywhere)

- **scDRS / seismicGWAS on the other 4 atlases**: Smillie (PRIMARY,
  account-gated), TAURUS, HCA Gut, Pan-GI (~1.6M cells — the memory run the
  VM exists for). None run.
- **All 5 concordance axes** (`06_concordance`, `08_cross_method`,
  `09_cross_gwas`, `10_broad_atlas_hca`, `11_broad_atlas_pangi`): code
  scaffolded, none run. These consume per-atlas results that mostly don't
  exist yet.
- **Test-retest** (seeds 1,2,3; seismic gate rho >= 0.999): not run.
- **Regime-2 / Brown's method** (`07_regime2_meta`): needs the per-null-draw
  cell-type-Z permutation feather (`null_aggregations_<tier>.feather`,
  ~5 hrs compute). NOT generated. `aggregate_null_draws.py` exists.
- **MAGMA gene-property sanity track** (`05_magma_geneprop`): not confirmed run.
- **Manuscript**: intro was drafted per M4 plan; rest unwritten.

---

## Honest completion estimate

~35-45% of v1. The 13-step code scaffold is complete and infra is solid,
but EXECUTED analysis is one atlas (Garrido) deep on one method (scDRS) plus
controls. Seismic is unsaved (must redo). The bulk — 4 atlases x 2 methods,
all concordance, test-retest, Brown's, writing — remains.

Per the project's own risk register: the binding constraint is FINISHING,
not more analysis or novelty. Highest-value next work is scaling the
established Garrido pattern to the remaining atlases, not deepening Garrido.

---

## Open threads / gotchas (don't lose these)

1. **Unpushed commits.** 3 commits sit on the VM only (scdrs-pin,
   portability fixes, MHC 2x2). Terminal won't accept a paste of a GitHub
   PAT (browser-SSH clipboard issue). Push from a machine where paste works
   (laptop). Until pushed, the DECISIONS record + fixes are VM-local.

2. **Scattered results / no manifest was the meta-blocker.** Results are
   gitignored, so nothing synced. Now that VM is canonical and a MANIFEST
   exists, keep the manifest current and consider tracking the small
   group-analysis TSVs (tiny) so status is never archaeology again.

3. **Garrido cell-count delta: 29,675 (VM) vs 30,068 (laptop), ~1.3%.**
   Same GEO tar, different pandas/anndata versions or join handling. Soft
   tripwire only, but reconcile before any headline Garrido number — the two
   machines are not byte-identical.

4. **Covariate gap on Garrido.** Locked cov set is
   log_n_genes, log_n_counts, donor, sample, sex. On Garrido: sample==donor
   (1:1, redundant) and SEX IS UNAVAILABLE (not in GEO RAW.tar or annotation
   CSV). Today's VM runs used log_n_genes, log_n_counts, donor only. This is
   NOT the locked headline config. Source per-donor sex externally (Garrido
   2023 supp) before any reportable Garrido number, or document the omission.

5. **Smillie (PRIMARY atlas) is account-gated** (Single Cell Portal login +
   data-use agreement). Start that registration BEFORE next session so it
   isn't the blocker — it's the long pole to a headline Smillie x de Lange run.

6. **scdrs 1.0.2 vs modern stack — 5 breaks** (documented in DECISIONS
   2026-07-07): hyphen->underscore flags, species positionals required, NO
   seed flag (seed=42 unenforceable via CLI, Python API only), get_dummies
   bool-dummy cast error (cov file must be pre-numericized to float64), CLI
   doesn't mkdir out-folder. `03_scdrs_compute.slurm` will NOT run
   as-committed on this stack. Every VM scDRS run needs a pre-numericized
   cov file.

7. **README is stale** — still lists Kong/Mennillo as core atlases; they were
   replaced by Garrido (corr 2/7) and TAURUS (corr 16). See README patch.

---

## Suggested next-session order (dependency-ordered)

1. Push the 3 commits (from laptop). Resolve seismicGWAS install on VM.
2. Transfer/regenerate atlas .h5ad + covariate files onto the canonical VM.
3. Smillie: complete SCP account, stage, run scDRS both MHC variants
   (replicate the MHC finding on the PRIMARY atlas — highest-value science).
4. Re-run Garrido seismicGWAS with outputs SAVED. Then Smillie seismic.
5. Roll remaining atlases (TAURUS, HCA, Pan-GI) through both methods.
6. Test-retest, then the 5 concordance axes.
7. Brown's null-draw feather -> regime-2.
8. Manuscript.
