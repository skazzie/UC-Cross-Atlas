# scripts/ — Hummingbird quickstart

Everything you need to go from a fresh Hummingbird account to a running
M1 + M3 pipeline. Per-template detail lives in `scripts/slurm/README.md`.

## TL;DR

```bash
# 1. SSH to Hummingbird, clone the repo
ssh <username>@hummingbird.ucsc.edu
git clone git@github.com:skazzie/UC-Cross-Atlas.git
cd UC-Cross-Atlas/uc-cross-atlas

# 2. Configure paths and SLURM defaults
$EDITOR scripts/config.sh
#  - UCC_PARTITION         your general partition (sinfo)
#  - UCC_PARTITION_HIGHMEM your high-memory partition (for Pan-GI)
#  - UCC_ACCOUNT           your project/account code if required
#  - UCC_EMAIL             where SLURM should mail
#  - UCC_SCRATCH           where atlases + results live

# 3. Set up the conda env + R packages (~10-20 min, login node)
bash scripts/setup_env.sh

# 4. Auto-fetch what's freely downloadable; manually fetch the rest
bash scripts/download_refs.sh
#  Manual: Smillie / Kong / Mennillo / HCA Gut / Pan-GI / Liu / SCZ
#  See the printed instructions and the M1 spot-checks in PLAN.md.

# 5. M1 MAGMA pipeline (3 GWAS x ~30 min each)
sbatch --export=ALL,GWAS=delange scripts/slurm/01_magma.slurm
sbatch --export=ALL,GWAS=liu     scripts/slurm/01_magma.slurm
sbatch --export=ALL,GWAS=scz     scripts/slurm/01_magma.slurm

# 6. Once M1 spot-checks pass, M3 regime-1 jobs (see slurm/README.md)
```

## What lives where on Hummingbird

| Path | What | Why |
|------|------|-----|
| `$HOME/uc-cross-atlas/` | the repo (small) | code lives here |
| `$UCC_SCRATCH/data/gwas/` | GWAS sumstats (~GB each) | scratch, not $HOME |
| `$UCC_SCRATCH/data/atlases/` | h5ad + sce.rds files (tens of GB) | scratch |
| `$UCC_SCRATCH/data/atlases/donor_metadata/` | donor-attribution CSVs | small, but kept with atlases |
| `$UCC_SCRATCH/data/reference/` | 1000G EUR, MAGMA binary, NCBI37 gene loc, CL ontology | scratch |
| `$UCC_SCRATCH/results/` | all generated outputs | scratch |
| `$UCC_SCRATCH/logs/` | optional log dir if you redirect SLURM `--output=` | scratch |

`config.sh` creates these directories on first source.

## What's locked vs what you have to fill in

Locked in `config.sh` defaults (override by editing or by exporting env
vars before submit):

- Conda env name: `uc-cross-atlas`
- Scratch path: `/hb/scratch/$USER/uc-cross-atlas`

You must fill in (or override at submit time):

- `UCC_PARTITION` and `UCC_PARTITION_HIGHMEM` — get from `sinfo` /
  Hummingbird docs.
- `UCC_ACCOUNT` — only if your cluster's QOS requires it.
- `UCC_EMAIL` — for SLURM notifications.
- The `#SBATCH --partition=128x24` and `#SBATCH --mail-user=...` lines
  in each `.slurm` file (one-time edit), OR override at submit:

  ```bash
  sbatch -p $UCC_PARTITION --mail-user=$UCC_EMAIL ... script.slurm
  ```

You also have to fill in the Liu 2023 sample size (`N_FIXED`) and the
Trubetskoy SCZ sample size in `scripts/slurm/01_magma.slurm` — verify
both at M1 from the source papers.

## Smoke testing before submitting batch jobs

After `setup_env.sh` succeeds, this should already be passing:

```bash
conda activate uc-cross-atlas
pytest code/06_concordance/test_metrics.py -q
python -c "import scdrs; print(scdrs.__version__)"
```

If the conda activation isn't sticking in your shell, source the conda
hook explicitly:

```bash
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate uc-cross-atlas
```

## When something is wrong

- `module load miniconda3` fails — run `module avail` and edit
  `scripts/setup_env.sh` + the `module load` line at the top of each
  `.slurm` to match your cluster's module name.
- `conda env create` is slow — that's normal first time (~10 min).
- A `wget` in `download_refs.sh` 404s — URLs at CTGLab and GWAS Catalog
  occasionally change; fall through to manual download per the printed
  instructions.
- `sbatch` errors about partition/account — `sinfo` and `sacctmgr show
  assoc user=$USER` show what's available; update `config.sh`.

## What this scaffolding does NOT do

- **Run anything for you.** You execute jobs; this is templating.
- **Download datasets that require accounts.** Single Cell Portal, PGC,
  cellxgene Discover all need you to log in interactively.
- **Decide if results pass the M3 gate.** Test-retest pass thresholds
  (scDRS ρ ≥ 0.9; seismicGWAS ρ ≥ 0.999) are evaluated by you against
  the outputs.
- **Substitute for `module avail`.** Hummingbird-specific module names
  vary; verify before your first run.

See `scripts/slurm/README.md` for per-template details and the full
M1 → M3 submission sequence.
