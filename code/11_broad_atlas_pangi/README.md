# 11_broad_atlas_pangi — Pan-GI comparator (integration-pipeline robustness)

Concordance axis #5b from PLAN.md §"Five concordance axes": Pan-GI
(Oliver 2024) as an integration-pipeline-robustness comparator, with
dual analysis (with/without donor-overlap exclusion).

**Critical framing:** Pan-GI integrated Smillie 2019 and anchored Kong
2023, so 2/3 of UC trio donors are inside Pan-GI. Pan-GI is **not** an
independent broad-atlas replication of the trio.

**Scope:** scDRS + seismicGWAS, 2 GWAS each = 4 runs total. Plus
donor-overlap-excluded re-analysis.

## Counts policy

Same as `10_broad_atlas_hca`: raw count layer with published cell-type
labels, no re-clustering.

## Pan-GI dual analysis

1. **With donor overlap retained:** report concordance with each UC atlas
   — sanity check that mostly-the-same-data produces high concordance.
2. **With donor-overlap exclusion (HEADLINE):** subset Pan-GI to donors
   *not* from Smillie 2019 or Kong 2023, using the donor-attribution
   metadata in `data/atlases/donor_metadata/pangi_donor_metadata.csv`
   (built in M1). Re-compute concordance.

Donor-overlap-excluded Pan-GI is **still an scVI-integrated multi-study
atlas** — just without our trio's donors. It is *not* "Pan-GI without
integration." Methods text states this explicitly.

## Output

- `results/broad_atlas/pangi_concordance_full.tsv` — donor-overlap retained.
- `results/broad_atlas/pangi_concordance_excluded.tsv` — donor-overlap excluded (headline).
- `results/broad_atlas/pangi_donor_audit.tsv` — cross-atlas donor audit
  across all 5 atlases.

## v1 sensitivities (DECISIONS.md 2026-05-20 (3/7))

Pan-GI Extended+ contains Elmentaite2021 (= HCA Gut, ~398k cells) and
Kong2023 (CD-only, ~235k cells); HCA Gut is therefore *nested* within
Pan-GI rather than independent. Smillie2019 was not found in the
`study` column but is scanned for empirically rather than assumed absent.

Two paired sensitivity runs are required for the Pan-GI x UC GWAS
analyses:

1. **HCA Gut overlap test** — re-run with Elmentaite2021 cells removed.
   Loader: `code/02_atlas_prep/load_pangi.load_pangi_no_elmentaite()`.
   Compares the resulting cell-type prioritization against the headline
   Pan-GI run to assess whether shared donors with HCA Gut drive the
   signal.

2. **Smillie overlap empirical scan** — re-run with any cells whose
   `donorID_unified` matches the Smillie 2019 donor-ID pattern
   (`^(N|UC)\d+$`) removed. Expected to be a no-op (0 cells matched in
   inspection); the loader still produces a report so the empirical
   overlap is documented rather than asserted.
   Loader: `code/02_atlas_prep/load_pangi.load_pangi_no_smillie()`.

## v1 filter chain

Defined in `load_pangi.py`:

```
disease in {normal, ulcerative colitis, inflammatory bowel disease}
organ_unified in {ascending colon, caecum, colon, descending colon,
                  rectum, sigmoid colon, transverse colon}
sample_type != "Organ_donor_resection"
```

Expected post-filter: ~150-200k cells (tractable on the 128x24 partition
with `--mem=48G --time=06:00:00`).

## Slice choice

Pan-GI v1 = **Extended+ - 18485 genes** slice (1,596,200 cells), the
only Pan-GI slice with all lineages (epithelial + immune + stromal +
endothelial + neural). Other slices (e.g. "Extended - Large Intestine")
are lineage-restricted and cannot support cross-lineage prioritization.

## Driver script

`code/11_broad_atlas_pangi/run_pangi_comparison.py`. Full + two sensitivities
+ donor audit. CLI:

```bash
python code/11_broad_atlas_pangi/run_pangi_comparison.py \
    --pangi-results-base results/pangi \
    --pangi-no-elmentaite-results-base results/pangi_no_elmentaite \
    --pangi-no-smillie-results-base results/pangi_no_smillie \
    --uc-atlases smillie garrido_trigo mennillo \
    --gwas delange liu \
    --methods scdrs seismic \
    --tiers broad fine \
    --scdrs-dir results/scdrs \
    --seismic-dir results/seismic \
    --out-full results/broad_atlas/pangi_concordance_full.tsv \
    --out-no-elmentaite results/broad_atlas/pangi_concordance_no_elmentaite.tsv \
    --out-no-smillie results/broad_atlas/pangi_concordance_no_smillie.tsv \
    --out-donor-audit results/broad_atlas/pangi_donor_audit.tsv \
    --pangi-h5ad data/atlases/pangi.h5ad
```

Donor audit reads h5ad files with `backed='r'`. Atlases not on disk are
skipped with a warning rather than failing. Runs on login node.
