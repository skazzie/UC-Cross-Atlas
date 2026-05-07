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
