# Variant effect prediction and gene constraint (SNVs & SVs)

Functional annotation and constraint analysis of biallelic SNVs and structural variants (SVs,
≤100 kb) across **human, chimpanzee, and bonobo**, using species-matched T2T references
(CHM13/ht2t, mPanTro3, mPanPan1). Variants are annotated with Ensembl VEP, binned into impact
classes, combined with gene-level constraint (sHet), and tested for site-frequency-spectrum
differences and GO enrichment of constrained / high-impact genes shared across species.

See **[`methods_detailed_description.md`](methods_detailed_description.md)** for the full
methods (variant input → VEP → SFS → HIGH-impact SV-vs-SNV enrichment → sHet → statistics).
**[`VEP_Readme.txt`](VEP_Readme.txt)** holds the raw command log (minimap2/paftools, bcftools
`+liftover`, liftoff, VEP/snpEff invocations).

## Top-level scripts

| Script | Purpose |
|---|---|
| [`vep_shet_combined.R`](vep_shet_combined.R) | Main analysis: merge VEP impact calls with sHet, build the combined SNV+SV tables |
| [`vep_shet_functios.R`](vep_shet_functios.R) | Shared helper functions used across the sHet/VEP scripts |
| [`plot_sv_snv_impact_long.R`](plot_sv_snv_impact_long.R) | Impact-stratified site frequency spectra (SNVs vs SVs) per species |

## `snvs_vep/` — SNV annotation

VEP of phased SNVs and per-gene sHet aggregation.

- [`vep_snvs.R`](snvs_vep/vep_snvs.R), [`vep_shet.R`](snvs_vep/vep_shet.R),
  [`vep_shet_by_subspecies.R`](snvs_vep/vep_shet_by_subspecies.R),
  [`get_vcf_long_chimp.R`](snvs_vep/get_vcf_long_chimp.R)
- `long_only/` — long-read impact outputs (`snv_impact_long.tsv.gz` is large and kept local).
- [`GO_analysis_shared_species/`](snvs_vep/GO_analysis_shared_species) — GO enrichment of
  constrained SNV genes shared across species (gene lists for the cross-species and
  human–chimp shared sets); see its [README](snvs_vep/GO_analysis_shared_species/README.md).

## `svs_vep/` — SV annotation

VEP of structural variants and per-gene sHet aggregation.

- [`vep_shet_SVs.R`](svs_vep/vep_shet_SVs.R),
  [`vep_shet_by_subspecies.R`](svs_vep/vep_shet_by_subspecies.R)
- [`GO_analysis_shared_species/`](svs_vep/GO_analysis_shared_species) — GO enrichment of
  high-impact SV genes shared across species, with enrichment `scripts/`
  ([`run_GO_enrichment.R`](svs_vep/GO_analysis_shared_species/scripts/run_GO_enrichment.R),
  `_build_loc_map.R`), gene lists and backgrounds; see its
  [README](svs_vep/GO_analysis_shared_species/README.md).

## Supporting data

| Folder | Contents |
|---|---|
| [`annotations/`](annotations) | Per-species gene models / BEDs (`ht2t_gene.bed.gz`, `mPanTro3_gene.bed.gz`, `mPanPan1_gene.bed.gz`, …) |
| [`input_shet_tables/`](input_shet_tables) | Constraint inputs — Pritchard sHet table and gnomAD v4.1 constraint metrics |
| [`output_tables/`](output_tables) | Combined VEP+sHet outputs (`vep_shet_snvs+svs_long/wide_data.csv`, strong/extreme subsets, high-impact outliers) |
