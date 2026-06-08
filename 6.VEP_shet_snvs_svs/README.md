# Variant effect prediction and gene constraint (SNVs & SVs)

VEP annotation + gene-constraint (sHet) analysis of SNVs and SVs (≤100 kb) in human, chimpanzee,
and bonobo against species-matched T2T refs (ht2t, mPanTro3, mPanPan1): impact classes, SFS,
HIGH-impact SV-vs-SNV enrichment, and GO enrichment of shared constrained/high-impact genes.
Full methods in [`methods_detailed_description.md`](methods_detailed_description.md); command log
in [`VEP_Readme.txt`](VEP_Readme.txt).

#### Top-level scripts
`vep_shet_combined.R` (merge VEP impact with sHet, build combined tables), `vep_shet_functios.R`
(helpers), `plot_sv_snv_impact_long.R` (impact-stratified SFS, SNVs vs SVs).

#### `snvs_vep/`
VEP of phased SNVs + per-gene sHet (`vep_snvs.R`, `vep_shet.R`, `vep_shet_by_subspecies.R`,
`get_vcf_long_chimp.R`); `long_only/` outputs (large `snv_impact_long.tsv.gz` kept local);
[`GO_analysis_shared_species/`](snvs_vep/GO_analysis_shared_species) — GO enrichment of shared
constrained SNV genes (see its README).

#### `svs_vep/`
VEP of SVs + per-gene sHet (`vep_shet_SVs.R`, `vep_shet_by_subspecies.R`);
[`GO_analysis_shared_species/`](svs_vep/GO_analysis_shared_species) — GO enrichment of shared
high-impact SV genes, with `scripts/`, gene lists and backgrounds (see its README).

#### Supporting data
`annotations/` (per-species gene BEDs), `input_shet_tables/` (Pritchard sHet + gnomAD v4.1
constraint), `output_tables/` (combined VEP+sHet long/wide tables, strong/extreme subsets,
high-impact outliers).
