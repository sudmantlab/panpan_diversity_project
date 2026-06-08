# Ancestral recombination graphs (ARGs)

SINGER-based ARG reconstruction for chimpanzees, bonobos, and humans, using per-window
coalescence times (TMRCA, T<sub>pooled</sub>, T<sub>within</sub>) to scan for population-specific
selection (refs: mPanTro3, mPanPan1, ht2t). Full methods in
[`methods_singer_ARG.md`](methods_singer_ARG.md).

#### `run_singer_tskit_scripts/`
Per-species SINGER + `tskit` pipeline ([setup](run_singer_tskit_scripts/README.txt)):
`haplo2diplo*.py` (diploid VCFs) → `*_singer_gwide.py` (SINGER on 5 Mb chunks) →
`merge_ARGs.py` (stitch chromosome tree sequences) → `windowed_tmrca_stats*.py` (per-1 kb
TMRCA/T<sub>pooled</sub>/T<sub>within</sub>). Human panels in `humans/` (`hprc_hgsvc/`,
`humans579_dantu/`, and superpop variants).

#### `output_processing_scripts/`
Turns windowed coalescence tables into scans, outliers, tables and figures: per-species scans
(`singer_chimp_pops_faceted.R`, `singer_bonobo.R`, `singer_human*.R`), cross-species intersection
(`singer_interesect_species.R`), table builders (`make_*_tables.R`, `apply_gene_rename.py`).
`bonobos/`, `chimps/`, `humans/` hold inputs, outlier `output_tables/` (by percentile) and plots.

#### Other
`Deng_singer_analysis/` (overlap with Deng *et al.* regions); `test_singer_hla/` (HLA test region
validating the pipeline).
