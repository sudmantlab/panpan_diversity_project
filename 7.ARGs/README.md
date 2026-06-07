# Ancestral recombination graphs (ARGs)

Reconstruction of local genealogies for **chimpanzees, bonobos, and humans** using Bayesian
inference of ancestral recombination graphs with **SINGER**, and use of per-window coalescence
times (TMRCA, T<sub>pooled</sub>, T<sub>within</sub>) to scan for population-specific positive
selection — broadly following the framework of Deng *et al.* Chimpanzees are analysed against
mPanTro3, bonobos against mPanPan1, and humans (HPRC + HGSVC; 1000G 579-sample panel) against
T2T-CHM13v2.0.

See **[`methods_singer_ARG.md`](methods_singer_ARG.md)** for the full methods (VCF preparation
→ genome partitioning → SINGER inference → ARG merging → per-window coalescence statistics →
two-stage filtering → statistics).

## `run_singer_tskit_scripts/` — inference pipeline

Reference/feature setup and the per-species SINGER + `tskit` pipeline.
[Setup commands](run_singer_tskit_scripts/README.txt) cover reference sizes, accessibility
masks, and centromere extraction.

Per-species pipeline (in `bonobo/`, `chimps/`, `humans/…`):

| Step | Script |
|---|---|
| Per-haplotype → diploid VCFs | `haplo2diplo.py` (and `…_human579.py`, `…_unphased.py` variants) |
| Run SINGER genome-wide on 5 Mb chunks | `*_singer_gwide.py` |
| Stitch posterior chunks into chromosome tree sequences | `merge_ARGs.py` |
| Per-1 kb-window TMRCA / T<sub>pooled</sub> / T<sub>within</sub> | `windowed_tmrca_stats.py`, `windowed_tmrca_stats_byPOP.py` |

The human side has several panels: [`hprc_hgsvc/`](run_singer_tskit_scripts/humans/hprc_hgsvc),
[`humans579_dantu/`](run_singer_tskit_scripts/humans/humans579_dantu), and the superpopulation
variants ([`humans579_dantu_superpops/`](run_singer_tskit_scripts/humans/humans579_dantu_superpops),
[`humans579_dantu_hgsvc_superpops/`](run_singer_tskit_scripts/humans/humans579_dantu_hgsvc_superpops)).
The [`chimps/`](run_singer_tskit_scripts/chimps) folder additionally has regional-tree
extraction and visualization (`extract_regional_trees.py`, `ts_to_newick.py`,
`relabel_rescale_nwk.py`, `plot_ts_to_svg.py`, `plot_regional_svg_labeled.py`).

## `output_processing_scripts/` — downstream analysis

R/Python that turns the windowed coalescence tables into filtered scans, outlier calls, tables,
and figures.

- Per-species scans: [`singer_chimp_pops_faceted.R`](output_processing_scripts/singer_chimp_pops_faceted.R),
  [`singer_bonobo.R`](output_processing_scripts/singer_bonobo.R),
  [`singer_human.R`](output_processing_scripts/singer_human.R),
  [`singer_human_579_plusDantu.R`](output_processing_scripts/singer_human_579_plusDantu.R),
  [`singer_human_hdgp1k.R`](output_processing_scripts/singer_human_hdgp1k.R).
- Cross-species intersection: [`singer_interesect_species.R`](output_processing_scripts/singer_interesect_species.R).
- Table builders: [`make_wide_tables.R`](output_processing_scripts/make_wide_tables.R),
  [`make_species_only_tables.R`](output_processing_scripts/make_species_only_tables.R),
  [`make_natgenetics_tables.R`](output_processing_scripts/make_natgenetics_tables.R),
  [`apply_gene_rename.py`](output_processing_scripts/apply_gene_rename.py).
- `bonobos/`, `chimps/`, `humans/` hold per-species inputs, outlier `output_tables/`
  (by percentile/quantile, e.g. 99.90 / 99.99, plus population-agnostic), `plots/` and zoom-ins.

## Other

- [`Deng_singer_analysis/`](Deng_singer_analysis) — overlap/comparison of the scan with regions
  from Deng *et al.* (`gencode47_rename.R`, gene-name and region-overlap tables).
- [`test_singer_hla/`](test_singer_hla) — HLA test region used to validate the pipeline
  (SINGER tutorial, regional-tree and windowed-TMRCA scripts).

> Note: some large outputs are kept local and git-ignored — `output_supp_tables_EXTRA/`,
> `test_singer_hla/tmrca_data.csv`, and the >100 MB tables under
> `output_processing_scripts/humans/input_tables/`.
