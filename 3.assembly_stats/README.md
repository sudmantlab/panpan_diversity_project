# Assembly statistics

Quality and completeness assessment of the PANPAN genome assemblies (hifiasm, hifiasm-HiC,
and verkko haplotypes) — contiguity, base accuracy, gene completeness, telomere-to-telomere
status, and reference alignments. Both haplotypes (H1/H2) are evaluated where available.

Reference T2T assemblies used throughout: **mPanPan1** (bonobo) and **mPanTro3** (chimpanzee).

## Scripts

| Script | Purpose |
|---|---|
| [`assemblystats.py`](assemblystats.py) | Computes per-assembly contiguity statistics — Nx / Lx and the area-under-the-Nx-curve (auN) from contig lengths |
| [`get_assembly_stats.py`](get_assembly_stats.py) | Gathers per-sample stats against the reference panel (hg38, T2T-CHM13, mPanTro3, mPanPan1, mPonAbe1) and writes the unified table |
| [`plotstats.R`](plotstats.R) | Plots genome-wide stats (NGx curves, auN H1 vs H2, coverage) from `genome_stats_summary.tsv` + metadata |
| [`plot_compleasm.R`](plot_compleasm.R) | Plots gene-completeness (compleasm) comparisons between haplotypes |

## Summary tables

| File | Contents |
|---|---|
| [`genome_stats_summary.tsv`](genome_stats_summary.tsv) | Per-individual summary: N50 (H1/H2), QV, k-mer + depth coverage, assembly method, metadata |
| [`compleasm/`](compleasm) | compleasm (BUSCO-style) gene-completeness counts per assembly/haplotype (`compleasm.csv`, single-copy H1/H2) |
| [`input_tables/`](input_tables) | Metadata and reference tables feeding the scripts (`PanPan_Metadata.tsv`, reference lists, unified stats) |

## Plots

- [`stats_plots/`](stats_plots) — NGx curves, auN H1-vs-H2 comparisons, coverage plots, and
  compleasm H1-vs-H2 figures.
- [`pafr_all2Ref_plots/`](pafr_all2Ref_plots) — whole-genome **pafr** alignment dot-plots of
  all assemblies against the reference (separate bonobo / chimpanzee / verkko-chimpanzee PDFs).

## Telomeres & gaps

[`telomeres_gaps/`](telomeres_gaps) — telomere-to-telomere (T2T) status, telomere counts, and
assembly gaps, derived from RagTag scaffolding of each assembly against mPanPan1 / mPanTro3.

- Scripts: [`plot_telomeres.R`](telomeres_gaps/plot_telomeres.R),
  [`plot_telomeres_new.R`](telomeres_gaps/plot_telomeres_new.R),
  [`plot_telomeres_new_yesXY.R`](telomeres_gaps/plot_telomeres_new_yesXY.R).
- `new_inputs/` — per-assembly RagTag telomere/gap/summary CSVs.
- `new_outputs/` — counts and figures (nT2T and nTelomeres per species, gaps, stacked plots),
  with autosome-only and all-chromosome variants.
- The `yesXY/` subfolders repeat the analysis **including the sex chromosomes**; the default
  outputs are autosome-focused.
