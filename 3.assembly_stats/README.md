# Assembly statistics

Quality/completeness of the PANPAN assemblies (hifiasm, hifiasm-HiC, verkko; H1/H2) —
contiguity, QV, gene completeness, T2T status, and reference alignments. Refs: mPanPan1, mPanTro3.

#### Scripts
`assemblystats.py` (Nx/Lx/auN), `get_assembly_stats.py` (stats vs reference panel),
`plotstats.R` (NGx/auN/coverage), `plot_compleasm.R` (completeness H1 vs H2).

#### Tables
`genome_stats_summary.tsv` (N50/QV/coverage per individual), `compleasm/` (BUSCO-style
completeness counts), `input_tables/` (metadata + reference lists).

#### Plots
`stats_plots/` (NGx, auN, coverage, compleasm); `pafr_all2Ref_plots/` (whole-genome alignment
dot-plots vs reference — bonobo/chimp/verkko).

#### `telomeres_gaps/`
T2T status, telomere counts, and gaps from RagTag scaffolding vs mPanPan1/mPanTro3. Scripts
`plot_telomeres*.R`; `new_inputs/` (RagTag CSVs), `new_outputs/` (counts + figures); `yesXY/`
repeats with sex chromosomes (default outputs are autosome-focused).
