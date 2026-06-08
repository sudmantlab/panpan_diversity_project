# SNV diversity

SNV calling in human, chimpanzee, and bonobo from long-read (assembly-based) and short-read
alignments, with downstream π, PCA, and long-vs-short comparisons (refs: mPanTro3, mPanPan1, ht2t).

#### `alignments/`
Alignments feeding SNV calling: `all2all.py` (assembly all-to-all), `assemblies2ref/`
(`all2ref*.py` + per-species configs), `reads2ref/` (`winnowmap_*.py`).

#### `calling_snv_long_from_assemblies/`
SNVs from assembly alignments (`minimap2snvs.py`), merged to diploid genotypes (`haplo2diplo*.py`).

#### `calling_snv_shortlong_from_reads/`
Joint SNV calling across BAMs (`allbams2snvs_*.py`); π/diversity via **pixy**
(`code_diversity_stats/`); analysis & PCA (`snv_diversity*.R`, `snv_pca.R`,
`stats_longread_vs_shortread.R`); PCA inputs (`input_pca_eigen/`); long-vs-short windowed π in
[`snv_pi_windows/`](calling_snv_shortlong_from_reads/snv_pi_windows) (see its README/METHODS).

#### `chimp_bonobo_range_maps/`
Sampling/range maps (`plot_dataset.R`, `plot_tree.R`, `mapchimps/mapRanges.R`).
