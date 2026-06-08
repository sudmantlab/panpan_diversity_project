# SNV diversity

Read-based and Assembly-based mapping, SNV calling, diversity and population structure analysis,
and long-vs-short comparisons

#### `alignments/`
Alignments for SNV calling: assembly all-to-all (`all2all.py`), assemblies-to-ref
(`all2ref*.py`, per-species configs), reads-to-ref (`winnowmap_*.py`).

#### `calling_snv_long_from_assemblies/`
SNVs called from assembly alignments (`minimap2snvs.py`), merged to diploid genotypes
(`haplo2diplo*.py`).

#### `calling_snv_shortlong_from_reads/`
SNV calling from BAMs (`allbams2snvs_*.py`); diversity stats code (`code_diversity_stats/`) and
plots (`snv_diversity*.R`, `snv_pca.R`, `stats_longread_vs_shortread.R`); PCA inputs
(`input_pca_eigen/`); long-vs-short windowed π in
[`snv_pi_windows/`](calling_snv_shortlong_from_reads/snv_pi_windows).

#### `chimp_bonobo_range_maps/`
Sampling/range maps (`plot_dataset.R`, `plot_tree.R`, `mapchimps/mapRanges.R`).
