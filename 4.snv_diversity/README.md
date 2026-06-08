# SNV diversity

Calling and analysis of single-nucleotide variants (SNVs) in **human, chimpanzee, and bonobo**
from both **long-read (assembly-based)** and **short-read** alignments, and downstream
nucleotide-diversity (π), PCA, and long-vs-short-read comparisons. Variants are called against
the species-matched references (mPanTro3, mPanPan1, T2T-CHM13/ht2t).

## `alignments/` — alignment

Read/assembly alignments used as input to SNV calling.

- [`all2all.py`](alignments/all2all.py) — all-to-all assembly alignment.
- [`assemblies2ref/`](alignments/assemblies2ref) — align assemblies to reference
  ([`all2ref.py`](alignments/assemblies2ref/all2ref.py),
  [`all2ref_primaries.py`](alignments/assemblies2ref/all2ref_primaries.py),
  [`all2ref_verkko.py`](alignments/assemblies2ref/all2ref_verkko.py)) with per-species configs
  (`config_bonobos.yaml`, `config_chimps.yaml`, `config_humans.yaml`).
- [`reads2ref/`](alignments/reads2ref) — align reads to reference
  ([`winnowmap_HPRC.py`](alignments/reads2ref/winnowmap_HPRC.py),
  [`winnowmap_PANPAN.py`](alignments/reads2ref/winnowmap_PANPAN.py)).

## `calling_snv_long_from_assemblies/` — SNVs from assemblies

SNV calling from primary assemblies aligned to the reference, then conversion to diploid genotypes.

- [`minimap2snvs.py`](calling_snv_long_from_assemblies/minimap2snvs.py) — call SNVs from
  assembly alignments.
- `haplo2diplo.py` (and `haplo2diplo_human579+Dantu.py`, `haplo2diplo_unphased.py`) — merge
  per-haplotype calls into single-column phased diploid genotypes.

## `calling_snv_shortlong_from_reads/` — SNVs from reads + diversity

Joint SNV calling from short- and long-read alignments, plus nucleotide-diversity analysis.

- [`allbams2snvs_chimps.py`](calling_snv_shortlong_from_reads/allbams2snvs_chimps.py),
  [`allbams2snvs_bonobos.py`](calling_snv_shortlong_from_reads/allbams2snvs_bonobos.py) — joint
  SNV calling across BAMs.
- [`code_diversity_stats/`](calling_snv_shortlong_from_reads/code_diversity_stats) — per-window
  π / diversity from all-sites VCFs with **pixy** (`diversity_stats_fromVCFs.py`,
  `diversity_stats_HSA_fromVCFs.py`, `calculate_pi_ratio.py`, `hgdp_tgp_1KG_pi.py`,
  `pixy_snakefile.py`, configs).
- Analysis & plots: [`snv_diversity.R`](calling_snv_shortlong_from_reads/snv_diversity.R),
  [`snv_diversity_plots.R`](calling_snv_shortlong_from_reads/snv_diversity_plots.R),
  [`snv_pca.R`](calling_snv_shortlong_from_reads/snv_pca.R),
  [`stats_longread_vs_shortread.R`](calling_snv_shortlong_from_reads/stats_longread_vs_shortread.R).
- `input_pca_eigen/` — PCA inputs (VCFs from long-read assemblies, short reads, and combined
  short+long reads).
- [`snv_pi_windows/`](calling_snv_shortlong_from_reads/snv_pi_windows) — windowed π comparison
  between long- and short-read calls; see its
  [README](calling_snv_shortlong_from_reads/snv_pi_windows/README.md) and
  [METHODS](calling_snv_shortlong_from_reads/snv_pi_windows/METHODS_pi_windowed.md).

## `chimp_bonobo_range_maps/` — sampling maps

Geographic range maps and sampling context for the chimpanzee subspecies and bonobo.

- [`plot_dataset.R`](chimp_bonobo_range_maps/plot_dataset.R),
  [`plot_tree.R`](chimp_bonobo_range_maps/plot_tree.R),
  [`mapchimps/mapRanges.R`](chimp_bonobo_range_maps/mapchimps/mapRanges.R).
