# SV diversity

Structural-variant diversity across human, chimpanzee, and bonobo — alignment, tandem-repeat
analyses, and long/short structural-variant sets.

#### `alignments/`
Read/assembly alignment for SV calling: `winnowmap_*.py` (followed by Sniffles2),
`minimap2svs.py` (followed by SVIM-asm; see also
[`all2ref.py`](https://github.com/sudmantlab/panpan_diversity_project/blob/main/4.snv_diversity/alignments/assemblies2ref/all2ref.py)
for the alignments used for SyRI).

#### `TRs/`
Tandem-repeat (TR) variation in humans and chimpanzees — heterozygosity and
pathogenic-expansion analyses (scripts + data); see its [README](TRs/README.md).

#### `longSVs/` · `shortSVs/`
Long- and short-read structural-variant call sets *(in progress)*.
