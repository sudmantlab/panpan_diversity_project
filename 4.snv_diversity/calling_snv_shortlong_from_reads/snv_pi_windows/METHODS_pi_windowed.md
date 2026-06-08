# Windowed π comparison between long- and short-read sequencing

This document describes the methodology used to compare per-window nucleotide
diversity (π) between long-read and short-read sequencing in chimpanzees,
including annotation overlap and outlier characterisation.

### 1. Per-window π estimation

- **Tool:** [pixy](https://pixy.readthedocs.io/) (Korunes & Samuk 2021).
- **Input:** all-sites VCFs called from long-read (assembly-based) and
  short-read alignments respectively, against the **mPanTro3** reference.
- **Window size:** **10,000 bp non-overlapping** (`--window_size 10000`).
- **Statistic:** within-population π, computed per window as
  `count_diffs / count_comparisons` after dropping invariant and
  fully-missing sites. Pixy reports `avg_pi`, `no_sites`, `count_diffs`,
  `count_comparisons`, and `count_missing` per (population, chromosome,
  window).
- Window coordinates were verified to be 1-based, inclusive, and
  non-overlapping (e.g. `1–10000`, `10001–20000`, ...).

### 2. Combining long and short outputs

- The per-population pixy outputs from long- and short-read pipelines were
  concatenated with an extra `dataset` column (`long` or `short`) and
  written to `combined_within_pops_pi.txt.gz`.

### 3. Filters applied for all downstream analyses

- **Populations:** kept Central, Eastern, Western. Dropped:
  - **Hybrid** — admixed individuals.
  - **NC** — short-read-only population, no long-read counterpart.
- **Chromosomes:** **autosomes only** (chr1–23 in the mPanTro3 naming;
  chrX and chrY were excluded with `^chr[XY]_`).
- **Per-window NA π:** windows where no sites were callable in a given
  (population, dataset) were dropped from comparisons that need both sides.

### 4. Genomic feature annotations

The following BED tracks were intersected with the 10 kb window grid using
`GenomicRanges::overlapsAny()` (any overlap > 0 bp counts as a hit):

| Track                                          | File                                             |
|------------------------------------------------|--------------------------------------------------|
| Centromere (`Cen`), Telomere (`Telo`)          | `mPanTro3_genomefeatures.bed` (col 4 = label)    |
| Centromeric satellite (`CenSat`)               | `mPanTro3_CenSat.bed`                             |
| Segmental duplications (`SegDup`)              | `mPanTro3_sedefSegDups.bed`                       |
| Tandem-repeat catalogue (`TR`)                 | `pantro_catalog.no_overlaps_simp.bed`             |
| Short-read mask (`SRA`)                        | `mPanTro3_SR_mask.bed`                            |

All intervals were `reduce()`-merged (overlapping/adjacent intervals
collapsed) before computing per-window flags. Annotations were restricted
to the chromosome set used in the π analysis (autosomes hap1 / hap2 as
applicable).

### 5. Region labelling (priority-based)

Each window was assigned a single region label by the priority

```
Centromere > Telomere > CenSat > SegDup > TR > Other (incl. SRA)
```

i.e. a window inside both `CenSat` and `SegDup` is labelled `CenSat`. The
"Other" bucket bundles the short-read mask (`SRA`) plus genuinely
unannotated euchromatin, because at the per-window distributional level
SRA-flagged windows are indistinguishable from unannotated windows
(median log₁₀(π) within ±0.05 of each other and KS D < 0.15).

For density plots only, a coarser labelling
`Centromere > Telomere > CenSat > Other (incl. SegDup, TR, SRA)`
was used to highlight the regions with a clear long-read advantage.

### 6. π_long / π_short ratio (per window)

- For each window per population we paired the long and short π values.
- **Both-zero windows dropped** (no information).
- **Short-side zeros / NAs** (long has a value, short has 0 or no callable
  sites) were replaced by ε = 1 / median(`count_comparisons_short`) ≈
  1.7 × 10⁻⁷ — the per-comparison rate corresponding to "one mismatch
  in a typical short-read window". This avoids ratio blow-up at zero
  without dominating real signal.
- Per-window ratio: `avg_pi_long / pi_short_adj`.
- log₁₀(ratio) is used for densities and direction of effect.

### 7. Outlier definition

- "Long-favoured outliers": **top 1% of windows by π_long/π_short ratio,
  per population**. Each population has ~2,800 such outlier windows.
- ~600 windows are top-1% outliers in **all 3** populations (vs. ~28
  expected by chance — ~21× enrichment for cross-population
  reproducibility).

### 8. Genome-wide fold-change estimates (per region, per population)

For each region and each population, weighted π was computed as

```
pi_pop_region = sum(count_diffs) / sum(count_comparisons)
```

across all windows assigned to that region in that population. The
**fold change** is

```
fold_pop_region = pi_long / pi_short
```

The **per-population fold values were then averaged (mean and median)
across the three chimpanzee populations** to give a sample-size-independent
summary. (A naive pool-then-divide approach is biased because the
short-read dataset has more individuals → more comparisons → π_short
estimate dominated by short-read sample size.)

Pop-pooled fold summary saved to `fold_change_per_region_pop_pooled.tsv`.

### 9. Distribution comparisons

- Per-region log₁₀(π) was compared between long and short reads with:
  - **Kolmogorov–Smirnov** (`ks.test`) — distribution shape.
  - **Wilcoxon rank-sum** (`wilcox.test`) — median shift.
- Tests run pop-by-pop and pop-pooled.
  - Pooled results: `ks_wilcox_long_vs_short_log10pi_pooled.tsv`
  - Per-pop results: `ks_wilcox_long_vs_short_log10pi.tsv`
- Caveat: with n ≈ 250,000 windows per pop in the "Other" category,
  even tiny distributional differences yield p < 1 × 10⁻³⁰⁰. Effect
  size (KS D and median log₁₀ shift) is the meaningful interpretation.


The full pipeline producing the plots in `final_plots/` is in
`make_final_plots.R`.

### 10. Outputs

| File                                                  | Content |
|--------------------------------------------------------|---------|
| `combined_within_pops_pi.txt.gz`                      | concatenated long+short pixy output, with `dataset` column |
| `fold_change_per_region_pop_pooled.tsv`               | per-region long/short fold (3 methods) |
| `ks_wilcox_long_vs_short_log10pi_pooled.tsv`          | distribution tests pop-pooled |
| `ks_wilcox_long_vs_short_log10pi.tsv`                 | distribution tests per pop |
| `final_plots/manhattan_pi_ratio_capped.pdf`           | Manhattan of π_long/π_short, top-1% outliers annotated |
| `final_plots/manhattan_pi_ratio_facet_chrom.pdf`      | same, faceted by chromosome × population |
| `final_plots/manhattan_pi_ratio_facet_chrom_clean.pdf`| same, no annotation overlay |
| `final_plots/density_log10_ratio_by_region.pdf`       | density of log₁₀(π_long/π_short), per pop, by region |
| `final_plots/density_pi_by_dataset_region.pdf`        | density of π, long vs short, per pop, by region |

### 11. stats results

#### Per-region fold (long π / short π), per-population mean

| Region            | mean fold | median fold | n windows / pop |
|-------------------|----------:|------------:|----------------:|
| Centromere        |    ~6.0×  |     ~5.8×   |             218 |
| Telomere          |    ~2.4×  |     ~2.7×   |               5 |
| CenSat            |    ~2.2×  |     ~2.3×   |          10,609 |
| SegDup            |    ~1.0×  |     ~1.0×   |          12,222 |
| TR                |    ~0.89× |     ~0.85×  |         245,965 |
| Other (SRA-only)  |    ~3.1×  |     ~3.3×   |              33 |

#### Outlier composition (top-1% long-favoured windows, mean across pops)

| Region     | % of outliers | background % | enrichment |
|------------|--------------:|-------------:|-----------:|
| CenSat     |        ~71%   |        ~12%  |    ~5.9×   |
| SegDup     |        ~15%   |         ~7%  |    ~2.1×   |
| Centromere |         ~5%   |        ~1.7% |    ~3.0×   |
| Telomere   |         ~0.1% |       ~0.01% |    ~10×    |
| Other      |         ~8%   |         —    |    —       |

#### Distribution tests, log₁₀(π) long vs short, pooled across 3 pops

| Region                       | n long / n short | median Δ log₁₀ | fold | KS D | KS p     |
|------------------------------|-----------------:|---------------:|-----:|-----:|---------:|
| Centromere                   |      617 / 612   |       +0.78    | 6.0× | 0.52 | < 1e-300 |
| Telomere                     |       12 / 13    |       +0.25    | 1.8× | 0.50 | 0.07     |
| CenSat                       |  30,363 / 30,345 |       −0.05    | 0.89×| 0.13 | < 1e-300 |
| Other (SegDup + TR + SRA)    | 771,330 / 773,289|       −0.09    | 0.81×| 0.14 | < 1e-300 |
