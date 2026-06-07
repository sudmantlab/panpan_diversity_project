# Methods detailed description

### Variant input

Biallelic SNVs and structural variants (SVs) were obtained from the long-read
diploid phased call  SNV sets and truvari-merged SV call sets ≤100 kb:

- **Human (HPRC year-1):** `hprc_mapped2ht2t.sorted.filtered.vcf.gz`
  (SNVs) and `hprc.concat.vcf.gz` (SVs), with CHM13 / T2T as reference.
- **Chimpanzee (n = 24):** `pantros_mapped2mPanTro3.sorted.filtered.vcf.gz`
  (SNVs) and `panpan-pt.concat.vcf.gz` (SVs), with mPanTro3 (curated v.20231122)
  as reference.
- **Bonobo (n = 5):** `panpa_mapped2mPanPan1.sorted.filtered.vcf.gz` (SNVs) and
  `panpan-pp.concat.vcf.gz` (SVs), with mPanPan1 as reference.


### Variant Effect Prediction (VEP)

SNVs and SVs were annotated with Ensembl VEP (v110) using species-matched
reference FASTAs and gene annotations:

- **Human:** CHM13/T2T (`GCF_009914755.1` GFF), reference FASTA `ht2t.fasta`.
- **Chimpanzee:** `mPanTro3_modified_chrom_alias.gtf`, reference FASTA
  `mPanTro3.pri.cur.20231122.fasta`.
- **Bonobo:** `mPanPan1_modified_chrom_alias.gtf`, reference FASTA
  `mPanPan1.fasta`.

Annotation files were sorted, bgzip-compressed and tabix-indexed
(`tabix -p gff -C`). VEP was invoked with `--phased`, `--fork 32`,
`--force_overwrite` and a custom GFF/GTF overlap track
(`--custom file=<species.gff/gtf>,format=gff/gtf,type=overlap`). For each
variant the most severe predicted consequence was retained and binned into one
of four impact classes: **HIGH**, **MODERATE**, **LOW**, **MODIFIER**. For SFS
analyses these were collapsed into **moderate/high** (HIGH + MODERATE) and
**modifier/low** (LOW + MODIFIER).

### Site frequency spectra (SFS)

For each species we computed the unfolded SFS of biallelic SNVs and SVs by
counting non-reference alleles across all phased diploid genomes (max
non-reference allele count = 2n − 1: 9 for bonobo, 47 for chimpanzee, 94 for
human; fixed sites with all-alternative alleles excluded). Counts were
converted to proportions within each `species × impact-class` stratum
(stratified by collapsed impact class moderate/high vs modifier/low) and
plotted with `ggplot2`/`ggh4x` (square-root y-axis to expose the rare end of
the spectrum). See `plot_sv_snv_impact_long.R`.

### Per-variant HIGH-impact enrichment of SVs vs SNVs

To compare the per-variant rate of predicted HIGH-impact effects between SNVs
and SVs, we calculated, separately for each species, two ratios from the same
VEP-annotated counts.

1. **Overall variant abundance** — total VEP-annotated SNVs ÷ total
   VEP-annotated SVs, giving the per-species excess of SNVs over SVs.
2. **HIGH-impact enrichment** —

   ```
   fold enrichment = (N_HIGH-SV / N_total-SV) / (N_HIGH-SNV / N_total-SNV)
   ```

   i.e. the proportion of all VEP-annotated SVs classified as HIGH impact
   divided by the corresponding proportion among SNVs.

Per-species values:

| Species     | Total SNVs  | Total SVs | SNV / SV | HIGH-SV % | HIGH-SNV % | SV/SNV HIGH fold |
|-------------|-------------|-----------|----------|-----------|------------|------------------|
| Bonobo      | 7,925,833   | 20,618    | ≈384×    | 3.69%     | 0.0142%    | ≈260×            |
| Chimpanzee  | 29,638,083  | 75,799    | ≈391×    | 2.75%     | 0.0164%    | ≈168×            |
| Human       | 28,235,840  | 62,915    | ≈449×    | 3.78%     | 0.0186%    | ≈203×            |

Across species, biallelic SNVs outnumber SVs by ~380–450-fold, while SVs are
~170–260-fold more likely than SNVs to carry a HIGH-impact predicted effect.

### Gene-level constraint (sHet)

Per-gene constraint was taken from the posterior-mean sHet estimates of
Zeng et al. (2024) GeneBayes
(`Pritchard_Table_genebayes_post_mean.xlsx`, written out as
`Pritchard_Table.csv`). Human Ensembl gene IDs were mapped to gene symbols via
g:Profiler (`gProfiler_hsapiens_6-29-2025_9-12-33 PM.csv`).

For each species, gene symbols were taken from the corresponding annotation
(CHM13 `ht2t_gene.bed.gz`, `mPanTro3_gene.bed.gz`, `mPanPan1_gene.bed.gz`),
restricted to protein-coding genes and filtered to remove uncharacterised
LOC-prefixed loci (entries with no description, "uncharacterized" /
"putative(ly) uncharacterized" descriptions, or non-protein-coding annotations)
and spliceosomal / small-RNA loci. Per gene × species, the maximum VEP impact
across all overlapping variants was retained and re-coded as
**Absent/Modifier/Low** (impact ≤ 2), **Moderate** (3) or **High** (4). The
resulting long- and wide-format tables (`vep_shet_snvs+svs_long_data.csv`,
`vep_shet_snvs+svs_wide_data.csv`) form the basis of all sHet plots.

Genes were further classified into selection regimes following Zeng et al.:
**nearly neutral** (sHet < 10⁻⁴), **weak selection** (10⁻⁴ – 10⁻³),
**strong selection** (10⁻³ – 10⁻²) and **extreme selection** (≥ 10⁻²); these
appear as background shading in all sHet plots.

### Statistics

For each variant type (SNVs / SVs), the distribution of per-gene sHet among
genes carrying HIGH-impact variants was compared between species using
two-sided Wilcoxon rank-sum tests (`wilcox.test` in base R; chimpanzee vs
human and bonobo vs human, run independently).

| Comparison                         | Variant | W       | P-value  |
|------------------------------------|---------|---------|----------|
| Chimpanzee vs Human (HIGH impact)  | SNVs    | 2,248,898 | 5.6 × 10⁻⁴ |
| Chimpanzee vs Human (HIGH impact)  | SVs     | 108,072 | 0.044    |
| Bonobo vs Human (HIGH impact)      | SNVs    | 500,415 | 0.29     |
| Bonobo vs Human (HIGH impact)      | SVs     | 25,979  | 0.17     |

