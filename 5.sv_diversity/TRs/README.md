# Tandem Repeat (TR) analyses

This directory contains files and code used for analyses and figures related to tandem repeat (TR) variation in humans and chimpanzees.

### files:
- `human_chimp_exp_lens.txt`: This file contains mean allele length and copy number of expansion disorder TR loci in humans and chimpanzees.
- `human_chimp_exp_sumstats.txt`: This file contains per-sample allele length and copy number of expansion disorder TR loci in humans and chimpanzees.

### scripts:
- `scripts/filter_merged_vcf.sh`: Filters multisample VCF input.
- `scrips/get_heterozygosity.py`: Calculates observed and expected TR heterozygosity from a VCF input.
- `scripts/plot_heterozygosity.R`: Generates boxplots comparing TR heterozygosity between humans and chimpanzees.
- `plot_pathogenic_trs.R`: Generates **1)** a line plot showing human and chimp mean allele length for expansion disorder TRs, alongside the pathogenic repeat length threshold in humans.
and **2)** allele length distributions for humans and chimpanzees at four example pathogenic TR loci.
