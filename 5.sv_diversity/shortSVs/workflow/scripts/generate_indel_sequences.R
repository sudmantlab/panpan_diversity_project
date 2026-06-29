#!/usr/bin/env Rscript

## Load required packages
library(tidyverse)
library(vcfR)
library(janitor)
library(plyranges)
## Read in the arguments
args = commandArgs(trailingOnly=TRUE)
indir <- args[1]
dataset <- args[2]
tandem_repeat_path <- args[3]
threads <- args[4] %>% as.integer()

vcf_raw <- read.vcfR(str_c(indir, "/", dataset, ".sniffles_svim-asm.truvari.sorted.merged.vcf.gz"))
tandem_repeat <- read_tsv(tandem_repeat_path, col_names = c("seqnames", "start", "end")) %>% as_granges()
sample_size_total <- (ncol(vcf_raw@gt)-1)/2
vcf <-vcfR2tidy(vcf_raw, single_frame = TRUE, format_fields = c("GT", "ID"))
vcf_dat <- vcf$dat %>%
  clean_names() %>%
  mutate(svtype = factor(svtype, levels = c("INS", "DEL", "DUP", "INV")))
vcf_long <- vcf_dat %>%
  dplyr::select(-gt_gt_alleles) %>%
  tidyfast::dt_separate(indiv, into=c("sample_id", "caller"), sep="-", remove = TRUE) %>%
  dplyr::select(id, chrom, pos, ref, alt, svlen, svtype, end, caller, sample_id, gt_id, gt_gt) %>%
  group_by(id, chrom, pos, ref, alt, svlen, svtype, end, caller) %>%
  mutate(n_id=sum(!is.na(unique(gt_id)))) %>%
  ungroup() %>%
  mutate(allele_count=case_when(is.na(gt_gt) | gt_gt=="0/0" ~ 0L,
                                gt_gt=="0/1" | gt_gt=="1/0" ~ 1L,
                                gt_gt=="1/1" ~ 2L,)) %>%
  ## INV called by svimasm doesn't have an svlen field, so I'll need to calculate it
  mutate(svlen=ifelse(is.na(svlen), str_count(alt), svlen))
sv_id_table <- vcf_long %>%
  distinct(id, caller, gt_id) %>%
  filter(!is.na(gt_id))
vcf_summary <- vcf_long %>%
  filter(n_id>0) %>%
  group_by(id, chrom, pos, ref, alt, svlen, svtype, end, caller, n_id) %>%
  summarize(total_allele_count=sum(allele_count), allele_count_vector=str_c(allele_count, collapse = "")) %>%
  ungroup()
sv_intersecting_tr <- vcf_summary %>%
  distinct(id, svtype, svlen=abs(svlen), chrom, pos, end) %>%
  ## when end is NA, fill it with pos if it's an insertion, or pos+svlen if it's a deletion or inversion
  mutate(end=ifelse(is.na(end) & svtype=="INS", pos, end),
         end=ifelse(is.na(end), pos+svlen, end)) %>%
  mutate(seqnames=chrom, start=(pos+end)/2, end=start) %>%
  as_granges() %>%
  filter_by_overlaps(tandem_repeat) %>%
  as_tibble()
vcf_summary_wide <- vcf_summary %>%
  pivot_wider(names_from = caller, values_from = c("n_id", "total_allele_count", "allele_count_vector")) %>%
  mutate(svimasm=!is.na(n_id_svimasm), sniffles=!is.na(n_id_sniffles)) %>%
  mutate(method_support = case_when(svimasm&sniffles ~ "consensus",
                                    sniffles ~ "sniffles_only",
                                    svimasm ~ "svimasm_only")) %>%
  mutate(intersecting_tr=(id %in% sv_intersecting_tr$id))
sv_id_table %>%
  write_tsv(str_c(indir, "/", dataset, ".sniffles_svim-asm.truvari.sorted.merged.id.tsv"))
vcf_summary_wide %>%
  write_tsv(str_c(indir, "/", dataset, ".sniffles_svim-asm.truvari.sorted.merged.summary.tsv"))
vcf_summary_wide %>%
  filter(svtype %in% c("INS", "DEL")) %>%
  mutate(alt=ifelse(svtype=="INS", alt, ref)) %>% 
  select(id, alt) %>%
  mutate(id=str_c(">", id), n=row_number()) %>%
  pivot_longer(1:2) %>%
  arrange(n, desc(name)) %>%
  dplyr::select(value) %>%
  write_tsv(str_c(indir, "/", dataset, ".sniffles_svim-asm.truvari.sorted.merged.indel.fasta"), col_names = FALSE)

