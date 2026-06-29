#!/usr/bin/env Rscript

## Load required packages
library(tidyverse)
library(cowplot)
library(knitr)
library(vcfR)
library(janitor)
library(plyranges)
library(parallel)
## Read in the arguments
args = commandArgs(trailingOnly=TRUE)
indir <- args[1]
dataset <- args[2]
threads <- args[3] %>% as.integer()

vcf_raw <- read.vcfR(str_c(indir, "/", dataset, ".vcf"))
#vcf_field_names(vcf_raw, tag = "FORMAT")
## add sniffles ID to info field
#vcf_raw@fix[,8] <- str_c(vcf_raw@fix[,8], ";SNIFFLESID=",vcf_raw@fix[,3])
## add sniffles ID to meta
#vcf_raw@meta <- c(vcf_raw@meta, "##INFO=<ID=SNIFFLESID,Number=1,Type=String,Description=\"Sniffles ID\">")

## add sniffles ID to gt
sample_size_total=ncol(vcf_raw@gt)-1
variant_id <- vcf_raw@fix[,3]
gt_mat <- vcf_raw@gt[, 2:(sample_size_total+1)]
gt_mat_new <- mclapply(1:sample_size_total, 
                       function(i) {str_replace(gt_mat[, i], "[^:]*$", variant_id)}, 
                       mc.cores = threads) %>%
  do.call(cbind, .)
vcf_raw@gt[, 2:(sample_size_total+1)] <- gt_mat_new
## add sniffles suffix to sample ID
colnames(vcf_raw@gt)[-1] <- str_c(colnames(vcf_raw@gt)[-1], "-sniffles")

vcf <-vcfR2tidy(vcf_raw, single_frame = TRUE, format_fields = c("GT"))
vcf_dat <- vcf$dat %>%
  clean_names() 
#write_tsv(vcf_dat, "../sniffles/mPanTro3/panpan-pt.tsv.gz")
vcf_filtered <- vcf_dat %>%
  tidyfast::dt_separate(col=coverage, into = c("V1", "V2", "depth", "V4", "V5"), sep = ",") %>%
  mutate(svlen = abs(svlen), 
         depth=parse_integer(depth),
         seqlen = case_when(svtype == "INS" ~ str_count(alt),
                            svtype == "DEL" ~ str_count(ref),
                            TRUE ~ svlen),
         length_ratio=svlen/seqlen) %>%
  filter(! svtype %in% c("BND"), 
         ## remove the ones with imprecise positions
         precise,
         filter=="PASS",
         ## remove the individual genotypes with missing data
         !is.na(gt_gt),
         ## remove svs that are too long or too short
         svlen >= 50, svlen <= 100000,
         ## remove insertion and deletion calls without an alternate allele sequence
         ! alt %in% c("<INS>", "<DEL>"),
         ## remove insertion and deletion calls with vastly different svlen and seqlen
         length_ratio > 0.95, length_ratio < (1/0.95)) %>%
  mutate(gt = case_when(gt_gt == "0/0" ~ 0L,
                        gt_gt == "0/1" ~ 1L,
                        gt_gt == "1/1" ~ 2L,),
         svtype = factor(svtype, levels = c("INS", "DEL", "DUP", "INV"))) %>%
  group_by(chrom, pos, end, id, qual, svlen, svtype, ref, alt, depth) %>%
  summarize(sample_size=dplyr::n(), allele_count=sum(gt), allele_frequency=allele_count/sample_size/2, het_frequency=sum(gt==1)/sample_size) %>%
  ungroup() %>%
  ## remove sites with fewer than 50% sample size in each species
  #filter(sample_size >= sample_size_total/2) %>%
  ## remove sites with higher than 90% heterozygote frequency in each species
  #filter(het_frequency<=0.9) %>%
  ## remove sites with extremely high sequencing depth
  filter((depth <= quantile(depth, 0.95) | svtype == "DUP")) 

filtered_indel_index <- (vcf_raw@fix[,3] %in% (vcf_filtered %>% filter(svtype %in% c("INS", "DEL")) %>% pull(id)))
vcf_indel_raw <- vcf_raw
vcf_indel_raw@fix <- vcf_raw@fix[filtered_indel_index,]
vcf_indel_raw@gt <- vcf_raw@gt[filtered_indel_index,]
write.vcf(vcf_indel_raw, str_c(indir, "/", dataset, ".filtered.indel.vcf.gz"))

filtered_dupinv_index <- (vcf_raw@fix[,3] %in% (vcf_filtered %>% filter(svtype %in% c("DUP", "INV")) %>% pull(id)))
vcf_dupinv_raw <- vcf_raw
vcf_dupinv_raw@fix <- vcf_raw@fix[filtered_dupinv_index,]
vcf_dupinv_raw@gt <- vcf_raw@gt[filtered_dupinv_index,]
write.vcf(vcf_dupinv_raw, str_c(indir, "/", dataset, ".filtered.dupinv.vcf.gz"))

#write_tsv(vcf_filtered, "../sniffles/mPanTro3/panpan-pt.filtered.tsv.gz")
## Output fasta formatted insertion and deletion sequences
#write_indel_fasta(vcf_filtered, "../sniffles/mPanTro3/panpan-pt.indel.fasta")