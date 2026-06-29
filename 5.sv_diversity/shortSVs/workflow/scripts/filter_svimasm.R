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

vcf_raw <- read.vcfR(str_c(indir, "/", dataset, ".truvari.merged.missing2ref.vcf.gz"))
#vcf_field_names(vcf_raw, tag = "FORMAT")
## make unique variant ID 
vcf_raw@fix[,3] <- str_c("svim_asm.", 1:nrow(vcf_raw@fix))

## add svimasm ID to info field
#vcf_raw@fix[,8] <- str_c(vcf_raw@fix[,8], ";SVIMASMID=",vcf_raw@fix[,3])
## add svimasm ID to meta
#vcf_raw@meta <- c(vcf_raw@meta, "##INFO=<ID=SVIMASMID,Number=1,Type=String,Description=\"Svim-asm ID\">")

## rename DUP:TANDEM as DUP
vcf_raw@fix[,8] <- str_replace(vcf_raw@fix[,8], "DUP:TANDEM", "DUP")

## add svimasm ID to gt
sample_size_total=ncol(vcf_raw@gt)-1
variant_id <- vcf_raw@fix[,3]
gt_mat <- vcf_raw@gt[, 2:(sample_size_total+1)]
gt_mat_new <- mclapply(1:sample_size_total, 
                       function(i) {str_c(gt_mat[, i], ":", variant_id)}, 
                       mc.cores = threads) %>%
  do.call(cbind, .)
vcf_raw@gt[, 2:(sample_size_total+1)] <- gt_mat_new
## add ID to GT format
vcf_raw@gt[, 1]  <- str_c(vcf_raw@gt[, 1], ":ID")
## add individual variant ID to meta
vcf_raw@meta <- c(vcf_raw@meta, '##FORMAT=<ID=ID,Number=1,Type=String,Description="Individual sample SV ID for multi-sample output">')
## remove INFO fields added by truvari
vcf_raw@meta <- (!str_detect(vcf_raw@meta, "Collapse|Consolidated")) %>%
  vcf_raw@meta[.]
vcf_raw@fix[,8] <- vcf_raw@fix[,8] %>%
  str_remove(";NumCollapsed.*")
## add svimasm suffix to sample ID
colnames(vcf_raw@gt)[-1] <- str_c(colnames(vcf_raw@gt)[-1], "-svimasm")

vcf <-vcfR2tidy(vcf_raw, single_frame = TRUE, format_fields = c("GT"))
vcf_dat <- vcf$dat %>%
  clean_names() 
#write_tsv(vcf_dat, "../sniffles/mPanTro3/panpan-pt.tsv.gz")
vcf_filtered <- vcf_dat %>%
  filter(! svtype %in% c("BND"), 
         ## remove the ones with imprecise positions
         filter=="PASS") %>%
  ## inversions don't have a svlen field, so I had to manually fill in
  mutate(svlen=ifelse(svtype=="INV", str_count(alt), abs(svlen))) %>%
  ## remove svs that are too long or too short
  filter(svlen >= 50, svlen <= 100000) %>%
  mutate(gt = case_when(gt_gt == "0/0" ~ 0L,
                        gt_gt == "0/1" ~ 1L,
                        gt_gt == "1/0" ~ 1L,
                        gt_gt == "1/1" ~ 2L,),
         svtype = factor(svtype, levels = c("INS", "DEL", "DUP", "INV")),
         ## many of the end coordinates got lost after truvari merge
         ## this may remove too many deletions that intersect with tandem repeats. need to rescue some of these
         end = case_when(svtype=="INS"~pos,
                         svtype=="DEL"~pos+svlen,
                         svtype=="DUP"~end,
                         svtype=="INV"~pos+svlen,)) %>%
  group_by(chrom, pos, end, id, svlen, svtype, ref, alt) %>%
  summarize(allele_count=sum(gt), allele_frequency=allele_count/sample_size_total/2, het_frequency=sum(gt==1)/sample_size_total) %>%
  ungroup()
## remove sites with higher than 90% heterozygote frequency in each species
#filter(het_frequency<=0.9)

filtered_indel_index <- (vcf_raw@fix[,3] %in% (vcf_filtered %>% filter(svtype %in% c("INS", "DEL")) %>% pull(id)))
vcf_indel_raw <- vcf_raw
vcf_indel_raw@fix <- vcf_raw@fix[filtered_indel_index,]
vcf_indel_raw@gt <- vcf_raw@gt[filtered_indel_index,]
write.vcf(vcf_indel_raw, str_c(indir, "/", dataset, ".truvari.merged.missing2ref.filtered.indel.vcf.gz"))

filtered_dupinv_index <- (vcf_raw@fix[,3] %in% (vcf_filtered %>% filter(svtype %in% c("DUP", "INV")) %>% pull(id)))
vcf_dupinv_raw <- vcf_raw
vcf_dupinv_raw@fix <- vcf_raw@fix[filtered_dupinv_index,]
vcf_dupinv_raw@gt <- vcf_raw@gt[filtered_dupinv_index,]
write.vcf(vcf_dupinv_raw, str_c(indir, "/", dataset, ".truvari.merged.missing2ref.filtered.dupinv.vcf.gz"))
