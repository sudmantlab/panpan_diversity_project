#!/usr/bin/env Rscript

## Load required packages
library(tidyverse)
library(plyranges)
## Read in the arguments
args = commandArgs(trailingOnly=TRUE)
indir <- args[1]
dataset <- args[2]

vcf_summary_wide <- read_tsv(str_c(indir, "/", dataset, ".sniffles_svim-asm.truvari.sorted.merged.summary.tsv")) %>%
  mutate(svtype = factor(svtype, levels = c("INS", "DEL", "DUP", "INV")))
sample_size_total <- vcf_summary_wide$allele_count_vector_sniffles[1] %>% str_count()
sv_id_table <- read_tsv(str_c(indir, "/", dataset, ".sniffles_svim-asm.truvari.sorted.merged.id.tsv"))
sv_categories <-  read_table(str_c(indir, "/", dataset, ".sniffles_svim-asm.truvari.sorted.merged.indel.fasta.out"), skip=3, col_names = FALSE) %>%
  transmute(motif=X10, type=X11) %>%
  count(motif, type) %>%
  group_by(motif) %>%
  slice_max(order_by=n, n=1, with_ties = FALSE) %>%
  dplyr::select(-n)
gff_original <- read_tsv(str_c(indir, "/", dataset, ".sniffles_svim-asm.truvari.sorted.merged.indel.fasta.out.gff"), comment = "#", col_names = FALSE) %>%
  separate(X9, into = c("X9", "X10", "X11", "X12"), sep = "[ \"\t]+") %>%
  mutate(motif=str_remove(X10, "Motif:"), divergence=X6) %>%
  left_join(sv_categories, by="motif")
gff_top_type <- gff_original %>%
  transmute(seqnames=X1, start=X4, end=X5, divergence, type) %>%
  mutate(type=case_when(type %in% c("Low_complexity", "Simple_repeat") ~ "Simple/tandem repeat",
                        str_detect(type, "^Satellite") ~ "Simple/tandem repeat",
                        TRUE ~ type)) %>%
  as_granges() %>%
  group_by(type) %>%
  reduce_ranges(divergence=mean(divergence)) %>%
  ungroup() %>%
  as_tibble() %>% 
  group_by(seqnames, type) %>%
  summarise(length=sum(width), divergence=sum(divergence*width)/sum(width)) %>%
  ungroup() %>%
  group_by(seqnames) %>%
  slice_max(length, n=1, with_ties = FALSE) %>%
  dplyr::rename(id=seqnames, length_top_type=length) %>%
  ungroup()
gff_all_types <- gff_original %>%
  transmute(seqnames=X1, start=X4, end=X5, divergence, type) %>%
  as_granges() %>%
  reduce_ranges() %>%
  ungroup() %>%
  as_tibble() %>%
  group_by(seqnames) %>%
  summarise(length_annotated=sum(width)) %>%
  ungroup() %>%
  dplyr::rename(id=seqnames)
gff <- left_join(gff_top_type, gff_all_types)
annotated_indel_complete <- vcf_summary_wide %>%
  mutate(svlen=abs(svlen),
         ## note that recording the length of the deletion/insertion sequence was useful because they are sometimes different from the svlen (e.g. DUPs for svimasm, variable breakpoints when merged by sniffles)
         seqlen=ifelse(svtype=="DEL", str_count(ref), str_count(alt))) %>%
  filter(svtype %in% c("INS", "DEL"), seqlen >=50, seqlen <= 100000) %>%
  dplyr::select(chrom, pos, id, svlen, seqlen, svtype, method_support, intersecting_tr) %>%
  distinct() %>%
  left_join(gff, by="id") %>%
  mutate(type=ifelse(is.na(type), "Non-repeat", type)) 
annotated_indel <- annotated_indel_complete %>%
  mutate(p=length_top_type/seqlen, p_annotated=length_annotated/seqlen) %>%
  mutate(type_simplified=case_when(p_annotated < 0.8 ~ "Non-repeat",
                        p/p_annotated < 0.8 ~ "Non-repeat",
                        type %in% c("SINE/Alu", "LINE/L1", "Retroposon/SVA", "Non-repeat", "Simple/tandem repeat") ~ type,
                        str_detect(type, "^SINE|^LINE|^LTR") ~ "Other Type I TE",
                        str_detect(type, "^DNA|^RC") ~ "Type II TE",
                        str_detect(type, "RNA$") ~ "RNA array",
                        type == "unknown" ~ "Unknown", 
                        TRUE ~ type)) %>%
  mutate(type_simplified=fct_relevel(type_simplified, c("SINE/Alu", "LINE/L1", "Retroposon/SVA", "Other Type I TE", "Type II TE", "Simple/tandem repeat", "RNA array", "Unknown", "Non-repeat")))
final_id_set <- annotated_indel %>%
  filter((! intersecting_tr) | (! type_simplified %in% c("Simple/tandem repeat", "Unknown", "Non-repeat"))) %>%
  bind_rows(vcf_summary_wide %>% filter(svtype %in% c("DUP", "INV"))) %>%
  dplyr::select(id, svtype, method_support)
vcf_summary_wide_final <- vcf_summary_wide %>%
  filter(id %in% final_id_set$id)
final_id_table <- sv_id_table %>%
  inner_join(final_id_set, by="id")

final_id_set %>%
  write_tsv(str_c(indir, "/", dataset, ".sniffles_svim-asm.truvari.sorted.merged.final_id_set.tsv"))
final_id_table %>%
  write_tsv(str_c(indir, "/", dataset, ".sniffles_svim-asm.truvari.sorted.merged.final_id_table.tsv"))
annotated_indel %>%
  write_tsv(str_c(indir, "/", dataset, ".sniffles_svim-asm.truvari.sorted.merged.annotated_indel.tsv"))
final_id_table %>%
  filter(method_support == "consensus", caller == "svimasm") %>%
  dplyr::select(gt_id) %>%
  write_tsv(str_c(indir, "/", dataset, ".svimasm.consensus.id"), col_names=FALSE)
final_id_table %>%
  filter(method_support == "consensus", caller == "sniffles") %>%
  dplyr::select(gt_id) %>%
  write_tsv(str_c(indir, "/", dataset, ".sniffles.consensus.id"), col_names=FALSE)
final_id_table %>%
  filter(method_support == "svimasm_only", caller == "svimasm") %>%
  dplyr::select(gt_id) %>%
  write_tsv(str_c(indir, "/", dataset, ".svimasm.only.id"), col_names=FALSE)
final_id_table %>%
  filter(method_support == "sniffles_only", caller == "sniffles") %>%
  dplyr::select(gt_id) %>%
  write_tsv(str_c(indir, "/", dataset, ".sniffles.only.id"), col_names=FALSE)
