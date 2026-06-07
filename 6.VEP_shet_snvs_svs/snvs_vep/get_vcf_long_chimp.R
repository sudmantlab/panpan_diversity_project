# Set the working directory
setwd("/Users/joanocha/Google Drive/My Drive/POSTDOC/PANPAN/analysis/Figure3_VEP_snvs_svs/snvs_vep/")
# Load all necessary libraries
# The 'tidyverse' package includes dplyr, stringr, readr, and others.
library(vcfR)
library(tidyverse)
library(janitor)
library(data.table)
library(dplyr)

# --- Process Metadata ---
panpan_sample_id <- read_tsv("/global/scratch/users/nicolas931010/sv_detection/docs/sample_table_combined.tsv") %>%
  pull(sample_id) 
panpan_metadata <- read_tsv("/global/scratch/users/nicolas931010/sv_detection/docs/PanPan_Metadata_updated.txt") %>%
  janitor::clean_names() %>%
  distinct(genome_name, species, species_code, common_name, new_name) %>%
  mutate(species_code = ifelse(species_code == "unknown" & genome_name %in% c("AG05253_1"), "Ppan", species_code),
         species_code = ifelse(species_code == "unknown" & genome_name %in% c("AG18356_6", "AG16618_12"), "Ptv x Ptt", species_code),
         species_code = ifelse(species_code == "unknown", "Ptv", species_code)) %>%
  mutate(common_name = ifelse(is.na(common_name) & genome_name %in% c("AG05253_1"), "Bonobo", common_name),
         common_name = ifelse(is.na(common_name) & genome_name %in% c("AG18356_6", "AG16618_12"), "Western x Central hybrid", common_name),
         common_name = ifelse(is.na(common_name), "Western chimpanzee", common_name)) %>%
  filter(genome_name %in% panpan_sample_id) %>%
  distinct(genome_name, common_name)

# --- Read and Tidy VCF Data ---
vcf_raw_chimps <- read.vcfR("/global/scratch/users/joana_rocha/PANPAN/PANPAN_snvs/lifted_pantros_mapped2mPanTro3.sorted.filtered.vcf.gz")
vcf_raw_bonobos <- read.vcfR("/global/scratch/users/joana_rocha/PANPAN/PANPAN_snvs/panpa_mapped2mPanPan1.sorted.filtered.vcf.gz")
vcf_raw_humans <- read.vcfR("/global/scratch/users/joana_rocha/PANPAN/PANPAN_snvs/hprc_mapped2ht2t.sorted.filtered.vcf.gz")


## Handle of 3 species for VEP Join ######
####chimps
vcf <- vcfR2tidy(vcf_raw_chimps, single_frame = TRUE, format_fields = c("GT"))
vcf_dat <- vcf$dat %>%
  clean_names()
vcf_long_chimp_tmp <- vcf_dat %>% # Create the unique variant identifier and name it 'id'.This combines chromosome, position, and alleles to create a robust ID that will match the VEP output format. It replaces the old 'id' column.
  mutate(id = str_c(chrom, "_", pos, "_", ref, "/", alt)) %>% # Calculate allele count, now correctly handling both unphased (e.g., "0/1") and phased (e.g., "0|1", "1|0") genotypes.
  mutate(allele_count = case_when(
    gt_gt %in% c("0/0", "0|0") ~ 0L,
    gt_gt %in% c("0/1", "1/0", "0|1", "1|0") ~ 1L,
    gt_gt %in% c("1/1", "1|1") ~ 2L,
    gt_gt %in% c("./.", ".|.") ~ 0L,
    is.na(gt_gt) ~ 0L  # If gt_gt is NA, allele_count is 0
  ),
  gt_gt = if_else(is.na(gt_gt), "0|0", gt_gt)
  ) %>%
  dplyr::select(id, chrom, pos, ref, alt, indiv, gt_gt, allele_count)
sample_size <- (ncol(vcf_raw_chimps@gt)-1)
vcf_annotated_chimp <- vcf_long_chimp_tmp %>%
  group_by(id, chrom, pos, ref, alt) %>%
  summarize(total_allele_count=sum(allele_count), allele_count_vector=str_c(allele_count, collapse = "")) %>%
  ungroup() %>%
  mutate(allele_frequency=total_allele_count/sample_size/2,
         maf=ifelse(allele_frequency<0.5, allele_frequency, 1-allele_frequency))
write_tsv(vcf_annotated_chimp, "vcf_annotated_chimp.tsv")
#vcf_annotated_chimp  <- vcf_summary_chimp

###bonobos
vcf_bonobo <- vcfR2tidy(vcf_raw_bonobos, single_frame = TRUE, format_fields = c("GT"))
vcf_dat_bonobo <- vcf_bonobo$dat %>%
  clean_names()
vcf_long_bonobos_tmp <- vcf_dat_bonobo %>%
  mutate(id = str_c(chrom, "_", pos, "_", ref, "/", alt)) %>%
  mutate(allele_count = case_when(
      gt_gt %in% c("0/0", "0|0") ~ 0L,
      gt_gt %in% c("0/1", "1/0", "0|1", "1|0") ~ 1L,
      gt_gt %in% c("1/1", "1|1") ~ 2L,
      gt_gt %in% c("./.", ".|.") ~ 0L,
      is.na(gt_gt) ~ 0L  # If gt_gt is NA, allele_count is 0
    ),
    gt_gt = if_else(is.na(gt_gt), "0|0", gt_gt)
  ) %>%
  dplyr::select(id, chrom, pos, ref, alt, indiv, gt_gt, allele_count)

sample_size <- (ncol(vcf_raw_bonobos@gt)-1)
vcf_annotated_bonobo <- vcf_long_bonobos_tmp %>%
  group_by(id, chrom, pos, ref, alt) %>%
  summarize(total_allele_count=sum(allele_count), allele_count_vector=str_c(allele_count, collapse = "")) %>%
  ungroup() %>%
  mutate(allele_frequency=total_allele_count/sample_size/2,
         maf=ifelse(allele_frequency<0.5, allele_frequency, 1-allele_frequency))
#vcf_annotated_bonobo <- vcf_summary_bonobos
write_tsv(vcf_annotated_bonobo, "vcf_annotated_bonobo.tsv")

###humans
vcf_human <- vcfR2tidy(vcf_raw_humans, single_frame = TRUE, format_fields = c("GT"))
vcf_dat_human <- vcf_human$dat %>%
  clean_names()
vcf_long_humans_tmp <- vcf_dat_human %>% # this combines chromosome, position, and alleles to create a robust ID that will match the VEP output format. It replaces the old 'id' column.
  mutate(id = str_c(chrom, "_", pos, "_", ref, "/", alt)) %>% # Calculate allele count, now correctly handling both unphased (e.g., "0/1") and phased (e.g., "0|1", "1|0") genotypes.
  mutate(allele_count = case_when(
    gt_gt %in% c("0/0", "0|0") ~ 0L,
    gt_gt %in% c("0/1", "1/0", "0|1", "1|0") ~ 1L,
    gt_gt %in% c("1/1", "1|1") ~ 2L,
    gt_gt %in% c("./.", ".|.") ~ 0L,
    is.na(gt_gt) ~ 0L  # If gt_gt is NA, allele_count is 0
  ),
  gt_gt = if_else(is.na(gt_gt), "0|0", gt_gt)
  ) %>%
  dplyr::select(id, chrom, pos, ref, alt, indiv, gt_gt, allele_count)
sample_size <- (ncol(vcf_raw_humans@gt)-1)
vcf_annotated_human <- vcf_long_humans_tmp %>%
  group_by(id, chrom, pos, ref, alt) %>%
  summarize(total_allele_count=sum(allele_count), allele_count_vector=str_c(allele_count, collapse = "")) %>%
  ungroup() %>%
  mutate(allele_frequency=total_allele_count/sample_size/2,
         maf=ifelse(allele_frequency<0.5, allele_frequency, 1-allele_frequency))
write_tsv(vcf_annotated_human, "vcf_annotated_human.tsv")


### all humans
vcf_raw_allhumans <- read.vcfR("/global/scratch/users/joana_rocha/PANPAN/PANPAN_snvs/hprc_hgsvc_mapped2ht2t.sorted.filtered.vcf.gz")
fix_df <- as_tibble(vcf_raw_allhumans@fix)
gt_df  <- as_tibble(vcf_raw_allhumans@gt)
chunk_size <- 1000000 # Process 1 million variants at a time
num_variants <- nrow(fix_df)
starts <- seq(1, num_variants, by = chunk_size)
results_list <- vector("list", length(starts))
sample_size <- (ncol(vcf_raw_allhumans@gt) - 1)
for (i in seq_along(starts)) {
  start_row <- starts[i]
  end_row <- min(start_row + chunk_size - 1, num_variants)
  cat(sprintf("Processing chunk %d of %d: variants %d to %d\n", i, length(starts), start_row, end_row))
  fix_chunk <- fix_df[start_row:end_row, ]
  gt_chunk  <- gt_df[start_row:end_row, ]
  chunk_long <- bind_cols(fix_chunk, gt_chunk) %>%
    pivot_longer(
      cols = -(CHROM:FORMAT),
      names_to = "indiv",
      values_to = "gt_gt"
    )
  chunk_summary <- chunk_long %>%
    mutate(
      # Correctly handle NAs, just like in your first script
      allele_count = case_when(
        gt_gt %in% c("0/0", "0|0") ~ 0L,
        gt_gt %in% c("0/1", "1/0", "0|1", "1|0") ~ 1L,
        gt_gt %in% c("1/1", "1|1") ~ 2L,
        gt_gt %in% c("./.", ".|.") ~ 0L,
        is.na(gt_gt) ~ 0L  # Treat NA as homozygous reference
      )
    ) %>%
    group_by(CHROM, POS, REF, ALT) %>%
    summarize(
      # Now that NAs are 0, we don't need na.rm = TRUE
      total_allele_count = sum(allele_count),
      # This also gets simpler because allele_count has no NAs
      allele_count_vector = str_c(allele_count, collapse = ""),
      .groups = 'drop'
    )
  
  results_list[[i]] <- chunk_summary
}
vcf_annotated_human <- bind_rows(results_list)
vcf_annotated_human <- vcf_annotated_human %>%
  mutate(
    id = str_c(CHROM, "_", POS, "_", REF, "/", ALT),
    allele_frequency = total_allele_count / (sample_size * 2),
    maf = ifelse(allele_frequency < 0.5, allele_frequency, 1 - allele_frequency)
  ) %>%
  select(
    id, 
    chrom = CHROM, 
    pos = POS, 
    ref = REF, 
    alt = ALT, 
    total_allele_count, 
    allele_count_vector, 
    allele_frequency, 
    maf
  )
write_tsv(vcf_annotated_human, "vcf_annotated_allhuman.tsv")

head(vcf_annotated_bonobo)
head(vcf_annotated_chimp)
head(vcf_annotated_human)


vcf_annotated_combined <- bind_rows(vcf_annotated_chimp %>% mutate(species="chimpanzee"),
                                    vcf_annotated_bonobo %>% mutate(species="bonobo"),
                                    vcf_annotated_human %>% mutate(species="human"))
write_tsv(vcf_annotated_combined, "vcf_annotated_combined_allhuman.tsv")

# ---Object needed for VEP within chimp analysis ---
vcf_for_vep_join <- vcf_annotated_combined %>%
  dplyr::select(id, chrom, pos, total_allele_count, allele_frequency, species)
write_tsv(vcf_for_vep_join, "vcf_for_vep_join_allhuman.tsv")

# ---Join Data and Finalize Table to get vcf_long_chimp and long_chimp_by_species needed for VEP within chimp analysis ---
vcf_long_chimp <- vcf_long_chimp_tmp %>%
  left_join(panpan_metadata, by = c("indiv" = "genome_name")) %>%
  mutate(common_name = str_remove(common_name, " chimpanzee"),
         common_name = str_remove(common_name, " hybrid"),
         # Re-level the factor for correct ordering in plots/summaries.
         common_name = fct_relevel(common_name, c("Western", "Western x Central", "Central", "Eastern")))
write_tsv(vcf_long_chimp, "mPanTro3/vcf_long_chimp.tsv")

#vcf_long_chimp <- read_tsv("mPanTro3/vcf_long_chimp.tsv.gz") %>%
#  mutate(common_name=fct_relevel(common_name, c("Western", "Western x Central", "Central", "Eastern")))

vcf_long_chimp_by_species <- vcf_long_chimp %>%
  mutate(common_name=fct_relevel(common_name, c("Western", "Western x Central", "Central", "Eastern")))  %>%
  group_by(id, chrom, pos, ref, alt, common_name) %>%
  summarise(allele_count=sum(allele_count), sample_size=dplyr::n()) %>%
  ungroup()
write_tsv(vcf_long_chimp_by_species, "mPanTro3/vcf_long_chimp_by_species.tsv")



