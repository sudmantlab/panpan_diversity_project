setwd("/Users/joanocha/Google Drive/My Drive/POSTDOC/PANPAN/analysis/Figure4_SINGER/chimps/")
library(dplyr)
library(ggplot2)
library(scales)
library(readr)
library(dplyr)
library(ggrepel)
library(tidyr)
library(patchwork)
library(valr)
library(readr)
library(glue)
library(tidyverse)
library(tidytext)
library(data.table)
library(stringr)
library(ggrastr)
library(ggnewscale) 
library(UpSetR)
library(cowplot)
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")

#bonobo "#FF8C69"
#chimp "#8A865D"
# human "#5C608A"

# Load and transform data
coalescence_df<- read_csv("genome-wide_metrics_annotatedbyPop_CenSat_removed.csv")%>%
  rename(chrom = chromosome)

annotation_df <- fread("mPanTro3_gene.bed.gz") %>%
  rename(genes = name) %>%
  select(genes, gene_biotype, description) %>%
  distinct()

full_annotation_df <- fread("mPanTro3_gene.bed.gz") 

gene_map <- read_csv("chimps_outliers_renamed.csv")
View(gene_map)

conflicts_prefer(dplyr::first)
coalescence_df <- coalescence_df %>%
  separate_rows(genes, sep = ",") %>%
  left_join(annotation_df, by = "genes") %>%
  group_by(chrom, start, end, population) %>%
  summarise(
    avg_tmrca = first(avg_tmrca),
    T_pooled = first(T_pooled),
    T_within = first(T_within),
    Tpooled_Twithin_ratio = first(Tpooled_Twithin_ratio),
    genes = paste(genes, collapse = ","),
    gene_biotypes = paste(na.omit(gene_biotype), collapse = ";"),
    descriptions = paste(na.omit(description), collapse = ";")
  ) %>%
  ungroup()

coalescence_df$avg_tmrca_myr <- coalescence_df$avg_tmrca * 25 / 1e6
coalescence_df$T_within_myr <- coalescence_df$T_within * 25 / 1e6
coalescence_df$T_pooled_myr <- coalescence_df$T_pooled * 25 / 1e6


features_df <- read_tsv(
  "mPanTro3_genomefeatures.bed",
  col_names = c("chrom", "start", "end", "feature"),
  col_types = "cddc"
)

tr_df <- read_tsv(
  "pantro_catalog.no_overlaps_simp.bed",
  col_names = c("chrom", "start", "end", "motif"),
  col_types = "cddc"
)

sedef_df <- read_tsv(
  "mPanTro3_sedefSegDups.bed",
  col_names = c("chrom", "start", "end", "coordinates"),
  col_types = "cddc",
  col_select = 1:4
)

censat_df <- read_tsv(
  "mPanTro3_CenSat.bed",
  col_names = c("chrom", "start", "end", "censat"),
  col_types = "cddc",
  col_select = 1:4
)

rmsk_df <- read_tsv(
  "mPanTro3_RepeatMasker.bed",
  col_names = c("chrom", "start", "end", "feature"),
  col_types = "cddc",
  col_select = 1:4
)

SRA_df <- read_tsv(
  "mPanTro3_SR_mask.bed",
  col_names = c("chrom", "start", "end"),
  col_types = "cddc",
  col_select = 1:3
)

# --- 2. Create Unique Windows & Calculate Overlaps ---
unique_windows_df <- coalescence_df %>%
  select(chrom, start, end) %>%
  distinct() # Create a dataframe with only unique windows to avoid multiplying overlaps.

all_features_df <- bind_rows(
  features_df,
  mutate(tr_df, feature = "TR"),
  mutate(sedef_df, feature = "SEDEF_SD"),
  mutate(censat_df, feature = "CENSAT"),
  mutate(SRA_df, feature = "SRA"),
  rmsk_df
) %>%
  select(chrom, start, end, feature)

feature_annotations <- bed_map(
  unique_windows_df, # Use unique windows
  all_features_df,
  feature_annotation = paste(unique(sort(feature)), collapse = ",")
) %>%
  select(chrom, start, end, feature_annotation)


features_overlap_summary <- bed_intersect(unique_windows_df, features_df) %>%
  group_by(chrom, start.x, end.x) %>%
  summarise(features_overlap_bp = sum(.overlap), .groups = "drop") %>%
  rename(start = start.x, end = end.x)


censat_overlap_summary <- bed_intersect(unique_windows_df, censat_df) %>%
  group_by(chrom, start.x, end.x) %>%
  summarise(censat_overlap_bp = sum(.overlap), .groups = "drop") %>%
  rename(start = start.x, end = end.x)

tr_overlap_summary <- bed_intersect(unique_windows_df, tr_df) %>% # Use unique windows
  group_by(chrom, start.x, end.x) %>%
  summarise(tr_overlap_bp = sum(.overlap), .groups = "drop") %>%
  rename(start = start.x, end = end.x)

sedef_overlap_summary <- bed_intersect(unique_windows_df, sedef_df) %>% # Use unique windows
  group_by(chrom, start.x, end.x) %>%
  summarise(sedef_overlap_bp = sum(.overlap), .groups = "drop") %>%
  rename(start = start.x, end = end.x)

rmsk_overlap_summary <- bed_intersect(unique_windows_df, rmsk_df) %>%
  group_by(chrom, start.x, end.x) %>%
  summarise(rmsk_overlap_bp = sum(.overlap), .groups = "drop") %>%
  rename(start = start.x, end = end.x)

SRA_overlap_summary <- bed_intersect(unique_windows_df, SRA_df) %>%
  group_by(chrom, start.x, end.x) %>%
  summarise(SRA_overlap_bp = sum(.overlap), .groups = "drop") %>%
  rename(start = start.x, end = end.x)

# --- 3. Join Results Back to the ORIGINAL Dataframe and preprocess data and FILTER OUT CHR X/Y---
coalescence_df_with_features <- coalescence_df %>%
  left_join(feature_annotations, by = c("chrom", "start", "end")) %>%
  left_join(tr_overlap_summary, by = c("chrom", "start", "end")) %>%
  left_join(sedef_overlap_summary, by = c("chrom", "start", "end")) %>%
  left_join(censat_overlap_summary, by = c("chrom", "start", "end")) %>%
  left_join(rmsk_overlap_summary, by = c("chrom", "start", "end")) %>%
  left_join(features_overlap_summary, by = c("chrom", "start", "end")) %>%
  left_join(SRA_overlap_summary, by = c("chrom", "start", "end")) %>%
  mutate(
    feature_annotation = if_else(is.na(feature_annotation), ".", feature_annotation),
    tr_overlap_bp = coalesce(tr_overlap_bp, 0L),
    sedef_overlap_bp = coalesce(sedef_overlap_bp, 0L),
    censat_overlap_bp = coalesce(censat_overlap_bp, 0L),
    rmsk_overlap_bp = coalesce(rmsk_overlap_bp, 0L),
    SRA_overlap_bp = coalesce(SRA_overlap_bp, 0L),
    features_overlap_bp = coalesce(features_overlap_bp, 0L),
    window_size = end - start,
    tr_overlap_percent = (tr_overlap_bp / window_size) * 100,
    sedef_overlap_percent = (sedef_overlap_bp / window_size) * 100,
    censat_overlap_percent = (censat_overlap_bp / window_size) * 100,
    rmsk_overlap_percent = (rmsk_overlap_bp / window_size) * 100,
    features_overlap_percent = (features_overlap_bp / window_size) * 100,
    SRA_overlap_percent = (SRA_overlap_bp / window_size) * 100
  ) %>%
  select(
    chrom, start, end, population, avg_tmrca_myr, T_pooled_myr, T_within_myr, Tpooled_Twithin_ratio, genes, feature_annotation, gene_biotypes, descriptions,
    tr_overlap_percent, 
    sedef_overlap_percent,
    censat_overlap_percent,
    rmsk_overlap_percent,
    features_overlap_percent,
    SRA_overlap_percent
  ) %>%
  rename(chromosome = chrom) %>%
  mutate(
    chr_clean = sub("_.*", "", chromosome)
  ) %>%
  filter(!chr_clean %in% c("chrX", "chrY")) %>% 
  mutate(
    chr_clean = factor(chr_clean, levels = paste0("chr", 1:23)),
    midpoint = (start + end) / 2,
    position_mb = midpoint / 1e6,
    chr_num = as.integer(gsub("chr", "", chr_clean))
  ) %>%
  arrange(chr_clean, start)



############ SUMMARY STATS ON TRs and other features #####################
unique_windows_summary <- coalescence_df_with_features %>%
  group_by(chromosome, start, end, genes, feature_annotation) %>%
  summarise(
    # These overlap values are the same for all rows of a unique window
    tr_overlap_percent = first(tr_overlap_percent),
    sedef_overlap_percent = first(sedef_overlap_percent),
    censat_overlap_percent = first(censat_overlap_percent),
    rmsk_overlap_percent = first(rmsk_overlap_percent),
    SRA_overlap_percent = first(SRA_overlap_percent),
    T_within_myr = mean(T_within_myr, na.rm = TRUE),
    Tpooled_Twithin_ratio = mean(Tpooled_Twithin_ratio, na.rm = TRUE),
  ) %>%
  ungroup() # Don't forget to ungroup
#View(unique_windows_summary)

features_summary<-unique_windows_summary %>%
  pivot_longer(
    cols = c(tr_overlap_percent, SRA_overlap_percent),
    names_to = "category",
    values_to = "percent"
  ) %>%
  filter(percent > 0) %>%
  ggplot(aes(x = percent, fill = category)) +
  geom_histogram(bins = 50, color = "white", show.legend = FALSE) +
  facet_wrap(~ category, scales = "free") +
  labs(
    title = "Distribution of Overlap with TMRCA Windows",
    x = "Overlap Percentage with 1kb Window (%)",
    y = "Number of Windows in Chimpanzees"
  ) +
  theme_minimal()
ggsave("features_summary.pdf", features_summary, width=6, height=3)  



############## FILTER FOR COMPLEX NOISY REGIONS ###########

# Calculate chromosome offsets with chr_len included
chrom_offsets <- coalescence_df_with_features %>%
  group_by(chr_clean) %>%
  summarize(chr_len = max(end), .groups = "drop") %>%
  arrange(chr_clean) %>%
  # The as.numeric() here prevents the integer overflow
  mutate(offset = cumsum(lag(as.numeric(chr_len), default = 0))) %>%
  select(chr_clean, chr_len, offset)

coalescence_df_with_features<- coalescence_df_with_features %>%
  left_join(chrom_offsets, by = "chr_clean") %>%
  mutate(genome_pos = start + offset,
         genome_pos_mb = genome_pos / 1e6)



unique(coalescence_df_with_features$feature_annotation)
#pattern_to_exclude <- "SAT|CEN"
coalescence_df_with_features2 <- coalescence_df_with_features %>%
  filter(
    censat_overlap_percent == 0,
    tr_overlap_percent <= 5,
    SRA_overlap_percent > 20,
    #sedef_overlap_percent == 0,
    !str_detect(feature_annotation, "Satellite"),
    !str_detect(feature_annotation, "Low_complexity"),
    !str_detect(feature_annotation, "Cen"),
    #!str_detect(feature_annotation, "TR"),
    !str_detect(feature_annotation, "CENSAT"),
    !str_detect(feature_annotation, "Gap"),
    #str_detect(feature_annotation, "SRA"),
    ) #%>%
  #mutate(
    #gene_feature = glue("{genes} ({feature_annotation}) (SRA: {round(SRA_overlap_percent, 2)}%)")
 # )




# Check the unique annotations in your new dataframe
unique(coalescence_df_with_features2$feature_annotation)

### plot distribution
unique_tmrca_summary <- coalescence_df_with_features2 %>%
  group_by(chromosome, start, end, genes) %>%
  summarise(
    avg_tmrca_myr = mean(avg_tmrca_myr, na.rm = TRUE),
    T_pooled_myr = mean(T_pooled_myr, na.rm = TRUE),
    T_within_myr = mean(T_within_myr, na.rm = TRUE),
    Tpooled_Twithin_ratio = mean(Tpooled_Twithin_ratio, na.rm = TRUE),
    .groups = "drop" # Using .groups is modern practice
  )

summary_lines <- unique_tmrca_summary %>%
  pivot_longer(
    cols = c(T_pooled_myr, avg_tmrca_myr),
    names_to = "statistic",
    values_to = "myr"
  ) %>%
  group_by(statistic) %>%
  summarise(
    mean_val = mean(myr, na.rm = TRUE)
  )

# --- 3. Build the plot ---
tmrca_plot <- unique_tmrca_summary %>%
  pivot_longer(
    cols = c(T_pooled_myr, avg_tmrca_myr),
    names_to = "statistic",
    values_to = "myr"
  ) %>%
  ggplot(aes(x = myr, fill = statistic)) + 
  geom_histogram(bins = 50, alpha = 0.6, position = "identity", color = "white") +
  geom_vline(
    data = summary_lines, 
    aes(xintercept = mean_val, color = statistic), 
    linetype = "dashed", 
    linewidth = 1,
    show.legend = FALSE # Hide duplicate legend for the lines
  ) +
  geom_text(
    data = summary_lines,
    aes(
      x = mean_val, 
      label = paste0("mean = ", round(mean_val, 2)), 
      color = statistic
    ),
    y = Inf,        # Position text at the top
    vjust = 1.5,      # Nudge down from the top
    hjust = "left",   # Align text
    nudge_x = 0.1,    # Nudge right from the line
    size = 3.5,
    show.legend = FALSE # Hide duplicate legend for the text
  ) +
  scale_x_continuous(breaks = seq(0, 20, by = 2)) +
  labs(
    x = "Chimpanzee TMRCA in 1kb windows (Mya)",
    y = "Number of windows",
    fill = "Statistic" # Optional: clean up legend title
  ) +
  theme_cowplot() +
  theme(legend.position = "bottom") +
  coord_cartesian(xlim = c(NA, 20), clip = "off")
tmrca_plot
ggsave("tmrca_summary.pdf", tmrca_plot, width = 5, height = 3)


tmrca_plot2 <-  unique_tmrca_summary %>%
  pivot_longer(
    cols = Tpooled_Twithin_ratio,
    names_to = "statistic",
    values_to = "myr"
  ) %>%
  ggplot(aes(x = myr, fill = statistic)) + 
  scale_x_continuous(breaks = seq(0, 20, by = 2)) +
  geom_histogram(bins = 50000, alpha = 0.6, position = "identity", color = "white", fill="darkgreen") +
  labs(
    x = "Chimpanzee Tpoooled/Twithin",
    y = "Number of windows",
    fill = "Statistic" # Optional: clean up legend title
  ) +
  theme_cowplot() +
  theme(legend.position = "bottom") + 
  coord_cartesian(xlim = c(NA, 20), clip = "off") 
tmrca_plot2
ggsave("tmrca_summary2.pdf", tmrca_plot2, width = 5, height = 3)


#################### GET TABLES OF TOP GENES ##############


#### ALL CHIMPS ##############
ranked_T_pooled_myr <- coalescence_df_with_features %>%
  distinct(chromosome, start, end, .keep_all = TRUE) %>%
  filter(T_pooled_myr > 6) %>%
  arrange(desc(T_pooled_myr)) %>%
  select(chromosome, start, end, T_pooled_myr, genes, descriptions,  gene_biotypes, SRA_overlap_percent, tr_overlap_percent, censat_overlap_percent, rmsk_overlap_percent, feature_annotation)
write_csv(ranked_T_pooled_myr, "ranked_T_pooled_myr.csv")
wide_ranked_genes_T_pooled_myr <- ranked_T_pooled_myr %>%
  separate_rows(genes, sep = "[,;]") %>%
  mutate(genes = trimws(genes)) %>%
  filter(genes != "" & !is.na(genes) & genes != ".") %>%
  group_by(genes) %>%
  slice_max(order_by = T_pooled_myr, n = 1, with_ties = FALSE) %>%
  filter(!str_detect(descriptions, "uncharacterized"))  %>%
  ungroup() %>%
  arrange(desc(T_pooled_myr))
View(wide_ranked_genes_T_pooled_myr)
write_csv(wide_ranked_genes_T_pooled_myr, "unfiltered_wide_ranked_genes_T_pooled_myr.csv")


get_pop_significant_genes <- function(data, metric_var, threshold_quantile) {
  
  # --- 1. Calculate the threshold for EACH population ---
  pop_thresholds <- data %>%
    group_by(population) %>%
    summarise(
      threshold = quantile({{ metric_var }}, probs = threshold_quantile, na.rm = TRUE),
      .groups = "drop"
    )
  
  # --- 2. Find significant genes using the per-population thresholds ---
  significant_genes <- data %>%
    left_join(pop_thresholds, by = "population") %>%
    filter({{ metric_var }} > threshold, genes != ".") %>%
    
    # ** THE FIX: Process each row individually to handle mismatches **
    rowwise() %>%
    mutate(
      # a. Split each column's string into a list of items
      genes_list = str_split(genes, "[,;]"),
      biotypes_list = str_split(gene_biotypes, "[,;]"),
      desc_list = str_split(descriptions, "[,;]"),
      
      # b. Use the number of genes as the "master" length
      master_len = length(genes_list[[1]]),
      
      # c. Pad the shorter lists with "NA" to match the master length
      biotypes_padded = list(c(biotypes_list[[1]], rep("NA", master_len - length(biotypes_list[[1]])))),
      desc_padded = list(c(desc_list[[1]], rep("NA", master_len - length(desc_list[[1]])))),
      
      # d. Paste the now-parallel lists back into corrected strings
      genes = paste(genes_list[[1]], collapse = ","),
      gene_biotypes = paste(biotypes_padded[[1]], collapse = ";"),
      descriptions = paste(desc_padded[[1]], collapse = ";")
    ) %>%
    ungroup() %>% # End the row-wise operation
    
    # This step will now work correctly with the padded strings
    separate_rows(genes, gene_biotypes, descriptions, sep = "[,;]") %>%
    
    mutate(genes = trimws(genes)) %>%
    filter(genes != "", genes != ".") %>%
    group_by(genes, population) %>%
    slice_max({{ metric_var }}, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(
      chromosome, 
      start, 
      end,
      population, 
      avg_tmrca_myr,
      T_within_myr,
      Tpooled_Twithin_ratio,
      genes,
      gene_biotypes, 
      descriptions,
      tr_overlap_percent,
      SRA_overlap_percent
    )
  
  return(significant_genes)
}

get_pop_significant_genes <- function(data, metric_var, threshold_value) {
  
  significant_genes <- data %>%
    filter({{ metric_var }} > threshold_value, genes != ".") %>%
    
    rowwise() %>%
    mutate(
      genes_list = str_split(genes, "[,;]"),
      biotypes_list = str_split(gene_biotypes, "[,;]"),
      desc_list = str_split(descriptions, "[,;]"),
      master_len = length(genes_list[[1]]),
      biotypes_padded = list(c(biotypes_list[[1]], rep("NA", master_len - length(biotypes_list[[1]])))),
      desc_padded = list(c(desc_list[[1]], rep("NA", master_len - length(desc_list[[1]])))),
      genes = paste(genes_list[[1]], collapse = ","),
      gene_biotypes = paste(biotypes_padded[[1]], collapse = ";"),
      descriptions = paste(desc_padded[[1]], collapse = ";")
    ) %>%
    ungroup() %>%
    
    separate_rows(genes, gene_biotypes, descriptions, sep = "[,;]") %>%
    mutate(genes = trimws(genes)) %>%
    filter(genes != "", genes != ".") %>%
    group_by(genes, population) %>%
    slice_max({{ metric_var }}, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(
      chromosome, start, end, population, genes,
      gene_biotypes, descriptions, T_within_myr, Tpooled_Twithin_ratio
    )
  
  return(significant_genes)
}


quantile_to_use <- 0.9990
significant_genes_Twithin <- get_pop_significant_genes(
  data = coalescence_df_with_features2 %>%
  filter(!population %in% c("hybrid", "Western x Central hybrid")), 
  metric_var = T_within_myr,
  threshold_quantile = quantile_to_use
) %>%
  arrange(desc(T_within_myr))


significant_genes_Twithin <- get_pop_significant_genes(
  data = coalescence_df_with_features2 %>%
    filter(!population %in% c("hybrid", "Western x Central hybrid")), 
  metric_var = T_within_myr,
  threshold_quantile = 6
) %>%
  arrange(desc(T_within_myr))



significant_genes_avgTMRCA <- get_pop_significant_genes(
  data = coalescence_df_with_features %>%
    filter(!population %in% c("hybrid", "Western x Central hybrid")), 
  metric_var = avg_tmrca_myr,
  threshold_quantile = quantile_to_use
) %>%
  arrange(desc(avg_tmrca_myr))


significant_genes_ratio <- get_pop_significant_genes(
  data = coalescence_df_with_features %>%
    filter(!population %in% c("hybrid", "Western x Central hybrid")), 
  metric_var = Tpooled_Twithin_ratio,
  threshold_quantile = quantile_to_use
) %>%
  arrange(desc(Tpooled_Twithin_ratio))


# Convert the 'significant_genes_Twithin' dataframe into a list format where names are populations and values are vectors of gene names.
upset_list_Twithin <- split(significant_genes_Twithin$genes, significant_genes_Twithin$population)
#upset(
#  fromList(upset_list_Twithin), 
#  nsets = 5, # Show the 5 populations
#  nintersects = 20, # Show the 20 most frequent intersections
#  order.by = "freq", # Order the intersections by frequency
#  text.scale = 1.5,
##  point.size = 3,
#  line.size = 1
#)
# Reshape the significant genes table from long to wide format
wide_significant_genes_Twithin <- significant_genes_Twithin %>%
  # First, for each gene, count how many populations it appears in
  group_by(genes) %>%
  mutate(
    num_populations = n(),
    is_shared = ifelse(n() > 1, "Yes", "No") # Add a simple Yes/No column
  ) %>%
  ungroup() %>%
  pivot_wider(
    id_cols = c(genes, descriptions, num_populations, is_shared), 
    names_from = population,
    values_from = T_within_myr 
  ) %>%
  # Arrange by the most shared genes first, then by the number of populations
  arrange(desc(is_shared), desc(num_populations))

# View the top of the new wide-format table
View(wide_significant_genes_Twithin)
write_csv(wide_significant_genes_Twithin, "unfiltered_chimps_nohybrid_wide_significant_genes_Twithin.csv")


wide_significant_genes_ratio <- significant_genes_ratio %>%
  group_by(genes) %>%
  mutate(
    num_populations = n(),
    is_shared = ifelse(n() > 1, "Yes", "No") 
  ) %>%
  ungroup() %>%
  pivot_wider(
    id_cols = c(genes, descriptions, num_populations, is_shared), 
    names_from = population,
    values_from = Tpooled_Twithin_ratio 
  ) %>%
  arrange(desc(is_shared), desc(num_populations))
View(wide_significant_genes_ratio)
write_csv(wide_significant_genes_ratio, "unfiltered_chimps_nohybrid_wide_significant_genes_ratio.csv")


plot_data_ranked_ratio <- coalescence_df_with_features2 %>% 
  filter(!population %in% c("hybrid", "Western x Central hybrid")) %>%
  arrange(desc(Tpooled_Twithin_ratio)) %>%
  filter(Tpooled_Twithin_ratio > 6)  %>%
  select(chromosome, start, end, population, Tpooled_Twithin_ratio, genes, gene_biotypes, descriptions)
write_csv(plot_data_ranked_ratio, "ranked_ratio.csv")

plot_data_ranked_tmrca <- coalescence_df_with_features2 %>%
  filter(!population %in% c("hybrid", "Western x Central hybrid")) %>%
  arrange(desc(T_within_myr)) %>%
  filter(T_within_myr > 6) %>% 
  select(chromosome, start, end, population, T_within_myr, genes, gene_biotypes, descriptions)
write_csv(plot_data_ranked_tmrca, "ranked_tmrca.csv")



#################### GET TABLES OF TOP GENES PER SUBSPECIES ##############
FIXED_THRESHOLD <- 6
PERCENTILE_THRESHOLD <- 0.9999

pop_specific_thresholds <- coalescence_df_with_features2 %>%
  filter(population %in% c("Central", "Eastern", "Western")) %>%
  group_by(population) %>%
  summarise(
    Twithin_perc_thresh = quantile(T_within_myr, probs = PERCENTILE_THRESHOLD, na.rm = TRUE),
    Tratio_perc_thresh = quantile(Tpooled_Twithin_ratio, probs = PERCENTILE_THRESHOLD, na.rm = TRUE),
    .groups = "drop"
  )

#population Twithin_perc_thresh Tratio_perc_thresh
# Central                  12.2             6617  
# Eastern                  10.2              232  
# Western                   8.28               2.35

# Outliers for T_within_myr (> 99.99th percentile)
outliers_Twithin_percentile <- coalescence_df_with_features2 %>%
  filter(population %in% c("Central", "Eastern", "Western")) %>%
  left_join(pop_specific_thresholds, by = "population") %>%
  filter(T_within_myr > Twithin_perc_thresh) %>%
  select(-ends_with("_perc_thresh")) %>% # Removes the temporary threshold columns
  arrange(population, desc(T_within_myr))

# Outliers for Tpooled_Twithin_ratio (> 99.99th percentile)
outliers_Tratio_percentile <- coalescence_df_with_features2 %>%
  filter(population %in% c("Central", "Eastern", "Western")) %>%
  left_join(pop_specific_thresholds, by = "population") %>%
  filter(Tpooled_Twithin_ratio > Tratio_perc_thresh) %>%
  select(-ends_with("_perc_thresh")) %>% # Removes the temporary threshold columns
  arrange(population, desc(Tpooled_Twithin_ratio))

# Outliers for T_within_myr (> 6)
outliers_Twithin_fixed <- coalescence_df_with_features2 %>%
  filter(population %in% c("Central", "Eastern", "Western")) %>%
  filter(T_within_myr > FIXED_THRESHOLD) %>%
  arrange(population, desc(T_within_myr))

# Outliers for Tpooled_Twithin_ratio (> 6)
outliers_Tratio_fixed <- coalescence_df_with_features2 %>%
  filter(population %in% c("Central", "Eastern", "Western")) %>%
  filter(Tpooled_Twithin_ratio > FIXED_THRESHOLD) %>%
  arrange(population, desc(Tpooled_Twithin_ratio))


create_wide_gene_summary <- function(outlier_df, metric_col) {
  metric_col <- enquo(metric_col)
  gene_summary <- outlier_df %>%
    separate_rows(genes, sep = "[,;]") %>%
    mutate(genes = trimws(genes)) %>%
    filter(genes != "" & !is.na(genes) & genes != ".") %>%
    group_by(genes, population) %>%
    summarise(
      max_metric = max(!!metric_col, na.rm = TRUE),
      descriptions = first(descriptions), 
      .groups = "drop"
    )
  wide_summary <- gene_summary %>%
    pivot_wider(
      id_cols = c(genes, descriptions), # These will be the unique row identifiers
      names_from = population,
      values_from = max_metric
    ) %>%
    rowwise() %>%
    mutate(
      num_populations = sum(!is.na(c_across(c(Central, Eastern, Western)))),
      is_shared = ifelse(num_populations > 1, "Yes", "No")
    ) %>%
    ungroup() %>%
    arrange(desc(is_shared), desc(num_populations)) %>%
    select(genes, descriptions, num_populations, is_shared, everything())
  return(wide_summary)
}

wide_genes_Twithin_percentile <- create_wide_gene_summary(outliers_Twithin_percentile, T_within_myr)
wide_genes_Tratio_percentile <- create_wide_gene_summary(outliers_Tratio_percentile, Tpooled_Twithin_ratio)
wide_genes_Twithin_fixed <- create_wide_gene_summary(outliers_Twithin_fixed, T_within_myr)
wide_genes_Tratio_fixed <- create_wide_gene_summary(outliers_Tratio_fixed, Tpooled_Twithin_ratio)

write_csv(wide_genes_Twithin_percentile, "outliers_by_population_wide/wide_genes_Twithin_percentile.csv")
write_csv(wide_genes_Tratio_percentile, "outliers_by_population_wide/wide_genes_Tratio_percentile.csv")
write_csv(wide_genes_Twithin_fixed, "outliers_by_population_wide/wide_genes_Twithin_fixed.csv")
write_csv(wide_genes_Tratio_fixed, "outliers_by_population_wide/wide_genes_Tratio_fixed.csv")

write_csv(wide_genes_Twithin_fixed, "chimp_wide_genes_Twithin_6MYA.csv")
write_csv(outliers_Twithin_fixed, "chimp_outliers_Twithin_fixed.csv" )
########### PLOTTTING MANHATTANS #############################
# Process remaining data
coalescence_df_with_features2  <- coalescence_df_with_features2 %>%
  filter(!population %in% c("hybrid", "Western x Central hybrid")) %>%
  arrange(chr_clean, start) %>%
  mutate(row_index = row_number())

# Recalculate chromosome label positions
axisdf <- coalescence_df_with_features2 %>%
  group_by(chr_clean) %>%
  summarize(center = mean(row_index))

# Custom population colors
pop_colors <- c(
  "Western" = "#9dced9",
  #"hybrid" = "#bbc671",
  "Eastern" = "#ffb35a",
  "Central" = "#4c5d4c"
)


# 1. Chromosome-faceted plots with FIXED Y-AXIS
chromosome_plot_pooled <- function(data, y_var, threshold_quantile, title, y_label, gene_map_df = NULL) {
  threshold <- quantile(data[[y_var]], probs = threshold_quantile, na.rm = TRUE)
  plot_data <- data %>%
    mutate(
      position_mb = start / 1e6,
      color_condition = ifelse(!!sym(y_var) > threshold, "Above", "Below")
    )
  gene_labels <- plot_data %>%
    filter(color_condition == "Above", !is.na(genes)) %>%
    mutate(gene_original = strsplit(as.character(genes), ",")) %>%
    unnest(gene_original) %>%
    mutate(gene_original = trimws(gene_original)) %>%
    filter(gene_original != "")
  
  if (!is.null(gene_map_df)) {
    clean_gene_map_df <- gene_map_df %>% distinct(genes, .keep_all = TRUE)
    gene_labels <- gene_labels %>%
      left_join(clean_gene_map_df, by = c("gene_original" = "genes")) %>%
      mutate(gene = coalesce(genes_renamed, gene_original))
  } else {
    gene_labels <- gene_labels %>% mutate(gene = gene_original)
  }
  gene_labels <- gene_labels %>%
    filter(!is.na(gene) & gene != "" & gene != "NA" & !grepl("^LOC", gene)) %>%
    group_by(gene, chr_clean) %>%
    slice_max(order_by = !!sym(y_var), n = 1, with_ties = FALSE) %>%
    ungroup()
  
  ggplot(plot_data, aes(x = position_mb, y = !!sym(y_var))) +
    
    ggrastr::geom_point_rast(aes(color = color_condition), size = 0.8, alpha = 0.7, dpi = 300) +
    
    scale_color_manual(values = c("Above" = "#8A865D", "Below" = "grey70"), guide = "none") +
    

    geom_hline(yintercept = threshold, linetype = "dashed", color = "red", linewidth = 0.5) +
    
    geom_text_repel(data = gene_labels, aes(label = gene),
                    size = 2, color = "black", min.segment.length = 0.1,
                    segment.size = 0.2, box.padding = 0.3, max.overlaps = Inf) +
    
    facet_wrap(~chr_clean, scales = "free_x", ncol = 3) +
    theme_classic() +
    theme(
      legend.position = "none",
      panel.grid.major.y = element_line(linewidth = 0.2, color = "grey90"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
      axis.text.y = element_text(size = 8),
      strip.text = element_text(size = 8, face = "bold"),
      strip.background = element_rect(fill = "grey85", color = NA),
      panel.spacing = unit(0.5, "lines")
    ) +
    labs(x = "Position (Mb)", y = y_label, title = title)
}

# Call the new population-agnostic function
pooled_plot <- chromosome_plot_pooled(
  data = coalescence_df_with_features2,
  y_var = "T_pooled_myr",
  threshold_quantile = 0.9990,
  title = "Pooled TMRCA 99.90th",
  y_label = "Pairwise TMRCA (Mya)",
  gene_map_df = gene_map
)
ggsave("chimp_manhattan_pooled_99.90th_filtered.pdf", plot = pooled_plot, width = 10, height = 12)


chromosome_plot_per_pop_thresholds <- function(data, y_var, threshold_quantile, title, y_label, 
                                               gene_map_df = NULL, 
                                               pop_colors = c("Western" = "#9dced9", "Eastern" = "#ffb35a", "Central" = "#4c5d4c")) {
  
  # --- 1. Filter data ---
  data_filtered <- data %>%
    filter(!population %in% c("hybrid", "Western x Central hybrid"))
  #  Calculate a threshold for EACH population
  pop_thresholds <- data_filtered %>%
    group_by(population) %>%
    summarise(threshold = quantile(!!sym(y_var), probs = threshold_quantile, na.rm = TRUE), .groups = "drop")
  
  # --- 2. Prepare base data for plotting ---
  plot_data <- data_filtered %>%
    mutate(position_mb = start / 1e6) %>%
    left_join(pop_thresholds, by = "population")
  
  # --- 3. Prepare gene labels ---
  # The filter now compares each point to its own population's threshold
  gene_labels <- plot_data %>%
    filter(!!sym(y_var) > threshold, !is.na(genes)) %>%
    # ... (rest of gene label logic is the same) ...
    mutate(gene_original = strsplit(as.character(genes), ",")) %>%
    unnest(gene_original) %>%
    mutate(gene_original = trimws(gene_original)) %>%
    filter(gene_original != "")
  
  if (!is.null(gene_map_df)) {
    clean_gene_map_df <- gene_map_df %>% distinct(genes, .keep_all = TRUE)
    gene_labels <- gene_labels %>%
      left_join(clean_gene_map_df, by = c("gene_original" = "genes")) %>%
      mutate(gene = coalesce(genes_renamed, gene_original))
  } else {
    gene_labels <- gene_labels %>% mutate(gene = gene_original)
  }
  
  gene_labels <- gene_labels %>%
    filter(!is.na(gene) & gene != "" & gene != "NA" & !grepl("^LOC", gene)) %>%
    group_by(gene, chr_clean) %>%
    slice_max(order_by = !!sym(y_var), n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # --- 4. Generate the plot ---
  ggplot(plot_data, aes(x = position_mb, y = !!sym(y_var))) +
    
    # Points below their pop-specific threshold
    ggrastr::geom_point_rast(data = . %>% filter(!!sym(y_var) <= threshold),
                             color = "grey70", size = 0.5, alpha = 0.7, dpi = 300) +
    
    # Points above their pop-specific threshold
    ggrastr::geom_point_rast(data = . %>% filter(!!sym(y_var) > threshold),
                             aes(color = population), size = 1.5, alpha = 0.9, dpi = 300) +
    
    scale_color_manual(name = "Population", values = pop_colors) +
    
    # **KEY CHANGE**: Draw a separate dashed line for each population's threshold
    geom_hline(data = pop_thresholds, 
               aes(yintercept = threshold, color = population), 
               linetype = "dashed", linewidth = 0.6) +
    
    geom_text_repel(data = gene_labels, aes(label = gene),
                    size = 2, color = "black", min.segment.length = 0.1,
                    segment.size = 0.2, box.padding = 0.3, max.overlaps = Inf) +
    
    facet_wrap(~chr_clean, scales = "free_x", ncol = 3) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_line(linewidth = 0.2, color = "grey90"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
      axis.text.y = element_text(size = 8),
      strip.text = element_text(size = 8, face = "bold"),
      strip.background = element_rect(fill = "grey85", color = NA),
      panel.spacing = unit(0.5, "lines")
    ) +
    labs(x = "Position (Mb)", y = y_label, title = title)
}


chromosome_plot_vary_thresh<-chromosome_plot_per_pop_thresholds(
  data = coalescence_df_with_features2,
  y_var = "T_within_myr",
  threshold_quantile = 0.9990,
  title = "Within TMRCA 99.90th", 
  y_label = "Pairwise TMRCA (Mya)", 
  gene_map_df = gene_map, pop_colors = pop_colors  
)
#ggsave("chimp_manhattan_by_chromosome_plot_vary_thresh_99.99th_filtered.png", plot = chromosome_plot_vary_thresh, width = 10, height = 12, dpi = 300)
ggsave("chimp_manhattan_by_chromosome_plot_vary_thresh_99.90th_filtered.pdf", plot = chromosome_plot_vary_thresh, width = 10, height = 12)

#2.Population-faceted plots with free y axis

create_advanced_manhattan_plot <- function(data, gene_map_df, y_var, threshold_quantile, title, y_label,
                                           label_colors = c("TRUE" = "red", "FALSE" = "black")) {
  
  pop_thresholds <- data %>%
    group_by(population) %>%
    summarise(
      threshold = quantile(!!sym(y_var), probs = threshold_quantile, na.rm = TRUE),
      .groups = "drop"
    )
  
  significant_genes_unnested <- data %>%
    left_join(pop_thresholds, by = "population") %>%
    filter(!!sym(y_var) > threshold, !is.na(genes)) %>%
    mutate(gene_original = strsplit(as.character(genes), ",")) %>%
    unnest(gene_original) %>%
    mutate(gene_original = trimws(gene_original)) %>%
    filter(gene_original != "")
  all_significant_genes <- significant_genes_unnested %>%
    # Join your mapping table. The key from the main data is 'gene_original',
    # and the key from your mapping table is assumed to be 'genes'.
    left_join(gene_map_df, by = c("gene_original" = "genes")) %>%
    # Use coalesce to create the definitive 'gene' label column.
    # It takes 'genes_renamed' if available, otherwise it falls back to 'gene_original'.
    mutate(gene = coalesce(genes_renamed, gene_original)) %>%
    # Final cleanup of the chosen gene names
    filter(gene != "", !str_detect(gene, "^LOC|^LINC|NA"))
  
  # c. Find which genes are unique (appear in only one population)
  # This section and the next remain the same, as they operate on the new 'gene' column
  unique_genes <- all_significant_genes %>%
    group_by(gene) %>%
    summarise(pop_count = n_distinct(population), .groups = "drop") %>%
    filter(pop_count == 1)
  
  # d. Prepare the final gene labels data, flagging unique genes
  gene_labels <- all_significant_genes %>%
    mutate(is_unique = gene %in% unique_genes$gene) %>%
    group_by(gene, population) %>%
    arrange(desc(!!sym(y_var))) %>%
    slice(1) %>%
    ungroup()
  
  # Calculate chromosome center points for x-axis labels
  chrom_centers <- chrom_offsets %>%
    mutate(center = offset + chr_len / 2,
           center_mb = center / 1e6)
  
  # --- 3. Point Coloring Logic ---
  data <- data %>%
    filter(!population %in% c("hybrid", "Western x Central hybrid")) %>%
    mutate(
      chr_parity = ifelse(chr_num %% 2 == 0, "even", "odd"),
      point_color = case_when(
        chr_parity == "even"  ~ "grey70",
        population == "Western" ~ "#9dced9",
        population == "Eastern" ~ "#ffb35a",
        population == "Central" ~ "#4c5d4c",
        TRUE                  ~ "black"
      )
    )
  
  # --- 4. Generate the Plot ---
  # The plotting code remains identical
  ggplot(data, aes(x = genome_pos_mb, y = !!sym(y_var))) +
    geom_vline(
      xintercept = c(0, chrom_offsets$offset + chrom_offsets$chr_len) / 1e6,
      color = "grey80",
      linewidth = 0.3
    ) +
    ggrastr::geom_point_rast(aes(color = point_color), size = 0.8, alpha = 0.7, show.legend = FALSE) +
    geom_hline(
      data = pop_thresholds,
      aes(yintercept = threshold),
      linetype = "dashed",
      color = "red",
      linewidth = 0.5
    ) +
    scale_x_continuous(
      breaks = chrom_centers$center_mb,
      labels = chrom_centers$chr_clean
    ) +
    scale_color_identity() +
    new_scale_color() +
    geom_text_repel(
      data = gene_labels,
      aes(label = gene, color = is_unique),
      size = 2, min.segment.length = 0.1,
      segment.size = 0.2, box.padding = 0.3, max.overlaps = Inf
    ) +
    scale_color_manual(
      name = "",
      values = label_colors, 
      labels = c("TRUE" = "Unique to Population", "FALSE" = "Shared")
    ) +
    facet_wrap(~population, ncol = 1, scales = "free_y") +
    theme_classic() +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_line(linewidth = 0.2, color = "grey90"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
      axis.text.y = element_text(size = 8),
      strip.text = element_text(size = 8, face = "bold"),
      strip.background = element_rect(fill = "grey85", color = NA),
      panel.spacing = unit(1, "lines")
    ) +
    labs(x = "Genomic Position (Mb)", y = y_label, title = title)
}

# A custom function to format log-axis labels with rounded exponents
label_log_rounded <- function(x) {
  log_x <- log10(x)
  rounded_log_x <- round(log_x, 1)
  parse(text = paste0("10^", rounded_log_x))
}





quantile_to_use <- 0.9990
plot_avg_pop <- create_advanced_manhattan_plot(
  coalescence_df_with_features2, 
  gene_map,
  "T_within_myr", 
  quantile_to_use,
  "Genome-wide Pairwise Coalescent Time, TMRCA (99.90th Percentile)",
  "Pairwise TMRCA (Mya)"
 )
ggsave("chimps_genome_plot_avg_pairwise_by_pop_99.90th.pdf", plot_avg_pop, width = 10, height=6)


quantile_to_use <- 0.9999
plot_avg_pop <- create_advanced_manhattan_plot(
  coalescence_df_with_features2, 
  gene_map,
  "T_within_myr", 
  quantile_to_use,
  "Genome-wide Pairwise Coalescent Time, TMRCA (99.99th Percentile)",
  "Pairwise TMRCA (Mya)"
)
plot_avg_pop



single_pop_data <- coalescence_df_with_features2 %>%
  filter(population == "Central") # - it doesn't matter
quantile_to_use <- 0.9990
plot_avg_pop2 <- create_advanced_manhattan_plot(
  single_pop_data, 
  gene_map,
  "T_pooled_myr", 
  quantile_to_use,
  "Genome-wide Pairwise Coalescent Time, TMRCA (99.90th Percentile)",
  "Pairwise TMRCA (Mya)"
)
ggsave("all_chimps_pairwise_tmrca_gwide_99.90th.pdf", plot_avg_pop2, width = 10, height=3)



plot_ratio_pop <- create_advanced_manhattan_plot(
  data = coalescence_df_with_features2, 
  gene_map_df = gene_map, # <-- Add the gene map here
  y_var = "Tpooled_Twithin_ratio", 
  threshold_quantile = quantile_to_use,
  title = "Genome-wide Tpooled/Twithin (99.99th Percentile)",
  y_label = "log10(Tpooled/Twithin)",
  # Pass the inverted color scheme: Shared = red, Unique = black
  label_colors = c("TRUE" = "black", "FALSE" = "red") 
)

# Add the log scale with the fixed upper limit and rounded labels
plot_ratio_pop_log_formatted <- plot_ratio_pop + 
  scale_y_log10(
    limits = c(10^0, NA), # Set top limit only
    labels = label_log_rounded # Use the custom rounding function
  )


#ggsave("chimps_genome_plot_avg_pairwise_by_pop.png", plot_avg_pop, width = 10, height = 6, dpi = 300)
#ggsave("chimps_genome_plot_avg_tmrca_by_pop.png", plot_avg_pop2, width = 10, height = 6, dpi = 300)
#ggsave("chimps_genome_plot_tpooled_ratio_by_pop.png", plot_ratio_pop_log_formatted, width = 10, height = 6, dpi = 300)
ggsave("chimps_genome_plot_avg_pairwise_by_pop.pdf", plot_avg_pop, width = 10, height=6)
ggsave("chimps_genome_plot_tpooled_ratio_by_pop.pdf", plot_ratio_pop_log_formatted, width = 10, height = 6)


# Create and save chromosome-faceted plots
plot_avg_points <- create_faceted_point_plot(
  coalescence_df_with_features2, "T_within_myr", 10,
  "Average Pairwise Coalescence Time by Chromosome",
  "Coalescence Time (myr)"
)
plot_ratio_points <- create_faceted_point_plot(
  coalescence_df_with_features2, "Tpooled_Twithin_ratio", 10,
  "Tpooled/Twithin Ratio by Chromosome",
  "Tpooled/Twithin Ratio (myr)"
)

ggsave("pointplot_avg_pairwise_fixedY.png", plot_avg_points, width = 9, height = 9, dpi = 300)
ggsave("pointplot_avg_pairwise_fixedY.pdf", plot_avg_points, width = 9, height = 9)
ggsave("pointplot_tpooled_ratio_freeY.png", plot_ratio_points, width = 9, height = 9, dpi = 300)
ggsave("pointplot_tpooled_ratio_freeY.pdf", plot_ratio_points, width = 9, height = 9)


tidy_mismatched_genes <- function(df) {
  results_list <- vector("list", nrow(df))
  for (i in 1:nrow(df)) {
    row_data <- df[i, ]
    genes_vec <- str_split(row_data$genes, "[,;]")[[1]]
    desc_vec <- str_split(row_data$descriptions, "[,;]")[[1]]
    master_len <- length(genes_vec)
    length(desc_vec) <- master_len
    results_list[[i]] <- tibble(
      row_data %>% select(-genes, -descriptions),
      gene = trimws(genes_vec),
      description = trimws(desc_vec)
    )
  }
  bind_rows(results_list)
}

plot_zoomed_tmrca <- function(data,
                              gene_map_df,
                              target_chr,
                              target_start,
                              target_end,
                              y_var,
                              y_label,
                              plot_title,
                              label_threshold = 0,
                              label_all_hits = FALSE) { # <-- NEW ARGUMENT
  
  region_df <- data %>%
    filter(
      chromosome == target_chr,
      start <= target_end,
      end >= target_start
    ) %>%
    filter(!population %in% c("hybrid", "Western x Central hybrid"))
  
  if (nrow(region_df) == 0) {
    # ... (error handling is the same) ...
  }
  
  gene_labels <- region_df %>%
    filter(!!sym(y_var) >= label_threshold) %>%
    filter(genes != "." & !is.na(genes)) %>%
    tidy_mismatched_genes() %>%
    left_join(gene_map_df, by = c("gene" = "genes")) %>%
    mutate(display_label = coalesce(genes_renamed, gene)) %>%
    filter(
      !is.na(display_label),
      display_label != "",
      display_label != "NA",
      !str_detect(display_label, "^LOC")
    )
  
  # **KEY CHANGE**: Conditionally apply the filtering for top hits
  if (!label_all_hits) {
    gene_labels <- gene_labels %>%
      group_by(display_label) %>%
      slice_max(order_by = !!sym(y_var), n = 1, with_ties = FALSE) %>%
      ungroup()
  }
  p <- ggplot(region_df, aes(x = start, y = !!sym(y_var), color = population)) +
    geom_hline(yintercept = label_threshold, linetype = "dashed", color = "red", linewidth = 0.4) +
    geom_line(linewidth = 0.8) +
    geom_text_repel(
      data = gene_labels,
      aes(label = display_label),
      nudge_y = 1.5,
      size = 3.5,
      segment.color = "grey50",
      color = "black",
      box.padding = 0.5,
      max.overlaps = Inf
    ) +
    scale_color_manual(
      name = "Population",
      values = c("Western" = "#9dced9", "Eastern" = "#ffb35a", "Central" = "#4c5d4c")
    ) +
    scale_x_continuous(
      labels = scales::label_number(scale = 1e-6, suffix = " Mb"),
      limits = c(target_start, target_end)
    ) +
    labs(
      title = plot_title,
      x = paste("Genomic Position (", target_chr, ")"),
      y = y_label
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    )
  
  return(p)
}

create_overlay_plot <- function(data,
                                gene_map_df, 
                                target_chr,
                                target_start,
                                target_end,
                                plot_title,
                                label_threshold,
                                label_all_hits = FALSE) { # <-- NEW ARGUMENT
  
  # Base plot now receives the new argument
  base_plot <- plot_zoomed_tmrca(
    data = data,
    gene_map_df = gene_map_df,
    target_chr = target_chr,
    target_start = target_start,
    target_end = target_end,
    y_var = "T_within_myr",
    y_label = "Pairwise TMRCA (Mya)",
    plot_title = plot_title,
    label_threshold = label_threshold,
    label_all_hits = label_all_hits # <-- Pass the argument down
  )
  
  # Prepare data for the overlay line (unchanged)
  avg_line_data <- data %>%
    filter(
      chromosome == target_chr,
      start <= target_end,
      end >= target_start
    )
  
  # Prepare labels for the overlay line (Average Pairwise TMRCA)
  avg_tmrca_labels <- data %>%
    filter(
      chromosome == target_chr,
      start <= target_end,
      end >= target_start
    ) %>%
    filter(T_pooled_myr >= label_threshold, !is.na(genes), genes != ".") %>%
    tidy_mismatched_genes() %>%
    left_join(gene_map_df, by = c("gene" = "genes")) %>%
    mutate(display_label = coalesce(genes_renamed, gene)) %>%
    filter(
      !is.na(display_label),
      display_label != "",
      display_label != "NA",
      !str_detect(display_label, "^LOC")
    )
  if (!label_all_hits) {
    avg_tmrca_labels <- avg_tmrca_labels %>%
      group_by(display_label) %>%
      slice_max(order_by = T_pooled_myr, n = 1, with_ties = FALSE) %>% 
      ungroup()
  }
  
  # Add the overlay layers to the base plot (unchanged)
  final_plot <- base_plot +
    geom_line(
      data = avg_line_data,
      aes(y = T_pooled_myr, linetype = "avg_tmrca_line"), # Note: Corrected to use avg_tmrca_myr
      color = "darkgrey",
      linewidth = 0.2
    ) +
    scale_linetype_manual(
      name = NULL,
      values = c("avg_tmrca_line" = "solid"),
      labels = c("avg_tmrca_line" = "Average Pairwise TMRCA")
    ) +
    geom_text_repel(
      data = avg_tmrca_labels,
      aes(y = T_pooled_myr, label = display_label), # Note: Corrected to use Tpooled
      color = "darkgrey",
      size = 2.5,
      nudge_y = 0.5,
      segment.color = "grey50",
      box.padding = 0.5,
      max.overlaps = Inf
    )
  
  return(final_plot)
}


create_overlay_plot2 <- function(data,
                                 gene_map_df, 
                                 target_chr,
                                 target_start,
                                 target_end,
                                 plot_title,
                                 label_threshold,
                                 label_all_hits = FALSE) { # <-- NEW ARGUMENT
  
  # Base plot uses the correct plot_zoomed_tmrca function
  base_plot <- plot_zoomed_tmrca(
    data = data,
    gene_map_df = gene_map_df,
    target_chr = target_chr,
    target_start = target_start,
    target_end = target_end,
    y_var = "T_within_myr",
    y_label = "Pairwise TMRCA (Mya)",
    plot_title = plot_title,
    label_threshold = label_threshold,
    label_all_hits = label_all_hits # <-- Pass the argument down
  )
  return(base_plot)
}

GYP_test <-coalescence_df_with_features  %>%
  filter(chromosome == "chr3_hap1_hsa4") %>%
  #filter(str_detect(descriptions, "MNS blood group")) %>%
  filter(start > 145031136) %>%
  filter(start < 146530513) %>%
  arrange(desc(T_within_myr))
View(GYP_test)
write_csv(GYP_test, "GYPtest.csv")


my_chr <- "chr3_hap1_hsa4"
my_start <- 145100000
my_end <- 146200000
mytitle<-"Pairwise & Average TMRCA GYP region"
my_label_threshold <-0

GYP <- create_overlay_plot2(
  data = coalescence_df_with_features,
  gene_map_df = gene_map,
  target_chr = my_chr,
  target_start = my_start,
  target_end = my_end,
  plot_title = mytitle,
  label_threshold = my_label_threshold,
  label_all_hits = FALSE# <-- Set the new option to TRUE for this plot only
)

print(GYP)

GYP_annotated <- GYP  +
  annotate(
    "rect",
    xmin = 145422074,
    xmax = 145932842,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.2,
    fill = "lightgrey"
  )
print(GYP_annotated)
ggsave("GYP_annotated.pdf", GYP_annotated, width=5, height=4)



HLA_test <-coalescence_df_with_features  %>%
  filter(chromosome == "chr5_hap1_hsa6") %>%
  #filter(str_detect(descriptions, "MNS blood group")) %>%
  filter(start > 35000000) %>%
  filter(start < 40000000) %>%
  arrange(desc(T_within_myr))
View(HLA_test)
write_csv(HLA_test, "HLAtest.csv")

my_chr <- "chr5_hap1_hsa6"
my_start <- 35000000
my_end <- 40000000
mytitle<-"Pairwise & Average TMRCA HLA region"
my_label_threshold <- 10  # Only plot TMRCA values >= 6 Myr

HLA <- create_overlay_plot2(
  data = coalescence_df_with_features,
  gene_map_df = gene_map,
  target_chr = my_chr,
  target_start = my_start,
  target_end = my_end,
  plot_title = mytitle,
  label_threshold = my_label_threshold
)

print(HLA)
ggsave("HLA_justpairwise.pdf", HLA, width=5, height=4)




LILS_test<-coalescence_df_with_features  %>%
  filter(chromosome == "chr20_hap1_hsa19") %>%
  filter(start > 57019858) %>%
  filter(start < 59005408) %>%
  arrange(desc(T_within_myr)) %>%
  View()
write_csv(LILS_test, "LILtest.csv")

my_chr <- "chr20_hap1_hsa19"
my_start <- 57019858
my_end <- 59005408
mytitle<-"TMRCA LIL,KIR,NLRP region"
my_label_threshold <- 10  # Only plot TMRCA values >= 6 Myr

LIL <- create_overlay_plot2(
  data = coalescence_df_with_features,
  gene_map_df = gene_map,
  target_chr = my_chr,
  target_start = my_start,
  target_end = my_end,
  plot_title = mytitle,
  label_threshold = my_label_threshold
)

print(LIL)
ggsave("LIL_KIR_NLRP_justpairwise.pdf", LIL, width=5, height=4)



my_chr <- "chr9_hap1_hsa11"
my_start <- 9300000
my_end <- 9900000
mytitle <- "HBE1 region"
my_label_threshold <- 8

HBE <- create_overlay_plot2(
  data = coalescence_df_with_features,
  gene_map_df = gene_map,
  target_chr = my_chr,
  target_start = my_start,
  target_end = my_end,
  plot_title = mytitle,
  label_threshold = my_label_threshold,
  label_all_hits = TRUE# <-- Set the new option to TRUE for this plot only
)
print(HBE)
HBE_annotated <- HBE +
  annotate(
    "rect",
    xmin = 9399914,
    xmax = 9654572,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.2,
    fill = "lightgrey"
  )
print(HBE_annotated)

ggsave("HBE.pdf", HBE_annotated, width=5, height=4)


my_chr <- "chr18_hap1_hsa16"
my_start <- 2420000
my_end <- 2810000
mytitle <- "Hemoglobin Subunit Alpha region (HBA1/2/3−like)"
my_label_threshold <- 6

HBA <- create_overlay_plot2(
  data = coalescence_df_with_features,
  gene_map_df = gene_map,
  target_chr = my_chr,
  target_start = my_start,
  target_end = my_end,
  plot_title = mytitle,
  label_threshold = my_label_threshold,
  label_all_hits = FALSE# <-- Set the new option to TRUE for this plot only
)
print(HBA)
HBA_annotated <- HBA +
  annotate(
    "rect",
    xmin = 2503753,
    xmax = 2514204,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.2,
    fill = "lightgrey"
  )
print(HBA_annotated)



ggsave("HBA_annotated.pdf", HBA_annotated, width=5, height=4)



my_chr <- "chr19_hap1_hsa17"
my_start <- 95500000
my_end <- 96500000
mytitle<-"TMRCA FOXK2 and FNK3 region"
my_label_threshold <- 10
FOXK2 <- create_overlay_plot2(
  data = coalescence_df_with_features2,
  gene_map_df = gene_map,
  target_chr = my_chr,
  target_start = my_start,
  target_end = my_end,
  plot_title = mytitle,
  label_threshold = my_label_threshold
)
print(FOXK2)
ggsave("FOXK2.pdf", FOXK2, width=5, height=4)



########### RATIO PLOTS #########

label_log_rounded <- function(x) {
  log_x <- log10(x)
  # Round the exponent for cleaner labels, e.g., to 1 decimal place
  rounded_log_x <- round(log_x, 1)
  parse(text = paste0("10^", rounded_log_x))
}

create_overlay_plot3 <- function(data,
                                 gene_map_df, 
                                 target_chr,
                                 target_start,
                                 target_end,
                                 plot_title,
                                 label_threshold,
                                 label_all_hits = FALSE) { # <-- NEW ARGUMENT
  
  base_plot <- plot_zoomed_tmrca(
    data = data,
    gene_map_df = gene_map_df,
    target_chr = target_chr,
    target_start = target_start,
    target_end = target_end,
    y_var = "Tpooled_Twithin_ratio",
    y_label = "log10(Tpooled/Twithin)",
    plot_title = plot_title,
    label_threshold = label_threshold,
    label_all_hits = label_all_hits # <-- Pass the argument down
  )
  
  # You might want to add the log scale here or outside the function
  # For example:
  # base_plot <- base_plot + scale_y_log10(labels = label_log_rounded)
  
  return(base_plot)
}


#### CENTRAL OUTLIERS #####
my_chr <- "chr6_hap1_hsa7"
my_start <- 109000000
my_end <- 116000000
mytitle<-"PMS2 and PRKRIP1 region"
my_label_threshold <- 10^3.9 # Only plot TMRCA values >= 6 Myr

PRKRIP1  <- create_overlay_plot3(data = coalescence_df_with_features,gene_map_df = gene_map, target_chr = my_chr, target_start = my_start,target_end = my_end, plot_title = mytitle, label_threshold = my_label_threshold) + 
  scale_y_log10(labels = label_log_rounded)

print(PRKRIP1)
ggsave("outliers_zoom/PRKRIP1_v2.pdf", PRKRIP1, width=4, height=4)


my_chr <- "chr20_hap1_hsa19"
my_start <- 2400000
my_end <- 3400000
mytitle<-"TPGS1 region"
my_label_threshold <- 10^3.9 # Only plot TMRCA values >= 6 Myr

TPGS1  <- create_overlay_plot3(data = coalescence_df_with_features,gene_map_df = gene_map, target_chr = my_chr, target_start = my_start,target_end = my_end, plot_title = mytitle, label_threshold = my_label_threshold) + 
  scale_y_log10(labels = label_log_rounded)
print(TPGS1)
ggsave("outliers_zoom/TPGS1.pdf", TPGS1, width=4, height=4)


my_chr <- "chr13_hap1_hsa2b"
my_start <- 14100000
my_end <- 14700000
mytitle<-"SLC35F5 region"
my_label_threshold <- 10^3.7 # Only plot TMRCA values >= 6 Myr

SLC35F5 <- create_overlay_plot3(data = coalescence_df_with_features,gene_map_df = gene_map, target_chr = my_chr, target_start = my_start,target_end = my_end, plot_title = mytitle, label_threshold = my_label_threshold) + 
  scale_y_log10(labels = label_log_rounded)
print(SLC35F5)
ggsave("outliers_zoom/SLC35F5.pdf", SLC35F5, width=4, height=4)


### EASTERN OUTLIERS"

my_chr <- "chr5_hap1_hsa6"
my_start <- 73200000
my_end <- 	74200000
mytitle<-"EYS region"
my_label_threshold <- 10^2.4 # Only plot TMRCA values >= 6 Myr
EYS<- create_overlay_plot3(data = coalescence_df_with_features,gene_map_df = gene_map, target_chr = my_chr, target_start = my_start,target_end = my_end, plot_title = mytitle, label_threshold = my_label_threshold) + 
  scale_y_log10(labels = label_log_rounded)
print(EYS)
ggsave("outliers_zoom/EYS.pdf", EYS, width=4, height=4)



my_chr <- "chr7_hap1_hsa8"
my_start <- 58500000
my_end <- 	60200000
mytitle<-"SNTG1 region"
my_label_threshold <- 10^2.4 # Only plot TMRCA values >= 6 Myr
SNTG1<- create_overlay_plot3(data = coalescence_df_with_features,gene_map_df = gene_map, target_chr = my_chr, target_start = my_start,target_end = my_end, plot_title = mytitle, label_threshold = my_label_threshold) + 
  scale_y_log10(labels = label_log_rounded)
print(SNTG1)
ggsave("outliers_zoom/SNTG1.pdf", SNTG1, width=4, height=4)




my_chr <- "chr15_hap1_hsa14"
my_start <- 103500000
my_end <- 104500000
mytitle<-"IGHV4 region"
my_label_threshold <- 10^2.4 # Only plot TMRCA values >= 6 Myr

IGHV4  <- create_overlay_plot3(data = coalescence_df_with_features,gene_map_df = gene_map, target_chr = my_chr, target_start = my_start,target_end = my_end, plot_title = mytitle, label_threshold = my_label_threshold) + 
  scale_y_log10(labels = label_log_rounded)
print(IGHV4)
ggsave("outliers_zoom/IGHV4.pdf", IGHV4, width=4, height=4)


my_chr <- "chr18_hap1_hsa16"
my_start <-87500000
my_end <-	88500000
mytitle<-"SPATAL region"
my_label_threshold <- 10^2.3 # Only plot TMRCA values >= 6 Myr
SPATAL <- create_overlay_plot3(data = coalescence_df_with_features,gene_map_df = gene_map, target_chr = my_chr, target_start = my_start,target_end = my_end, plot_title = mytitle, label_threshold = my_label_threshold) + 
  scale_y_log10(labels = label_log_rounded)
print(SPATAL)
ggsave("outliers_zoom/SPATAL.pdf", SPATAL, width=4, height=4)


### "VERUS OUTLIERS"
my_chr <-"chr19_hap1_hsa17"
my_start <- 96300000
my_end <- 96600000
mytitle<-"B3GNTL1 region"
my_label_threshold <- 10^1 # Only plot TMRCA values >= 6 Myr
B3GNTL1_verus  <- create_overlay_plot3(data = coalescence_df_with_features,gene_map_df = gene_map, target_chr = my_chr, target_start = my_start,target_end = my_end, plot_title = mytitle, label_threshold = my_label_threshold) + 
  scale_y_log10(labels = label_log_rounded)
print(B3GNTL1_verus)
ggsave("outliers_zoom/B3GNTL1_verus.pdf", B3GNTL1_verus, width=4, height=4)


my_chr <-"chr14_hap1_hsa13"
my_start <-119400000
my_end <- 119900000
mytitle<-"RASA3 region"
my_label_threshold <- 10^1  # Only plot TMRCA values >= 6 Myr
RASA3_verus  <- create_overlay_plot3(data = coalescence_df_with_features2,gene_map_df = gene_map, target_chr = my_chr, target_start = my_start,target_end = my_end, plot_title = mytitle, label_threshold = my_label_threshold) + 
  scale_y_log10(labels = label_log_rounded)
print(RASA3_verus)
ggsave("outliers_zoom/RASA3_verus.pdf", RASA3_verus, width=4, height=4)

my_chr <- "chr5_hap1_hsa6"
my_start <- 38000000
my_end <- 39000000
mytitle<-"HLA region"
my_label_threshold <- 10^1  # Only plot TMRCA values >= 6 Myr
HLA_verus  <- create_overlay_plot3(data = coalescence_df_with_features,gene_map_df = gene_map, target_chr = my_chr, target_start = my_start,target_end = my_end, plot_title = mytitle, label_threshold = my_label_threshold) + 
  scale_y_log10(labels = label_log_rounded)
  
print(HLA_verus)
ggsave("outliers_zoom/HLA_verus.pdf", HLA_verus, width=4, height=4)