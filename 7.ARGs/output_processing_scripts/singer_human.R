setwd("/Users/joanocha/Google Drive/My Drive/POSTDOC/PANPAN/analysis/Figure4_SINGER/")
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
conflict_prefer("mutate", "dplyr")




coalescence_df <- fread("humans/hprc_hgsvc3_hgdp1k_results/genome-wide_metrics_annotatedbyPop.csv") %>%
  rename(chrom = chromosome)

annotation_df <- fread("humans/ht2t_gene.bed.gz") %>%
  rename(genes = name) %>%
  select(genes,gene_biotype, description) %>%
  distinct()

full_annotation_df <- fread("humans/ht2t_gene.bed.gz") 

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

coalescence_df$avg_tmrca_myr <- coalescence_df$avg_tmrca * 28 / 1e6
coalescence_df$T_pooled_myr <- coalescence_df$T_pooled * 28 / 1e6
coalescence_df$T_within_myr <- coalescence_df$T_within * 28 / 1e6


features_df <- read_tsv(
  "humans/ht2t_genomefeatures.bed",
  col_names = c("chrom", "start", "end", "feature"),
  col_types = "cddc"
  )
  
tr_df <- read_tsv(
  "humans/homo_catalog.no_overlaps_simp.bed",
  col_names = c("chrom", "start", "end", "motif"),
  col_types = "cddc"
  )
  
sedef_df <- read_tsv(
  "humans/ht2t_sedefSegDups.bed",
  col_names = c("chrom", "start", "end", "coordinates"),
  col_types = "cddc",
  col_select = 1:4
  )

rmsk_df <- read_tsv(
  "humans/ht2t_RepeatMasker.bed",
  col_names = c("chrom", "start", "end", "feature"),
  col_types = "cddc",
  col_select = 1:4
)


SRA_df <- read_tsv(
  "humans/ht2t_SR_mask.bed",
  col_names = c("chrom", "start", "end"),
  col_types = "cddc",
  col_select = 1:3
)

unique_windows_df <- coalescence_df %>%
  select(chrom, start, end) %>%
  distinct() # Create a dataframe with only unique windows to avoid multiplying overlaps.

all_features_df <- bind_rows(
  features_df,
  mutate(tr_df, feature = "TR"),
  mutate(sedef_df, feature = "SEDEF_SD"),
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
  left_join(rmsk_overlap_summary, by = c("chrom", "start", "end")) %>%
  left_join(features_overlap_summary, by = c("chrom", "start", "end")) %>%
  left_join(SRA_overlap_summary, by = c("chrom", "start", "end")) %>%
  mutate(
    feature_annotation = if_else(is.na(feature_annotation), ".", feature_annotation),
    tr_overlap_bp = coalesce(tr_overlap_bp, 0L),
    sedef_overlap_bp = coalesce(sedef_overlap_bp, 0L),
    rmsk_overlap_bp = coalesce(rmsk_overlap_bp, 0L),
    SRA_overlap_bp = coalesce(SRA_overlap_bp, 0L),
    features_overlap_bp = coalesce(features_overlap_bp, 0L),
    window_size = end - start,
    tr_overlap_percent = (tr_overlap_bp / window_size) * 100,
    sedef_overlap_percent = (sedef_overlap_bp / window_size) * 100,
    rmsk_overlap_percent = (rmsk_overlap_bp / window_size) * 100,
    features_overlap_percent = (features_overlap_bp / window_size) * 100,
    SRA_overlap_percent = (SRA_overlap_bp / window_size) * 100
  ) %>%
  select(
    chrom, start, end, population, avg_tmrca_myr, T_pooled_myr, T_within_myr, Tpooled_Twithin_ratio, genes, feature_annotation, gene_biotypes, descriptions,
    tr_overlap_percent, 
    sedef_overlap_percent,
    rmsk_overlap_percent,
    features_overlap_percent,
    SRA_overlap_percent
  ) %>%
  rename(chromosome = chrom)  %>%
  mutate(
    chr_clean = sub("_.*", "", chromosome)
  ) %>%
    filter(!chr_clean %in% c("chrX", "chrY")) %>% 
    mutate(
      chr_clean = factor(chr_clean, levels = paste0("chr", 1:22)),
      midpoint = (start + end) / 2,
      position_mb = midpoint / 1e6,
      chr_num = as.integer(gsub("chr", "", chr_clean))
    ) %>%
    arrange(chr_clean, start)
  

############ SUMMARY STATS ON TRs and other features #####################
unique_windows_summary <- coalescence_df_with_features %>%
    group_by(chromosome, start, end, genes, feature_annotation) %>%
    summarise(
      tr_overlap_percent = first(tr_overlap_percent),
      sedef_overlap_percent = first(sedef_overlap_percent),
      SRA_overlap_percent = first(SRA_overlap_percent),
      rmsk_overlap_percent = first(rmsk_overlap_percent),
      T_within_myr = mean(T_within_myr, na.rm = TRUE),
      Tpooled_Twithin_ratio = mean(Tpooled_Twithin_ratio, na.rm = TRUE),
    ) %>%
    ungroup() # Don't forget to ungroup

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
    x = "Overlap Percentage with 1kb Window (%)",
    y = "Number of Windows in Humans"
    ) +
    theme_bw()
ggsave("humans/features_summary.pdf", features_summary, width=6, height=3)  

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

  # Custom population colors  
  pop_colors <- c(
    "GWD" = "#CC9900",
    "ESN" = "#FFCC33",
    "MKK" = "#CC9900",
    "PJL" = "#CC0066",
    "MSL" = "#CC9933",
    "ASW" = "#FF6600",
    "PEL" = "#FF3333",
    "CHS" = "#66FF33",
    "CLM" = "#993333",
    "YRI" = "#CC9900",
    "ACB" = "#FF9933",
    "KHV" = "#99CC33",
    "ITU" = "#CC3366",
    "ASK" = "#999999",
    "PUR" = "#993300",
    "LWK" = "#996600",
    "STU" = "#FF33FF",
    "BEB" = "#660066",
    "JPT" = "#006600",
    "FIN" = "#66CCCC"
  )
  

  
################ RANKING ##########################
create_ranked_gene_plot <- function(data, y_variable, y_axis_label) {
  
  # Use the {{ }} operator to pass the column name to dplyr and ggplot
  # This is the modern way to program with these packages.
  
  plot_data_faceted <- data %>%
    # Reorder the gene labels based on the chosen y-variable
    mutate(
      gene_feature = fct_reorder(gene_feature, {{ y_variable }}),
      start_ranked = reorder_within(midpoint, -{{ y_variable }}, population)
    )
  
  # Filter out hybrid cases
  plot_data_filtered <- plot_data_faceted %>%
    filter(population != "hybrid")
  
  # --- Create the plot ---
  ggplot(
    plot_data_filtered,
    aes(
      x = start_ranked,
      y = {{ y_variable }}, # Use the y_variable here
      color = population
    )
  ) +
    geom_point(size = 2, alpha = 0.8, show.legend = FALSE) +
    geom_text_repel(
      aes(label = gene_feature),
      size = 2,
      box.padding = 0.5,
      max.overlaps = Inf
    ) +
    scale_x_reordered() +
    scale_color_manual(name = "Population", values = pop_colors) +
    facet_wrap(~ population, scales = "free_x", ncol = 1) +
    labs(
      title = "Top Genes Ranked by Highest Coalescence Metric",
      x = "Genomic Position (Ordered by Rank within Population)",
      y = y_axis_label # Use the custom y-axis label
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold")
    )
}
#data_for_plot_a <- coalescence_df_with_features %>% filter(T_within_myr > 10) %>%
#  mutate(
#    gene_feature = glue("{genes} ({feature_annotation}) (TR: {round(tr_overlap_percent, 2)}%)")
#  ) 
#data_for_plot_b <- coalescence_df_with_features %>% filter(Tpooled_Twithin_ratio > 10) %>%
#  mutate(
#    gene_feature = glue("{genes} ({feature_annotation}) (TR: {round(tr_overlap_percent, 2)}%)")
#  )
# Call the function
#plot_data_faceted_a <-create_ranked_gene_plot(
#  data = data_for_plot_a, 
#  y_variable = T_within_myr, 
#  y_axis_label = "Coalescence Time (Myr)"
#)
#plot_data_faceted_b <-create_ranked_gene_plot(
#  data = data_for_plot_b, 
#  y_variable = Tpooled_Twithin_ratio, 
#  y_axis_label = "T_pooled / T_within Ratio"
#)


### filtering TRs out

GYP_test <-coalescence_df_with_features  %>%
  filter(chromosome == "chr4") %>%
  filter(str_detect(descriptions, "MNS blood group")) %>%
  arrange(desc(T_within_myr))
View(GYP_test)
write_csv(GYP_test, "GYPtest.csv")


############## FILTER FOR COMPLEX NOISY REGIONS ###########
unique(coalescence_df_with_features$feature_annotation)
coalescence_df_with_features2 <- coalescence_df_with_features %>%
  filter(
    tr_overlap_percent <= 5,
    SRA_overlap_percent > 20,
    #sedef_overlap_percent == 0,
    !str_detect(feature_annotation, "Satellite"),
    !str_detect(feature_annotation, "Low_complexity"),
    !str_detect(feature_annotation, "Cen"),
    !str_detect(feature_annotation, "Gap"),
    !str_detect(feature_annotation, "CEN"),
    !str_detect(feature_annotation, "SAT"),
    !str_detect(feature_annotation, "rDNA"),
  ) 

# Check the unique annotations in your new dataframe
unique(coalescence_df_with_features2$feature_annotation)


# CASES KEPT FOR NOW include combinations of windows overlapped by
#"TEL"                  
#"SD",
# SEDEF_SD
#But not by SAT/SD, SAT or CEN, rDNA or rDNA-MODEL 

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
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  labs(
    x = "Human TMRCA in 1kb windows (Mya)",
    y = "Number of windows",
    fill = "Statistic" # Optional: clean up legend title
  ) +
  theme_cowplot() +
  theme(legend.position = "bottom") +
  coord_cartesian(xlim = c(NA, 14), clip = "off")
tmrca_plot
ggsave("humans/tmrca_summary.pdf", tmrca_plot, width = 5, height = 3)


tmrca_plot2 <-  unique_tmrca_summary %>%
  pivot_longer(
    cols = Tpooled_Twithin_ratio,
    names_to = "statistic",
    values_to = "myr"
  ) %>%
  ggplot(aes(x = myr, fill = statistic)) + 
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity", color = "white", fill="#5C608A") +
  labs(
    x = "Human Tpoooled/Twithin",
    y = "Number of windows",
    fill = "Statistic" # Optional: clean up legend title
  ) +
  theme_cowplot() +
  theme(legend.position = "bottom") #+ 
  coord_cartesian(xlim = c(NA, 14), clip = "off") 
tmrca_plot2
ggsave("humans/tmrca_summary2.pdf", tmrca_plot2, width = 5, height = 3)










#################### GET TABLES OF TOP GENES ##############
#### ALL Humans ##############
ranked_T_pooled_myr <- coalescence_df_with_features %>%
  distinct(chromosome, start, end, .keep_all = TRUE) %>%
  filter(T_pooled_myr > 6) %>%
  arrange(desc(T_pooled_myr)) %>%
  select(chromosome, start, end, T_pooled_myr, genes, descriptions,  gene_biotypes, SRA_overlap_percent, tr_overlap_percent,  rmsk_overlap_percent, feature_annotation)
#write_csv(ranked_T_pooled_myr, "ranked_T_pooled_myr.csv")
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
#write_csv(wide_ranked_genes_T_pooled_myr, "humans/humans_wide_ranked_genes_T_pooled_6myr.csv")



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
      genes,
      gene_biotypes, 
      descriptions,
      T_within_myr,
      Tpooled_Twithin_ratio
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
significant_genes_Twithin <- get_pop_significant_genes(
  data = coalescence_df_with_features2, # Use the corrected data frame
  metric_var = T_within_myr,
  threshold_value = 6
) %>%
  arrange(desc(T_within_myr))

quantile_to_use <- 0.9999
significant_genes_ratio <- get_pop_significant_genes(
  data = coalescence_df_with_features2, # Use the corrected data frame
  metric_var = Tpooled_Twithin_ratio,
  threshold_quantile = quantile_to_use
) %>%
  arrange(desc(Tpooled_Twithin_ratio))



#View(significant_genes_Twithin)
#View(significant_genes_ratio)
write_csv(significant_genes_Twithin, "humans/humans_high_Twithin_significant_genes_6MY.csv")


# Convert the 'significant_genes_Twithin' dataframe into a list format where names are populations and values are vectors of gene names.
upset_list_Twithin <- split(significant_genes_Twithin$genes, significant_genes_Twithin$population)
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
write_csv(wide_significant_genes_Twithin, "humans/humans_wide_6myt_genes_Twithin.csv")



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
write_csv(wide_significant_genes_ratio, "humans/humans_wide_99.00th_genes_ratio.csv")


plot_data_ranked_ratio <- coalescence_df_with_features2 %>%
  arrange(desc(Tpooled_Twithin_ratio)) %>%
  filter(Tpooled_Twithin_ratio > 10)  %>%
  select(chromosome, start, end, population, Tpooled_Twithin_ratio, genes, gene_biotypes, descriptions, feature_annotation, tr_overlap_percent, SRA_overlap_percent)
write_csv(plot_data_ranked_ratio, "humans/ranked_ratio.csv")


plot_data_ranked_tmrca <- coalescence_df_with_features2 %>%
  arrange(desc(T_within_myr)) %>%
  filter(T_within_myr > 10) %>% 
  select(chromosome, start, end, population, T_within_myr, genes, gene_biotypes, descriptions, feature_annotation, tr_overlap_percent, SRA_overlap_percent)
write_csv(plot_data_ranked_tmrca, "humans/ranked_tmrca.csv")

########### PLOTTTING MANHATTANS ALL TOGETHER#############################
# Process remaining data
coalescence_df_with_features2  <- coalescence_df_with_features2 %>%
  arrange(chr_clean, start) %>%
  mutate(row_index = row_number())

# Recalculate chromosome label positions
axisdf <- coalescence_df_with_features2 %>%
  group_by(chr_clean) %>%
  summarize(center = mean(row_index))


# Custom population colors
pop_colors <- c(
  "AFR" = "#FFCC33",
  "AMR" = "#FF3333",
  "EAS" = "#006600",
  "SAS" = "#FF33FF",
  "EUR" = "#66CCCC"
  
)
#threshold<-10

########### PLOTTTING MANHATTANS BY POP#############################
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
    filter(gene_original != "") %>%
    filter(gene_original !=str_detect(genes, "LOC|LINC"))
  
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
    
    scale_color_manual(values = c("Above" = "#5C608A", "Below" = "grey70"), guide = "none") +
    
    geom_hline(yintercept = threshold, linetype = "dashed", color = "red", linewidth = 0.5) +
    
    geom_text_repel(data = gene_labels, aes(label = gene),
                    size = 2, color = "black", min.segment.length = 0.1,
                    segment.size = 0.2, box.padding = 0.3, max.overlaps = Inf) +
    
    facet_wrap(~chr_clean, scales = "free_x", ncol = 4) +
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
  threshold_quantile = 0.9900,
  title = "Pooled TMRCA 99.00th",
  y_label = "Pairwise TMRCA (Mya)",
  gene_map_df = NULL
)
ggsave("humans/humans_manhattan_pooled_99.00th_filtered.pdf", plot = pooled_plot, width = 12, height = 12)

chromosome_plot_per_pop_thresholds <- function(data, y_var, threshold_quantile, title, y_label, 
                                               gene_map_df = NULL, 
                                               pop_colors =  c(
                                                 "AFR" = "#FFCC33",
                                                 "AMR" = "#FF3333",
                                                 "EAS" = "#006600",
                                                 "SAS" = "#FF33FF",
                                                 "EUR" = "#66CCCC"
                                                 
                                               )) {
  
  # --- 1. Filter data ---

  pop_thresholds <-  data %>%
    group_by(population) %>%
    summarise(threshold = quantile(!!sym(y_var), probs = threshold_quantile, na.rm = TRUE), .groups = "drop")
  
  # --- 2. Prepare base data for plotting ---
  plot_data <-  data %>%
    mutate(position_mb = start / 1e6) %>%
    left_join(pop_thresholds, by = "population")
  
  gene_labels <- plot_data %>%
    filter(!!sym(y_var) > threshold, !is.na(genes)) %>%
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
  gene_map_df = NULL, 
  pop_colors = pop_colors  
)
ggsave("humans/human_manhattan_by_chromosome_plot_vary_thresh_99.90th_filtered.pdf", plot = chromosome_plot_vary_thresh, width = 10, height = 12)







create_population_genome_plot <- function(data, y_var, threshold_quantile, title, y_label, 
                                          label_colors = c("TRUE" = "red", "FALSE" = "black")) {
  
  pop_thresholds <- data %>%
    group_by(population) %>%
    summarise(
      threshold = quantile(!!sym(y_var), probs = threshold_quantile, na.rm = TRUE),
      .groups = "drop"
    )
  
  all_significant_genes <- data %>%
    left_join(pop_thresholds, by = "population") %>%
    filter(!!sym(y_var) > threshold, !is.na(genes)) %>%
    mutate(gene = strsplit(as.character(genes), ",")) %>%
    unnest(gene) %>%
    mutate(gene = trimws(gene)) %>%
    filter(gene != "", !str_detect(gene, "LOC|LINC"))
  
  unique_genes <- all_significant_genes %>%
    group_by(gene) %>%
    summarise(pop_count = n_distinct(population)) %>%
    filter(pop_count == 1)
  
  gene_labels <- all_significant_genes %>%
    mutate(is_unique = gene %in% unique_genes$gene) %>%
    group_by(gene, population) %>%
    arrange(desc(!!sym(y_var))) %>%
    slice(1) %>%
    ungroup()
  
  chrom_centers <- chrom_offsets %>%
    mutate(center = offset + chr_len / 2,
           center_mb = center / 1e6)
  
  data <- data %>%
    mutate(
      chr_parity = ifelse(chr_num %% 2 == 0, "even", "odd"),
      point_color = case_when(
        chr_parity == "even"  ~ "grey70",
        population == "AFR"   ~ "#FFCC33",
        population == "AMR"   ~ "#FF3333",
        population == "EAS"   ~ "#CCFF66", 
        population == "SAS"   ~ "#FF33FF",
        population == "EUR"   ~ "#66CCCC",
        TRUE                  ~ "black"
      )
    )
  
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
      labels = c("TRUE" = "unique", "FALSE" = "shared")
    ) +
    
    facet_wrap(~population, ncol = 1, scales = "free_y") +
    theme_classic() +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_line(linewidth = 0.2, color = "grey90"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
      axis.text.y = element_text(size = 7),
      strip.text = element_text(size = 8, face = "bold"),
      strip.background = element_rect(fill = "grey85", color = NA),
      panel.spacing = unit(1, "lines")
    ) +
    labs(x = "Genomic Position (Mb)", y = y_label, title = title)
}


# A custom function to format log-axis labels with rounded exponents
label_log_rounded <- function(x) {
  # Take the log10 of the break value
  log_x <- log10(x)
  # Round the exponent to two decimal places
  rounded_log_x <- round(log_x, 2)
  parse(text = paste0("10^", rounded_log_x))
}

# Define the quantile for the threshold
quantile_to_use <- 0.9990

# --- Plot for T_within_myr (uses the default label colors) ---
plot_avg_pop <- create_population_genome_plot(
  coalescence_df_with_features2, 
  "T_within_myr", 
  quantile_to_use,
  "Genome-wide Pairwise Coalescent Time, TMRCA (99.90th Percentile)",
  "Pairwise TMRCA (Mya)"
)

# --- Plot for Tpooled_Twithin_ratio (with inverted label colors) ---
plot_ratio_pop <- create_population_genome_plot(
  coalescence_df_with_features2, 
  "Tpooled_Twithin_ratio", 
  quantile_to_use,
  "Genome-wide Tpooled/Twithin (99.90th Percentile)",
  "log10(Tpooled/Twithin)",
  # Pass the inverted color scheme: Shared = red, Unique = black
  label_colors = c("TRUE" = "black", "FALSE" = "red") 
)

# Add the log scale with the fixed upper limit and rounded labels
plot_ratio_pop_log_formatted <- plot_ratio_pop + 
  scale_y_log10(
    #limits = c(NA, 10^2), # Set top limit only
    labels = label_log_rounded # Use the custom rounding function
  )

#ggsave("humans/humans_genome_plot_avg_pairwise_by_pop.png", plot_avg_pop, width = 10, height = 8, dpi = 300)
#ggsave("humans/humans_genome_plot_tpooled_ratio_by_pop.png", plot_ratio_pop_log_formatted, width = 10, height = 8, dpi = 300)

ggsave("humans/humans_genome_plot_avg_pairwise_by_pop_99.90th.pdf", plot_avg_pop, width = 10, height=8)
ggsave("humans/humans_genome_plot_tpooled_ratio_by_pop_90.90.pdf", plot_ratio_pop_log_formatted, width = 10, height = 8)



# --- ZOOMED  ---
plot_zoomed_tmrca <- function(data,
                              target_chr,
                              target_start,
                              target_end,
                              y_var,
                              y_label,
                              plot_title,
                              label_threshold = 0) {
  
  # Filter data for the region
  region_df <- data %>%
    filter(
      chromosome == target_chr,
      start <= target_end,
      end >= target_start
    )
  
  if (nrow(region_df) == 0) {
    warning(paste("No data found for the specified region:", target_chr, target_start, "-", target_end))
    return(ggplot() + labs(title = plot_title, subtitle = "No data available") + theme_minimal())
  }
  
  # Prepare gene labels with internal data cleaning
  gene_labels <- region_df %>%
    filter(!!sym(y_var) >= label_threshold, genes != "." & !is.na(genes)) %>%
    
    # ** THE FIX: Process each row to handle and fix mismatched lists **
    rowwise() %>%
    mutate(
      # a. Split each column's string into a list of items
      genes_list = list(str_split(genes, "[,;]")[[1]]),
      desc_list = list(str_split(descriptions, "[,;]")[[1]]),
      bio_list = list(str_split(gene_biotypes, "[,;]")[[1]]),
      
      # b. Use the number of genes as the "master" length
      master_len = length(genes_list),
      
      # c. Pad the shorter lists with "NA" to match the master length
      desc_padded = list(c(desc_list, rep("NA", master_len - length(desc_list)))),
      bio_padded = list(c(bio_list, rep("NA", master_len - length(bio_list)))),
      
      # d. Paste the now-parallel lists back into corrected strings using a single separator
      genes = paste(genes_list, collapse = ","),
      descriptions = paste(desc_padded, collapse = ","),
      gene_biotypes = paste(bio_padded, collapse = ",")
    ) %>%
    ungroup() %>% # End the row-wise operation
    
    # This will now work correctly with the padded strings
    separate_rows(genes, descriptions, gene_biotypes, sep = ",") %>%
    
    mutate(
      gene = trimws(genes),
      description = trimws(descriptions)
    ) %>%
    filter(gene != "") %>%
    mutate(display_label = gene) %>%
    group_by(gene) %>%
    slice_max(order_by = !!sym(y_var), n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # Create the plot
  p <- ggplot(region_df, aes(x = start, y = !!sym(y_var), color = population)) +
    geom_hline(yintercept = label_threshold, linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_line(linewidth = 0.4, alpha = 0.6) +
    geom_text_repel(
      data = gene_labels,
      aes(label = display_label),
      nudge_y = 1.5, size = 2, segment.color = "grey50", color = "black",
      box.padding = 0.5, max.overlaps = Inf
    ) +
    scale_color_manual(
      name = "Population",
      values = c("AFR"="#FFCC33", "AMR"="#FF3333", "EAS"="#CCFF66", "EUR"="#66CCCC", "SAS"="#FF33FF")
    ) +
    scale_x_continuous(
      labels = scales::label_number(scale = 1e-6, suffix = " Mb"),
      limits = c(target_start, target_end)
    ) +
    labs(
      title = plot_title,
      x = paste("Genomic Position on", target_chr),
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

# --- Call The Function For Both Metrics ---
# This code will now run without errors.
############ HLA ########################
my_chr <- "chr6"
my_start <- 29507887
my_end <- 33110601
my_label_threshold <- 8

tmrca_within_plot <- plot_zoomed_tmrca(
  data = coalescence_df_with_features,
  target_chr = my_chr,
  target_start = my_start,
  target_end = my_end,
  y_var = "T_within_myr",
  y_label = "Pairwise TMRCA (Mya)",
  plot_title = "HLA Region",
  label_threshold = my_label_threshold
)
#avg_tmrca_myr
#tmrca_within_plot2 <- plot_zoomed_tmrca(
#  data = coalescence_df_with_features,
#  target_chr = my_chr,
#  target_start = my_start,
#  target_end = my_end,
#  y_var = "avg_tmrca_myr",
#  y_label = "Average TMRCA (Mya)",
#  plot_title = "",
#  label_threshold = my_label_threshold
#)
tmrca_ratio_plot <- plot_zoomed_tmrca(
  data = coalescence_df_with_features,
  target_chr = my_chr,
  target_start = my_start,
  target_end = my_end,
  y_var = "Tpooled_Twithin_ratio", # Changed to the ratio variable
  y_label = "Tpooled/Twithin",   # Updated y-axis label
  plot_title = "HLA region", # Updated title
  label_threshold = my_label_threshold
)
HLA<-tmrca_within_plot + tmrca_ratio_plot +
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom')
HLA
ggsave("humans/HLA_humans.pdf", HLA, width=8, height=4)



############ GYP ########################
my_chr <- "chr4"
my_start <- 147000000
my_end <- 147700000
my_label_threshold <- 4

tmrca_within_plot <- plot_zoomed_tmrca(
  data = coalescence_df_with_features,
  target_chr = my_chr,
  target_start = my_start,
  target_end = my_end,
  y_var = "T_within_myr",
  y_label = "Pairwise TMRCA (Mya)",
  plot_title = "GYP Region",
  label_threshold = my_label_threshold
)
tmrca_ratio_plot <- plot_zoomed_tmrca(
  data = coalescence_df_with_features,
  target_chr = my_chr,
  target_start = my_start,
  target_end = my_end,
  y_var = "Tpooled_Twithin_ratio", # Changed to the ratio variable
  y_label = "Tpooled/Twithin",   # Updated y-axis label
  plot_title = "GY Pregion", # Updated title
  label_threshold = my_label_threshold
)
GYP<-tmrca_within_plot + tmrca_ratio_plot +
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom')
GYP
ggsave("humans/GYP_humans.pdf", GYP, width=7, height=4)



############ HBE ########################
my_chr <- "chr11"
my_start <- 5250000
my_end <- 5350000
my_label_threshold <- 1

tmrca_within_plot <- plot_zoomed_tmrca(
  data = coalescence_df_with_features,
  target_chr = my_chr,
  target_start = my_start,
  target_end = my_end,
  y_var = "T_within_myr",
  y_label = "Pairwise TMRCA (Mya)",
  plot_title = "HBE1 TMRCA",
  label_threshold = my_label_threshold
)
tmrca_ratio_plot <- plot_zoomed_tmrca(
  data = coalescence_df_with_features,
  target_chr = my_chr,
  target_start = my_start,
  target_end = my_end,
  y_var = "Tpooled_Twithin_ratio", # Changed to the ratio variable
  y_label = "Tpooled/Twithin",   # Updated y-axis label
  plot_title = "HBE1 Ratio", # Updated title
  label_threshold = my_label_threshold
)
HBE1<-tmrca_within_plot + tmrca_ratio_plot +
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom')
HBE1
ggsave("humans/HBE1_humans.pdf", HBE1, width=7, height=4)


############ HBA ########################
my_chr <- "chr16"
my_start <- 155000
my_end <- 175000
my_label_threshold <- 1

tmrca_within_plot <- plot_zoomed_tmrca(
  data = coalescence_df_with_features,
  target_chr = my_chr,
  target_start = my_start,
  target_end = my_end,
  y_var = "T_within_myr",
  y_label = "Pairwise TMRCA (Mya)",
  plot_title = "HBA TMRCA",
  label_threshold = my_label_threshold
)
tmrca_ratio_plot <- plot_zoomed_tmrca(
  data = coalescence_df_with_features,
  target_chr = my_chr,
  target_start = my_start,
  target_end = my_end,
  y_var = "Tpooled_Twithin_ratio", # Changed to the ratio variable
  y_label = "Tpooled/Twithin",   # Updated y-axis label
  plot_title = "HBA Ratio", # Updated title
  label_threshold = my_label_threshold
)
HBA<-tmrca_within_plot + tmrca_ratio_plot +
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom')
HBA
ggsave("humans/HBA_humans.pdf", HBA, width=7, height=4)



############ FOXK2 ########################
my_chr <- "chr17"
my_start <- 83300000
my_end <- 83700000

my_label_threshold <- 1

tmrca_within_plot <- plot_zoomed_tmrca(
  data = coalescence_df_with_features,
  target_chr = my_chr,
  target_start = my_start,
  target_end = my_end,
  y_var = "T_within_myr",
  y_label = "Pairwise TMRCA (Mya)",
  plot_title = "FOXK2 TMRCA",
  label_threshold = my_label_threshold
)
tmrca_ratio_plot <- plot_zoomed_tmrca(
  data = coalescence_df_with_features,
  target_chr = my_chr,
  target_start = my_start,
  target_end = my_end,
  y_var = "Tpooled_Twithin_ratio", # Changed to the ratio variable
  y_label = "Tpooled/Twithin",   # Updated y-axis label
  plot_title = "FOXK2  Ratio", # Updated title
  label_threshold = my_label_threshold
)
FOXK2 <-tmrca_within_plot + tmrca_ratio_plot +
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom')
FOXK2 
ggsave("humans/FOXK2_humans.pdf",FOXK2 , width=7, height=4)


############ ATP9B ########################
my_chr <- "chr18"
my_start <- 79300000
my_end <- 79600000


my_label_threshold <- 3.2

tmrca_within_plot <- plot_zoomed_tmrca(
  data = coalescence_df_with_features,
  target_chr = my_chr,
  target_start = my_start,
  target_end = my_end,
  y_var = "T_within_myr",
  y_label = "Pairwise TMRCA (Mya)",
  plot_title = "ATP9B TMRCA",
  label_threshold = my_label_threshold
)
tmrca_ratio_plot <- plot_zoomed_tmrca(
  data = coalescence_df_with_features,
  target_chr = my_chr,
  target_start = my_start,
  target_end = my_end,
  y_var = "Tpooled_Twithin_ratio", # Changed to the ratio variable
  y_label = "Tpooled/Twithin",   # Updated y-axis label
  plot_title = "ATP9B  Ratio", # Updated title
  label_threshold = my_label_threshold
)
ATP9B <-tmrca_within_plot + tmrca_ratio_plot +
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom')
ATP9B
ggsave("humans/ATP9B_humans.pdf",ATP9B , width=7, height=4)




############ KIRs ########################
my_chr <- "chr19"
my_start <- 57700000
my_end <- 59161007


my_label_threshold <- 8
tmrca_within_plot <- plot_zoomed_tmrca(
  data = coalescence_df_with_features,
  target_chr = my_chr,
  target_start = my_start,
  target_end = my_end,
  y_var = "T_within_myr",
  y_label = "Pairwise TMRCA (Mya)",
  plot_title = "KIR and NLRP Region",
  label_threshold = my_label_threshold
)
tmrca_within_plot
#avg_tmrca_myr
#tmrca_within_plot2 <- plot_zoomed_tmrca(
#  data = coalescence_df_with_features,
#  target_chr = my_chr,
#  target_start = my_start,
#  target_end = my_end,
#  y_var = "avg_tmrca_myr",
#  y_label = "Average TMRCA (Mya)",
#  plot_title = "",
#  label_threshold = my_label_threshold
#)
my_label_threshold <- 6
tmrca_ratio_plot <- plot_zoomed_tmrca(
  data = coalescence_df_with_features,
  target_chr = my_chr,
  target_start = my_start,
  target_end = my_end,
  y_var = "Tpooled_Twithin_ratio", # Changed to the ratio variable
  y_label = "Tpooled/Twithin",   # Updated y-axis label
  plot_title = "KIR and NLRP Region", # Updated title
  label_threshold = my_label_threshold
)
tmrca_ratio_plot
KIR<-tmrca_within_plot + tmrca_ratio_plot +
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom')
KIR
ggsave("humans/KIR_humans.pdf", KIR, width=8, height=4)




############ KANSL########################
my_chr <- "chr17"
my_start <- 47055163
my_end <- 	47087028


my_label_threshold <- 2
tmrca_within_plot <- plot_zoomed_tmrca(
  data = coalescence_df_with_features,
  target_chr = my_chr,
  target_start = my_start,
  target_end = my_end,
  y_var = "T_within_myr",
  y_label = "Pairwise TMRCA (Mya)",
  plot_title = "KANSL1",
  label_threshold = my_label_threshold
)
tmrca_within_plot
#avg_tmrca_myr
#tmrca_within_plot2 <- plot_zoomed_tmrca(
#  data = coalescence_df_with_features,
#  target_chr = my_chr,
#  target_start = my_start,
#  target_end = my_end,
#  y_var = "avg_tmrca_myr",
#  y_label = "Average TMRCA (Mya)",
#  plot_title = "",
#  label_threshold = my_label_threshold
#)
my_label_threshold <- 2
tmrca_ratio_plot <- plot_zoomed_tmrca(
  data = coalescence_df_with_features,
  target_chr = my_chr,
  target_start = my_start,
  target_end = my_end,
  y_var = "Tpooled_Twithin_ratio", # Changed to the ratio variable
  y_label = "Tpooled/Twithin",   # Updated y-axis label
  plot_title = "KANSL1", # Updated title
  label_threshold = my_label_threshold
)
tmrca_ratio_plot
KANSL<-tmrca_within_plot + tmrca_ratio_plot +
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom')
KANSL
ggsave("humans/KANSL_humans.pdf", KANSL, width=8, height=4)