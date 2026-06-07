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

##NaN for avg_tmrca: 
#This happens when `res["tmrca_count"]` is zero, meaning no tree overlapped the window in any ARG sample.

coalescence_df <- read_csv("bonobos/genome-wide_metrics_annotated_noCensat.csv", quote = "\"") %>%
  rename(chrom = chromosome)


annotation_df <- fread("bonobos/mPanPan1_gene.bed.gz") %>%
  rename(genes = name) %>%
  select(genes, gene_biotype, description) %>%
  distinct()

#Perform the left join to add the annotation data.
coalescence_df <- coalescence_df %>%
  separate_rows(genes, sep = ",") %>%
  left_join(annotation_df, by = "genes") %>%
  group_by(chrom, start, end) %>%
  summarise(
    avg_tmrca = first(avg_tmrca),
    avg_pairwise_coalescence_time  = first(avg_pairwise_coalescence_time),
    genes = paste(genes, collapse = ","),
    gene_biotypes = paste(na.omit(gene_biotype), collapse = ";"),
    descriptions = paste(na.omit(description), collapse = ";")
  ) %>%
  ungroup()

# Convert generations to millions of years
coalescence_df$avg_tmrca_myr <- coalescence_df$avg_tmrca * 25 / 1e6
coalescence_df$avg_pairwise_myr <- coalescence_df$avg_pairwise_coalescence_time * 25 / 1e6


features_df <- read_tsv(
  "bonobos/mPanPan1_genomefeatures.bed",
  col_names = c("chrom", "start", "end", "feature"),
  col_types = "cddc", # Specify: character, double, double, character
  col_select = 1:4
)
features_df <- semi_join(features_df, coalescence_df, by = "chrom")


tr_df <- read_tsv(
  "bonobos/panpan_catalog.no_overlaps_simp.bed",
  col_names = c("chrom", "start", "end", "motif"),
  col_types = "cddc", # Specify: character, double, double, character
) %>%
  mutate(feature = "TR")

sedef_df <- read_tsv(
  "bonobos/mPanPan1_sedefSegDups.bed",
  col_names = c("chrom", "start", "end", "coordinates"),
  col_types = "cddc", # Specify: character, double, double, character
  col_select = 1:4 # Select first 4 columns
) %>%
  mutate(feature = "SEDEF_SD") # Add the feature name
sedef_df <- semi_join(sedef_df, coalescence_df, by = "chrom")

censat_df <- read_tsv(
  "bonobos/mPanPan1_CenSat.bed",
  col_names = c("chrom", "start", "end", "censat"),
  col_types = "cddc",
  col_select = 1:4
)

rmsk_df <- read_tsv(
  "bonobos/mPanPan1_RepeatMasker.bed",
  col_names = c("chrom", "start", "end", "feature"),
  col_types = "cddc",
  col_select = 1:4
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

# --- 3. Join Results Back to the ORIGINAL Dataframe and preprocess data and FILTER OUT CHR X/Y---
coalescence_df_with_features <- coalescence_df %>%
  left_join(feature_annotations, by = c("chrom", "start", "end")) %>%
  left_join(tr_overlap_summary, by = c("chrom", "start", "end")) %>%
  left_join(sedef_overlap_summary, by = c("chrom", "start", "end")) %>%
  left_join(censat_overlap_summary, by = c("chrom", "start", "end")) %>%
  left_join(rmsk_overlap_summary, by = c("chrom", "start", "end")) %>%
  left_join(features_overlap_summary, by = c("chrom", "start", "end")) %>%
  mutate(
    feature_annotation = if_else(is.na(feature_annotation), ".", feature_annotation),
    tr_overlap_bp = coalesce(tr_overlap_bp, 0L),
    sedef_overlap_bp = coalesce(sedef_overlap_bp, 0L),
    censat_overlap_bp = coalesce(censat_overlap_bp, 0L),
    rmsk_overlap_bp = coalesce(rmsk_overlap_bp, 0L),
    features_overlap_bp = coalesce(features_overlap_bp, 0L),
    window_size = end - start,
    tr_overlap_percent = (tr_overlap_bp / window_size) * 100,
    sedef_overlap_percent = (sedef_overlap_bp / window_size) * 100,
    censat_overlap_percent = (censat_overlap_bp / window_size) * 100,
    rmsk_overlap_percent = (rmsk_overlap_bp / window_size) * 100,
    features_overlap_percent = (features_overlap_bp / window_size) * 100
  ) %>%
  select(
    chrom, start, end,
    avg_tmrca_myr, avg_pairwise_myr, 
    genes, feature_annotation, gene_biotypes, descriptions,
    tr_overlap_percent, 
    sedef_overlap_percent,
    censat_overlap_percent,
    rmsk_overlap_percent,
    features_overlap_percent
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



# Track initial number of windows
#original_window_count <- nrow(coalescence_df_with_features)
# Identify and track all problematic windows
#problematic_windows <- coalescence_df_with_features  %>%
#  mutate(
#    issue_type = case_when(
#      chr_clean %in% c("chrX", "chrY") ~ "X/Y_chromosome",
#      is.na(avg_tmrca) | is.nan(avg_tmrca) ~ "NA/Nan_TMRCA",
#      avg_tmrca_myr < 0 ~ "Negative_TMRCA",
#      avg_pairwise_myr < 0 ~ "Negative_pairwise",
#      TRUE ~ "Valid"
#    )
#  ) %>%
#  filter(issue_type != "Valid")

# Calculate proportions
#removed_count <- nrow(problematic_windows)
#proportion_removed <- removed_count / original_window_count

# Breakdown of issues
#issue_summary <- problematic_windows %>%
#  count(issue_type) %>%
#  mutate(proportion = n / original_window_count)


# Remove problematic windows and chromosomes
#coalescence_df_with_features  <- coalescence_df_with_features  %>%
#  filter(
#    !chr_clean %in% c("chrX", "chrY"),
#    !(is.na(avg_tmrca) | is.nan(avg_tmrca)),
#    avg_tmrca_myr >= 0,
#    avg_pairwise_myr >= 0
#  )


############ SUMMARY STATS ON TRs and other features #####################
unique_windows_summary <- coalescence_df_with_features %>%
  group_by(chromosome, start, end, genes, feature_annotation) %>%
  summarise(
    # These overlap values are the same for all rows of a unique window
    tr_overlap_percent = first(tr_overlap_percent),
    sedef_overlap_percent = first(sedef_overlap_percent),
    censat_overlap_percent = first(censat_overlap_percent),
    rmsk_overlap_percent = first(rmsk_overlap_percent),
    avg_pairwise_myr = mean(avg_pairwise_myr, na.rm = TRUE),
  ) %>%
  ungroup() # Don't forget to ungroup

features_summary<-unique_windows_summary %>%
  pivot_longer(
    cols = c(tr_overlap_percent, censat_overlap_percent),
    names_to = "category",
    values_to = "percent"
  ) %>%
  filter(percent > 0) %>%
  ggplot(aes(x = percent, fill = category)) +
  geom_histogram(bins = 50, color = "white", show.legend = FALSE) +
  facet_wrap(~ category, scales = "free") +
  labs(
    x = "Overlap Percentage with 1kb Window (%)",
    y = "Number of Windows - Bonobos"
  ) +
  theme_minimal()
ggsave("bonobos/features_summary.pdf", features_summary, width=3, height=3) 

# Calculate chromosome offsets with chr_len included
chrom_offsets <- coalescence_df_with_features %>%
  group_by(chr_clean) %>%
  summarize(chr_len = max(end), .groups = "drop") %>%
  arrange(chr_clean) %>%
  # The as.numeric() here prevents the integer overflow
  mutate(offset = cumsum(lag(as.numeric(chr_len), default = 0))) %>%
  select(chr_clean, chr_len, offset)

# Add genome positions to main data frame
coalescence_df_with_features<- coalescence_df_with_features %>%
  left_join(chrom_offsets, by = "chr_clean") %>%
  mutate(genome_pos = start + offset,
         genome_pos_mb = genome_pos / 1e6)
coalescence_df_with_features  <- coalescence_df_with_features %>%
  mutate(midpoint = (start + end) / 2) %>%
  arrange(chr_clean, start) %>%
  mutate(row_index = row_number())
axisdf <- coalescence_df_with_features  %>%
  group_by(chr_clean) %>%
  summarize(center = mean(row_index))


full_annotation_df <- fread("bonobos/mPanPAn1_gene.bed.gz") 
GYP_test <-coalescence_df_with_features  %>%
  filter(chromosome == "chr3_pat_hsa4") %>%
  #filter(str_detect(descriptions, "MNS blood group")) %>%
  filter(start >= 141449208) %>%
  filter(start <= 142620159) %>%
  arrange(desc(avg_pairwise_myr))
View(GYP_test)
write_csv(GYP_test, "bonobos/GYPtest.csv")


###### RANKED GENES PLOT
#significant_genes2 <- coalescence_df_with_features %>%
#  filter(avg_pairwise_myr > threshold, genes != ".") %>%
#  separate_rows(genes, sep = ";") %>%
#  filter(genes != "", genes != ".") %>%
#  group_by(genes) %>%
#  slice_max(avg_pairwise_myr, n = 1, with_ties = FALSE) %>%
#  ungroup() %>%
#  arrange(desc(avg_pairwise_myr)) %>%
#  mutate(
#    rank = row_number(),
#    genomic_pos = paste0(chr_clean, ":", round(start/1e6, 2), "Mb"),
#    gene_label = paste0(genes, " (", feature_annotation, ")") 
#  )
#ranked_plot<-ggplot(significant_genes2 %>% filter(avg_pairwise_myr > 10),  # Filter main data too
#              aes(x = reorder(genomic_pos, -avg_pairwise_myr), 
#                  y = avg_pairwise_myr)
#              ) +
#  geom_point(color = "salmon", size = 1) +
#  geom_text_repel(
#    aes(label = gene_label),
#    size = 2,
#    box.padding = 0.5,
#    max.overlaps = Inf
#  )  +
#  scale_x_discrete(name = "Genomic Position (Chromosome:Mb)") +
#  scale_y_continuous(name = "Average Pairwise Coalescence Time (Myr)") +
#  labs(title = "Bonobo Top Genes ranked by Highest TMRCA",
#       subtitle = paste("Threshold =", round(threshold, 2), "Myr")) +
#  theme_classic() +
#  theme(
#    axis.text.x = element_blank(),
#    panel.grid.major.x = element_blank()
#  )
#ranked_plot
#ggsave("bonobos/bonobos_dcreasing_avgpawise_tmrca.png", ranked_plot, dpi=300, width=12, height=6)


### MANHATTAN FILTERED! ######
############## FILTER FOR COMPLEX NOISY REGIONS ###########
unique(coalescence_df_with_features$feature_annotation)
#pattern_to_exclude <- "SAT|CEN"
coalescence_df_with_features2 <- coalescence_df_with_features %>%
  filter(
    #avg_pairwise_myr >= 0,
    censat_overlap_percent == 0,
    tr_overlap_percent <= 5,
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


unique(coalescence_df_with_features2$feature_annotation)
# CASES KEPT FOR NOW include combinations of windows overlapped by
#"Telo"                  
# SEDEF_SD
coalescence_df_with_features2  <- coalescence_df_with_features2 %>%
mutate(
  population = "Bonobo")



#################### GET TABLES OF TOP GENES ##############

get_pop_significant_genes <- function(data, metric_var, threshold_quantile) {
  pop_thresholds <- data %>%
    group_by(population) %>%
    summarise(
      threshold = quantile({{ metric_var }}, probs = threshold_quantile, na.rm = TRUE),
      .groups = "drop"
    )
  significant_genes <- data %>%
    left_join(pop_thresholds, by = "population") %>%
    filter({{ metric_var }} > threshold, genes != ".") %>%
    separate_rows(genes, sep = ";") %>%
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
      feature_annotation,
      gene_biotypes,
      descriptions,
      avg_tmrca_myr, 
      avg_pairwise_myr,
      tr_overlap_percent,
      sedef_overlap_percent,
      censat_overlap_percent,
      rmsk_overlap_percent,
      midpoint,
      row_index
    )
  
  return(significant_genes)
}

# Define the quantile to use for significance
quantile_to_use <- 0.9999


# Get significant genes based on T_within_myr for each population
significant_genes_avg_pairwise_myr <- get_pop_significant_genes(
  data = coalescence_df_with_features2,
  metric_var = avg_pairwise_myr,
  threshold_quantile = quantile_to_use # Corrected argument name
) %>%
  arrange(desc(avg_pairwise_myr))

# Get significant genes based on Tpooled_Twithin_ratio for each population
significant_genes_avg_tmrca_myr<- get_pop_significant_genes(
  data = coalescence_df_with_features2,
  metric_var = avg_tmrca_myr,
  threshold_quantile = quantile_to_use # Corrected argument name
) %>%
  arrange(desc(avg_tmrca_myr))

View(significant_genes_avg_tmrca_myr)
View(significant_genes_avg_pairwise_myr)

write_csv(significant_genes_avg_pairwise_myr,"bonobos/bonobos_high_avg_pairwise_TMRCA_significant_genes.csv")
write_csv(significant_genes_avg_tmrca_myr, "bonobos/bonobos_high_avg_TMRCA_significant_genes.csv")


########### PLOTTTING MANHATTANS #############################
# Process remaining data
coalescence_df_with_features2  <- coalescence_df_with_features2 %>%
  arrange(chr_clean, start) %>%
  mutate(row_index = row_number())

# Recalculate chromosome label positions
axisdf <- coalescence_df_with_features2 %>%
  group_by(chr_clean) %>%
  summarize(center = mean(row_index))




plot_data_ranked_avg_tmrca_myr<- coalescence_df_with_features2 %>%
  arrange(desc(avg_tmrca_myr)) %>%
  filter(avg_tmrca_myr > 6)  %>%
  select(chromosome, start, end, population, avg_tmrca_myr, genes, gene_biotypes, descriptions, feature_annotation)
write_csv(plot_data_ranked_avg_tmrca_myr, "bonobos/ranked_avg_tmrca_myr.csv")


plot_data_ranked_avg_pairwise_myr <- coalescence_df_with_features2 %>%
  arrange(desc(avg_pairwise_myr)) %>%
  filter(avg_pairwise_myr > 6) %>% 
  select(chromosome, start, end, population, avg_pairwise_myr, genes, gene_biotypes, descriptions, feature_annotation)
write_csv(plot_data_ranked_avg_pairwise_myr , "bonobos/ranked_avg_pairwise_myr.csv")

# Custom population colors
pop_colors <- c(
  "Bonobo" = "salmon"
)


# 1. Chromosome-faceted plots with FIXED Y-AXIS
pop_colors <- c(
  "Bonobo" = "salmon"
)

# 1. Chromosome-faceted plots with FIXED Y-AXIS and conditional coloring
create_faceted_point_plot <- function(data, y_var, threshold, title, y_label) {
  data <- data %>%
    mutate(color_condition = ifelse(!!sym(y_var) > threshold, "Above", "Below"))
  
  gene_labels <- data %>%
    filter(!!sym(y_var) > threshold, !is.na(genes)) %>%
    mutate(gene = strsplit(as.character(genes), ",")) %>%
    unnest(gene) %>%
    mutate(gene = trimws(gene)) %>%
    filter(gene != "") %>%
    group_by(gene, chr_clean) %>%
    arrange(desc(!!sym(y_var))) %>%
    slice(1) %>%
    ungroup()
  
  ggplot(data, aes(x = position_mb, y = !!sym(y_var))) +
    
    # **KEY CHANGE**: Replaced geom_point() with geom_point_rast()
    ggrastr::geom_point_rast(aes(color = color_condition), size = 1, alpha = 0.7, dpi = 300) +
    
    geom_hline(yintercept = threshold, linetype = "dashed", color = "red", linewidth = 0.5) +
    geom_text_repel(
      data = gene_labels,
      aes(label = gene),
      size = 2,
      color = "black",
      min.segment.length = 0.1,
      segment.size = 0.2,
      box.padding = 0.3,
      max.overlaps = Inf,
      force = 1,
      force_pull = 0.5
    ) +
    scale_color_manual(values = c("Above" = "salmon", "Below" = "lightgrey"),
                       guide = "none") +
    facet_wrap(~chr_clean, scales = "free_x", nrow = 6) +
    theme_classic() +
    theme(
      legend.position = "none",
      panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
      axis.text.y = element_text(size = 8),
      strip.text = element_text(size = 8, face = "bold"),
      strip.background = element_rect(fill = "grey85", color = NA),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.4, "cm"),
      panel.spacing = unit(0.5, "lines")
    ) +
    labs(x = "Position (Mb)", y = y_label, title = title) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
}



create_population_genome_plot <- function(data, y_var, threshold_quantile, title, y_label, y_max = NULL) {
  
  # --- 1. Calculate the threshold for EACH population ---
  pop_thresholds <- data %>%
    group_by(population) %>%
    summarise(
      threshold = quantile(!!sym(y_var), probs = threshold_quantile, na.rm = TRUE),
      .groups = "drop"
    )
  
  # --- 2. (REVISED) Logic to label each gene only once at its max value ---
  gene_labels <- data %>%
    left_join(pop_thresholds, by = "population") %>%
    filter(!!sym(y_var) > threshold) %>%
    mutate(gene = strsplit(as.character(genes), ",")) %>%
    unnest(gene) %>%
    mutate(gene = trimws(gene)) %>%
    filter(!is.na(gene) & gene != "" & gene != ".") %>%
    group_by(gene, population) %>%
    slice_max(order_by = !!sym(y_var), n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # Calculate chromosome center points for x-axis labels
  chrom_centers <- chrom_offsets %>%
    mutate(center = offset + chr_len / 2,
           center_mb = center / 1e6)
  
  # --- 3. Create the point_color column using your original case_when ---
  data <- data %>%
    mutate(
      chr_parity = ifelse(chr_num %% 2 == 0, "even", "odd"),
      point_color = case_when(
        chr_parity == "even"  ~ "grey70",
        population == "Bonobo" ~  "salmon",
        TRUE                  ~ "black"
      )
    )
  
  # --- 4. Generate the plot ---
  p <- ggplot(data, aes(x = genome_pos_mb, y = !!sym(y_var))) +
    geom_vline(
      xintercept = c(0, chrom_offsets$offset + chrom_offsets$chr_len) / 1e6, 
      color = "grey80", 
      linewidth = 0.3
    ) +
    
    # **KEY CHANGE**: Replaced geom_point() with geom_point_rast()
    ggrastr::geom_point_rast(aes(color = point_color), size = 0.8, alpha = 0.7, show.legend = FALSE, dpi = 300) +
    
    geom_hline(
      data = pop_thresholds,
      aes(yintercept = threshold), 
      linetype = "dashed", 
      color = "red", 
      linewidth = 0.5
    ) +
    geom_text_repel(
      data = gene_labels,
      aes(label = gene),
      size = 2, color = "black", min.segment.length = 0.1,
      segment.size = 0.2, box.padding = 0.3, max.overlaps = Inf # Allow all labels
    ) +
    scale_x_continuous(
      breaks = chrom_centers$center_mb,
      labels = chrom_centers$chr_clean
    ) +
    scale_color_identity() +  
    facet_wrap(~population, ncol = 1, scales = "free_y") +  
    theme_classic() +
    theme(
      legend.position = "none",
      panel.grid.major.y = element_line(linewidth = 0.2, color = "grey90"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
      axis.text.y = element_text(size = 8),
      strip.text = element_text(size = 8, face = "bold"),
      strip.background = element_rect(fill = "grey85", color = NA),
      panel.spacing = unit(1, "lines")
    ) +
    labs(x = "Genomic Position (Mb)", y = y_label, title = title)
  
  # --- 5. Conditionally apply the y-axis limit ---
  if (!is.null(y_max)) {
    p <- p + coord_cartesian(ylim = c(NA, y_max))
  }
  
  return(p)
}
# Define the quantile for the threshold
quantile_to_use <- 0.9999
plot_avg_pop_bonobo <- create_population_genome_plot(
  coalescence_df_with_features2, 
  "avg_pairwise_myr", 
  quantile_to_use, # Pass the quantile here
  "",
  "Pairwise TMRCA (Mya)"
)
plot_avg_pop_bonobo

ggsave("bonobos/all_bonobos_avg_pairwise_tmrca_gwide.png", plot_avg_pop_bonobo, width = 10, height= 3, dpi=300)
ggsave("bonobos/all_bonobos_avg_pairwise_tmrca_gwide.pdf", plot_avg_pop_bonobo, width = 10, height= 3)

plot_avg_points <- create_faceted_point_plot(
  coalescence_df_with_features2, "avg_pairwise_myr", 9,
  "",
  "Pairwise TMRCA (Mya)"
)
plot_avg_points
ggsave("bonobos/pointplot_avg_pairwise_fixedY.png", plot_avg_points, width = 7, height = 9, dpi = 300)
ggsave("bonobos/pointplot_avg_pairwise_fixedY.pdf", plot_avg_points, width = 7, height = 9)



# ZOOM
plot_zoomed_tmrca_single_pop <- function(data,
                                         target_chr,
                                         target_start,
                                         target_end,
                                         y_var,
                                         y_label,
                                         plot_title,
                                         label_threshold = 0) {
  
  # 1. Filter data for the region
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
  
  # 2. Prepare gene labels, with corrected data cleaning
  gene_labels <- region_df %>%
    filter(!!sym(y_var) >= label_threshold, genes != "." & !is.na(genes)) %>%
    rowwise() %>%
    mutate(
      genes_vec = str_split(genes, "[,;]")[[1]],
      desc_vec = str_split(descriptions, "[,;]")[[1]],
      bio_vec = str_split(gene_biotypes, "[,;]")[[1]],
      master_len = length(genes_vec),
      desc_padded = list(c(desc_vec, rep("NA", master_len - length(desc_vec)))),
      bio_padded = list(c(bio_vec, rep("NA", master_len - length(bio_vec)))),
      genes = paste(genes_vec, collapse = ","),
      descriptions = paste(desc_padded[[1]], collapse = ","),
      gene_biotypes = paste(bio_padded[[1]], collapse = ",")
    ) %>%
    ungroup() %>%
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
  
  # 3. Create the plot
  p <- ggplot(region_df, aes(x = start, y = !!sym(y_var))) +
    geom_hline(yintercept = label_threshold, linetype = "dashed", color = "red", linewidth = 0.8) +
    
    # **KEY CHANGE**: Set the line color directly, since there's only one
    geom_line(color = "salmon", linewidth = 0.4, alpha = 0.6) +
    
    geom_text_repel(
      data = gene_labels,
      aes(label = display_label),
      nudge_y = 1.5, size = 2, segment.color = "grey50", color = "black",
      box.padding = 0.5, max.overlaps = Inf
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
      # **KEY CHANGE**: Hide the legend, as it's not necessary
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    )
  
  return(p)
}



my_chr <- "chr3_pat_hsa4"
my_start <- 141449208
my_end <- 142620159
my_label_threshold <-6
tmrca_within_plot <- plot_zoomed_tmrca_single_pop(
  data = coalescence_df_with_features,
  target_chr = my_chr,
  target_start = my_start,
  target_end = my_end,
  y_var = "avg_pairwise_myr",
  y_label = "Pairwise TMRCA (Millions of Years)",
  plot_title = "Glycophorin Region",
  label_threshold = my_label_threshold
)
tmrca_within_plot