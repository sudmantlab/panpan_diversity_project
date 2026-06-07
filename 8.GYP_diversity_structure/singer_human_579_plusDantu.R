# Resolve the project root to whichever path actually exists on this machine.
# macOS Google Drive sometimes lives under ~/Google Drive (symlink) and
# sometimes under ~/Library/CloudStorage/GoogleDrive-<account>/My Drive.
.project_root_candidates <- c(
  "/Users/joanocha/Google Drive/My Drive/POSTDOC/PANPAN/analysis/Figure4_SINGER/singer_output_processing_scripts",
  "/Users/joanocha/Library/CloudStorage/GoogleDrive-joana.laranjeira.rocha@gmail.com/My Drive/POSTDOC/PANPAN/analysis/Figure4_SINGER/singer_output_processing_scripts"
)
.project_root <- .project_root_candidates[dir.exists(.project_root_candidates)][1]
if (is.na(.project_root)) stop("Could not locate project root from known candidates.")
setwd(.project_root)
message("Working directory set to: ", getwd())

# All PDFs from this run go here. Switch this single line to repoint everything.
.out_dir <- file.path(.project_root, "humans", "humans579_dantu_pops")
dir.create(.out_dir, recursive = TRUE, showWarnings = FALSE)
message("Outputs will be saved to: ", .out_dir)
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




coalescence_df <- fread("humans/humans579_dantu_pops/genome-wide_metrics_annotatedbyPop.csv") %>%
  rename(chrom = chromosome)

annotation_df <- fread("humans/ht2t_gene.bed.gz") %>%
  rename(genes = name) %>%
  select(genes,gene_biotype, description) %>%
  distinct()

full_annotation_df <- fread("humans/ht2t_gene.bed.gz")

# --- The per-population CSV does not carry gene annotations in its `genes`
# column. The per-superpop CSV uses identical windows and IS annotated, so we
# borrow its (chromosome, start, end) -> genes mapping rather than recomputing
# from the bed (much faster).
.superpop_csv <- "humans/humans579_dantu_superpops/genome-wide_metrics_annotatedbyPop.csv"
if (file.exists(.superpop_csv)) {
  superpop_window_anno <- fread(.superpop_csv, select = c("chromosome", "start", "end", "genes")) %>%
    distinct() %>%
    rename(chrom = chromosome)
  message("Loaded superpop window annotation: ", nrow(superpop_window_anno), " unique windows.")
} else {
  superpop_window_anno <- NULL
  message("Superpop annotation CSV not found at ", .superpop_csv,
          "; the `genes` column will remain whatever the input CSV provides.")
}

conflicts_prefer(dplyr::first)

# If the input CSV has no usable gene names but we loaded the superpop
# window-annotation above, fill `genes` from it now (same windows in both).
if (!is.null(superpop_window_anno) &&
    (!"genes" %in% names(coalescence_df) ||
     all(is.na(coalescence_df$genes) | coalescence_df$genes == ""))) {
  message("Filling `genes` from superpop window annotation.")
  coalescence_df <- coalescence_df %>%
    select(-any_of("genes")) %>%
    left_join(superpop_window_anno, by = c("chrom", "start", "end"))
}

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

  # Custom population colors - grouped by region.
  # Within each region, max lightness/saturation contrast for legibility.
  # AFR  -> yellow/orange/brown family
  # AMR  -> red family
  # EAS  -> green family
  # SAS  -> magenta/purple family
  # EUR  -> blue family
  # MENA -> tan/khaki/grey-teal (sits in WestEurasia facet alongside EUR)
  # ASK (Ashkenazi -> WestEurasia), ASL (African in St Louis -> Africa)
  pop_colors <- c(
    # AFR (9: includes ASL)
    "YRI" = "#FFD700",   # gold
    "ESN" = "#FFEC8B",   # light goldenrod
    "ACB" = "#FF7F00",   # bright orange
    "ASW" = "#D62728",   # vivid red-orange
    "GWD" = "#B8860B",   # dark goldenrod
    "MSL" = "#E6BE8A",   # light tan
    "LWK" = "#5D2F0E",   # very dark brown
    "MKK" = "#FF4500",   # orangered
    "ASL" = "#8B0000",   # dark red / maroon
    # AMR (4)
    "PEL" = "#E41A1C",   # red
    "CLM" = "#A50026",   # dark red
    "MXL" = "#FB6A4A",   # salmon
    "PUR" = "#67001F",   # very dark crimson
    # EAS (6)
    "CHB" = "#006400",   # dark green
    "CHS" = "#7CFC00",   # lawn green
    "JPT" = "#228B22",   # forest green
    "KHV" = "#98FB98",   # pale green
    "CDX" = "#3CB371",   # medium sea green
    "SouthKorea" = "#556B2F",  # dark olive green
    # SAS (6)
    "PJL" = "#FF00FF",   # bright magenta
    "BEB" = "#8B008B",   # dark magenta
    "GIH" = "#DDA0DD",   # plum
    "ITU" = "#FF1493",   # deep pink
    "STU" = "#9400D3",   # dark violet
    "India" = "#C71585", # medium violet red
    # EUR (4)
    "FIN" = "#000080",   # navy
    "GBR" = "#4169E1",   # royal blue
    "IBS" = "#1E90FF",   # dodger blue
    "TSI" = "#87CEEB",   # sky blue
    # MENA (8) - earth/teal tones to sit beside the blues in WestEurasia
    "Egypt" = "#A0522D",            # sienna
    "Morocco" = "#D2B48C",          # tan
    "Jordan" = "#BC8F8F",           # rosy brown
    "Syria" = "#8B7765",            # dark khaki
    "SaudiArabia" = "#A52A2A",      # brown
    "Oman" = "#708090",             # slate grey
    "Yemen" = "#2F4F4F",            # dark slate grey
    "UnitedArabEmirates" = "#5F9EA0", # cadet blue
    # Ashkenazi (WestEurasia)
    "ASK" = "#9370DB",   # medium purple
    # Superpopulations (kept so the same palette works for the *_superpops inputs)
    "Africa"      = "#FFCC33",
    "America"     = "#FF3333",
    "EastAsia"    = "#006600",
    "SouthAsia"   = "#FF33FF",
    "WestEurasia" = "#66CCCC",
    "AFR" = "#FFCC33",
    "AMR" = "#FF3333",
    "EAS" = "#006600",
    "SAS" = "#FF33FF",
    "EUR" = "#66CCCC"
  )

  # Population -> superpopulation (continent) lookup. Used to facet the
  # GYP plot by superpop. ASK = Ashkenazi (WestEurasia), ASL = African in
  # Saint Louis (Africa). MENA pops folded into WestEurasia.
  pop_to_superpop <- c(
    # Africa
    YRI = "Africa", ACB = "Africa", ASW = "Africa", ESN = "Africa",
    GWD = "Africa", LWK = "Africa", MKK = "Africa", MSL = "Africa",
    ASL = "Africa",
    # America
    PEL = "America", CLM = "America", MXL = "America", PUR = "America",
    # EastAsia
    CHB = "EastAsia", CHS = "EastAsia", JPT = "EastAsia", KHV = "EastAsia",
    CDX = "EastAsia", SouthKorea = "EastAsia",
    # SouthAsia
    PJL = "SouthAsia", BEB = "SouthAsia", GIH = "SouthAsia",
    ITU = "SouthAsia", STU = "SouthAsia", India = "SouthAsia",
    # WestEurasia (incl. MENA + Ashkenazi)
    FIN = "WestEurasia", GBR = "WestEurasia", IBS = "WestEurasia",
    TSI = "WestEurasia", Egypt = "WestEurasia", Morocco = "WestEurasia",
    Jordan = "WestEurasia", Syria = "WestEurasia",
    SaudiArabia = "WestEurasia", Oman = "WestEurasia",
    Yemen = "WestEurasia", UnitedArabEmirates = "WestEurasia",
    ASK = "WestEurasia"
  )
  


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
ggsave(file.path(.out_dir, "tmrca_chrm4_summary.pdf"), tmrca_plot, width = 5, height = 3)


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
ggsave(file.path(.out_dir, "tmrca_summary2.pdf"), tmrca_plot2, width = 5, height = 3)




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

quantile_to_use <- 0.9900
significant_genes_Twithin <- get_pop_significant_genes(
  data = coalescence_df_with_features2, # Use the corrected data frame
  metric_var = T_within_myr,
  threshold_quantile = quantile_to_use
) %>%
  arrange(desc(T_within_myr))

significant_genes_ratio <- get_pop_significant_genes(
  data = coalescence_df_with_features2, # Use the corrected data frame
  metric_var = Tpooled_Twithin_ratio,
  threshold_quantile = quantile_to_use
) %>%
  arrange(desc(Tpooled_Twithin_ratio))



#View(significant_genes_Twithin)
#View(significant_genes_ratio)
#write_csv(significant_genes_ratio, "humans/humans_high_Tratio_significant_genes.csv")


# Convert the 'significant_genes_Twithin' dataframe into a list format where names are populations and values are vectors of gene names.
upset_list_Twithin <- split(significant_genes_Twithin$genes, significant_genes_Twithin$population)
#upset(
#  fromList(upset_list_Twithin), 
#  nsets = 5, # Show the 5 populations
#  nintersects = 20, # Show the 20 most frequent intersections
#  order.by = "freq", # Order the intersections by frequency
#  text.scale = 1.5,
#  point.size = 3,
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
if (interactive()) View(wide_significant_genes_Twithin)
#write_csv(wide_significant_genes_Twithin, "humans/unfiltered_humans_wide_99.00th_genes_Twithin.csv")



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
if (interactive()) View(wide_significant_genes_ratio)
#write_csv(wide_significant_genes_ratio, "humans/humans_wide_99.00th_genes_ratio.csv")


plot_data_ranked_ratio <- coalescence_df_with_features2 %>%
  arrange(desc(Tpooled_Twithin_ratio)) %>%
  filter(Tpooled_Twithin_ratio > 10)  %>%
  select(chromosome, start, end, population, Tpooled_Twithin_ratio, genes, gene_biotypes, descriptions, feature_annotation, tr_overlap_percent, SRA_overlap_percent)
#write_csv(plot_data_ranked_ratio, "humans/ranked_ratio.csv")


plot_data_ranked_tmrca <- coalescence_df_with_features2 %>%
  arrange(desc(T_within_myr)) %>%
  filter(T_within_myr > 10) %>% 
  select(chromosome, start, end, population, T_within_myr, genes, gene_biotypes, descriptions, feature_annotation, tr_overlap_percent, SRA_overlap_percent)
#write_csv(plot_data_ranked_tmrca, "humans/ranked_tmrca.csv")

########### PLOTTTING MANHATTANS ALL TOGETHER#############################
# Process remaining data
coalescence_df_with_features2  <- coalescence_df_with_features2 %>%
  arrange(chr_clean, start) %>%
  mutate(row_index = row_number())

# Recalculate chromosome label positions
axisdf <- coalescence_df_with_features2 %>%
  group_by(chr_clean) %>%
  summarize(center = mean(row_index))


# pop_colors is defined once near the top and contains BOTH the 38 1KG+Dantu
# population codes AND the 5 superpopulation codes (Africa/America/EastAsia/
# SouthAsia/WestEurasia and AFR/AMR/EAS/SAS/EUR). The same palette therefore
# works whether the input CSV is annotated by population or by superpopulation —
# entries not present in the data are simply ignored by ggplot.
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
  threshold_quantile = 0.9995,
  title = "Pooled TMRCA 99.95th",
  y_label = "Pairwise TMRCA (Mya)",
  gene_map_df = NULL
)
pooled_plot
ggsave(file.path(.out_dir, "humans_manhattan_pooled_99.95th_filtered.pdf"), plot = pooled_plot, width = 12, height = 12)

chromosome_plot_per_pop_thresholds <- function(data, y_var, threshold_quantile, title, y_label, 
                                               gene_map_df = NULL, 
                                               pop_colors =  c(
                                                 "Africa" = "#FFCC33",
                                                 "America" = "#FF3333",
                                                 "EastAsia" = "#006600",
                                                 "SouthAsia" = "#FF33FF",
                                                 "WestEurasia" = "#66CCCC"
                                                 
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
  threshold_quantile = 0.9995,
  title = "Within TMRCA 99.95th", 
  y_label = "Pairwise TMRCA (Mya)", 
  gene_map_df = NULL, 
  pop_colors = pop_colors  
)
tryCatch(
  ggsave(file.path(.out_dir, "human_manhattan_by_chromosome_plot_vary_thresh_99.95th_filtered.pdf"),
         plot = chromosome_plot_vary_thresh, width = 6, height = 3),
  error = function(e) message(
    "Skipping per-pop manhattan ggsave (likely ggrastr 'Empty raster' from sparse populations): ",
    conditionMessage(e)
  )
)







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
        population == "Africa"   ~ "#FFCC33",
        population == "America"   ~ "#FF3333",
        population == "EastAsia"   ~ "#CCFF66", 
        population == "SouthAsia"   ~ "#FF33FF",
        population == "WestEurasia"   ~ "#66CCCC",
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
quantile_to_use <- 0.9995

# --- Plot for T_within_myr (uses the default label colors) ---
plot_avg_pop <- create_population_genome_plot(
  coalescence_df_with_features2, 
  "T_within_myr", 
  quantile_to_use,
  "Genome-wide Pairwise Coalescent Time, TMRCA (99.95th Percentile)",
  "Pairwise TMRCA (Mya)"
)

# --- Plot for Tpooled_Twithin_ratio (with inverted label colors) ---
plot_ratio_pop <- create_population_genome_plot(
  coalescence_df_with_features2, 
  "Tpooled_Twithin_ratio", 
  quantile_to_use,
  "Genome-wide Tpooled/Twithin (99.95th Percentile)",
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

#ggsave("humans/humans_genome_plot_avg_pairwise_by_pop.png", plot_avg_pop, width = 4, height = 4, dpi = 300)
#ggsave("humans/humans_genome_plot_tpooled_ratio_by_pop.png", plot_ratio_pop_log_formatted, width = 4, height = 4, dpi = 300)

#ggsave("humans/humans_genome_plot_avg_pairwise_by_pop_99.99th.pdf", plot_avg_pop, width = 4, height=4)
#ggsave("humans/humans_genome_plot_tpooled_ratio_by_pop_90.99.pdf", plot_ratio_pop_log_formatted, width = 4, height = 4)



# --- ZOOMED  ---
plot_zoomed_tmrca <- function(data,
                              target_chr,
                              target_start,
                              target_end,
                              y_var,
                              y_label,
                              plot_title,
                              label_threshold = 0,
                              upper_threshold = NULL,        # purple dashed line (e.g. 6 Mya)
                              upper_threshold_color = "#7B3F99",
                              pop_colors = NULL,
                              pop_to_superpop = NULL,
                              facet_by_superpop = FALSE,
                              superpop_subset = NULL,        # restrict to these superpops
                              facet_ncol = 2) {
  
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

  # Drop populations with no usable values for `y_var` in this region
  # (NA, non-finite, or all-zero windows are treated as "no estimate").
  pops_with_data <- region_df %>%
    group_by(population) %>%
    summarise(
      n_valid = sum(is.finite(.data[[y_var]]) & .data[[y_var]] > 0),
      .groups = "drop"
    ) %>%
    filter(n_valid > 0) %>%
    pull(population)

  dropped_pops <- setdiff(unique(as.character(region_df$population)), as.character(pops_with_data))
  if (length(dropped_pops) > 0) {
    message("plot_zoomed_tmrca: dropping populations with no valid ", y_var,
            " in region: ", paste(dropped_pops, collapse = ", "))
  }
  region_df <- region_df %>%
    filter(population %in% pops_with_data) %>%
    mutate(population = droplevels(factor(population)))

  if (nrow(region_df) == 0) {
    warning("No populations with valid values remain after filtering.")
    return(ggplot() + labs(title = plot_title, subtitle = "No data available") + theme_minimal())
  }

  # Optional: attach superpop and (if requested) prepare for facet_wrap.
  if (!is.null(pop_to_superpop)) {
    region_df <- region_df %>%
      mutate(superpop = unname(pop_to_superpop[as.character(population)]))
    unmapped <- region_df %>%
      filter(is.na(superpop)) %>%
      pull(population) %>%
      unique()
    if (length(unmapped) > 0) {
      message("plot_zoomed_tmrca: populations not in pop_to_superpop (dropped from facet): ",
              paste(unmapped, collapse = ", "))
      region_df <- region_df %>% filter(!is.na(superpop))
    }
    # Restrict to a chosen subset of superpops (e.g. Africa-only, or
    # everything except Africa) when caller requested.
    if (!is.null(superpop_subset)) {
      region_df <- region_df %>% filter(superpop %in% superpop_subset)
      if (nrow(region_df) == 0) {
        warning("No data remains after filtering to superpop_subset = ",
                paste(superpop_subset, collapse = ", "))
        return(ggplot() + labs(title = plot_title, subtitle = "No data available") + theme_minimal())
      }
    }
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
    mutate(display_label = gene)

  # When faceting by superpop, label the peak per (gene, superpop) so labels
  # appear in each facet where the gene actually peaks. Otherwise one label
  # per gene (global peak).
  if (facet_by_superpop && "superpop" %in% names(gene_labels)) {
    gene_labels <- gene_labels %>%
      group_by(gene, superpop) %>%
      slice_max(order_by = !!sym(y_var), n = 1, with_ties = FALSE) %>%
      ungroup()
  } else {
    gene_labels <- gene_labels %>%
      group_by(gene) %>%
      slice_max(order_by = !!sym(y_var), n = 1, with_ties = FALSE) %>%
      ungroup()
  }
  
  # Maximally-distinct qualitative palette, cycled fresh per facet.
  # ColorBrewer Set1 (9) extended with Set2/Set3 picks for cases with more
  # populations in a single facet (e.g. WestEurasia = 13 pops).
  distinct_pal <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
    "#FFD92F", "#A65628", "#F781BF", "#999999",
    "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
    "#A6D854", "#FFD92F", "#B3B3B3"
  )

  # ============================================================
  # Faceted-by-superpop path: one subplot per superpop, each with
  # its own legend (patchwork::wrap_plots), each using the same
  # max-contrast palette so within-facet pops are easy to tell apart.
  # ============================================================
  if (facet_by_superpop && "superpop" %in% names(region_df)) {
    sp_order <- intersect(
      c("Africa", "America", "WestEurasia", "SouthAsia", "EastAsia"),
      unique(as.character(region_df$superpop))
    )

    build_subplot <- function(sp) {
      sub_df <- region_df %>% filter(superpop == sp)
      sub_pops <- sort(unique(as.character(sub_df$population)))
      if (length(sub_pops) > length(distinct_pal)) {
        sub_pal <- setNames(scales::hue_pal()(length(sub_pops)), sub_pops)
      } else {
        sub_pal <- setNames(distinct_pal[seq_along(sub_pops)], sub_pops)
      }

      sub_genes <- if ("superpop" %in% names(gene_labels)) {
        gene_labels %>% filter(superpop == sp)
      } else {
        gene_labels
      }

      p_sub <- ggplot(sub_df, aes(x = start, y = !!sym(y_var), color = population)) +
        geom_hline(yintercept = label_threshold, linetype = "dashed",
                   color = "red", linewidth = 0.5)
      if (!is.null(upper_threshold)) {
        p_sub <- p_sub + geom_hline(yintercept = upper_threshold,
                                    linetype = "dashed",
                                    color = upper_threshold_color,
                                    linewidth = 0.5)
      }
      p_sub +
        geom_line(linewidth = 0.45, alpha = 0.85) +
        geom_text_repel(
          data = sub_genes,
          aes(label = display_label),
          nudge_y = 0.6, size = 2, segment.color = "grey50", color = "black",
          box.padding = 0.3, max.overlaps = Inf
        ) +
        scale_color_manual(name = NULL, values = sub_pal) +
        scale_x_continuous(
          labels = scales::label_number(scale = 1e-6, suffix = " Mb"),
          limits = c(target_start, target_end)
        ) +
        guides(color = guide_legend(nrow = 2, override.aes = list(linewidth = 1.4))) +
        labs(title = sp, x = NULL, y = y_label) +
        theme_classic() +
        theme(
          legend.position = "bottom",
          legend.key.size = unit(0.4, "cm"),
          legend.text = element_text(size = 7),
          legend.margin = margin(t = 0, b = 2),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 9)
        )
    }

    subplots <- lapply(sp_order, build_subplot)

    # For a single-superpop case let the panel be plot-wide; otherwise
    # honour the requested ncol.
    eff_ncol <- if (length(subplots) == 1) 1 else facet_ncol
    combined <- patchwork::wrap_plots(subplots, ncol = eff_ncol) +
      patchwork::plot_annotation(
        title = plot_title,
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
      )
    return(combined)
  }

  # ============================================================
  # Non-faceted (single-plot) path - uses supplied pop_colors palette
  # ============================================================
  pops_in_data <- sort(unique(as.character(region_df$population)))
  if (is.null(pop_colors) || !any(pops_in_data %in% names(pop_colors))) {
    color_scale <- scale_color_discrete(name = "Population")
  } else {
    missing_pops <- setdiff(pops_in_data, names(pop_colors))
    if (length(missing_pops) > 0) {
      extra <- setNames(scales::hue_pal()(length(missing_pops)), missing_pops)
      pop_colors <- c(pop_colors, extra)
    }
    color_scale <- scale_color_manual(name = "Population", values = pop_colors)
  }

  p <- ggplot(region_df, aes(x = start, y = !!sym(y_var), color = population)) +
    geom_hline(yintercept = label_threshold, linetype = "dashed", color = "red", linewidth = 0.8)
  if (!is.null(upper_threshold)) {
    p <- p + geom_hline(yintercept = upper_threshold, linetype = "dashed",
                        color = upper_threshold_color, linewidth = 0.8)
  }
  p <- p +
    geom_line(linewidth = 0.4, alpha = 0.6) +
    geom_text_repel(
      data = gene_labels,
      aes(label = display_label),
      nudge_y = 1.5, size = 2, segment.color = "grey50", color = "black",
      box.padding = 0.5, max.overlaps = Inf
    ) +
    color_scale +
    guides(color = guide_legend(ncol = 6, override.aes = list(linewidth = 1.2))) +
    scale_x_continuous(
      labels = scales::label_number(scale = 1e-6, suffix = " Mb"),
      limits = c(target_start, target_end)
    ) +
    labs(title = plot_title, x = paste("Genomic Position ", target_chr), y = y_label) +
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

############ GYP ########################
my_chr <- "chr4"
my_start <- 147000000
my_end <- 147700000
my_label_threshold <- 1.18

.gyp_rect <- annotate(
  "rect",
  xmin = 147188085,
  xmax = 147458090,
  ymin = -Inf,
  ymax = Inf,
  alpha = 0.2,
  fill = "lightgrey"
)

# `&` broadcasts a layer to every subplot when the object is a patchwork;
# fall back to `+` for the non-faceted ggplot case.
add_gyp_rect <- function(p) {
  if (inherits(p, "patchwork")) p & .gyp_rect else p + .gyp_rect
}

# ---- (1) Africa-only standalone (for main figure) ----
gyp_plot_africa <- plot_zoomed_tmrca(
  data = coalescence_df_with_features,
  target_chr = my_chr,
  target_start = my_start,
  target_end = my_end,
  y_var = "T_within_myr",
  y_label = "Pairwise TMRCA (Mya)",
  plot_title = "Glycophorin region in humans (Africa)",
  label_threshold = my_label_threshold,
  upper_threshold = 6.0,
  pop_colors = pop_colors,
  pop_to_superpop = pop_to_superpop,
  facet_by_superpop = TRUE,
  superpop_subset = "Africa"
)
GYP_africa_annotated <- add_gyp_rect(gyp_plot_africa)
.gyp_out_africa <- file.path(.out_dir, "GYP_humans_Africa.pdf")
ggsave(.gyp_out_africa, GYP_africa_annotated, width = 5, height = 4)
message("Saved: ", .gyp_out_africa)

# ---- (2) Non-Africa continents (for supplementary material) ----
gyp_plot_other <- plot_zoomed_tmrca(
  data = coalescence_df_with_features,
  target_chr = my_chr,
  target_start = my_start,
  target_end = my_end,
  y_var = "T_within_myr",
  y_label = "Pairwise TMRCA (Mya)",
  plot_title = "Glycophorin region in humans",
  label_threshold = my_label_threshold,
  upper_threshold = 6.0,
  pop_colors = pop_colors,
  pop_to_superpop = pop_to_superpop,
  facet_by_superpop = TRUE,
  superpop_subset = c("America", "WestEurasia", "SouthAsia", "EastAsia")
)
GYP_other_annotated <- add_gyp_rect(gyp_plot_other)
print(GYP_other_annotated)
.gyp_out_other <- file.path(.out_dir, "GYP_humans.pdf")
ggsave(.gyp_out_other, GYP_other_annotated, width = 14, height = 8)
message("Saved: ", .gyp_out_other)

#tmrca_ratio_plot <- plot_zoomed_tmrca(
#  data = coalescence_df_with_features,
#  target_chr = my_chr,
#  target_start = my_start,
#  target_end = my_end,
#  y_var = "Tpooled_Twithin_ratio", # Changed to the ratio variable
#  y_label = "Tpooled/Twithin",   # Updated y-axis label
#  plot_title = "GYP Pregion", # Updated title
#  label_threshold = my_label_threshold
#)
#GYP<-tmrca_within_plot + tmrca_ratio_plot +
#  plot_layout(guides = 'collect') & 
#  theme(legend.position = 'bottom')
#GYP





