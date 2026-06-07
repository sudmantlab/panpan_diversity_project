setwd("/Users/joanocha/Google Drive/My Drive/POSTDOC/PANPAN/analysis/Figure4_SINGER/")

library(dplyr)
library(ggplot2)
library(scales)
library(readr)
library(ggrepel)
library(tidyr)
library(patchwork)
library(valr)
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


library(RColorBrewer)


conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("mutate", "dplyr")
conflicts_prefer(dplyr::first)

# --- 1. DEFINE POPULATION MAPPING (Continent / Superpopulation) ---
# Add or modify populations here as needed
pop_map <- tibble(
  population = c(
    # AFR
    "Mbuti", "Biaka", "Mandenka", "Yoruba", "San", "BantuSouthAfrica", 
    "BantuKenya", "JuhoanNorth", "BantuHerero", "BantuTswana", "Esan", 
    "Gambian", "Mende", "Luhya", "ACB", "ASW", "MSL", "GWD", "YRI", "ESN", "LWK",
    # AMR
    "Colombian", "Mayan", "MXL", "PEL", "PUR", "CLM", "Karitiana", "Surui", "Pima", "Piapoco",
    # EAS
    "Dai", "Han", "Japanese", "Kinh", "CHS", "Daur", "Hezhen", "Lahu", "Miao", 
    "Mongola", "Naxi", "Oroqien", "She", "Tu", "Tujia", "Xibo", "Yi", "JPT", "CDX", "KHV", "CHB",
    "NorthernHan", "Oroqen", "Yakut", "Mongolian", "Cambodian",
    # EUR
    "Finnish", "Tuscan", "CEU", "GBR", "IBS", "Orcadian", "French", "Basque", 
    "Italian", "Sardinian", "Russian", "Adygei", "TSI", "FIN", "BergamoItalian",
    # SAS
    "Bengali", "Punjabi", "Gujarati", "ITU", "STU", "Brahui", "Balochi", 
    "Hazara", "Makrani", "Sindhi", "Pathan", "Burusho", "Kalash", "BEB", "PJL", "GIH", "Uygur",
    # MID
    "Bedouin", "Druze", "Mozabite", "Palestinian",
    # OCE
    "Papuan", "Bougainville"
  ),
  continent = c(
    rep("AFR", 21),
    rep("AMR", 10),
    rep("EAS", 26),
    rep("EUR", 15),
    rep("SAS", 17),
    rep("MID", 4),
    rep("OCE", 2)
  )
)

# --- 2. LOAD & PROCESS DATA ---

# Load main metrics
coalescence_df <- fread("humans/chr16-wide_metrics_annotatedbyPop.csv.gz") %>%
  rename(chrom = chromosome) %>%
  filter(chrom == "chr16")

coalescence_df <- coalescence_df %>%
  separate_rows(genes, sep = ",") %>%
  left_join(pop_map, by = "population")


#coalescence_df<- coalescence_df %>% 
  #filter(continent == "AFR") %>% 
 # group_by(chrom, start, end, population, continent) %>% 
  #summarise(
   # avg_tmrca = first(avg_tmrca),
    #T_within = first(T_within),
    #Tajimas_D = first(Tajimas_D),
    #genes = paste(unique(genes), collapse = ","), # Keep genes if needed later
    #.groups = "drop"
  #)

# Convert to Millions of Years (Myr)
coalescence_df$avg_tmrca_myr <- coalescence_df$avg_tmrca * 28 / 1e6
coalescence_df$T_within_myr <- coalescence_df$T_within * 28 / 1e6

# --- 3. SINGLE CONTINENT PLOTTING FUNCTION ---

plot_continent_specific <- function(data, target_continent, target_chr, target_start, target_end, 
                                    y_var, y_label, plot_title, label_threshold = 0) {
  
  # A. Filter for the specific continent and region
  region_df <- data %>%
    filter(
      continent == target_continent,
      chrom == target_chr, 
      start <= target_end, 
      end >= target_start
    )
  
  if (nrow(region_df) == 0) {
    message(paste("No data found for continent:", target_continent))
    return(NULL)
  }
  
  # B. Generate Dynamic Color Palette for this Continent's Populations
  pops_in_cont <- unique(region_df$population)
  n_pops <- length(pops_in_cont)
  
  # Use RColorBrewer to make enough distinct colors
  # If fewer than 3 pops, manual colors; otherwise dynamic palette
  if(n_pops < 3) {
    cont_palette <- brewer.pal(3, "Set1")[1:n_pops]
  } else {
    cont_palette <- colorRampPalette(brewer.pal(min(n_pops, 9), "Set1"))(n_pops)
  }
  names(cont_palette) <- sort(pops_in_cont)
  
  
  # C. Prepare Gene Labels (Top peak per gene)
  gene_labels <- region_df %>%
    filter(!!sym(y_var) >= label_threshold, genes != "." & !is.na(genes)) %>%
    separate_rows(genes, sep = "[,;]") %>%
    mutate(genes = trimws(genes)) %>%
    filter(genes != "", genes != "NA") %>% # ADDED check for literal "NA" string
    # Label the highest point for the gene across the whole continent's data
    group_by(genes) %>%
    slice_max(order_by = !!sym(y_var), n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # D. Create Plot
  p <- ggplot(region_df, aes(x = start, y = !!sym(y_var))) +
    
    # Threshold Line
    geom_hline(yintercept = label_threshold, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    
    # Lines: Colored by Population
    geom_line(aes(color = population, group = population), linewidth = 0.8, alpha = 0.8) +
    
    # Gene Labels
    geom_text_repel(
      data = gene_labels,
      aes(label = genes),
      size = 3, nudge_y = 0.5, min.segment.length = 0, color = "black",
      box.padding = 0.5, max.overlaps = Inf
    ) +
    
    # Scales
    scale_x_continuous(
      labels = scales::label_number(scale = 1e-6, suffix = " Mb"),
      limits = c(target_start, target_end)
    ) +
    
    # Apply specific palette
    scale_color_manual(values = cont_palette) +
    
    labs(
      title = paste0(plot_title, " (", target_continent, ")"),
      x = paste("Genomic Position on", target_chr),
      y = y_label,
      color = "Population"
    ) +
    
    theme_classic() +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 9),
      legend.title = element_text(face="bold"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11, face="bold"),
      plot.title = element_text(hjust = 0.5, size=12, face="bold")
    )
  
  return(p)
}


# --- EXECUTE LOOP FOR ALL CONTINENTS (GYP Region) ---

# Get list of unique continents (including "Unknown")
target_continents <- unique(coalescence_df$continent)

# Config for GYP
#gyp_chr <- "chr4"
#gyp_start <- 143666000
#gyp_end <- 144200000


# Config for HBA1,HBA2
gyp_chr <- "chr16"
gyp_start <- 100000 
gyp_end <- 200000

print(paste("Generating plots for continents:", paste(target_continents, collapse=", ")))

for (cont in target_continents) {
  
  # 1. T_within Plot
  p_twithin <- plot_continent_specific(
    data = coalescence_df,
    target_continent = cont,
    target_chr = gyp_chr,
    target_start = gyp_start,
    target_end = gyp_end,
    y_var = "T_within_myr",
    y_label = "Pairwise TMRCA (Mya)",
    plot_title = "HBA Region: Average Pairwise TMRCA",
    label_threshold = 2
  )
  
  if (!is.null(p_twithin)) {
    ggsave(paste0("humans/HBA_Twithin_", cont, ".pdf"), p_twithin, width = 8, height = 5)
  }
  
  # 2. Avg TMRCA Plot
  p_avg <- plot_continent_specific(
    data = coalescence_df,
    target_continent = cont,
    target_chr = gyp_chr,
    target_start = gyp_start,
    target_end = gyp_end,
    y_var = "avg_tmrca_myr",
    y_label = "TMRCA (Mya)",
    plot_title = "HBA Region: Root TMRCA",
    label_threshold = 2
  )
  
  if (!is.null(p_avg)) {
    ggsave(paste0("humans/HBA_TMRCA_", cont, ".pdf"), p_avg, width = 8, height = 5)
  }
  
  # 3. Tajima's D Plot
  p_tajima <- plot_continent_specific(
    data = coalescence_df,
    target_continent = cont,
    target_chr = gyp_chr,
    target_start = gyp_start,
    target_end = gyp_end,
    y_var = "Tajimas_D",
    y_label = "Tajima's D",
    plot_title = "HBA Region: Tajima's D",
    label_threshold = 0 # Highlight strongly negative values
  )
  
  if (!is.null(p_tajima)) {
    ggsave(paste0("humans/HBA_TajimasD_", cont, ".pdf"), p_tajima, width = 8, height = 5)
  }
}