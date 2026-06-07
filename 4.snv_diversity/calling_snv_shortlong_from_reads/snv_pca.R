setwd("/Users/joanocha/Library/CloudStorage/GoogleDrive-joana.laranjeira.rocha@gmail.com/My Drive/POSTDOC/PANPAN/analysis/Figure1_Diversity_Assembly_stats/PANPAN_snvs")
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyr)
library(tibble)
library(cowplot)
library(ggplot2)
library(dplyr)
library(data.table)
library(patchwork) # For side-by-side twin figures
library(tibble) # For deframe()


metadata <- fread("PANPAN_PradoMartine_deManuel_METADATA.tsv")

meta_dict <- bind_rows(
  metadata %>% select(IID = id, Population = CommonName, Color = colors, type3),
  metadata %>% select(IID = sample, Population = CommonName, Color = colors, type3)
) %>%
  distinct(IID, .keep_all = TRUE) %>%
  mutate(
    type3 = trimws(type3),
    Data_Type = case_when(
      grepl("hap1", IID) ~ "Long-reads haplotype 1",
      grepl("hap2", IID) ~ "Long-reads haplotype 2",
      type3 == "short-read" ~ "Short-reads",
      TRUE ~ "Long-reads"
    )
  )
color_map <- meta_dict %>%
  select(Population, Color) %>%
  distinct() %>%
  deframe()


# --- 2. DEFINE SHAPES ---
shape_map <- c(
  "Short-reads" = 1,                 # Open circle
  "Long-reads" = 16,                 # Filled circle
  "Long-reads haplotype 1" = 16,     
  "Long-reads haplotype 2" = 17      
)


# --- 3. THE SIMPLIFIED HELPER FUNCTION ---
load_pca_data <- function(vec_file, val_file, dataset_label, meta_dictionary) {
  if (!file.exists(vec_file) | !file.exists(val_file)) {
    warning(paste("File missing:", vec_file))
    return(NULL)
  }
  vec <- fread(vec_file)
  val <- fread(val_file)
  names(vec)[1:4] <- c("FID", "IID", "PC1", "PC2")
  vars <- val$V1 / sum(val$V1) * 100
  vec <- vec %>%
    left_join(meta_dictionary, by = "IID") %>%
    mutate(
      Dataset = dataset_label,
      PC1_var = round(vars[1], 1),
      PC2_var = round(vars[2], 1)
    )
  
  return(vec)
}


# --- 4. PLOTTING FUNCTION (Using theme_classic) ---

create_pca_plot <- function(data, title) {
  if(is.null(data)) return(NULL)
  
  pc1_label <- paste0("PC1 (", unique(data$PC1_var), "%)")
  pc2_label <- paste0("PC2 (", unique(data$PC2_var), "%)")
  
  ggplot(data, aes(x = PC1, y = PC2, color = Population, shape = Data_Type)) +
    geom_point(size = 1.5) +
    scale_color_manual(values = color_map) +
    scale_shape_manual(values = shape_map, name = "Data Type") +
    labs(title = title, x = pc1_label, y = pc2_label) +
    
    # 1. Use theme_bw() which natively draws the outer box
    theme_bw() +  
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      axis.text = element_text(color = "black", size = 10),
      
      # 2. Remove the internal grid lines so the data pops
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # 3. Enforce a strong, black "closed box" border
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
}
# --- 4. LOAD THE DATA ---
# Bonobos
bonobo_dip <- load_pca_data("pca/vcfs_from_reads_shortlong-reads/panpaniscus_reads_mapped2mPanPan1_BIALLELIC_SNPS.eigenvec", 
                            "pca/vcfs_from_reads_shortlong-reads/panpaniscus_reads_mapped2mPanPan1_BIALLELIC_SNPS.eigenval", 
                            "Bonobo Diploid", meta_dict)
bonobo_hap <- load_pca_data("pca/long-reads/panpa_mapped2mPanPan1.eigenvec", 
                            "pca/long-reads/panpa_mapped2mPanPan1.eigenval", 
                            "Bonobo Haploid", meta_dict)

# Chimpanzees
chimp_dip  <- load_pca_data("pca/vcfs_from_reads_shortlong-reads/pantros_reads_mapped2mPanTro3_BIALLELIC_SNPS.eigenvec", 
                            "pca/vcfs_from_reads_shortlong-reads/pantros_reads_mapped2mPanTro3_BIALLELIC_SNPS.eigenval", 
                            "Chimp Diploid", meta_dict)

chimp_hap  <- load_pca_data("pca/long-reads/pantros_mapped2mPanTro3.eigenvec", 
                            "pca/long-reads/pantros_mapped2mPanTro3.eigenval", 
                            "Chimp Haploid", meta_dict)


# --- TWIN FIGURES ---
# 1. Bonobos
p_bonobo_dip <- create_pca_plot(bonobo_dip, "Bonobo Diploid (SNVs)")
p_bonobo_hap <- create_pca_plot(bonobo_hap, "Bonobo Haploid (SNVs)")

twin_bonobo <-  p_bonobo_hap + p_bonobo_dip + plot_layout(guides = "collect")
print(twin_bonobo)
ggsave('pca_bonobo.pdf', twin_bonobo, width=8, height=2.5)
# 2. Chimpanzees
p_chimp_dip <- create_pca_plot(chimp_dip, "Chimpanzee Diploid (SNVs)")
p_chimp_hap <- create_pca_plot(chimp_hap, "Chimpanzee Haploid (SNVs)")
twin_chimp <-  p_chimp_hap + p_chimp_dip + plot_layout(guides = "collect")
print(twin_chimp)
