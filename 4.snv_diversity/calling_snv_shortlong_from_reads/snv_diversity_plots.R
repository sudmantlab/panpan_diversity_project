suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(tidyr)
  library(tibble)
  library(patchwork)
})

dir.create("figures", showWarnings = FALSE)


#SNV DIVERSITY 

afr      <- fread("diversity_stats/combined_pi_stats_AFR.tsv")
non_afr  <- fread("diversity_stats/combined_pi_stats_nonAFR.tsv")
pantros  <- fread("diversity_stats/pantros_combined_pi_stats.tsv")
ppan     <- fread("diversity_stats/ppan_combined_pi_stats.tsv")
metadata <- fread("PANPAN_PradoMartine_deManuel_METADATA.tsv")

callable_human_g <- 2875101262
avg_afr     <- mean(afr$mean_pi * afr$n_sites / callable_human_g)
avg_non_afr <- mean(non_afr$mean_pi * non_afr$n_sites / callable_human_g)

apes_raw <- bind_rows(pantros, ppan)
meta_unique <- metadata %>% distinct(sample, .keep_all = TRUE)
apes_processed <- apes_raw %>% left_join(meta_unique, by = "sample") %>%
  mutate(type3 = trimws(type3))

pop_order <- c("Bonobo", "Central", "Eastern", "Nigeria-Cameroon", "Western", "Western x Central")
apes_processed$CommonName <- factor(apes_processed$CommonName, levels = pop_order)

color_map <- apes_processed %>%
  select(CommonName, colors) %>%
  distinct() %>%
  filter(!is.na(CommonName)) %>%
  deframe()

# ---- individual barplot ----
plot1_data <- apes_processed %>%
  filter(!is.na(CommonName)) %>%
  arrange(CommonName, mean_pi) %>%
  mutate(sample_factor = factor(sample, levels = unique(sample)))

plot1 <- ggplot(plot1_data, aes(x = sample_factor, y = mean_pi, fill = CommonName, alpha = type3)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = avg_afr,     color = "black", linetype = "dotted", size = 0.7) +
  geom_hline(yintercept = avg_non_afr, color = "black", linetype = "dotted", size = 0.7) +
  annotate("text", x = 5, y = avg_afr,     label = "Human AFR", vjust = -0.5, color = "black", fontface="bold") +
  annotate("text", x = 5, y = avg_non_afr, label = "non-AFR",   vjust =  1.5, color = "black", fontface="bold") +
  scale_alpha_manual(values = c("long-read" = 1, "short-read" = 0.5), name = "Read Type") +
  scale_fill_manual(values = color_map) +
  labs(x = "Individuals", y = "SNV heterozygosity",
       title = "Fig 1E - Individual Diversity (Ordered by Population & Pi)") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave("figures/Fig1E_Pi_Bars_Evolutionary_Order.pdf", plot1, width = 12, height = 5)
ggsave("figures/Fig1E_Pi_Bars_Evolutionary_Order.png", plot1, width = 12, height = 5, dpi = 200)

# population boxplot 
plot2_data <- apes_processed %>%
  filter(!is.na(CommonName), CommonName != "Western x Central")

plot2 <- ggplot(plot2_data, aes(x = CommonName, y = mean_pi, fill = CommonName)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1.5) +
  scale_fill_manual(values = color_map) +
  labs(x = "", y = "SNV heterozygosity", title = "ED Fig 1G - Population SNV heterozygosity") +
  theme_classic() +
  coord_flip() +
  theme(legend.position = "none",
        axis.text.y = element_text(angle = 0, hjust = 1))

ggsave("figures/EDFig1G_Pi_Population_Boxplots.pdf", plot2, width = 7, height = 5)
ggsave("figures/EDFig1G_Pi_Population_Boxplots.png", plot2, width = 7, height = 5, dpi = 200)

# faceted by pop, short vs long read
plot3_data <- apes_processed %>%
  filter(!is.na(CommonName), !CommonName %in% c("Nigeria-Cameroon", "Western x Central"))

plot3 <- ggplot(plot3_data, aes(x = type3, y = mean_pi)) +
  geom_boxplot(fill = "white", color = "gray30", alpha = 0.8, outlier.shape = NA) +
  geom_jitter(aes(color = CommonName), width = 0.15, alpha = 0.8, size = 1.5) +
  coord_flip() +
  facet_wrap(~CommonName, ncol = 1, strip.position = "left") +
  scale_color_manual(values = color_map, name = "Population") +
  labs(x = "", y = "SNV heterozygosity", title = "ED Fig 1H - SNV heterozygosity by read type") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9),
    legend.position = "right",
    strip.placement = "outside",
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text.y.left = element_text(angle = 90, size = 9, face = "plain")
  )

ggsave("figures/EDFig1H_Pi_Faceted_ReadType.pdf", plot3, width = 10, height = 5)
ggsave("figures/EDFig1H_Pi_Faceted_ReadType.png", plot3, width = 10, height = 5, dpi = 200)


# PCA 

meta_dict <- bind_rows(
  metadata %>% select(IID = id,     Population = CommonName, Color = colors, type3),
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

color_map_pca <- meta_dict %>%
  select(Population, Color) %>%
  filter(!is.na(Population)) %>%
  distinct() %>%
  deframe()

shape_map <- c(
  "Short-reads" = 1,
  "Long-reads"  = 16,
  "Long-reads haplotype 1" = 16,
  "Long-reads haplotype 2" = 17
)

load_pca_data <- function(vec_file, val_file, dataset_label, meta_dictionary) {
  if (!file.exists(vec_file) || !file.exists(val_file)) {
    warning(paste("File missing:", vec_file))
    return(NULL)
  }
  vec <- fread(vec_file)
  val <- fread(val_file)
  names(vec)[1:4] <- c("FID", "IID", "PC1", "PC2")
  vars <- val$V1 / sum(val$V1) * 100
  vec %>%
    left_join(meta_dictionary, by = "IID") %>%
    mutate(Dataset = dataset_label,
           PC1_var = round(vars[1], 1),
           PC2_var = round(vars[2], 1))
}

create_pca_plot <- function(data, title) {
  if (is.null(data)) return(NULL)
  pc1_label <- paste0("PC1 (", unique(data$PC1_var), "%)")
  pc2_label <- paste0("PC2 (", unique(data$PC2_var), "%)")
  ggplot(data, aes(x = PC1, y = PC2, color = Population, shape = Data_Type)) +
    geom_point(size = 1.8) +
    scale_color_manual(values = color_map_pca) +
    scale_shape_manual(values = shape_map, name = "Data Type") +
    labs(title = title, x = pc1_label, y = pc2_label) +
    theme_bw() +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      axis.text = element_text(color = "black", size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
}

bonobo_dip <- load_pca_data(
  "pca/vcfs_from_reads_shortlong-reads/panpaniscus_reads_mapped2mPanPan1_BIALLELIC_SNPS.eigenvec",
  "pca/vcfs_from_reads_shortlong-reads/panpaniscus_reads_mapped2mPanPan1_BIALLELIC_SNPS.eigenval",
  "Bonobo Diploid", meta_dict)
bonobo_hap <- load_pca_data(
  "pca/long-reads/panpa_mapped2mPanPan1.eigenvec",
  "pca/long-reads/panpa_mapped2mPanPan1.eigenval",
  "Bonobo Haploid", meta_dict)
chimp_dip  <- load_pca_data(
  "pca/vcfs_from_reads_shortlong-reads/pantros_reads_mapped2mPanTro3_BIALLELIC_SNPS.eigenvec",
  "pca/vcfs_from_reads_shortlong-reads/pantros_reads_mapped2mPanTro3_BIALLELIC_SNPS.eigenval",
  "Chimpanzee Diploid", meta_dict)

p_bonobo_hap <- create_pca_plot(bonobo_hap, "Bonobo Haploid (SNVs)")
p_bonobo_dip <- create_pca_plot(bonobo_dip, "Bonobo Diploid (SNVs)")
p_chimp_dip  <- create_pca_plot(chimp_dip,  "Chimpanzee Diploid (SNVs)")

twin_bonobo <- p_bonobo_hap + p_bonobo_dip + plot_layout(guides = "collect")
ggsave("figures/Fig1F_EDFig1I_pca_bonobo.pdf", twin_bonobo, width = 10, height = 3.5)
ggsave("figures/Fig1F_EDFig1I_pca_bonobo.png", twin_bonobo, width = 10, height = 3.5, dpi = 200)

if (!is.null(p_chimp_dip)) {
  ggsave("figures/Fig1F_pca_chimp_diploid.pdf", p_chimp_dip, width = 6, height = 3.5)
  ggsave("figures/Fig1F_pca_chimp_diploid.png", p_chimp_dip, width = 6, height = 3.5, dpi = 200)
}

cat("\nDONE. Outputs written to figures/:\n")
print(list.files("figures"))
