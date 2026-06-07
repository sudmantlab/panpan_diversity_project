setwd("/Users/joanocha/Desktop/SINGER")
library(ggplot2)
library(scales)
library(readr)
library(ggplot2)
library(dplyr)
library(ggrepel) 
library(tidyr)
library(ggplot2)
library(scales)
library(RColorBrewer)  # For better color palettes
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork) 

# Read both datasets
bonobo <- read.csv("bonobo_significant_99.95_pwise_coalescence_genes.csv")
chimp <- read.csv("chim_significant_99.95_pwise_coalescence_genes.csv")


common_genes <- intersect(bonobo$genes, chim$genes)
common_data <- bonobo %>%
  filter(genes %in% common_genes) %>%
  inner_join(chim %>% filter(genes %in% common_genes),
             by = "genes",
             suffix = c("_bonobo", "_chim"))
write.csv(common_data, "common_significant_genes_bonobo_chim.csv", row.names = FALSE)
cat("Common genes found:", paste(common_genes, collapse = ", "), "\n")
cat("Total common genes:", length(common_genes), "\n")
cat("Merged data saved to 'common_significant_genes_bonobo_chim.csv'")


# Find unique genes
bonobo_genes <- unique(bonobo$genes)
chim_genes <- unique(chim$genes)
common_genes <- intersect(bonobo_genes, chim_genes)
# Bonobo-specific genes
bonobo_only <- bonobo %>% 
  filter(!genes %in% common_genes) %>%
  distinct(genes, .keep_all = TRUE)  # Keep unique entries
# Chimp-specific genes
chim_only <- chim %>% 
  filter(!genes %in% common_genes) %>%
  distinct(genes, .keep_all = TRUE)
write.csv(bonobo_only, "bonobo_specific_genes.csv", row.names = FALSE)
write.csv(chim_only, "chim_specific_genes.csv", row.names = FALSE)
cat("Common genes:", length(common_genes), "\n")
cat("Bonobo-specific genes:", nrow(bonobo_only), "\n")
cat("Chimp-specific genes:", nrow(chim_only), "\n")



### plot ranked genes
all_data <- bind_rows(
  # Bonobo-specific genes
  bonobo %>%
    filter(!genes %in% common_genes) %>%
    mutate(group = "Bonobo-specific"),
  
  # Chimp-specific genes
  chim %>%
    filter(!genes %in% common_genes) %>%
    mutate(group = "Chimp-specific"),
  
  # Common genes - use highest value per gene
  bonobo %>%
    filter(genes %in% common_genes) %>%
    inner_join(chim %>% filter(genes %in% common_genes),
               by = "genes", suffix = c("_bonobo", "_chim")) %>%
    mutate(avg_pairwise_myr = pmax(avg_pairwise_myr_bonobo, avg_pairwise_myr_chim)) %>%
    select(genes, avg_pairwise_myr, chromosome = chromosome_bonobo, start = start_bonobo) %>%
    mutate(group = "Common")
) %>%
  # Create global ranking
  arrange(desc(avg_pairwise_myr)) %>%
  mutate(global_rank = row_number())

# Create color palette
group_colors <- c(
  "Common" = "#56B4E9",        # Blue
  "Bonobo-specific" = "orange", # 
  "Chimp-specific" = "#009E73"  
)

# Create integrated plot
integrated_plot <- ggplot(all_data, aes(x = global_rank, y = avg_pairwise_myr, color = group)) +
  geom_point(alpha = 0.8, size = 1) +
  geom_text_repel(
    aes(label = ifelse(global_rank <= 100, genes, "")),  # Label top 100 genes
    size = 2.5,
    max.overlaps = 80,
    show.legend = FALSE,
    segment.color = "grey50"
  ) +
  scale_color_manual(values = group_colors) +
  labs(
    x = "Global Rank (Highest to Lowest TMRCA) for top genes > 99.95%, only top 100 are labeled",
    y = "Average Pairwise Coalescence Time (Myr)",
    title = "Ranking of Average Pairwise Coalescence Time (Myr) for windows with overlapping genes with ",
    color = "",
    caption = paste("Total genes > 99.95%:", nrow(all_data), 
                    "| Common genes:", length(common_genes),
                    "| Bonobo-specific:", sum(all_data$group == "Bonobo-specific"),
                    "| Chimp-specific:", sum(all_data$group == "Chimp-specific"))
  ) +
  theme_classic() +
  theme(
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    axis.text = element_text(size = 10)
  ) +
  scale_x_continuous(expand = expansion(mult = 0.02)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
  guides(color = guide_legend(nrow = 1))



# Display plot
print(integrated_plot)
ggsave("integrated_gene_ranking.png", integrated_plot, width = 14, height = 12, dpi = 300)