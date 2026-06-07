setwd("/Users/joanocha/Desktop/SINGER")
library(ggplot2)
library(scales)
library(readr)
library(ggplot2)
library(dplyr)
library(ggrepel) 
library(tidyr)

coalescence_df<- read_csv("chimps_genome-wide_metrics_annotated.csv", quote = "\"")
coalescence_df$avg_tmrca_myr <- coalescence_df$avg_tmrca * 25 / 1e6
coalescence_df$avg_pairwise_myr <- coalescence_df$avg_pairwise_coalescence_time * 25 / 1e6
coalescence_df <- coalescence_df %>%
  mutate(
    chr_clean = sub("_.*", "", chromosome),  # Remove everything after first underscore
    chr_clean = factor(chr_clean, levels = c(paste0("chr", 1:23), "chrX", "chrY"))
  ) %>%
  # Calculate window midpoint
  mutate(midpoint = (start + end)/2) %>%
  # Sort by chromosome order and position
  arrange(chr_clean, start) %>%
  # Add sequential index column
  mutate(row_index = row_number())
# Calculate chromosome label positions
axisdf <- coalescence_df %>%
  group_by(chr_clean) %>%
  summarize(center = mean(row_index))  # Use row index for positioning

coalescence_df <- coalescence_df %>%
  arrange(chr_clean, midpoint) %>%
  mutate(row_index = row_number())

threshold <- quantile(coalescence_df$avg_pairwise_myr, 
                      probs = 0.9995, 
                      na.rm = TRUE)
print(paste("99.95% threshold:", threshold))

# Preprocess significant points: 
# 1. Filter significant windows
# 2. Split multi-gene entries into separate rows
# 3. Keep only the highest coalescence time entry per gene
significant_genes <- coalescence_df %>%
  filter(avg_pairwise_myr > threshold, genes != ".") %>%
  # Split genes separated by ";" into individual rows
  separate_rows(genes, sep = ";") %>%
  # Remove any remaining empty gene entries
  filter(genes != "", genes != ".") %>%
  # For each gene, keep only the window with highest coalescence time
  group_by(genes) %>%
  slice_max(avg_pairwise_myr, n = 1, with_ties = FALSE) %>%
  ungroup()

# Create plot with gene annotations
test <- ggplot(coalescence_df, aes(x = row_index, y = avg_pairwise_myr)) +
  geom_point(
    aes(color = factor(as.numeric(chr_clean) %% 2)),
    alpha = 0.7,
    size = 0.8
  ) +
  geom_hline(
    yintercept = threshold,
    color = "black",
    linetype = "dashed",
    linewidth = 0.7
  ) +
  # Annotate only selected genes
  geom_text_repel(
    data = significant_genes,
    aes(label = genes),
    size = 2.5,
    color = "black",
    max.overlaps = 60,
    min.segment.length = 0.2,
    box.padding = 0.3
  ) +
  scale_color_manual(values = c("lightgrey", "darkgreen")) +
  scale_x_continuous(
    name = "Chromosome",
    breaks = axisdf$center,
    labels = gsub("chr", "", axisdf$chr_clean),
    expand = expansion(mult = 0.01)
  ) +
  labs(
    y = "Average Pairwise Coalescence Time (Myr)",
    title = "Chimpanzees (annotated genes for windows > 99.95% of the genome-wide empirical distribution)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(hjust = 1, vjust = 1)
  )
test

significant_genes <-significant_genes %>%
  select(chromosome, chr_clean, start, end, genes, avg_pairwise_myr) %>%
  arrange(desc(avg_pairwise_myr))  # Sort by most extreme values
head(significant_genes, 30)
write_csv(significant_genes, "chim_significant_99.95_pwise_coalescence_genes.csv")
ggsave("all_Chimpanzees_avg_pairwise_tmrca_gwide.png", test, dpi = 300, width = 10, height= 3.5)


###### RANKED GENES PLOT
significant_genes <- coalescence_df %>%
  filter(avg_pairwise_myr > threshold, genes != ".") %>%
  separate_rows(genes, sep = ";") %>%
  filter(genes != "", genes != ".") %>%
  group_by(genes) %>%
  slice_max(avg_pairwise_myr, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(avg_pairwise_myr)) %>%
  mutate(
    rank = row_number(),  # Create ranking
    genomic_pos = paste0(chr_clean, ":", round(start/1e6, 2), "Mb")  # Clean position
  )
test3<-ggplot(significant_genes %>% filter(avg_pairwise_myr > 10),  # Filter main data too
              aes(x = reorder(genomic_pos, -avg_pairwise_myr), 
                  y = avg_pairwise_myr)
              ) +
  geom_point(color = "darkgreen", size = 1) +
  geom_text_repel(
    aes(label = genes),
    size = 2,
    box.padding = 0.5,
    max.overlaps = Inf
  )  +
  scale_x_discrete(name = "Genomic Position (Chromosome:Mb)") +
  scale_y_continuous(name = "Average Pairwise Coalescence Time (Myr)") +
  labs(title = "Chimpanzees Top Genes ranked by Highest TMRCA",
       subtitle = paste("Threshold =", round(threshold, 2), "Myr (99.95% percentile)")) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank()
  )
ggsave("Chimpanzees_dcreasing_avgpawise_tmrca.png", test3, dpi=300, width=6, height=10)



#### CHROMOSOME-WIDE 
test2 <- ggplot(coalescence_df, aes(x = midpoint, y = avg_pairwise_myr)) +
  # Points with alternating chromosome colors
  geom_point(
    aes(color = factor(as.numeric(chr_clean) %% 2)),
    alpha = 0.7,
    size = 0.8
  ) +
  # Genome-wide threshold line
  geom_hline(
    yintercept = threshold,
    color = "black",
    linetype = "dashed",
    linewidth = 0.7
  ) +
  # Add gene annotations if they exist
  {
    if (nrow(significant_genes) > 0) {
      geom_text_repel(
        data = significant_genes,
        aes(label = genes),
        size = 2.5,
        color = "black",
        max.overlaps = 60,
        min.segment.length = 0.2,
        box.padding = 0.3
      )
    }
  } +
  # Color scheme
  scale_color_manual(values = c("lightgrey", "darkgreen")) +
  # Position in Mb
  scale_x_continuous(
    name = "Position (Mb)",
    labels = ~ . / 1000000  # Convert bp to Mb
  ) +
  # Facet by chromosome
  facet_wrap(
    ~ chr_clean, 
    scales = "free_x",
    ncol = 5
  ) +
  # Labels and title
  labs(
    y = "Average Pairwise Coalescence Time (Myr)",
    title = "Chimpanzees Coalescence Times by Chromosome",
    caption = "Genome-wide threshold shown"
  ) +
  # Theme settings
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12)
  )

print(test2)
ggsave("all_Chimpanzees_avg_pairwise_tmrca_chromowide.png", test2, dpi = 300, width = 12, height= 14)