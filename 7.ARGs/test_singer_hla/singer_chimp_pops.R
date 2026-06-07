setwd("/Users/joanocha/Desktop/SINGER")
library(ggplot2)
library(scales)
library(readr)
library(dplyr)
library(ggrepel)
library(tidyr)

coalescence_df <- read_csv("chimpanzee_genome-wide_metrics_annotatedbyPop.csv")
coalescence_df$T_pooled_myr <- as.numeric(coalescence_df$T_pooled * 25 / 1e6) 
coalescence_df$T_within_myr <- as.numeric(coalescence_df$T_within * 25 / 1e6)
coalescence_df$avg_tmrca_myr <- as.numeric(coalescence_df$avg_tmrca) * 25 / 1e6

# Track initial number of windows
original_window_count <- nrow(coalescence_df)

# Create chromosome cleaning upfront for filtering
coalescence_df <- coalescence_df %>%
  mutate(
    chr_clean = sub("_.*", "", chromosome),
    chr_clean = factor(chr_clean, levels = c(paste0("chr", 1:23), "chrX", "chrY"))
  )

# Identify and track all problematic windows
problematic_windows <- coalescence_df %>%
  filter(population == "Western") %>%
  mutate(
    issue_type = case_when(
      chr_clean %in% c("chrX", "chrY") ~ "X/Y_chromosome",
      is.na(avg_tmrca) | is.nan(avg_tmrca) ~ "NA/Nan_TMRCA",
      avg_tmrca_myr < 0 ~ "Negative_TMRCA",
      T_within_myr < 0 ~ "Negative_pairwise",
      TRUE ~ "Valid"
    )
  ) %>%
  filter(issue_type != "Valid")

# Assuming you have a gene annotation data frame with columns:
# First, ensure genes column is cleaned and prepared
problematic_windows_summary <- coalescence_df %>%
  filter(population == "Western") %>%
  mutate(
    issue_type = case_when(
      chr_clean %in% c("chrX", "chrY") ~ "X/Y_chromosome",
      is.na(avg_tmrca) | is.nan(avg_tmrca) ~ "NA/Nan_TMRCA",
      avg_tmrca_myr < 0 ~ "Negative_TMRCA",
      T_within_myr < 0 ~ "Negative_pairwise",
      TRUE ~ "Valid"
    )
  ) %>%
  filter(issue_type != "Valid") %>%
  # Remove sex chromosomes
  filter(!(chr_clean %in% c("chrX", "chrY"))) %>%
  # Handle NA values in genes column
  mutate(genes = replace_na(genes, "")) %>%
  # Group by original chromosome
  group_by(chromosome) %>%
  summarize(
    chr_clean = first(chr_clean),
    region_start = min(start),
    region_end = max(end),
    total_problematic_windows = n(),
    unique_issue_types = paste(sort(unique(issue_type)), collapse = ", "),
    min_tmrca = min(avg_tmrca_myr, na.rm = TRUE),
    max_tmrca = max(avg_tmrca_myr, na.rm = TRUE),
    # Extract gene information from existing genes column
    gene_count = {
      all_genes <- paste(genes, collapse = ",")
      split_genes <- unlist(strsplit(all_genes, ","))
      clean_genes <- trimws(split_genes)
      clean_genes <- clean_genes[clean_genes != ""]
      length(unique(clean_genes))
    },
    gene_names = {
      all_genes <- paste(genes, collapse = ",")
      split_genes <- unlist(strsplit(all_genes, ","))
      clean_genes <- trimws(split_genes)
      clean_genes <- clean_genes[clean_genes != ""]
      unique_genes <- unique(clean_genes)
      if (length(unique_genes) == 0) "" else paste(unique_genes, collapse = ", ")
    },
    .groups = "drop"
  )

#write_csv(problematic_windows_summary, "problematic_windows_summary.csv")




# Calculate proportions
removed_count <- nrow(problematic_windows)
proportion_removed <- removed_count / original_window_count

# Breakdown of issues
issue_summary <- problematic_windows %>%
  count(issue_type) %>%
  mutate(proportion = n / original_window_count)

# Print removal statistics
print(paste("Original windows:", original_window_count))
print(paste("Total removed windows:", removed_count))
print(paste("Proportion removed:", round(proportion_removed, 4)))
print("Issue breakdown:")
print(issue_summary)

# Remove problematic windows and chromosomes
coalescence_df <- coalescence_df %>%
  filter(
    !chr_clean %in% c("chrX", "chrY"),
    !(is.na(avg_tmrca) | is.nan(avg_tmrca)),
    avg_tmrca_myr >= 0,
  )

# Process remaining data
coalescence_df <- coalescence_df %>%
  mutate(midpoint = (start + end) / 2) %>%
  arrange(chr_clean, start) %>%
  mutate(row_index = row_number())

# Recalculate chromosome label positions
axisdf <- coalescence_df %>%
  group_by(chr_clean) %>%
  summarize(center = mean(row_index))


# Calculate threshold (99.99%)
threshold <- quantile(coalescence_df$T_within_myr, probs = 0.9999, na.rm = TRUE)
print(paste("99.99% threshold:", threshold))

# Custom population colors
pop_colors <- c(
  "Western" = "#9dced9",
  "hybrid" = "#bbc671",
  "Eastern" = "#ffb35a",
  "Central" = "#4c5d4c"
)

# Create color groups for plotting
coalescence_df <- coalescence_df %>%
  mutate(
    color_group = if_else(T_within_myr > threshold, population, "Below threshold")
  )

all_colors <- c("Below threshold" = "grey60", pop_colors)

# Prepare gene labels - get top window per gene above threshold
gene_labels <- coalescence_df %>%
  filter(T_within_myr > threshold, !is.na(genes)) %>%
  # Split multiple genes in same window
  mutate(gene = strsplit(as.character(genes), ",")) %>%
  unnest(gene) %>%
  mutate(gene = trimws(gene)) %>%
  filter(gene != "") %>%
  # For each gene, keep only the window with highest TMRCA
  group_by(gene) %>%
  arrange(desc(T_within_myr)) %>%
  slice(1) %>%
  ungroup()

manhattan_plot <- ggplot(coalescence_df, aes(x = row_index, y = T_within_myr)) +
  # All points colored by group
  geom_point(aes(color = color_group), size = 1, alpha = 0.7) +
  # Threshold line
  geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
  # Gene labels for top windows
  geom_text_repel(
    data = gene_labels,
    aes(label = gene),
    size = 3,
    color = "black",
    min.segment.length = 0.1,
    segment.size = 0.3,
    box.padding = 0.5,
    max.overlaps = 50,
    force = 2,
    force_pull = 1
  ) +
  # Chromosome labels
  scale_x_continuous(
    name = "Chromosome",
    breaks = axisdf$center,
    labels = axisdf$chr_clean,
    expand = c(0.01, 0.01)
  ) +
  # Custom color scale
  scale_color_manual(
    name = "Population",
    values = all_colors,
    breaks = c("Below threshold", names(pop_colors)),
    labels = c("Below threshold", names(pop_colors))
  ) +
  # Theme adjustments
  theme_bw() +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  ) +
  labs(
    x = "Chromosome",
    y = "Average pairwise coalescence time (myr)",
    title = "Genome-wide Coalescence Time Distribution"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

ggsave("manhattan_plot_labeled.png", plot = manhattan_plot, 
       width = 14, height = 7, dpi = 300)


# Calculate NEW threshold for Tpooled_Twithin_ratio (99.99%)
threshold_tpooled <- quantile(coalescence_df$Tpooled_Twithin_ratio, 
                              probs = 0.9999, 
                              na.rm = TRUE)
print(paste("99.99% threshold for Tpooled/Twithin ratio:", threshold_tpooled))
coalescence_df <- coalescence_df %>%
  mutate(
    color_group = if_else(Tpooled_Twithin_ratio > threshold_tpooled, population, "Below threshold")
  )
all_colors <- c("Below threshold" = "grey60", pop_colors)
gene_labels_tpooled <- coalescence_df %>%
  filter(Tpooled_Twithin_ratio > threshold_tpooled, !is.na(genes)) %>%
  mutate(gene = strsplit(as.character(genes), ",")) %>%
  unnest(gene) %>%
  mutate(gene = trimws(gene)) %>%
  filter(gene != "") %>%
  group_by(gene) %>%
  arrange(desc(Tpooled_Twithin_ratio)) %>%
  slice(1) %>%
  ungroup()

#Manhattan plot for Tpooled/Twithin ratio
manhattan_plot_tpooled <- ggplot(coalescence_df, aes(x = row_index, y = Tpooled_Twithin_ratio)) +
  geom_point(aes(color = color_group), size = 1, alpha = 0.7) +
  geom_hline(yintercept = threshold_tpooled, linetype = "dashed", color = "red") +
  geom_text_repel(
    data = gene_labels_tpooled,
    aes(label = gene),
    size = 3,
    color = "black",
    min.segment.length = 0.1,
    segment.size = 0.3,
    box.padding = 0.5,
    max.overlaps = 50,
    force = 2,
    force_pull = 1
  ) +
  scale_x_continuous(
    name = "Chromosome",
    breaks = axisdf$center,
    labels = axisdf$chr_clean,
    expand = c(0.01, 0.01)
  ) +
  scale_color_manual(
    name = "Population",
    values = all_colors,
    breaks = c("Below threshold", names(pop_colors)),
    labels = c("Below threshold", names(pop_colors))
  ) +
  # Theme adjustments
  theme_bw() +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  ) +
  labs(
    x = "Chromosome",
    y = "Tpooled/Twithin ratio",
    title = "Genome-wide Tpooled/Twithin Ratio Distribution"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
print(manhattan_plot_tpooled)
ggsave("manhattan_plot_Tpooled_ratio.png", plot = manhattan_plot_tpooled, 
       width = 14, height = 7, dpi = 300)

