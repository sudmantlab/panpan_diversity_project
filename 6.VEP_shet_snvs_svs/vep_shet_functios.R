setwd("/Users/joanocha/Google Drive/My Drive/POSTDOC/PANPAN/analysis/Figure3_VEP_snvs_svs/")
# --- 1. SETUP: Load Libraries and Define Paths ---
# Load all necessary libraries once
library(readr)
library(tidyverse)
library(ggrepel)
library(ggnewscale)

# Define base path to avoid repeating it
base_path <- "/Users/joanocha/Google Drive/My Drive/POSTDOC/PANPAN/analysis/Figure3_VEP_snvs_svs"

# --- 2. LOAD COMMON DATA (Done once) ---
# Load Pritchard, gProfiler, and annotation files that are shared between analyses
gprofiler_data <- read_csv(file.path(base_path, "gProfiler_hsapiens_6-29-2025_9-12-33 PM.csv"))
pritchard_file_data <- read_csv(file.path(base_path, "Pritchard_Table.csv"))

# Function to load and filter annotation data to reduce repetition
load_annotation <- function(file_path) {
  read_delim(file_path) %>%
    distinct(name, .keep_all = TRUE) %>%
    filter(gene_biotype == "protein_coding",
           !(str_detect(name, "^LOC") & is.na(description)),
           !(str_detect(name, "^LOC") & str_detect(description, "^putatively uncharacterized")),
           !(str_detect(name, "^LOC") & str_detect(description, "^putative uncharacterized")),
           !(str_detect(name, "^LOC") & str_detect(description, "^uncharacterized")),
           !(str_detect(description, "non-protein coding"))
    )
}

human_unique_annotation_data <- load_annotation(file.path(base_path, "ht2t_gene.bed.gz"))
chimp_unique_annotation_data <- load_annotation(file.path(base_path, "mPanTro3_gene.bed.gz"))
bonobo_unique_annotation_data <- load_annotation(file.path(base_path, "mPanPan1_gene.bed.gz"))

# Process Pritchard data once
gprofiler_name_map <- select(gprofiler_data, initial_alias, name)
pritchard_merged <- left_join(pritchard_file_data, gprofiler_name_map, by = c("ensg" = "initial_alias"))

# --- 3. PROCESS AND COMBINE VARIANT DATA ---
# Create a function to process gene impact data for a given variant type
process_variant_data <- function(variant_type_path, variant_type_name) {
  gene_impact_data <- read_tsv(file.path(base_path, variant_type_path)) %>%
    filter(!(str_detect(name, "^LOC") & str_detect(description, "spliceosomal RNA|small.*RNA")))
  
  gene_impact_data %>%
    pivot_longer(
      cols = c("human", "chimpanzee", "bonobo"),
      names_to = "species",
      values_to = "impact_level"
    ) %>%
    mutate(variant_type = variant_type_name) # Add the new variant_type column
}

# Process SV and SNV data and combine them
combined_impact_data <- bind_rows(
  process_variant_data("svs_vep/gene_impact_wide.tsv", "SVs"),
  process_variant_data("snvs_vep/gene_impact_wide.tsv", "SNVs")
)

# --- 4. MERGE AND FINALIZE DATASET ---
filter_pritchard <- function(pritchard_df, annotation_df, species_name) {
  pritchard_df %>%
    filter(!is.na(post_mean) & !is.na(name)) %>%
    filter(name %in% annotation_df$name) %>%
    mutate(species = species_name)
}

pritchard_merged_filtered <- bind_rows(
  filter_pritchard(pritchard_merged, human_unique_annotation_data, "human"),
  filter_pritchard(pritchard_merged, chimp_unique_annotation_data, "chimpanzee"),
  filter_pritchard(pritchard_merged, bonobo_unique_annotation_data, "bonobo")
) %>%
  mutate(
    selection_regime = case_when(
      post_mean < 1e-4 ~ "Nearly neutral",
      post_mean < 1e-3 ~ "Weak selection",
      post_mean < 1e-2 ~ "Strong selection",
      TRUE ~ "Extreme selection"
    ),
    selection_regime = factor(selection_regime, levels = c("Nearly neutral", "Weak selection", "Strong selection", "Extreme selection"))
  )

# Final merged dataset
final_data <- left_join(pritchard_merged_filtered, combined_impact_data, by = c("name", "species")) %>%
  mutate(
    impact_level = ifelse(is.na(impact_level), 0, impact_level),
    impact_level_factor = case_when(
      impact_level %in% c(0, 1, 2) ~ "Absent/Modifier/Low",
      impact_level == 3 ~ "Moderate",
      impact_level == 4 ~ "High"
    ),
    impact_level_factor = factor(impact_level_factor, levels = c("Absent/Modifier/Low", "Moderate", "High")),
    species = tools::toTitleCase(species)
  )

regime_rects <- data.frame(
  selection_regime = factor(c("Nearly neutral", "Weak selection", "Strong selection", "Extreme selection"), levels = c("Nearly neutral", "Weak selection", "Strong selection", "Extreme selection")),
  ymin = c(0, 1e-4, 1e-3, 1e-2),
  ymax = c(1e-4, 1e-3, 1e-2, Inf)
)


summary_stats <- final_data %>%
  group_by(species, variant_type, impact_level_factor) %>%
  summarise(
    n_genes = n(),
    mean_shet = mean(post_mean, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    combined_label = str_glue("mean sHet={round(mean_shet, 2)}\n(n={n_genes})")
  )

# --- 5. PLOTTING ---

box_plot_combined <- ggplot(final_data, aes(x = impact_level_factor)) +
  # Add the background rectangles for selection regimes
  geom_rect(
    data = regime_rects,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = selection_regime),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  scale_fill_manual(
    name = "Selection Regime",
    values = c("Nearly neutral" = "#a1d99b", "Weak selection" = "#fee08b", "Strong selection" = "#fdae61", "Extreme selection" = "#f46d43")
  ) +
  new_scale_fill() +
  # Add the boxplots, defining the y aesthetic here
  geom_boxplot(aes(y = post_mean, fill = impact_level_factor), outlier.shape = NA) +
  
  # Add the combined mean and count labels below the boxplots
  geom_text(
    data = summary_stats,
    aes(label = combined_label, y = 10^-4.5),
    size = 2.5,
    vjust = 0.5,
    lineheight = .9
  ) +
  
  # Add a small black dot for the mean, defining the y aesthetic here
  stat_summary(
    aes(y = post_mean),
    fun = mean,
    geom = "point",
    shape = 16,
    size = 1.5,
    color = "black"
  ) +
  
  scale_fill_manual(
    name = "VEP Impact",
    values = c("Absent/Modifier/Low" = "lightgrey", "Moderate" = "darkorange", "High" = "darkred")
  ) +
  scale_y_log10(
    name = "Gene Constraint Score log10(sHet)",
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    limits = c(10^-4.6, 1.5)
  ) +
  facet_grid(variant_type ~ species, scales = "free_x", space = "free_x") +
  labs(title = "sHet by VEP Impact", x = "") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    strip.background = element_rect(color="black"),
    # Add this line to remove x-axis text labels
    axis.text.x = element_blank()
  )

# Display the plot
print(box_plot_combined)
ggsave("box_plot_combined.pdf", box_plot_combined, dpi=300, width =12, height=5)

## PLOT 2: Density Plot
density_plot_combined <- ggplot(final_data, aes(x = post_mean)) +
  geom_rect(data = regime_rects, aes(xmin = ymin, xmax = ymax, ymin = -Inf, ymax = Inf, fill = selection_regime), inherit.aes = FALSE, alpha = 0.2) +
  scale_fill_manual(name = "Selection Regime", values = c("Nearly neutral" = "#a1d99b", "Weak selection" = "#fee08b", "Strong selection" = "#fdae61", "Extreme selection" = "#f46d43")) +
  new_scale_fill() +
  geom_density(aes(fill = impact_level_factor), alpha = 0.6) +
  scale_fill_manual(name = "VEP Impact", values = c("Absent/Modifier/Low" = "lightgrey", "Moderate" = "darkorange", "High" = "darkred")) +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  facet_grid(variant_type ~ species) + # Facet by both variant type and species
  labs(title = "Density of sHet by VEP Impact", x = "Gene Constraint Score log10(sHet)", y = "Density") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right", strip.background = element_rect(color="black"))

print(density_plot_combined)
ggsave("density_plot_combined.pdf", density_plot_combined, dpi=300, width =12, height=5)

# --- PLOT 3: Ranked Line Plot (High Impact Only, Top/Bottom 10) ---

# Filter for HIGH impact genes only
high_impact_genes <- final_data %>%
  filter(impact_level_factor == "High") %>%
  group_by(species, variant_type) %>%
  arrange(desc(post_mean)) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# Create dataframes for the labels: top 10 and bottom 10
top_10_labels <- high_impact_genes %>%
  filter(rank <= 10)

bottom_10_labels <- high_impact_genes %>%
  group_by(species, variant_type) %>%
  slice_tail(n = 10) %>%
  ungroup()

# Create the plot
ranked_line_plot_combined <- ggplot(high_impact_genes, aes(x = rank, y = post_mean)) +
  geom_rect(data = regime_rects, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = selection_regime), inherit.aes = FALSE, alpha = 0.2) +
  geom_line(color = "gray40") +
  # Add points and labels for the TOP 10 genes
  geom_point(data = top_10_labels, aes(color = species), size = 1) +
  geom_text_repel(
    data = top_10_labels,
    aes(label = name, color = species),
    size = 2,
    fontface = "bold",
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = 'grey50'
  ) +
  # Add points and labels for the BOTTOM 10 genes
  geom_point(data = bottom_10_labels, aes(color = species), size = 1) +
  geom_text_repel(
    data = bottom_10_labels,
    aes(label = name, color = species),
    size = 2,
    fontface = "bold",
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = 'grey50'
  ) +
  facet_grid(variant_type ~ species, scales = "free_x") +
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_fill_manual(name = "Selection Regime", values = c("Nearly neutral" = "#a1d99b", "Weak selection" = "#fee08b", "Strong selection" = "#fdae61", "Extreme selection" = "#f46d43")) +
  labs(
    title = "Constraint (sHet) of High Impact Genes",
    x = "Gene Rank",
    y = "Gene Constraint Score log10(sHet)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(color="black")
  )

# Display the plot
print(ranked_line_plot_combined)
ggsave("ranked_line_plot_combined.pdf", ranked_line_plot_combined, dpi=300, width =12, height=6)




############# STATS ###########################
# --- 1. Isolate the data for comparison ---
# High-impact SNVs in Chimps and Humans
high_impact_snvs <- final_data %>%
  filter(variant_type == "SNVs" & impact_level_factor == "High") %>%
  filter(species %in% c("Chimpanzee", "Human"))

high_impact_svs <- final_data %>%
  filter(variant_type == "SVs" & impact_level_factor == "High") %>%
  filter(species %in% c("Chimpanzee", "Human"))


# --- 2. Perform the Wilcoxon Rank-Sum Tests ---
wilcox_test_snv <- wilcox.test(post_mean ~ species, data = high_impact_snvs)
wilcox_test_sv <- wilcox.test(post_mean ~ species, data = high_impact_svs)

print("--- Wilcoxon Test for HIGH impact SNVs (Chimp vs. Human) ---")
print(wilcox_test_snv)

print("--- Wilcoxon Test for HIGH impact SVs (Chimp vs. Human) ---")
print(wilcox_test_sv)


# --- 3. Perform the T-tests ---

t_test_snv <- t.test(post_mean ~ species, data = high_impact_snvs)
t_test_sv <- t.test(post_mean ~ species, data = high_impact_svs)
print("--- T-test for HIGH impact SNVs (Chimp vs. Human) ---")
print(t_test_snv)
print("--- T-test for HIGH impact SVs (Chimp vs. Human) ---")
print(t_test_sv)