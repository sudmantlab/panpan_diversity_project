setwd("/Users/joanocha/Google Drive/My Drive/POSTDOC/PANPAN/analysis/Figure3_VEP_snvs_svs/svs_vep/")
library(readr)
library(readxl)
library(readr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(dunn.test)
library(ggridges)
library(ggpubr)
library(ggrepel)
library(ggnewscale)
# --- Load the gene impact wide file ---
gene_impact_file <- "gene_impact_wide_by_subspecies.tsv"
gene_impact_data <- read_tsv(gene_impact_file)
gene_impact_data <- gene_impact_data %>%
  rename(name = gene) %>%
  filter(!(str_detect(name, "^LOC") & str_detect(description, "spliceosomal RNA"))) %>%
  filter(!(str_detect(name, "^LOC") & str_detect(description, "small") & str_detect(description, "RNA"))
  )

gprofiler_file <- "/Users/joanocha/Google Drive/My Drive/POSTDOC/PANPAN/analysis/Figure3_VEP_snvs_svs/gProfiler_hsapiens_6-29-2025_9-12-33 PM.csv"
gprofiler_data <- read_csv(gprofiler_file)

pritchard_file <- "/Users/joanocha/Google Drive/My Drive/POSTDOC/PANPAN/analysis/Figure3_VEP_snvs_svs/Pritchard_Table.csv"
pritchard_file <- read_csv(pritchard_file)

gprofiler_name_map <- select(gprofiler_data, initial_alias, name)
pritchard_merged <- left_join(pritchard_file, gprofiler_name_map, by = c("ensg" = "initial_alias"))


chimp_annotation_file<-read_delim("/Users/joanocha/Google Drive/My Drive/POSTDOC/PANPAN/analysis/Figure3_VEP_snvs_svs/mPanTro3_gene.bed.gz")
chimp_unique_annotation_data <- distinct(chimp_annotation_file, name, .keep_all = TRUE) %>%
  filter(gene_biotype == "protein_coding",
         !(str_detect(name, "^LOC") & is.na(description)),
         !(str_detect(name, "^LOC") & str_detect(description, "^putatively uncharacterized")),
         !(str_detect(name, "^LOC") & str_detect(description, "^putative uncharacterized")),
         !(str_detect(name, "^LOC") & str_detect(description, "^uncharacterized")),
         !(str_detect(description, "non-protein coding"))
  )


chimp_pritchard_merged_filtered <- pritchard_merged %>%
  filter(!is.na(post_mean) & !is.na(name)) %>%
  filter(name %in% chimp_unique_annotation_data$name) %>%
  mutate(species = "chimpanzee") 

pritchard_merged_filtered <- chimp_pritchard_merged_filtered %>%
  mutate(
    selection_regime = case_when(
      post_mean < 1e-4 ~ "Nearly neutral",
      post_mean >= 1e-4 & post_mean < 1e-3 ~ "Weak selection",
      post_mean >= 1e-3 & post_mean < 1e-2 ~ "Strong selection",
      post_mean >= 1e-2 ~ "Extreme selection",
      TRUE ~ "Other"
    ),
    # Convert to a factor to control plotting order
    selection_regime = factor(selection_regime, levels = c("Nearly neutral", "Weak selection", "Strong selection", "Extreme selection"))
  )

absent_genes_in_shet<-filter(gene_impact_data, !name%in%pritchard_merged_filtered$name) 


# --- MERGE GENE IMPACT WITH PRITCHARD DATA ---


all_subspecies <- c("Western", "Western x Central", "Central", "Eastern")

pritchard_post_mean <- select(pritchard_merged_filtered, name, post_mean, selection_regime, species)

gene_impact_data_long <- gene_impact_data %>%
  pivot_longer(
    cols = all_of(all_subspecies), # Use all_of() for safety
    names_to = "subspecies",
    values_to = "impact_level"
  )


# grid of all genes from Pritchard data and all subspecies
gene_subspecies_grid <- expand_grid(
  name = unique(pritchard_post_mean$name),
  subspecies = all_subspecies
)

# Join the Pritchard selection data to this grid. This adds post_mean and selection_regime to every gene-subspecies row.
gene_impact_with_post_mean <- gene_subspecies_grid %>%
  left_join(pritchard_post_mean, by = "name") %>%
  # and join the actual impact data. Where a gene-subspecies combo has no impact data, `impact_level` will be NA.
  left_join(gene_impact_data_long, by = c("name", "subspecies")) %>%
  # Replace the NA impact levels with 0, representing "Absent".
  mutate(
    impact_level = replace_na(impact_level, 0)
  )


# Define the descriptive labels and the desired order
impact_labels <- c(
  `0` = "absent",
  `1` = "modifier",
  `2` = "low",
  `3` = "moderate",
  `4` = "high"
)
impact_levels <- c("absent", "modifier", "low", "moderate", "high")

# Calculate the counts
impact_summary <- gene_impact_with_post_mean %>%
  mutate(
    impact = recode(impact_level, !!!impact_labels)
  ) %>%
  mutate(
    impact = factor(impact, levels = impact_levels)
  ) %>%
  count(subspecies, impact, name = "number_of_cases")
print(impact_summary)

#analysis_data <- gene_impact_with_post_mean %>%
#  filter(!is.na(post_mean))

## fill NA fields with 0 


# Convert numeric impact levels into a new two-level factor
gene_impact_with_post_mean<- gene_impact_with_post_mean %>%
  mutate(
    impact_level_factor = case_when(
      impact_level %in% c(0, 1, 2) ~ "Absent/Modifier/Low",
      impact_level == 3 ~ "Moderate",
      impact_level == 4 ~ "High",
      #impact_level %in% c(0, 1) ~ "Not impactful",
      #impact_level == 1 ~ "Modifier",
      #impact_level %in% c(2, 3, 4) ~ "Impactful",
      TRUE ~ NA_character_ # Should not happen due to prior filter
    ),
    subspecies = tools::toTitleCase(subspecies)
  ) %>%
  # Convert to factor to control plotting order
  # mutate(impact_level_factor = factor(impact_level_factor, levels = c("Absent", "Not impactful", "Impactful")))
  mutate(impact_level_factor = factor(impact_level_factor, levels = c("Absent/Modifier/Low", "Moderate", "High")))


absent_gene_impact_with_post_mean <- gene_impact_with_post_mean %>%
  filter(impact_level == 0)
write_csv(absent_gene_impact_with_post_mean, "snvs_vep/absent_gene_impact_with_sHet_chimp_subspecies.csv")

gene_impact_with_post_mean <- gene_impact_with_post_mean %>%
  #filter(impact_level_factor != "Absent")  %>%
  filter(subspecies != "Western x Central") 




# Create a data frame for the background rectangles representing selection regimes
regime_rects <- data.frame(
  selection_regime = factor(c("Nearly neutral", "Weak selection", "Strong selection", "Extreme selection"),
                            levels = c("Nearly neutral", "Weak selection", "Strong selection", "Extreme selection")),
  ymin = c(0, 1e-4, 1e-3, 1e-2),
  ymax = c(1e-4, 1e-3, 1e-2, Inf)
)
# Create the plot
box_plot_with_regimes <- ggplot(gene_impact_with_post_mean, aes(x = impact_level_factor, y = post_mean)) +
  # Add the background rectangles for selection regimes
  geom_rect(
    data = regime_rects,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = selection_regime),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  # Define the fill colors for the background rectangles
  scale_fill_manual(
    name = "Selection Regime",
    values = c("Nearly neutral" = "#a1d99b", "Weak selection" = "#fee08b", "Strong selection" = "#fdae61", "Extreme selection" = "#f46d43")
  ) +
  # Introduce a new fill scale for the boxplots
  new_scale_fill() +
  # Add the boxplots
  geom_boxplot(aes(fill = impact_level_factor)) +
  # Define the fill colors for the boxplots
  scale_fill_manual(
    name = "VEP Impact",
    values = c("Absent/Modifier/Low" = "white" , "Moderate" = "lightgray", "High" = "darkgray")
    #values = c("Absent" = "white" , "Modifier" = "lightgray", "Impactful" = "darkgray")
  ) +
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  facet_wrap(~subspecies, ncol = 4) +
  labs(
    title = "sHet by VEP SVs - CHIMPANZEE SUBSPECIES",
    x = "",
    y = "Gene Constraint Score log10(sHet)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

print(box_plot_with_regimes)



density_plot_with_regimes <- ggplot(gene_impact_with_post_mean, aes(x = post_mean)) +
  # Add the background rectangles for selection regimes (x-axis version)
  geom_rect(
    data = regime_rects,
    aes(xmin = ymin, xmax = ymax, ymin = -Inf, ymax = Inf, fill = selection_regime),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  # Define the fill colors for the background rectangles
  scale_fill_manual(
    name = "Selection Regime",
    values = c("Nearly neutral" = "#a1d99b", "Weak selection" = "#fee08b", "Strong selection" = "#fdae61", "Extreme selection" = "#f46d43")
  ) +
  # Introduce a new fill scale for the density curves
  new_scale_fill() +
  # Add the density curves
  geom_density(aes(fill = impact_level_factor), alpha = 0.6) +
  # Define the fill colors for the density curves
  scale_fill_manual(
    name = "VEP",
    values = c("Absent/Modifier/Low" = "white" , "Moderate" = "lightgray", "High" = "darkgray")
    #values = c("Absent" = "white" , "Modifier" = "lightgray", "Impactful" = "black")
  ) +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  facet_wrap(~subspecies, ncol = 4) +
  labs(
    title = "Density of sHet by VEP for SVs - CHIMPANZEE SUBSPECIES",
    x = "Gene Constraint Score log10(sHet)",
    y = "Density SVs"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

print(density_plot_with_regimes)


impactful_genes <- present_gene_impact_with_post_mean %>%
  filter(impact_level_factor == "Impactful")

# Rank genes within each species from highest to lowest post_mean
ranked_impactful_genes <- impactful_genes %>%
  group_by(subspecies) %>%
  arrange(desc(post_mean)) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# Create dataframes for the labels: top 10 and bottom 10
# Top 50 (most constrained)
top_5_labels <- ranked_impactful_genes %>%
  filter(rank <= 5)

# Bottom 5 (least constrained)
bottom_5_labels <- ranked_impactful_genes %>%
  group_by(subspecies) %>%
  slice_tail(n = 5) %>% # Takes the last 10 rows for each species
  ungroup()



ranked_line_plot_combined <- ggplot(ranked_impactful_genes, aes(x = rank, y = post_mean)) +
  # Add the background rectangles for selection regimes
  geom_rect(
    data = regime_rects,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = selection_regime),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  # Add the line connecting all impactful genes
  geom_line(color = "gray40") +
  geom_point(data =top_5_labels, aes(color = subspecies), size = 1) +
  geom_text_repel(
    data =top_5_labels,
    aes(label = name, color = subspecies),
    size = 2,
    fontface = "bold",
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = 'grey50'
  ) +
  # Add points and labels for the BOTTOM 10 genes
  geom_point(data = bottom_5_labels, aes(color = subspecies), size = 1) +
  geom_text_repel(
    data = bottom_5_labels,
    aes(label = name, color = subspecies),
    size = 2,
    fontface = "bold",
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = 'grey50'
  ) +
  # Facet by species
  facet_wrap(~subspecies, scales = "free_x") +
  # Use a log scale for the y-axis for better visualization
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  # Define the fill colors for the background rectangles
  scale_fill_manual(
    name = "Selection Regime",
    values = c("Nearly neutral" = "#a1d99b", "Weak selection" = "#fee08b", "Strong selection" = "#fdae61", "Extreme selection" = "#f46d43")
  ) +
  # Add informative labels and a title
  labs(
    title = "Constraint (sHet) of Impactful Genes SVs - CHIMPANZEE SUBSPECIES",
    x = "Gene Rank",
    y = "Gene Constraint Score log10(sHet)"
  ) +
  # Use a clean theme
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    strip.text = element_text(face = "bold", size = 12)
  )

# Display the plot
print(ranked_line_plot_combined)