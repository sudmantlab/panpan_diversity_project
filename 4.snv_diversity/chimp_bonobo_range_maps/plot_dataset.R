# Load required libraries
library(tidyverse)

# Read the data
data <- read_tsv("PANPAN_HPRC_HGSVC_METADATA.txt")

# Create a function to get simplified dataset type
get_dataset_type2 <- function(dataset) {
  if (grepl("PANPAN", dataset)) return("PANPAN")
  return("HPRC+HGSVC3")
}

# Define colors for simplified dataset types
dataset_colors2 <- c(
  "PANPAN" = "darkgreen",
  "HPRC+HGSVC3" = "darkgrey"
)

# Prepare the data
plot_data <- data %>%
  mutate(dataset2 = sapply(dataset, get_dataset_type2),
         dataset_color = dataset_colors2[dataset2]) %>%
  group_by(p2, HPRC_meta.Population, dataset2, dataset_color, Continent.colors, HPRC_PANPAN.color) %>%
  summarise(count = n()) %>%
  ungroup()

# Define the correct order of continents
continent_order <- c("AFR", "AMR", "EAS", "EUR", "SAS")

# Reorder the data based on the defined continent order and ensure PANPAN samples are grouped
plot_data <- plot_data %>%
  mutate(p2 = factor(p2, levels = continent_order)) %>%
  arrange(p2, desc(dataset2), HPRC_meta.Population) %>%
  mutate(group = row_number())

# Calculate cumulative counts and angles for positioning
plot_data <- plot_data %>%
  mutate(
    end = cumsum(count),
    start = end - count,
    total = sum(count),
    start_angle = 2 * pi * start / total,
    end_angle = 2 * pi * end / total,
    mid_angle = (start_angle + end_angle) / 2
  )

# Function to create a data frame for each ring
create_ring_data <- function(data, r_inner, r_outer, explode_factor = 0) {
  data %>%
    rowwise() %>%
    mutate(
      explode_mid = if_else(dataset2 == "PANPAN", explode_factor, 0),
      x_mid = (r_inner + r_outer) / 2 * cos(mid_angle),
      y_mid = (r_inner + r_outer) / 2 * sin(mid_angle),
      x_explode = x_mid + explode_mid * cos(mid_angle),
      y_explode = y_mid + explode_mid * sin(mid_angle),
      x_diff = x_explode - x_mid,
      y_diff = y_explode - y_mid,
      x = list(c(r_inner * cos(seq(start_angle, end_angle, length.out = 50)),
                 r_outer * cos(seq(end_angle, start_angle, length.out = 50))) + x_diff),
      y = list(c(r_inner * sin(seq(start_angle, end_angle, length.out = 50)),
                 r_outer * sin(seq(end_angle, start_angle, length.out = 50))) + y_diff)
    ) %>%
    unnest(cols = c(x, y))
}

# Create data for each ring
dataset_ring <- create_ring_data(plot_data, 0.5, 1, explode_factor = 0.2) %>% mutate(ring = "dataset")
continent_ring <- create_ring_data(plot_data, 1.1, 1.6, explode_factor = 0.2) %>% mutate(ring = "continent")
population_ring <- create_ring_data(plot_data, 1.7, 2.2, explode_factor = 0.2) %>% mutate(ring = "population")

# Combine all rings
all_rings <- bind_rows(dataset_ring, continent_ring, population_ring)

# Create summary data for dataset and continent labels
dataset_labels <- plot_data %>%
  group_by(dataset2) %>%
  summarise(
    mid_angle = mean(mid_angle),
    count = sum(count)
  ) %>%
  mutate(
    x = 0.75 * cos(mid_angle) + if_else(dataset2 == "PANPAN", 0.2 * cos(mid_angle), 0),
    y = 0.75 * sin(mid_angle) + if_else(dataset2 == "PANPAN", 0.2 * sin(mid_angle), 0)
  )

continent_labels <- plot_data %>%
  group_by(p2) %>%
  summarise(
    mid_angle = mean(mid_angle),
    count = sum(count),
    is_panpan = any(dataset2 == "PANPAN")
  ) %>%
  mutate(
    x = 1.35 * cos(mid_angle) + if_else(is_panpan, 0.2 * cos(mid_angle), 0),
    y = 1.35 * sin(mid_angle) + if_else(is_panpan, 0.2 * sin(mid_angle), 0)
  )

# Create the plot
p <- ggplot(all_rings) +
  geom_polygon(aes(x = x, y = y, group = interaction(group, ring), 
                   fill = case_when(
                     ring == "dataset" ~ dataset_color,
                     ring == "continent" ~ Continent.colors,
                     ring == "population" ~ HPRC_PANPAN.color
                   )),
               alpha = 0.8) +
  # Dataset labels
  geom_text(data = dataset_labels,
            aes(x = x, y = y, label = dataset2),
            hjust = 0.5, vjust = 0.5, size = 4, color = "black") +
  # Continent labels
  geom_text(data = continent_labels,
            aes(x = x, y = y, label = p2),
            hjust = 0.5, vjust = 0.5, size = 4) +
  # Population labels (outside)
  geom_text(data = mutate(plot_data,
                          x = 2.3 * cos(mid_angle) + if_else(dataset2 == "PANPAN", 0.2 * cos(mid_angle), 0),
                          y = 2.3 * sin(mid_angle) + if_else(dataset2 == "PANPAN", 0.2 * sin(mid_angle), 0)),
            aes(x = x, y = y, 
                label = HPRC_meta.Population,
                angle = (mid_angle * 180 / pi + 90) %% 180 - 90),
            hjust = 0.5, vjust = 0.5, size = 3) +
  coord_fixed() +
  scale_fill_identity() +
  theme_void() +
  theme(legend.position = "none")

# Display the plot
print(p)

# Save the plot
ggsave("multi_level_population_plot_panpan_exploded.png", p, width = 12, height = 12, units = "in", dpi = 300)