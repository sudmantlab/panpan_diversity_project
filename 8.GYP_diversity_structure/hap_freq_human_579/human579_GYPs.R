library(tidyverse)
library(scales) 
library(malariaAtlas)
library(terra)
library(ggnewscale)
library(geodata)
library(scatterpie)
library(patchwork)
library(ggrepel)

# --- 1. DATA LOADING ---
struct_df <- read_tsv("gyp_human_structure_grouping.tsv")
meta_df <- read_tsv("human579.meta.tsv", 
                    col_names = c("sample", "sex_chrom", "id", "sex", "continent", "sub_region", "country", "data_source"))
new_meta <- read_tsv("METADATA_4099_p1p2.tsv")

structure_order_stacked <- c(
  "H3.1 (GYP E-B-A)", 
  "H2.1", "H2.2", "H3.2", "H3.3", "H4.1", "H4.2", "H7"
)

# Using your specific base_cols with lightgrey for Ref
base_cols <- c("#E41A1C", "#377EB8", 'lightgrey', "#984EA3", "#FF7F00", 
               "#FFFF33", "#A65628", "#F781BF", "#999999", "#000000","#4DAF4A")

haplo_colors_final <- c(
  "H3.1 (GYP E-B-A)"  = base_cols[3], # lightgrey
  "H2.1"              = base_cols[1], # Red
  "H2.2"              = base_cols[2], # Blue
  "H3.2"              = base_cols[4], # Purple
  "H3.3"              = base_cols[5], # Orange
  "H4.1"              = base_cols[6], # Yellow
  "H4.2"              = base_cols[7], # Brown
  "H7" = base_cols[8]  # Pink
)


df_unified <- struct_df %>%
  inner_join(meta_df, by = "sample") %>%
  mutate(location = if_else(country != ".", country, sub_region)) 

pop_lookup <- new_meta %>%
  select(
    hgdp_tgp_meta.Population,
    Latitude = hgdp_tgp_meta.Latitude,
    Longitude = hgdp_tgp_meta.Longitude
  ) %>%
  distinct(hgdp_tgp_meta.Population, .keep_all = TRUE)

df_master <- df_unified %>%
  left_join(pop_lookup, by = c("location" = "hgdp_tgp_meta.Population")) %>%
  mutate(hap_structure = case_when(
    group == 1 ~ "H2.1",
    group == 2 ~ "H2.2",
    group == 3 ~ "H3.1 (GYP E-B-A)",
    group == 4 ~ "H3.2",
    group == 5 ~ "H3.3",
    group == 6 ~ "H4.1",
    group == 7 ~ "H4.2",
    group == 8 ~ "H7 (Dantu Hybrid)",
    TRUE ~ as.character(group)
  )) %>%
  mutate(hap_structure = factor(hap_structure, levels = structure_order_stacked))

df_freq <- df_master %>%
  group_by(continent, location, hap_structure) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(continent, location) %>%
  mutate(frequency = count / sum(count)) %>%
  filter(continent != "Ref") %>%
  group_by(location) %>%
  mutate(pop_n = sum(count)) %>%
  ungroup() %>%
  group_by(continent) %>%
  mutate(cont_n = sum(count)) %>%
  ungroup() %>%
  mutate(
    location_n = paste0(location, " (n=", pop_n, ")"),
    continent_n = paste0(continent, " (n=", cont_n, ")")
  ) %>%
  arrange(continent, location) %>%
  mutate(location_n = factor(location_n, levels = unique(location_n)))

unified_stack_plot_annotated <- ggplot(df_freq, aes(x = location_n, y = frequency, fill = hap_structure)) +
  geom_col(position = "stack", color = "white", linewidth = 0.1) +
  facet_grid(~ continent_n, scales = "free_x") + 
  theme_classic() +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) + 
  scale_fill_manual(values = haplo_colors_final, breaks = structure_order_stacked) + 
  labs(title = "Composition of GYP Gene Structures per Location", x = "Population (Sample Size)", y = "Frequency", fill = "Structure Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(face = "bold", size = 10),
        panel.border = element_rect(fill = NA, color = "lightgray"),
        legend.position = "bottom", panel.spacing = unit(1, "lines"))

ggsave("Figure5_GYP_Stacked_Full.pdf", unified_stack_plot_annotated, width = 12, height = 5)


df_plot_dodged_filtered <- df_freq %>% filter(hap_structure != "H3.1 (GYP E-B-A)")

unified_dodge_zoomed <- ggplot(df_plot_dodged_filtered, aes(x = location_n, y = frequency, fill = hap_structure)) +
  geom_col(position = "dodge", color = "white", linewidth = 0.1) +
  facet_grid(~ continent_n, scales = "free_x") + 
  theme_classic() +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.1))) + 
  scale_fill_manual(values = haplo_colors_final, breaks = setdiff(structure_order_stacked, "H3.1 (GYP E-B-A)")) + 
  labs(title = "Relative Frequency of GYP Structural Variants", x = "Population (Sample Size)", y = "Frequency (Excluding Ref)", fill = "Variant Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(face = "bold", size = 10),
        panel.border = element_rect(fill = NA, color = "grey80"),
        legend.position = "bottom", panel.spacing = unit(1, "lines"))

ggsave("Figure5_GYP_Dodged_Zoomed.pdf", unified_dodge_zoomed, width = 12, height = 5)


