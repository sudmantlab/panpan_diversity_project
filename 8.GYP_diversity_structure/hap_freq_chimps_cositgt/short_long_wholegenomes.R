# ==========================================
library(data.table)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(scatterpie)
library(tidyr)
library(dplyr)
library(ggnewscale)
library(scales)
library(stringr)
library(readr)
library(ggrastr)
library(terra)
library(geodata)

structure_order <- c(
  "H2.1", "H3.1", "H3.2", "H4.1", "H4.2", "H4.3", 
  "H5.1", "H5.2", "H5.3", "H5.4", "H5.2+H5.4", "H6"
)


base_cols_chimp <- c(
  "#E41A1C", "lightgrey", "#377EB8", "#4DAF4A", "#984EA3", 
  "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", 
  "#66c2a5", "#000000"
)

# The missing object definition:
haplo_colors_final <- setNames(base_cols_chimp, structure_order)

# Subspecies / Map Background Colors
subspecies_colors <- c(
  "verus"          = '#9dced9', 
  "troglodytes"    = '#4c5d4c', 
  "ellioti"        = '#fe604c', 
  "schweinfurthii" = '#ffb35a',
  "paniscus"       = '#e9d7cb'
)

subspecies_labels <- c(
  "verus"          = "Western",
  "troglodytes"    = "Central",
  "ellioti"        = "Nigeria-Cameroon",
  "schweinfurthii" = "Eastern",
  "paniscus"       = "Bonobo"  
)

shp_pt  <- st_read("pt/data_0.shp", quiet = TRUE)  # Chimp Ranges
shp_ppa <- st_read("ppa/data_0.shp", quiet = TRUE) # Bonobo Ranges

# --- 1. DATA LOADING & METADATA PREP ---
chimp_struct_df   <- read_tsv("cositgt_tables_results/gyp_chimp_structure_grouping.tsv")
cositgt_genotypes <- read_tsv("cositgt_tables_results/cositgt_combined_genotypes.tsv")
cositgt_conversion <- read_tsv("cositgt_tables_results/cosigt_group_conversion.tsv")

meta_data <- fread("METADATA.tsv")
meta_unique <- meta_data %>%
  mutate(sample_match_id = sub("_[A-Za-z].*", "", genomeName)) %>%
  select(sample_match_id, genomeName, Species_mtDNA, SpeciesCode, CommonName) %>%
  distinct()

# --- 2. DEFINE MAPPING & COLORS ---
structure_mapping <- c(
  "1"  = "H2.1", "2"  = "H3.1", "3"  = "H3.2", "4"  = "H4.1", 
  "5"  = "H4.2", "6"  = "H4.3", "7"  = "H5.1", "8"  = "H5.2", 
  "9"  = "H5.3", "10" = "H5.4", "11" = "H6",   "8-10" = "H5.2+H5.4"
)

structure_order <- c("H2.1", "H3.1", "H3.2", "H4.1", "H4.2", "H4.3", 
                     "H5.1", "H5.2", "H5.3", "H5.4", "H5.2+H5.4", "H6")

# --- 3. HARMONIZE SHORT-READ (COSITGT) DATA ---
cositgt_long <- cositgt_genotypes %>%
  select(sample = `#sample.id`, cluster.1, cluster.2) %>%
  pivot_longer(cols = starts_with("cluster"), 
               names_to = "hap_index", 
               values_to = "haplotype.group") %>%
  left_join(cositgt_conversion, by = "haplotype.group") %>%
  mutate(
    hap_structure = structure_mapping[as.character(cosigt_group_new)],
    data_type = "short_read",
    sample_match_id = sample 
  )

# --- 4. PREPARE LONG-READ DATA ---
chimp_long_read <- chimp_struct_df %>%
  mutate(
    hap_structure = structure_mapping[as.character(group)],
    data_type = "long_read",
    sample_match_id = sub("\\..*", "", sample) 
  )

# --- 5. UNIFY & MERGE WITH METADATA ---
df_chimp_final <- bind_rows(chimp_long_read, cositgt_long) %>%
  mutate(final_group = coalesce(as.character(group), cosigt_group_new)) %>%
  left_join(meta_unique, by = "sample_match_id") %>%
  mutate(
    CommonName = ifelse(grepl("Nigeri", CommonName, ignore.case = TRUE), 
                        "Nigeria-Cameroon", CommonName),
    hap_structure = factor(hap_structure, levels = structure_order)
  ) %>%
  filter(sample_match_id != "mPanTro3")


# --- CALCULATE INDEPENDENT FREQUENCIES ---
chimp_freq_sep <- df_chimp_final %>%
  filter(!is.na(CommonName), !is.na(hap_structure)) %>%  # Group by data_type so frequencies are independent
  group_by(data_type, CommonName, hap_structure) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(data_type, CommonName) %>%
  mutate(
    frequency = count / sum(count),
    total_n = sum(count)
  ) %>%
  ungroup()

chimp_freq_sep <- chimp_freq_sep %>%
  mutate(subspecies_label = paste0(CommonName, "\n(n=", total_n, ")"))

chimp_bar_faceted <- ggplot(chimp_freq_sep, aes(x = subspecies_label, y = frequency, fill = hap_structure)) +
  geom_col(position = "stack", color = "white", linewidth = 0.1) +
  facet_grid(. ~ data_type, scales = "free_x", space = "free_x") +
  theme_classic() +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_fill_manual(values = base_cols_chimp, name = "GYP Structure") +
  labs(
    title = "GYP Haplotype Composition by Sequencing Technology",
    x = "Subspecies (Sample Size)",
    y = "Frequency"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

print(chimp_bar_faceted)

# Coordinates for Pies
pie_coords <- data.frame(
  CommonName = c("Western", "Nigeria-Cameroon", "Central", "Eastern", "Bonobo", "Western x Central hybrid"),
  Pie_Lon    = c(-12,  8.5,  10,  35,  22, -10), 
  Pie_Lat    = c( 12,   9,  -8,  -5,  -5,   0)  
)

pie_data_sep <- chimp_freq_sep %>%
  select(data_type, CommonName, hap_structure, count) %>%
  pivot_wider(names_from = hap_structure, values_from = count, values_fill = 0) %>%
  left_join(pie_coords, by = "CommonName") %>%
  mutate(radius = sqrt(rowSums(across(any_of(structure_order)))) * 0.6)

africa <- ne_countries(continent = "Africa", returnclass = "sf")

map_chimp_faceted <- ggplot() +
  rasterise(geom_sf(data = africa, fill = "#efefef", color = "white"), dpi = 300) +
  new_scale_fill() +
  geom_scatterpie(data = pie_data_sep, aes(x = Pie_Lon, y = Pie_Lat, r = radius),
                  cols = intersect(structure_order, names(pie_data_sep)), color = "white", linewidth = 0.2) +
  scale_fill_manual(values = haplo_colors_final, name = "GYP Structure") +
  facet_wrap(~data_type, ncol = 2) +
  coord_sf(xlim = c(-20, 45), ylim = c(-15, 20), expand = FALSE) +
  theme_void() +
  theme(strip.text = element_text(face = "bold", size = 12), legend.position = "bottom")

print(map_chimp_faceted)



# --- FACETED  ---
pie_data_sep <- chimp_freq_sep %>%
  select(data_type, CommonName, hap_structure, count) %>%
  pivot_wider(names_from = hap_structure, values_from = count, values_fill = 0) %>%
  left_join(pie_coords, by = "CommonName") %>%
  mutate(
    Total_N = rowSums(across(any_of(structure_order))),
    radius = sqrt(Total_N) * 0.5 
  )

pie_cols_present <- intersect(structure_order, names(pie_data_sep))

africa <- ne_countries(continent = "Africa", returnclass = "sf")

map_chimp_aesthetic <- ggplot() +
  rasterise(
    geom_sf(data = africa, color="#6F7378", fill = "#6F7378", linewidth = 0.1),
    dpi = 300
  ) +

  rasterise(
    geom_sf(data = shp_pt, aes(fill = SUBSPECIES), color = "black", alpha = 0.8, linewidth = 0.2),
    dpi = 300
  ) +
  rasterise(
    geom_sf(data = shp_ppa, aes(fill = "paniscus"), color = "black", alpha = 0.8, linewidth = 0.2),
    dpi = 300
  ) +
  scale_fill_manual(values = subspecies_colors, labels = subspecies_labels, name = "IUCN Range") +
  

  new_scale_fill() + 
  

  geom_scatterpie(
    data = pie_data_sep, 
    aes(x = Pie_Lon, y = Pie_Lat, r = radius),
    cols = pie_cols_present, 
    color = "black", 
    linewidth = 0.2
  ) +
  scale_fill_manual(values = haplo_colors_final, name = "Structural Haplogroup") +
  
 
  geom_text(
    data = pie_data_sep, 
    aes(x = Pie_Lon, y = Pie_Lat + radius + 1.8, label = str_wrap(CommonName, width = 15)), 
    fontface = "bold", 
    size = 3.5,
    lineheight = 0.9
  ) +
  
  facet_wrap(~data_type, ncol = 2) +

  coord_sf(xlim = c(-20, 45), ylim = c(-15, 20), expand = FALSE) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"), 
    strip.text = element_text(face = "bold", size = 12, margin = margin(b = 10)),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  ) +
  labs(title = "Chimpanzee GYP Haplotype Diversity")
map_chimp_aesthetic

ggsave(
  'Figure_Chimp_Map_Faceted.pdf', 
  map_chimp_aesthetic, 
  width = 14, 
  height = 7
)



# --- PREPARE THE UNIFIED FREQUENCY DATA ---
chimp_freq_bar <- df_chimp_unified %>%
  # Group by subspecies and the new unified structure
  group_by(CommonName, hap_structure) %>%
  summarise(count = n(), .groups = "drop") %>%
  # Calculate proportions per subspecies
  group_by(CommonName) %>%
  mutate(
    frequency = count / sum(count),
    total_n = sum(count),
    subspecies_label = paste0(CommonName, " (n=", total_n, ")")
  ) %>%
  ungroup()
subspecies_order <- c("Western", "Nigeria-Cameroon", "Central", "Eastern", "Bonobo", "Western x Central hybrid")
chimp_freq_bar <- chimp_freq_bar %>%
  mutate(subspecies_label = factor(subspecies_label, 
                                   levels = unique(subspecies_label[order(match(CommonName, subspecies_order))])))

# --- RE-ORDER LEVELS (H5.2+H5.4 First) ---
unified_structure_order <- c(
  "H5.2+H5.4", 
   "H2.1", "H3.1","H3.2", "H4.1", "H4.2", "H4.3", "H5.1", "H5.3", "H6"
)

subspecies_order <- c(
  "Western", 
  "Nigeria-Cameroon", 
  "Central", 
  "Eastern", 
  "Western x Central hybrid",
  "Bonobo" # Now last
)

chimp_freq_bar <- df_chimp_unified %>%
  group_by(CommonName, hap_structure) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(CommonName) %>%
  mutate(
    frequency = count / sum(count),
    total_n = sum(count),
    subspecies_label_str = paste0(CommonName, " (n=", total_n, ")")
  ) %>%
  ungroup() %>%
  mutate(
    hap_structure = factor(hap_structure, levels = unified_structure_order),
    subspecies_label = factor(subspecies_label_str, 
                              levels = unique(subspecies_label_str[order(match(CommonName, subspecies_order))]))
  )

#final plot
chimp_bar_final_ordered <- ggplot(chimp_freq_bar, aes(x = subspecies_label, y = frequency, fill = hap_structure)) +
  geom_col(position = "stack", color = "white", linewidth = 0.1) +
  
  theme_classic() +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  
  scale_fill_manual(values = unified_colors, breaks = unified_structure_order, name = "GYP Structure") +
  
  labs(
    x = "Subspecies (Sample Size)",
    y = "Relative Frequency",
    fill = "Structure Type"
  ) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.title.x = element_text(margin = margin(t = 10)),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right",
    aspect.ratio = 0.7 
  )
