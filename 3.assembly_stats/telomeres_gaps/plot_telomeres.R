setwd("/Users/joanocha/Library/CloudStorage/GoogleDrive-joana.laranjeira.rocha@gmail.com/My Drive/POSTDOC/PANPAN/analysis/Figure1_Diversity_Assembly_stats/PANPAN_assembly_stats/telomeres_gaps")
library(dplyr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(tibble)

metadata <- fread("/Users/joanocha/Google Drive/My Drive/POSTDOC/PANPAN/metadata_datasets/PANPAN_HPRC_HGSVC_METADATA.txt")

telomere_data <- read.csv("telomere_summary_panpan.csv", 
                          header = TRUE, 
                          stringsAsFactors = FALSE)
telomere_data <- telomere_data %>%
  mutate(num_telomeres = teloBegin + teloEnd)
head(telomere_data)

telomere_all <- read.csv("telomere_all_panpan.csv")  %>%
  mutate(T_status = case_when(
    endLength > 0 & startLength > 0 ~ "2T",
    endLength == 0 & startLength == 0 ~ "Missing",
    TRUE ~ "1T"
  ))
head(telomere_all)
fwrite(telomere_data, "telomere_stats_data_allchroms.csv")

gap_data <- read.csv("gaps.csv") %>%
  mutate(gapless_status = case_when(
    gaps == 0 ~ "gapless",
    TRUE ~ "with gaps"
  ))



### TELOMERE STATS #####################

### GAPS STATS ###

merged_data <- gap_data %>%
  left_join(metadata, by = "id") %>%
  filter(SpeciesCode != "Hsa ") %>%
  # Filter out T2T reference samples
  filter(!id %in% c("mPanPan1.hap1", "mPanPan1.hap2", 
                    "mPanTro3.hap1", "mPanTro3.hap2")) %>%
  # Create species grouping with Chimpanzee first
  mutate(
    Species_group = ifelse(grepl("Ppan", SpeciesCode), "Bonobo (n=5, 10 haplotypes)", "Chimpanzee (n=24, 48 haplotypes)"),
    Species_group = factor(Species_group, levels = c("Chimpanzee (n=24, 48 haplotypes)", "Bonobo (n=5, 10 haplotypes)"))  # Changed order here
  )

ggplot(merged_data, aes(x = gaps)) +
  geom_histogram(aes(fill = species), 
                 position = "identity", 
                 alpha = 0.7,
                 bins = 30,
                 color = "black") +  # Add black borders to bars
  facet_wrap(~species, ncol = 1, scales = "free") +  # Changed to scales = "free" for free x and y
  scale_fill_manual(values = c("bonobo" = "#e9d7cb", 
                               "chimpanzee" = "darkblue")) +
  scale_x_continuous(breaks = seq(0, max(merged_data$gaps), by = 5)) +
  theme_classic() +
  labs(x = "Number of gaps per chromosome (reference-scaffold) per haplotype",
       y = "Count",
       fill = "Species") +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Summary statistics
gap_stats <- merged_data %>%
  group_by(species) %>%
  summarise(
    mean_gaps = mean(gaps),
    median_gaps = median(gaps),
    max_gaps = max(gaps),
    min_gaps = min(gaps),
    sd_gaps = sd(gaps)
  )
print(gap_stats)


# Print summary statistics
summary_stats <- merged_data %>%
  group_by(Species_group, gapless_status) %>%
  summarise(count = n()) %>%
  spread(gapless_status, count)
print(summary_stats)

### TELOMERE + GAPS 
telomere_gaps_all <- left_join(gap_data, telomere_all,  by = c("id", "chrom"))

telomere_gaps_all <- telomere_gaps_all %>%
  mutate(is_T2T = (gapless_status == "gapless") & (T_status == "2T"))
head(telomere_gaps_all)

#t2t_counts <- table(telomere_gaps_all$species[telomere_gaps_all$is_T2T == TRUE])
#t2t_summary <- telomere_gaps_all %>%
#  filter(is_T2T == TRUE) %>%
#  group_by(species) %>%
#  summarise(t2t_count = n())

telomere_gaps_all <- telomere_gaps_all %>%
  mutate(T2T_status = case_when(
    gapless_status == "gapless" & T_status == "2T" ~ "T2T",
    TRUE ~ paste(gapless_status, T_status)
  ))
head(telomere_gaps_all)
merged_data <- telomere_gaps_all  %>%
  left_join(metadata, by = "id") %>%
  filter(SpeciesCode != "Hsa ") %>%
  # Filter out T2T reference samples
  filter(!id %in% c("mPanPan1.hap1", "mPanPan1.hap2", 
                    "mPanTro3.hap1", "mPanTro3.hap2")) %>%
  # Create species grouping with Chimpanzee first
  mutate(
    Species_group = ifelse(grepl("Ppan", SpeciesCode), "Bonobo (n=5, 10 haplotypes)", "Chimpanzee (n=24, 48 haplotypes)"),
    Species_group = factor(Species_group, levels = c("Chimpanzee (n=24, 48 haplotypes)", "Bonobo (n=5, 10 haplotypes)"))  # Changed order here
  )

merged_data_filtered <- telomere_gaps_all %>%
  left_join(metadata, by = "id") %>%
  
  # 1. Filter out human 
  filter(!grepl("Hsa", SpeciesCode)) %>%
  
  # 2. Filter out sex chromosomes using the 'chrom' column and your specific prefixes
  filter(!grepl("chrX_|chrY_", chrom)) %>%
  
  # Filter out T2T reference samples
  filter(!id %in% c("mPanPan1.hap1", "mPanPan1.hap2", 
                    "mPanTro3.hap1", "mPanTro3.hap2")) %>%
  
  # Create species grouping with Chimpanzee first
  mutate(
    Species_group = ifelse(grepl("Ppan", SpeciesCode), "Bonobo (n=5, 10 haplotypes)", "Chimpanzee (n=24, 48 haplotypes)"),
    Species_group = factor(Species_group, levels = c("Chimpanzee (n=24, 48 haplotypes)", "Bonobo (n=5, 10 haplotypes)")) 
  )

fwrite(merged_data_filtered, "telomere_gaps_autosomes.csv")

#t2t_counts <- table(merged_data_filtered$species[merged_data_filtered$is_T2T == TRUE])
#t2t_counts

T2T_only_filtered <- telomere_gaps_all %>%
  left_join(metadata, by = "id") %>%
  
  # 1. Filter out human 
  filter(!grepl("Hsa", SpeciesCode)) %>%
  
  # 2. Filter out sex chromosomes using the 'chrom' column
  filter(!grepl("chrX_|chrY_", chrom)) %>%
  
  # 3. Filter out T2T reference samples
  filter(!id %in% c("mPanPan1.hap1", "mPanPan1.hap2", 
                    "mPanTro3.hap1", "mPanTro3.hap2")) %>%
  
  # 4. KEEP ONLY T2T TRUE
  filter(is_T2T == TRUE) %>%
  
  filter(Assembler == 'verkko') %>%
  
  # Create species grouping with Chimpanzee first
  mutate(
    Species_group = ifelse(grepl("Ppan", SpeciesCode), "Bonobo (n=5, 10 haplotypes)", "Chimpanzee (n=24, 48 haplotypes)"),
    Species_group = factor(Species_group, levels = c("Chimpanzee (n=24, 48 haplotypes)", "Bonobo (n=5, 10 haplotypes)")) 
  )
t2t_counts <- T2T_only_filtered %>%
  count(id, name = "T2T_count") %>%
  arrange(desc(T2T_count)) # Sorts from highest to lowest count

t2t_counts
fwrite(T2T_only_filtered, 'T2T_only_filtered.csv')
fwrite(t2t_counts, 't2t_counts.csv')



# 1. Prepare Data
telomere_plot_data <- telomere_data %>%
  left_join(metadata, by = "id") %>%
  mutate(
    # Order populations logically
    p1 = factor(p1, levels = c("Bonobo", "Western chimpanzee", "Western x Central hybrid", 
                               "Central chimpanzee", "Eastern chimpanzee")),
    dataset = ifelse(id %in% c("mPanPan1.hap1", "mPanPan1.hap2", "mPanTro3.hap1", "mPanTro3.hap2"), 
                     "T2T reference", "PANPAN")
  )

# 2. Calculate summary stats per population (for the black segments)
telo_stats <- telomere_plot_data %>%
  group_by(p1) %>%
  summarise(
    median_val = median(num_telomeres),
    q1 = quantile(num_telomeres, 0.25),
    q3 = quantile(num_telomeres, 0.75)
  )

# 3. Plot: Number of Telomeres
ggplot(telomere_plot_data, aes(x = p1, y = num_telomeres)) +
  # Single violin per population
  geom_violin(fill = "white", color = "gray80") +
  # Custom segments for quartiles/median
  geom_segment(data = telo_stats,
               aes(x = as.numeric(p1) - 0.2, xend = as.numeric(p1) + 0.2,
                   y = q1, yend = q1), linetype = "dotted", color = "black") +
  geom_segment(data = telo_stats,
               aes(x = as.numeric(p1) - 0.2, xend = as.numeric(p1) + 0.2,
                   y = median_val, yend = median_val), size = 1, color = "black") +
  geom_segment(data = telo_stats,
               aes(x = as.numeric(p1) - 0.2, xend = as.numeric(p1) + 0.2,
                   y = q3, yend = q3), linetype = "dotted", color = "black") +
  # Points colored by population and shaped by dataset
  geom_jitter(aes(color = p1, shape = dataset), width = 0.15, size = 2.5, alpha = 0.8) +
  # Aesthetics
  scale_color_manual(values = deframe(color_mapping)) +
  scale_shape_manual(values = c("T2T reference" = 17, "PANPAN" = 16)) +
  theme_classic() +
  labs(x = "Population", y = "# Telomeres per haplotype", 
       color = "Population", shape = "Dataset") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  coord_cartesian(ylim = c(0, 50))


# 1. Calculate summary stats for T2T
t2t_stats <- telomere_plot_data %>%
  group_by(p1) %>%
  summarise(
    median_val = median(T2T),
    q1 = quantile(T2T, 0.25),
    q3 = quantile(T2T, 0.75)
  )





# --- 1. Data Preparation ---
# Define T2T reference IDs for highlighting
refs <- c("mPanPan1.hap1", "mPanPan1.hap2", "mPanTro3.hap1", "mPanTro3.hap2")

# Create the population-to-color mapping used in your previous plots
color_mapping <- metadata %>%
  select(p1, HPRC_PANPAN.color) %>%
  distinct() %>%
  filter(!is.na(p1))

# 2. PER-CHROMOSOME ASSEMBLY STATUS DATA
# Filtering for Autosomes and extracting native 1-23 indexing
status_data <- read.csv("gaps.csv") %>%
  inner_join(read.csv("telomere_all_panpan.csv"), by = c("id", "chrom")) %>%
  left_join(metadata, by = "id") %>%
  # Exclude X/Y, Human, and References
  filter(!grepl("X|Y", chrom),
         !SpeciesCode %in% c("Hsa", "Hsa "),
         !id %in% refs) %>%
  mutate(
    # Native Extraction: chr1_hap1_... -> 1, chr10_mat_... -> 10
    chr_label = gsub("^chr([0-9]+).*", "\\1", chrom),
    chr_label = factor(chr_label, levels = as.character(1:23)),
    
    # 6-Level Status Definition
    T_status = case_when(
      startLength > 0 & endLength > 0 ~ "2T",
      startLength == 0 & endLength == 0 ~ "0T",
      TRUE ~ "1T"
    ),
    Assembly_Status = factor(case_when(
      gaps == 0 & T_status == "2T" ~ "T2T",
      gaps == 0 ~ paste("gapless", T_status),
      TRUE ~ paste("with gaps", T_status)
    ), levels = c("T2T", "gapless 1T", "gapless 0T", "with gaps 2T", "with gaps 1T", "with gaps 0T")),
    
    Species = factor(ifelse(grepl("Ppan", SpeciesCode), "Bonobo", "Chimpanzee"), 
                     levels = c("Bonobo", "Chimpanzee"))
  )

# 3. HAPLOTYPE SUMMARY DATA (Autosomes Only)
summary_data <- read.csv("telomere_summary_panpan_autosomes.csv") %>%
  left_join(metadata, by = "id") %>%
  filter(!SpeciesCode %in% c("Hsa", "Hsa "), !is.na(SpeciesCode)) %>%
  mutate(
    num_telomeres = teloBegin + teloEnd,
    dataset = ifelse(id %in% refs, "T2T reference", "PANPAN"),
    Species = factor(ifelse(grepl("Ppan", SpeciesCode), "Bonobo", "Chimpanzee"), 
                     levels = c("Bonobo", "Chimpanzee")),
    p1 = factor(p1, levels = color_mapping$p1)
  )

# Calculate Stats for the segments
get_summary_stats <- function(df, col_name) {
  df %>%
    group_by(Species) %>%
    summarise(
      median_val = median(get(col_name)),
      q1 = quantile(get(col_name), 0.25),
      q3 = quantile(get(col_name), 0.75)
    )
}

telo_stats <- get_summary_stats(summary_data, "num_telomeres")
t2t_stats <- get_summary_stats(summary_data, "T2T")


# 2. STATUS PLOT (Stacked Bar)
library(dplyr)
library(ggplot2)

# --- 1. Fix the Data Prep (Using gap_data and telomere_all) ---
status_data <- gap_data %>%
  inner_join(telomere_all, by = c("id", "chrom")) %>%
  left_join(metadata, by = "id") %>%
  # Exclude X/Y, Human, and T2T References
  filter(!grepl("X|Y", chrom), 
         !SpeciesCode %in% c("Hsa", "Hsa "), 
         !id %in% refs) %>%
  mutate(
    # Native Extraction (chr1_hap1 -> 1)
    chr_label = gsub("^chr([0-9]+).*", "\\1", chrom),
    chr_label = factor(chr_label, levels = as.character(1:23)),
    
    # 6-Level Status Definition
    T_status = case_when(
      startLength > 0 & endLength > 0 ~ "2T",
      startLength == 0 & endLength == 0 ~ "0T",
      TRUE ~ "1T"
    ),
    Assembly_Status = factor(case_when(
      gaps == 0 & T_status == "2T" ~ "T2T",
      gaps == 0 ~ paste("gapless", T_status),
      TRUE ~ paste("with gaps", T_status)
    ), levels = c("T2T", "gapless 1T", "gapless 0T", "with gaps 2T", "with gaps 1T", "with gaps 0T")),
    
    # Species Order: Chimpanzee First (Top), Bonobo Second (Bottom)
    Species = factor(ifelse(grepl("Ppan", SpeciesCode), "Bonobo", "Chimpanzee"), 
                     levels = c("Chimpanzee", "Bonobo"))
  )

# --- 2. Create Dynamic Labels for the Facets ---
species_labels <- status_data %>%
  group_by(Species) %>%
  summarise(n = n_distinct(id)) %>%
  mutate(label = paste0(Species, " (n=", n, " haplotypes)"))

label_map <- setNames(species_labels$label, species_labels$Species)

status_data <- status_data %>%
  mutate(
    Species_labeled = factor(label_map[as.character(Species)], 
                             levels = label_map[c("Chimpanzee", "Bonobo")])
  )

# --- 3. Plot 1: All Assemblies ---
status_pal <- c("T2T"="#495d23", "gapless 1T"="#6a8f23", "gapless 0T"="#9ab973", 
                "with gaps 2T"="#8a4513", "with gaps 1T"="#cd8440", "with gaps 0T"="lightgrey")

stacked_all <- ggplot(status_data, aes(x = chr_label, fill = Assembly_Status)) +
  geom_bar(position = "stack") +
  facet_wrap(~Species_labeled, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = status_pal) +
  labs(x = "Chromosomes", 
       y = "# of haplotypes (All)", 
       fill = "Assembly status") +
  theme_classic() +
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 12, face = "bold"))

print(stacked_all)
ggsave("stacked.pdf", stacked_all, width = 6, height = 5, units = "in", dpi = 300)


species_labels <- status_data %>%
  group_by(Species) %>%
  summarise(n = n_distinct(id)) %>%
  mutate(label = paste0(Species, " (n=", n, " haplotypes)"))


# --- 1. Filter for Verkko (HiFi+HiC+ONT) ---
# Ensure your metadata actually has a column named exactly 'method'. 
# If it's called something else, replace 'method' below.
status_data_verkko <- status_data %>%
  filter(Assembler == "verkko")

# --- 2. Recalculate dynamic labels for the smaller subset ---
species_labels_verkko <- status_data_verkko %>%
  group_by(Species) %>%
  summarise(n = n_distinct(id)) %>%
  mutate(label = paste0(Species, " (n=", n, " haplotypes)"))

label_map_verkko <- setNames(species_labels_verkko$label, species_labels_verkko$Species)

status_data_verkko <- status_data_verkko %>%
  mutate(
    Species_labeled_verkko = factor(label_map_verkko[as.character(Species)], 
                                    levels = label_map_verkko[c("Chimpanzee", "Bonobo")])
  )

# --- 3. Plot 2: Verkko Only ---
stacked_verkko <- ggplot(status_data_verkko, aes(x = chr_label, fill = Assembly_Status)) +
  geom_bar(position = "stack") +
  # Use the newly generated Verkko-specific labels
  #facet_wrap(~Species_labeled_verkko, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = status_pal) +
  labs(x = "Chromosomes", 
       y = "# T2T chromosomes (HiFi+HiC+ONT)", 
       fill = "Assembly status") +
  theme_classic() +
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 12, face = "bold"))

print(stacked_verkko)
ggsave("stacked_verkko.pdf", stacked_verkko, width = 6, height = 3, units = "in", dpi = 300)




# PLOT 2: Number of Telomeres (Species comparison, colored by Population)
nTelos_per_species_autosomes<-ggplot(summary_data, aes(x = Species, y = num_telomeres)) +
  geom_violin(fill = "white", color = "gray85") +
  # Median and Quartile lines
  geom_segment(data = telo_stats, aes(x = as.numeric(Species)-0.2, xend = as.numeric(Species)+0.2, 
                                      y = q1, yend = q1), linetype = "dotted") +
  geom_segment(data = telo_stats, aes(x = as.numeric(Species)-0.2, xend = as.numeric(Species)+0.2, 
                                      y = median_val, yend = median_val), size = 1) +
  geom_segment(data = telo_stats, aes(x = as.numeric(Species)-0.2, xend = as.numeric(Species)+0.2, 
                                      y = q3, yend = q3), linetype = "dotted") +
  geom_jitter(aes(color = p1, shape = dataset), width = 0.1, size = 2, alpha = 0.8) +
  scale_color_manual(values = deframe(color_mapping)) +
  scale_shape_manual(values = c("T2T reference" = 17, "PANPAN" = 16)) +
  labs(y = "# Telomeres per haplotype (Autosomes)", x = "", color = "Population", shape = "Dataset") +
  theme_classic() +
  coord_cartesian(ylim = c(0, 50))
ggsave("nTeloneres_per_species_autosomes.pdf", nTelos_per_species_autosomes, width = 4, height = 3, units = "in", dpi = 300)

# PLOT 3: Number of T2T Chromosomes (Species comparison, colored by Population)
nT2T_per_species_autosomes<-ggplot(summary_data, aes(x = Species, y = T2T)) +
  geom_violin(fill = "white", color = "gray85") +
  # Median and Quartile lines
  geom_segment(data = t2t_stats, aes(x = as.numeric(Species)-0.2, xend = as.numeric(Species)+0.2, 
                                     y = q1, yend = q1), linetype = "dotted") +
  geom_segment(data = t2t_stats, aes(x = as.numeric(Species)-0.2, xend = as.numeric(Species)+0.2, 
                                     y = median_val, yend = median_val), size = 1) +
  geom_segment(data = t2t_stats, aes(x = as.numeric(Species)-0.2, xend = as.numeric(Species)+0.2, 
                                     y = q3, yend = q3), linetype = "dotted") +
  geom_jitter(aes(color = p1, shape = dataset), width = 0.1, size = 2, alpha = 0.8) +
  scale_color_manual(values = deframe(color_mapping)) +
  scale_shape_manual(values = c("T2T reference" = 17, "PANPAN" = 16)) +
  labs(y = "# T2T chromosomes per haplotype (Autosomes)", x = "", color = "Population", shape = "Dataset") +
  theme_classic() +
  coord_cartesian(ylim = c(0, 25))
#ggsave("nT2T_per_species_autosomes.pdf", nT2T_per_species_autosomes, width = 4, height = 3, units = "in", dpi = 300)