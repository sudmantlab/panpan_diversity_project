setwd("/Users/joanocha/Library/CloudStorage/GoogleDrive-joana.laranjeira.rocha@gmail.com/My Drive/POSTDOC/PANPAN/analysis/Figure1_Diversity_Assembly_stats/PANPAN_assembly_stats/telomeres_gaps")
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(tibble)

metadata <- fread("/Users/joanocha/Google Drive/My Drive/POSTDOC/PANPAN/metadata_datasets/PANPAN_HPRC_HGSVC_METADATA.txt")

out_dir <- "new_outputs"
dir.create(out_dir, showWarnings = FALSE)

refs <- c("mPanPan1.hap1", "mPanPan1.hap2", "mPanTro3.hap1", "mPanTro3.hap2")

# ---- Build old-style inputs from new_inputs ----

# 1. Per-chrom telomeres  (equivalent to telomere_all_panpan.csv)
telo_pan <- fread("new_inputs/ragtag_mPanPan1_all.csv")
telo_tro <- fread("new_inputs/ragtag_mPanTro3_all.csv")
telomere_all <- rbind(
  telo_pan %>% mutate(species = "bonobo"),
  telo_tro %>% mutate(species = "chimpanzee")
) %>%
  rename(id = sample) %>%
  # Defensively drop sex chromosomes
  filter(!grepl("chrX|chrY|_hsaX|_hsaY", chrom)) %>%
  mutate(T_status = case_when(
    endLength > 0 & startLength > 0 ~ "2T",
    endLength == 0 & startLength == 0 ~ "Missing",
    TRUE ~ "1T"
  ))

# 2. Per-chrom gaps (equivalent to gaps.csv)
gaps_pan <- fread("new_inputs/ragtag_mPanPan1_all~gaps.csv") %>%
  mutate(species = "bonobo")
gaps_tro <- fread("new_inputs/ragtag_mPanTro3_all~gaps.csv") %>%
  mutate(species = "chimpanzee")
gap_data <- rbind(gaps_pan, gaps_tro) %>%
  # Strip trailing ".chr" from sample column
  mutate(id = sub("\\.chr$", "", sample)) %>%
  rename(chrom = sequence, gaps = count) %>%
  select(species, id, chrom, gaps) %>%
  filter(!grepl("chrX|chrY|_hsaX|_hsaY", chrom)) %>%
  mutate(gapless_status = case_when(
    gaps == 0 ~ "gapless",
    TRUE ~ "with gaps"
  ))

# 3. Per-haplotype summary (equivalent to telomere_summary_panpan.csv)
sum_pan <- fread("new_inputs/ragtag_mPanPan1_all~summary.csv")
sum_tro <- fread("new_inputs/ragtag_mPanTro3_all~summary.csv")
telomere_data <- rbind(sum_pan, sum_tro) %>%
  rename(id = sample) %>%
  mutate(num_telomeres = teloBegin + teloEnd)

# 4. Autosomes summary — recompute from per-chrom telomere file (autosomes only)
summary_autosomes <- telomere_all %>%
  group_by(id) %>%
  summarise(
    teloBegin = sum(startLength > 0),
    teloEnd = sum(endLength > 0),
    T2T = sum(startLength > 0 & endLength > 0)
  ) %>%
  ungroup()

# Sanity: mirror old behaviour using provided summary for num_telomeres when possible.
# The new ~summary.csv treats teloBegin/teloEnd as COUNTS (see PR00251: 15/13/12), so
# using the same convention for autosomes is consistent.

fwrite(telomere_data, file.path(out_dir, "telomere_stats_data_allchroms.csv"))

### GAPS STATS ###
merged_data_gaps <- gap_data %>%
  left_join(metadata, by = "id") %>%
  filter(SpeciesCode != "Hsa ") %>%
  filter(!id %in% refs) %>%
  mutate(
    Species_group = ifelse(grepl("Ppan", SpeciesCode),
                           "Bonobo (n=5, 10 haplotypes)",
                           "Chimpanzee (n=24, 48 haplotypes)"),
    Species_group = factor(Species_group,
                           levels = c("Chimpanzee (n=24, 48 haplotypes)",
                                      "Bonobo (n=5, 10 haplotypes)"))
  )

### TELOMERE + GAPS merge ###
telomere_gaps_all <- left_join(gap_data, telomere_all, by = c("id", "chrom", "species")) %>%
  mutate(is_T2T = (gapless_status == "gapless") & (T_status == "2T")) %>%
  mutate(T2T_status = case_when(
    gapless_status == "gapless" & T_status == "2T" ~ "T2T",
    TRUE ~ paste(gapless_status, T_status)
  ))

merged_data_filtered <- telomere_gaps_all %>%
  left_join(metadata, by = "id") %>%
  filter(!grepl("Hsa", SpeciesCode)) %>%
  filter(!grepl("chrX_|chrY_|_hsaX|_hsaY", chrom)) %>%
  filter(!id %in% refs) %>%
  mutate(
    Species_group = ifelse(grepl("Ppan", SpeciesCode),
                           "Bonobo (n=5, 10 haplotypes)",
                           "Chimpanzee (n=24, 48 haplotypes)"),
    Species_group = factor(Species_group,
                           levels = c("Chimpanzee (n=24, 48 haplotypes)",
                                      "Bonobo (n=5, 10 haplotypes)"))
  )

fwrite(merged_data_filtered, file.path(out_dir, "telomere_gaps_autosomes.csv"))

T2T_only_filtered <- telomere_gaps_all %>%
  left_join(metadata, by = "id") %>%
  filter(!grepl("Hsa", SpeciesCode)) %>%
  filter(!grepl("chrX_|chrY_|_hsaX|_hsaY", chrom)) %>%
  filter(!id %in% refs) %>%
  filter(is_T2T == TRUE) %>%
  filter(Assembler == 'verkko') %>%
  mutate(
    Species_group = ifelse(grepl("Ppan", SpeciesCode),
                           "Bonobo (n=5, 10 haplotypes)",
                           "Chimpanzee (n=24, 48 haplotypes)"),
    Species_group = factor(Species_group,
                           levels = c("Chimpanzee (n=24, 48 haplotypes)",
                                      "Bonobo (n=5, 10 haplotypes)"))
  )
t2t_counts <- T2T_only_filtered %>%
  count(id, name = "T2T_count") %>%
  arrange(desc(T2T_count))

fwrite(T2T_only_filtered, file.path(out_dir, "T2T_only_filtered.csv"))
fwrite(t2t_counts, file.path(out_dir, "t2t_counts.csv"))

# Color mapping for populations
color_mapping <- metadata %>%
  select(p1, HPRC_PANPAN.color) %>%
  distinct() %>%
  filter(!is.na(p1))

# Haplotype summary data (autosomes) for violins
summary_data <- summary_autosomes %>%
  left_join(metadata, by = "id") %>%
  filter(!SpeciesCode %in% c("Hsa", "Hsa "), !is.na(SpeciesCode)) %>%
  mutate(
    num_telomeres = teloBegin + teloEnd,
    dataset = ifelse(id %in% refs, "T2T reference", "PANPAN"),
    Species = factor(ifelse(grepl("Ppan", SpeciesCode), "Bonobo", "Chimpanzee"),
                     levels = c("Bonobo", "Chimpanzee")),
    p1 = factor(p1, levels = color_mapping$p1)
  )

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

# --- STACKED BAR: assembly status per chromosome ---
status_data <- gap_data %>%
  inner_join(telomere_all, by = c("id", "chrom", "species")) %>%
  left_join(metadata, by = "id") %>%
  filter(!grepl("X|Y", chrom),
         !SpeciesCode %in% c("Hsa", "Hsa "),
         !id %in% refs) %>%
  mutate(
    chr_label = gsub("^chr([0-9]+).*", "\\1", chrom),
    chr_label = factor(chr_label, levels = as.character(1:23)),
    T_status = case_when(
      startLength > 0 & endLength > 0 ~ "2T",
      startLength == 0 & endLength == 0 ~ "0T",
      TRUE ~ "1T"
    ),
    Assembly_Status = factor(case_when(
      gaps == 0 & T_status == "2T" ~ "T2T",
      gaps == 0 ~ paste("gapless", T_status),
      TRUE ~ paste("with gaps", T_status)
    ), levels = c("T2T", "gapless 1T", "gapless 0T",
                  "with gaps 2T", "with gaps 1T", "with gaps 0T")),
    Species = factor(ifelse(grepl("Ppan", SpeciesCode), "Bonobo", "Chimpanzee"),
                     levels = c("Chimpanzee", "Bonobo"))
  )

species_labels <- status_data %>%
  group_by(Species) %>%
  summarise(n = n_distinct(id)) %>%
  mutate(label = paste0(Species, " (n=", n, " haplotypes)"))
label_map <- setNames(species_labels$label, species_labels$Species)
status_data <- status_data %>%
  mutate(Species_labeled = factor(label_map[as.character(Species)],
                                  levels = label_map[c("Chimpanzee", "Bonobo")]))

status_pal <- c("T2T"="#495d23", "gapless 1T"="#6a8f23", "gapless 0T"="#9ab973",
                "with gaps 2T"="#8a4513", "with gaps 1T"="#cd8440",
                "with gaps 0T"="lightgrey")

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
ggsave(file.path(out_dir, "stacked.pdf"), stacked_all,
       width = 9, height = 3.5, units = "in", dpi = 300)

# Verkko-only stacked bar
status_data_verkko <- status_data %>% filter(Assembler == "verkko")
if (nrow(status_data_verkko) > 0) {
  species_labels_verkko <- status_data_verkko %>%
    group_by(Species) %>%
    summarise(n = n_distinct(id)) %>%
    mutate(label = paste0(Species, " (n=", n, " haplotypes)"))
  label_map_verkko <- setNames(species_labels_verkko$label, species_labels_verkko$Species)
  status_data_verkko <- status_data_verkko %>%
    mutate(Species_labeled_verkko = factor(
      label_map_verkko[as.character(Species)],
      levels = label_map_verkko[intersect(c("Chimpanzee", "Bonobo"), names(label_map_verkko))]))

  n_verkko_hap <- n_distinct(status_data_verkko$id)
  stacked_verkko <- ggplot(status_data_verkko, aes(x = chr_label, fill = Assembly_Status)) +
    geom_bar(position = "stack") +
    scale_fill_manual(values = status_pal) +
    labs(x = "Chromosomes",
         y = "# of haplotypes (HiFi+HiC+ONT)",
         fill = "Assembly status",
         title = paste0("Verkko (n=", n_verkko_hap, " haplotypes)")) +
    theme_classic() +
    theme(plot.title = element_text(size = 12, face = "bold"))
  ggsave(file.path(out_dir, "stacked_verkko.pdf"), stacked_verkko,
         width = 6, height = 3, units = "in", dpi = 300)
}

# --- Violin plots (autosomes) ---
nTelos_per_species_autosomes <- ggplot(summary_data, aes(x = Species, y = num_telomeres)) +
  geom_violin(fill = "white", color = "gray85") +
  geom_segment(data = telo_stats, aes(x = as.numeric(Species)-0.2, xend = as.numeric(Species)+0.2,
                                      y = q1, yend = q1), linetype = "dotted") +
  geom_segment(data = telo_stats, aes(x = as.numeric(Species)-0.2, xend = as.numeric(Species)+0.2,
                                      y = median_val, yend = median_val), size = 1) +
  geom_segment(data = telo_stats, aes(x = as.numeric(Species)-0.2, xend = as.numeric(Species)+0.2,
                                      y = q3, yend = q3), linetype = "dotted") +
  geom_jitter(aes(color = p1, shape = dataset), width = 0.1, size = 2, alpha = 0.8) +
  scale_color_manual(values = deframe(color_mapping)) +
  scale_shape_manual(values = c("T2T reference" = 17, "PANPAN" = 16)) +
  labs(y = "# Telomeres per haplotype (Autosomes)", x = "",
       color = "Population", shape = "Dataset") +
  theme_classic() +
  coord_cartesian(ylim = c(0, 50))
ggsave(file.path(out_dir, "nTeloneres_per_species_autosomes.pdf"),
       nTelos_per_species_autosomes, width = 4, height = 3, units = "in", dpi = 300)

nT2T_per_species_autosomes <- ggplot(summary_data, aes(x = Species, y = T2T)) +
  geom_violin(fill = "white", color = "gray85") +
  geom_segment(data = t2t_stats, aes(x = as.numeric(Species)-0.2, xend = as.numeric(Species)+0.2,
                                     y = q1, yend = q1), linetype = "dotted") +
  geom_segment(data = t2t_stats, aes(x = as.numeric(Species)-0.2, xend = as.numeric(Species)+0.2,
                                     y = median_val, yend = median_val), size = 1) +
  geom_segment(data = t2t_stats, aes(x = as.numeric(Species)-0.2, xend = as.numeric(Species)+0.2,
                                     y = q3, yend = q3), linetype = "dotted") +
  geom_jitter(aes(color = p1, shape = dataset), width = 0.1, size = 2, alpha = 0.8) +
  scale_color_manual(values = deframe(color_mapping)) +
  scale_shape_manual(values = c("T2T reference" = 17, "PANPAN" = 16)) +
  labs(y = "# T2T chromosomes per haplotype (Autosomes)", x = "",
       color = "Population", shape = "Dataset") +
  theme_classic() +
  coord_cartesian(ylim = c(0, 25))
ggsave(file.path(out_dir, "nT2T_per_species_autosomes.pdf"),
       nT2T_per_species_autosomes, width = 4, height = 3, units = "in", dpi = 300)

# --- All-chroms versions (same as autosomes here since no X/Y in new inputs) ---
summary_all_data <- telomere_data %>%
  left_join(metadata, by = "id") %>%
  filter(!SpeciesCode %in% c("Hsa", "Hsa "), !is.na(SpeciesCode)) %>%
  mutate(
    dataset = ifelse(id %in% refs, "T2T reference", "PANPAN"),
    Species = factor(ifelse(grepl("Ppan", SpeciesCode), "Bonobo", "Chimpanzee"),
                     levels = c("Bonobo", "Chimpanzee")),
    p1 = factor(p1, levels = color_mapping$p1)
  )
telo_stats_all <- get_summary_stats(summary_all_data, "num_telomeres")
t2t_stats_all  <- get_summary_stats(summary_all_data, "T2T")

nTelos_per_species <- ggplot(summary_all_data, aes(x = Species, y = num_telomeres)) +
  geom_violin(fill = "white", color = "gray85") +
  geom_segment(data = telo_stats_all, aes(x = as.numeric(Species)-0.2, xend = as.numeric(Species)+0.2,
                                          y = q1, yend = q1), linetype = "dotted") +
  geom_segment(data = telo_stats_all, aes(x = as.numeric(Species)-0.2, xend = as.numeric(Species)+0.2,
                                          y = median_val, yend = median_val), size = 1) +
  geom_segment(data = telo_stats_all, aes(x = as.numeric(Species)-0.2, xend = as.numeric(Species)+0.2,
                                          y = q3, yend = q3), linetype = "dotted") +
  geom_jitter(aes(color = p1, shape = dataset), width = 0.1, size = 2, alpha = 0.8) +
  scale_color_manual(values = deframe(color_mapping)) +
  scale_shape_manual(values = c("T2T reference" = 17, "PANPAN" = 16)) +
  labs(y = "# Telomeres per haplotype", x = "",
       color = "Population", shape = "Dataset") +
  theme_classic() +
  coord_cartesian(ylim = c(0, 50))
ggsave(file.path(out_dir, "nTelomeres_per_species.pdf"),
       nTelos_per_species, width = 4, height = 3, units = "in", dpi = 300)

nT2T_per_species <- ggplot(summary_all_data, aes(x = Species, y = T2T)) +
  geom_violin(fill = "white", color = "gray85") +
  geom_segment(data = t2t_stats_all, aes(x = as.numeric(Species)-0.2, xend = as.numeric(Species)+0.2,
                                         y = q1, yend = q1), linetype = "dotted") +
  geom_segment(data = t2t_stats_all, aes(x = as.numeric(Species)-0.2, xend = as.numeric(Species)+0.2,
                                         y = median_val, yend = median_val), size = 1) +
  geom_segment(data = t2t_stats_all, aes(x = as.numeric(Species)-0.2, xend = as.numeric(Species)+0.2,
                                         y = q3, yend = q3), linetype = "dotted") +
  geom_jitter(aes(color = p1, shape = dataset), width = 0.1, size = 2, alpha = 0.8) +
  scale_color_manual(values = deframe(color_mapping)) +
  scale_shape_manual(values = c("T2T reference" = 17, "PANPAN" = 16)) +
  labs(y = "# T2T chromosomes per haplotype", x = "",
       color = "Population", shape = "Dataset") +
  theme_classic() +
  coord_cartesian(ylim = c(0, 25))
ggsave(file.path(out_dir, "nT2T_per_species.pdf"),
       nT2T_per_species, width = 4, height = 3, units = "in", dpi = 300)

# Write the rebuilt canonical CSVs into new_outputs
fwrite(telomere_data, file.path(out_dir, "telomere_summary_panpan.csv"))
fwrite(summary_autosomes, file.path(out_dir, "telomere_summary_panpan_autosomes.csv"))
fwrite(gap_data, file.path(out_dir, "gaps.csv"))
fwrite(telomere_all %>% select(chrom, endLength, startLength, id),
       file.path(out_dir, "telomere_all_panpan.csv"))

cat("Done. Outputs written to:", out_dir, "\n")
