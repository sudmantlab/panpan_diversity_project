library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(data.table)


df <- read.csv("compleasm/compleasm.csv", stringsAsFactors = FALSE)
sc <- df %>%
  filter(type == "S") %>%
  mutate(
    sample_id = str_extract(sample, "^[^.]+"),
    hap = case_when(
      str_detect(sample, "\\.hap1\\.") ~ "H1",
      str_detect(sample, "\\.hap2\\.") ~ "H2",
      TRUE ~ NA_character_
    ),
    PC_pct = PC * 100
  ) %>%
  filter(!is.na(hap))

#Prioritize assembler per sample: Verkko > Hifiasm+HiC > Hifiasm ---
dir_priority <- c(
  "Verkko-fasta_shortcut"      = 1,
  "Hifiasm-fasta-HiC_shortcut" = 2,
  "Hifiasm-fasta_shortcut"     = 3
)

sc <- sc %>%
  mutate(prio = dir_priority[directory]) %>%
  group_by(sample_id) %>%
  filter(prio == min(prio)) %>%
  ungroup()

# Pivot to wide: H1 vs H2 ---
wide <- sc %>%
  select(sample_id, directory, hap, PC_pct) %>%
  distinct() %>%
  group_by(sample_id, directory, hap) %>%
  summarise(PC_pct = mean(PC_pct), .groups = "drop") %>%
  pivot_wider(names_from = hap, values_from = PC_pct) %>%
  filter(!is.na(H1) & !is.na(H2)) %>%
  mutate(
    assembler = ifelse(directory == "Verkko-fasta_shortcut", "Verkko", "Hifiasm")
  )

df_meta <- fread("input_tables/PanPan_Metadata.tsv")
meta_unique <- df_meta %>%
  select(genomeName, CommonName, colors, sex) %>%
  distinct()

wide <- wide %>%
  left_join(meta_unique, by = c("sample_id" = "genomeName")) %>%
  filter(!is.na(CommonName), !is.na(sex))

# --- 5. Sex labels with counts ---
wide <- wide %>%
  group_by(sex) %>%
  mutate(sex_label = paste0(sex, " (n=", n(), ")")) %>%
  ungroup()

label_levels <- sort(unique(wide$sex_label))

# --- 6. Color mapping from metadata ---
color_mapping <- wide %>%
  select(CommonName, colors) %>%
  distinct() %>%
  filter(!is.na(CommonName))

named_colors <- setNames(color_mapping$colors, color_mapping$CommonName)

avg_verkko  <- mean(c(wide$H1[wide$assembler == "Verkko"],  wide$H2[wide$assembler == "Verkko"]),  na.rm = TRUE)
avg_hifiasm <- mean(c(wide$H1[wide$assembler == "Hifiasm"], wide$H2[wide$assembler == "Hifiasm"]), na.rm = TRUE)

label_verkko  <- paste0("Verkko: ",              round(avg_verkko,  2), "%")
label_hifiasm <- paste0("Hifiasm/Hifiasm-HiC: ", round(avg_hifiasm, 2), "%")

cat("Samples plotted:", nrow(wide), "\n")
print(wide)

# Shared plot base (points, scales, labels, theme) 
base_theme <- theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11, color = "black"),
    legend.background = element_rect(fill = "transparent", color = NA),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 13)
  )

base_layers <- list(
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "grey70", linewidth = 0.6),
  geom_point(aes(color = CommonName, shape = sex_label), size = 3.5, stroke = 1.2, alpha = 0.85),
  scale_shape_manual(values = setNames(c(4, 15), label_levels), name = "Sex"),
  scale_color_manual(values = named_colors, name = "Population"),
  xlab("H1 BUSCO single-copy genes (%)"),
  ylab("H2 BUSCO single-copy genes (%)"),
  base_theme
)

#  Extra version: colored dashed lines with numeric labels 
plot_compleasm_extra <- ggplot(data = wide, aes(x = H1, y = H2)) +
  # Verkko average (dark blue, dotdash) — single value applied to both axes
  geom_vline(xintercept = avg_verkko, linetype = "dotdash", color = "darkblue", linewidth = 0.7, alpha = 0.7) +
  geom_hline(yintercept = avg_verkko, linetype = "dotdash", color = "darkblue", linewidth = 0.7, alpha = 0.7) +
  annotate("text", x = -Inf, y = avg_verkko, label = label_verkko,
           vjust = -0.5, hjust = -0.05, color = "darkblue", size = 3.2, fontface = "italic") +

  geom_vline(xintercept = avg_hifiasm, linetype = "dashed", color = "darkred", linewidth = 0.7, alpha = 0.7) +
  geom_hline(yintercept = avg_hifiasm, linetype = "dashed", color = "darkred", linewidth = 0.7, alpha = 0.7) +
  annotate("text", x = -Inf, y = avg_hifiasm, label = label_hifiasm,
           vjust = -0.5, hjust = -0.05, color = "darkred", size = 3.2, fontface = "italic") +

  base_layers

print(plot_compleasm_extra)
ggsave("compleasm_H1_vs_H2_extra.pdf", plot_compleasm_extra, width = 6.5, height = 4, dpi = 300)
ggsave("compleasm_H1_vs_H2_extra.png", plot_compleasm_extra, width = 6.5, height = 4, dpi = 300)

#Simple version: single grey dashed H1 avg (vertical) and H2 avg (horizontal), no annotations ---
avg_h1 <- mean(wide$H1, na.rm = TRUE)
avg_h2 <- mean(wide$H2, na.rm = TRUE)

plot_compleasm_simple <- ggplot(data = wide, aes(x = H1, y = H2)) +
  geom_vline(xintercept = avg_h1, linetype = "dashed", color = "grey30", linewidth = 0.6, alpha = 0.8) +
  geom_hline(yintercept = avg_h2, linetype = "dashed", color = "grey30", linewidth = 0.6, alpha = 0.8) +
  base_layers

print(plot_compleasm_simple)
ggsave("compleasm_H1_vs_H2.pdf", plot_compleasm_simple, width = 6.5, height = 4, dpi = 300)
ggsave("compleasm_H1_vs_H2.png", plot_compleasm_simple, width = 6.5, height = 4, dpi = 300)

write.csv(wide, "compleasm_H1_H2_single_copy.csv", row.names = FALSE)

cat("Verkko avg (H1+H2 pooled):  ", round(avg_verkko, 3), "%\n")
cat("Hifiasm avg (H1+H2 pooled): ", round(avg_hifiasm, 3), "%\n")
