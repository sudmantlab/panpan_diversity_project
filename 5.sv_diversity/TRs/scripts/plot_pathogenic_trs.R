required_packages <- c("ggplot2", "dplyr", "ggrepel", "tidyr", "patchwork")

# Check and install missing packages
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

invisible(lapply(required_packages, install_if_missing))

library(ggplot2)
library(dplyr)
library(ggrepel)
library(tidyr)
library(patchwork)

df<-read.table("human_chimp_exp_lens.txt")

df <- df %>%
  rename(
    species = V1,
    TRID = V2,
    motif_seq = V3,
    length = V4,
    copy_number = V5,
    pathogenic = V6)

# Reshape to get human and chimp values side-by-side
df_lines <- df %>%
  select(species, TRID, copy_number) %>%
  pivot_wider(names_from = species, values_from = copy_number)

# Reorder TRID based on descending human (homo) copy number
df_lines <- df_lines %>%
  arrange(homo) %>%
  mutate(TRID = factor(TRID, levels = TRID))

# Apply same factor levels to the full dataset
df <- df %>%
  mutate(TRID = factor(TRID, levels = levels(df_lines$TRID)))

df_highlight <- df %>%
  filter(species == "homo") %>%
  mutate(TRID = factor(TRID, levels = levels(df_lines$TRID)))

p <- ggplot() +
  geom_segment(data = df_lines, aes(x = chimp, xend = homo, y = TRID, yend = TRID), color = "gray70") +
  geom_segment(data = df_highlight, aes(x = copy_number, xend = pathogenic, y = TRID, yend = TRID), linetype = "dashed", color = "red", linewidth = 0.5) +
  geom_point(data = df, aes(x = copy_number, y = TRID, color = species), size = 2, show.legend = FALSE) +
  geom_point(data = df_highlight, aes(x = pathogenic, y = TRID), shape = 8, color = "red", size = 2) +
  geom_text(data = df_highlight, aes(x = pathogenic, y = TRID, label = TRID), hjust = -0.30, size = 6/2.845, color = "black", family="Helvetica") +
  coord_cartesian(clip = "off") +
  scale_x_log10() +
  scale_color_manual(values = c("chimp" = "#8A865D", "homo" = "#5C608A")) +
  labs(x = "TR copy number", y = "Genes") +
theme(
  axis.line = element_line(linewidth = 1, colour = "black"),
  panel.grid.major.x = element_line(color = "gray95"), 
  panel.grid.minor.x = element_line(color = "gray95", linetype = "dotted"),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  plot.background = element_blank(),
  plot.margin = margin(5.5, 10, 5.5, 5.5),
  text = element_text(size = 8, family="Helvetica", color = "black"),
  axis.title.x = element_text(family="Helvetica", size = 8),
  axis.title.y = element_text(family="Helvetica", size = 8),
  axis.text.x = element_text(family="Helvetica", color = "black", size = 8),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  legend.position = "none")

ggsave("homo_pantro_pathogenic_tr_lens.pdf", plot = p, device = cairo_pdf, width = 6, height = 6)

# Plot allele length distribution for example genes

df <- read.table("human_chimp_exp_sumstats.txt")

df <- df %>%
  rename(
    TRID = V1,
    sample = V2,
    motif_len = V3,
    tr_len = V4,
    motif_seq = V5,
    copy_number = V6,
    spp = V7)

df <- df %>%
  filter(TRID %in% c("FXN", "HTT", "C9ORF72", "ATXN10"))

df_long <- df %>%
  group_by(TRID, spp, sample) %>%
  mutate(
    allele_num = row_number()) %>%
  ungroup()

df_order <- df_long %>%
  mutate(copy_number = as.numeric(copy_number)) %>%
  group_by(TRID) %>%
  arrange(
    spp == "homo", copy_number, .by_group = TRUE) %>%
  mutate(
    plot_id = row_number()) %>%
  ungroup() %>%
  mutate(
    spp = factor(spp, levels = c("homo", "chimp")),
    plot_id = factor(plot_id, levels = unique(plot_id)))

df_order <- df_order %>%
  mutate(TRID = factor(TRID, levels = c("FXN", "HTT", "C9ORF72", "ATXN10")))

annotations <- df_order %>%
  group_by(TRID) %>%
  summarise(y = max(as.numeric(plot_id))) %>%
  mutate(
    label = c(
      "Friedreich ataxia (GAA)70",
      "Huntington disease (CAG)36",
      "Amyotrophic lateral sclerosis (GGCCC)24",
      "Spinocerebellar ataxia (ATTCT)800"))

p <- ggplot(df_order,
            aes(y = plot_id, x = copy_number, fill = spp)) + 
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(data = annotations, aes(x = Inf, y = y, label = label), 
    inherit.aes = FALSE, hjust = 1, vjust = -0.5, family = "Helvetica", size = 3) +
  facet_wrap(~TRID, nrow = 2, ncol = 2, scales = "free") +
  scale_fill_manual(values = c("homo" = "#5C608A", "chimp" = "#8A865D")) +
  scale_x_continuous(name = "Copy number", expand = c(0, 0)) +
  scale_y_discrete(name = "Alleles", expand = c(0, 0)) +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 7, family = "Helvetica"),
    strip.placement = "outside", 
    text = element_text(family = "Helvetica", size = 6),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(size = 6),
    axis.title = element_text(size = 7))

ggsave("patch_pathogenic_trs.pdf", plot = p, device = cairo_pdf, width = 6, height = 6)
