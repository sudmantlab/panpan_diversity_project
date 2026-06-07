#!/usr/bin/env Rscript
# Generates all plots in final_plots/.
#
# Inputs:
#   combined_within_pops_pi.txt.gz                per-window pi for long + short
#   mPanTro3_genomefeatures.bed copy              Cen / Telo intervals
#   mPanTro3_CenSat.bed copy                      centromeric satellite intervals
#   mPanTro3_sedefSegDups.bed copy                segmental duplication intervals
#   pantro_catalog.no_overlaps_simp.bed copy      tandem-repeat catalogue (TR)
#
# Annotation categories considered (Manhattans + density):
#   Centromere (Cen)  -- from genomefeatures
#   Telomere   (Telo) -- from genomefeatures
#   CenSat            -- from CenSat.bed
#   SegDup            -- from sedefSegDups.bed
#   TR                -- from pantro_catalog.no_overlaps_simp.bed
#   "Other (incl. SRA)" -- everything else; this bucket includes the
#                          short-read mask (SR_mask.bed) plus single-copy
#                          euchromatin not flagged by any of the above.
#
# Outputs (final_plots/):
#   Manhattans:
#     manhattan_pi_ratio_capped.pdf
#     manhattan_pi_ratio_facet_chrom.pdf
#     manhattan_pi_ratio_facet_chrom_clean.pdf
#   Density:
#     density_log10_ratio_by_region.pdf
#     density_pi_by_dataset_region.pdf
#
# Filters:
#   - Drop Hybrid (admixed) and NC (short-only)
#   - Autosomes only (chrX/chrY excluded)
#   - For ratio: drop windows where both pi_long and pi_short are 0/NA;
#     where long has a value but short = 0/NA, replace short by epsilon
#     = 1 / median(count_comparisons_short)
#

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(ggplot2)
  library(scales)
  library(ggnewscale)
  library(ggrastr)
})

dir.create("final_plots", showWarnings = FALSE)

# --- Parameters --------------------------------------------------------
OUTLIER_PCT <- 0.01

# Manhattan colour/shape palette (5 named complex categories)
REGION_COLORS <- c("Centromere" = "#2ca02c", "Telomere" = "#d62728",
                   "CenSat"     = "#ff9933", "SegDup"   = "#9467bd",
                   "TR"         = "#17becf")
REGION_SHAPES <- c("Centromere" = 1, "Telomere" = 0,
                   "CenSat"     = 2, "SegDup"   = 5,
                   "TR"         = 6)


OTHER_LBL <- "Other (incl. SegDup, TR, SRA)"
REGION_COLORS_DENS <- c("Centromere" = unname(REGION_COLORS["Centromere"]),
                        "Telomere"   = unname(REGION_COLORS["Telomere"]),
                        "CenSat"     = unname(REGION_COLORS["CenSat"]),
                        structure("grey60", names = OTHER_LBL))
REGION_LEVELS_DENS <- c("Centromere", "Telomere", "CenSat", OTHER_LBL)

compact_theme <- theme_bw(base_size = 8) +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.spacing.x    = unit(0.05, "lines"),
        panel.spacing.y    = unit(0.1,  "lines"),
        plot.title         = element_text(size = 8, face = "bold"),
        plot.subtitle      = element_text(size = 6.5),
        plot.margin        = margin(2, 4, 2, 2),
        axis.title         = element_text(size = 7),
        axis.text          = element_text(size = 6),
        strip.background   = element_rect(fill = "grey95"),
        strip.text         = element_text(size = 7),
        strip.text.y       = element_text(angle = 0),
        legend.position    = "top",
        legend.box.spacing = unit(1, "pt"),
        legend.key.size    = unit(8, "pt"),
        legend.title       = element_text(size = 7),
        legend.text        = element_text(size = 6.5),
        legend.margin      = margin(0, 0, 0, 0))


panel_theme <- theme_bw(base_size = 7) +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major   = element_line(linewidth = 0.2, color = "grey92"),
        panel.spacing      = unit(0.05, "lines"),
        plot.title         = element_blank(),
        plot.subtitle      = element_blank(),
        plot.margin        = margin(1, 2, 1, 1),
        axis.title         = element_text(size = 6.5),
        axis.text          = element_text(size = 5.5),
        axis.ticks         = element_line(linewidth = 0.25),
        strip.background   = element_rect(fill = "grey95",
                                          color = NA),
        strip.text         = element_text(size = 6, margin = margin(1,1,1,1)),
        strip.text.y       = element_text(angle = 0),
        legend.position    = "top",
        legend.direction   = "horizontal",
        legend.box.spacing = unit(0, "pt"),
        legend.key.size    = unit(6, "pt"),
        legend.spacing.x   = unit(2, "pt"),
        legend.title       = element_text(size = 6),
        legend.text        = element_text(size = 5.5),
        legend.margin      = margin(0, 0, 0, 0))

save_pdf <- function(plot, filename, w, h) {
  ggsave(file.path("final_plots", filename), plot,
         width = w, height = h, units = "in", device = cairo_pdf)
}


d <- fread("combined_within_pops_pi.txt.gz")
d <- d[pop != "Hybrid" & pop != "NC"]
d <- d[!grepl("^chr[XY]_", chromosome)]


read_bed3 <- function(path, name_col = NA) {
  cols <- if (is.na(name_col)) 1:3 else 1:max(3, name_col)
  x <- fread(path, header = FALSE, sep = "\t", quote = "", fill = TRUE,
             select = cols)
  setnames(x, c("V1", "V2", "V3"), c("chrom", "start", "end"))
  if (!is.na(name_col)) setnames(x, paste0("V", name_col), "name")
  x
}

cat("Loading annotations...\n")
gfeat   <- read_bed3("mPanTro3_genomefeatures.bed copy",         name_col = 4)
censat  <- read_bed3("mPanTro3_CenSat.bed copy",                  name_col = 4)
segdup  <- read_bed3("mPanTro3_sedefSegDups.bed copy")
tr      <- read_bed3("pantro_catalog.no_overlaps_simp.bed copy")

chroms_used <- unique(d$chromosome)
gfeat   <- gfeat [chrom %in% chroms_used]
censat  <- censat[chrom %in% chroms_used]
segdup  <- segdup[chrom %in% chroms_used]
tr      <- tr    [chrom %in% chroms_used]

cat(sprintf("  gfeat=%d  censat=%d  segdup=%d  TR=%d\n",
            nrow(gfeat), nrow(censat), nrow(segdup), nrow(tr)))

to_gr <- function(x) GRanges(seqnames = x$chrom,
                              ranges   = IRanges(start = x$start + 1,
                                                 end   = x$end))

cen_gr    <- reduce(to_gr(gfeat[name == "Cen"]))
telo_gr   <- reduce(to_gr(gfeat[name == "Telo"]))
censat_gr <- reduce(to_gr(censat))
segdup_gr <- reduce(to_gr(segdup))
tr_gr     <- reduce(to_gr(tr))

cat("Computing per-window overlap flags...\n")
windows <- unique(d[, .(chromosome, window_pos_1, window_pos_2)])
win_gr  <- GRanges(windows$chromosome,
                   IRanges(windows$window_pos_1, windows$window_pos_2))
windows[, in_Cen    := overlapsAny(win_gr, cen_gr)]
windows[, in_Telo   := overlapsAny(win_gr, telo_gr)]
windows[, in_CenSat := overlapsAny(win_gr, censat_gr)]
windows[, in_SegDup := overlapsAny(win_gr, segdup_gr)]
windows[, in_TR     := overlapsAny(win_gr, tr_gr)]

cat(sprintf("Background fractions: Cen=%.4f Telo=%.4f CenSat=%.3f SegDup=%.3f TR=%.3f\n",
            mean(windows$in_Cen),  mean(windows$in_Telo),
            mean(windows$in_CenSat), mean(windows$in_SegDup),
            mean(windows$in_TR)))

d <- merge(d, windows, by = c("chromosome", "window_pos_1", "window_pos_2"))

window_region <- function(in_Cen, in_Telo, in_CenSat) {
  fcase(in_Cen,    "Centromere",
        in_Telo,   "Telomere",
        in_CenSat, "CenSat",
        default = OTHER_LBL)
}


chrom_order <- function(x) {
  raw <- sub("^chr([^_]+)_.*", "\\1", x)
  suppressWarnings(as.integer(raw))
}

d[, chrom_num := chrom_order(chromosome)]
d[, mid := (window_pos_1 + window_pos_2) / 2]

chr_lengths <- d[!is.na(avg_pi),
                 .(chr_max = max(window_pos_2)),
                 by = .(chromosome, chrom_num)]
setorder(chr_lengths, chrom_num)
chr_lengths[, offset := data.table::shift(cumsum(as.numeric(chr_max)),
                                           fill = 0)]

axis_dt <- chr_lengths[, .(chromosome, chrom_num, chr_max, offset,
                            center = offset + chr_max / 2,
                            chrom_label = sub("^chr([^_]+)_.*", "\\1",
                                              chromosome))]
setorder(axis_dt, chrom_num)
axis_dt[, chrom_label := factor(chrom_label, levels = chrom_label)]


#  pi_long / pi_short ratio  (raw, no logs)

wide <- dcast(d[, .(pop, chromosome, window_pos_1, window_pos_2,
                    dataset, avg_pi, no_sites, count_comparisons,
                    in_Cen, in_Telo, in_CenSat, in_SegDup, in_TR)],
              pop + chromosome + window_pos_1 + window_pos_2 +
              in_Cen + in_Telo + in_CenSat + in_SegDup + in_TR ~ dataset,
              value.var = c("avg_pi", "no_sites", "count_comparisons"))

EPS <- 1 / median(wide$count_comparisons_short[wide$count_comparisons_short > 0],
                  na.rm = TRUE)
cat(sprintf("Epsilon (1 / median count_comparisons_short) = %.3e\n", EPS))

wide <- wide[!is.na(avg_pi_long)]
both_zero <- (is.na(wide$avg_pi_long)  | wide$avg_pi_long  == 0) &
             (is.na(wide$avg_pi_short) | wide$avg_pi_short == 0)
wide <- wide[!both_zero]

wide[, pi_short_adj := fifelse(is.na(avg_pi_short) | avg_pi_short == 0,
                                EPS, avg_pi_short)]
wide[, ratio       := avg_pi_long / pi_short_adj]
wide[, log10_ratio := log10(ratio)]

wide[, chrom_num := chrom_order(chromosome)]
wide[, mid := (window_pos_1 + window_pos_2) / 2]
wide <- merge(wide, chr_lengths[, .(chromosome, offset)], by = "chromosome")
wide[, cum_pos := mid + offset]

# Top 1% ratio outliers per pop
thr_ratio <- wide[, .(thresh = quantile(ratio, 1 - OUTLIER_PCT)), by = pop]
wide <- merge(wide, thr_ratio, by = "pop", all.x = TRUE)
wide[, is_ratio_outlier := ratio >= thresh]

assign_outlier_region <- function(is_outlier, in_Cen, in_Telo,
                                   in_CenSat, in_SegDup, in_TR) {
  fcase(!is_outlier, NA_character_,
        in_Cen,      "Centromere",
        in_Telo,     "Telomere",
        in_CenSat,   "CenSat",
        in_SegDup,   "SegDup",
        in_TR,       "TR",
        default = "Other")
}
wide[, region := assign_outlier_region(is_ratio_outlier, in_Cen, in_Telo,
                                        in_CenSat, in_SegDup, in_TR)]

hilite_ratio <- wide[is_ratio_outlier == TRUE & region != "Other"]
hilite_ratio[, region := factor(region, levels = names(REGION_COLORS))]

y_cap <- quantile(wide$ratio, 0.995)


p_ratio <- ggplot() +
  rasterise(geom_point(data = wide,
                       aes(x = cum_pos, y = ratio,
                           color = factor(chrom_num %% 2)),
                       size = 0.18, alpha = 0.45, shape = 16),
            dpi = 300) +
  scale_color_manual(values = c("0" = "grey25", "1" = "grey60"),
                     guide = "none") +
  ggnewscale::new_scale_color() +
  geom_point(data = hilite_ratio,
             aes(x = cum_pos, y = ratio, color = region, shape = region),
             size = 1.6, stroke = 0.5, fill = NA) +
  scale_color_manual(values = REGION_COLORS, name = "Outlier", drop = FALSE) +
  scale_shape_manual(values = REGION_SHAPES, name = "Outlier", drop = FALSE) +
  scale_x_continuous(breaks = axis_dt$center, labels = axis_dt$chrom_label,
                     expand = expansion(mult = c(0.005, 0.005))) +
  facet_wrap(~ pop, ncol = 1, strip.position = "right") +
  coord_cartesian(ylim = c(0, y_cap)) +
  labs(title = expression(pi[long] / pi[short] ~
                            "(10 kb windows, raw, top 1% outliers)"),
       subtitle = paste0("y capped at 99.5th percentile (",
                         round(y_cap, 1), ")"),
       x = "Chromosome", y = expression(pi[long] / pi[short])) +
  compact_theme
save_pdf(p_ratio, "manhattan_pi_ratio_capped.pdf", 8, 3.5)


wide[, chrom_label := factor(sub("^chr([^_]+)_.*", "\\1", chromosome),
                              levels = levels(axis_dt$chrom_label))]
hilite_ratio[, chrom_label := factor(sub("^chr([^_]+)_.*", "\\1", chromosome),
                                       levels = levels(axis_dt$chrom_label))]

p_facet <- ggplot() +
  rasterise(geom_point(data = wide,
                       aes(x = mid / 1e6, y = ratio),
                       size = 0.18, alpha = 0.45, shape = 16,
                       color = "grey40"),
            dpi = 300) +
  geom_point(data = hilite_ratio,
             aes(x = mid / 1e6, y = ratio, color = region, shape = region),
             size = 1.3, stroke = 0.4, fill = NA) +
  scale_color_manual(values = REGION_COLORS, name = "Outlier", drop = FALSE) +
  scale_shape_manual(values = REGION_SHAPES, name = "Outlier", drop = FALSE) +
  facet_grid(pop ~ chrom_label, scales = "free_x", space = "free_x",
             switch = "y") +
  coord_cartesian(ylim = c(0, y_cap)) +
  labs(title = expression(pi[long] / pi[short] ~
                            "(10 kb, pop x chromosome, top 1% outliers)"),
       x = "Position (Mb)", y = expression(pi[long] / pi[short])) +
  compact_theme +
  theme(axis.text.x = element_text(size = 4.5, angle = 45, hjust = 1),
        strip.text  = element_text(size = 5.5))
save_pdf(p_facet, "manhattan_pi_ratio_facet_chrom.pdf", 13, 3)


p_facet_clean <- ggplot(wide) +
  rasterise(geom_point(aes(x = mid / 1e6, y = ratio),
                       size = 0.18, alpha = 0.45, shape = 16,
                       color = "grey40"),
            dpi = 300) +
  facet_grid(pop ~ chrom_label, scales = "free_x", space = "free_x",
             switch = "y") +
  coord_cartesian(ylim = c(0, y_cap)) +
  labs(title = expression(pi[long] / pi[short] ~
                            "(10 kb, pop x chromosome, clean)"),
       x = "Position (Mb)", y = expression(pi[long] / pi[short])) +
  compact_theme +
  theme(axis.text.x = element_text(size = 4.5, angle = 45, hjust = 1),
        strip.text  = element_text(size = 5.5))
save_pdf(p_facet_clean, "manhattan_pi_ratio_facet_chrom_clean.pdf", 13, 3)


#  Density plots

wide[, win_region := factor(
  window_region(in_Cen, in_Telo, in_CenSat),
  levels = REGION_LEVELS_DENS)]

# log10(ratio) density (per pop)
p_dens_log10 <- ggplot(wide,
                        aes(x = log10_ratio, fill = win_region,
                            color = win_region)) +
  geom_density(alpha = 0.35, linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black",
             linewidth = 0.25) +
  scale_fill_manual(values = REGION_COLORS_DENS, name = NULL) +
  scale_color_manual(values = REGION_COLORS_DENS, name = NULL) +
  guides(color = guide_legend(nrow = 2),
         fill  = guide_legend(nrow = 2)) +
  facet_wrap(~ pop, ncol = 3) +
  labs(x = expression(log[10](pi[long] / pi[short])), y = "Density") +
  panel_theme
save_pdf(p_dens_log10, "density_log10_ratio_by_region.pdf", 3, 1.4)

# pi distribution per dataset, coloured by region
d_dist <- d[!is.na(avg_pi) & avg_pi > 0,
            .(dataset, pop, avg_pi,
              in_Cen, in_Telo, in_CenSat, in_SegDup, in_TR)]
d_dist[, win_region := factor(
  window_region(in_Cen, in_Telo, in_CenSat),
  levels = REGION_LEVELS_DENS)]
d_dist[, dataset := factor(dataset, levels = c("long", "short"))]

p_pi_dens <- ggplot(d_dist,
                    aes(x = avg_pi, color = win_region,
                        fill = win_region)) +
  geom_density(alpha = 0.25, linewidth = 0.3) +
  scale_x_log10(labels = label_number(accuracy = 0.0001),
                breaks = c(1e-4, 1e-2)) +
  scale_color_manual(values = REGION_COLORS_DENS, name = NULL) +
  scale_fill_manual(values  = REGION_COLORS_DENS, name = NULL) +
  guides(color = guide_legend(nrow = 2),
         fill  = guide_legend(nrow = 2)) +
  facet_grid(dataset ~ pop) +
  labs(x = expression(pi ~ "(log"[10] * ")"), y = "Density") +
  panel_theme
save_pdf(p_pi_dens, "density_pi_by_dataset_region.pdf", 3, 1.8)

cat("\n=== Done. PDFs written to final_plots/ ===\n")
