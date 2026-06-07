required_packages <- c("ggplot2", "dplyr", "tidyr", "ggsignif", "stringr")

# Check and install missing packages
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

invisible(lapply(required_packages, install_if_missing))

library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(ggsignif)

df <- read.table("homo_pantro_heteroz.txt")

df <- df %>%
  rename(
    spp = V1,
    TRID = V2,
    Ho = V3,
    He = V4,
    pathogenic_status = V5,
    motif_len = V6,
    motif_seq = V7,
    feature = V8)

df <- df %>%
  mutate(spp_pat_status = interaction(spp, pathogenic_status, sep = "_"))

mycolors=c('homo'='#5C608A','chimp'='#8A865D')

# Plot expected heterozygosity relative to the pathogenic status of TR

test1 <- wilcox.test(He ~ pathogenic, data = df)
test2 <- wilcox.test(He ~ spp, data = filter(df, pathogenic == "pathogenic"))
test3 <- wilcox.test(He ~ spp, data = filter(df, pathogenic == "non_pathogenic"))
max_y <- max(df$He, na.rm = TRUE)

pvals <- data.frame(
  group1 = c("non_patHogenic", "chimp_pathogenic", "chimp_non_pathogenic"),
  group2 = c("pathogenic", "homo_pathogenic", "homo_non_pathogenic"),
  y.position = c(max_y * 1.05, max_y * 1.15, max_y * 1.25),
  p.adj = c(test1$p.value, test2$p.value, test3$p.value),
  p.label = c(
    scales::pvalue(test1$p.value),
    scales::pvalue(test2$p.value),
    scales::pvalue(test3$p.value)))

p <- ggplot(df, aes(x = factor(pathogenic), y = He, fill = spp, colour= spp)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, position = position_dodge(width = 0.75), show.legend=F) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, position = position_dodge(width = 0.75), show.legend=F) +
  # pathogenic vs non-pathogenic (across x-axis groups)
  geom_signif(
    comparisons = list(c("non_pathogenic", "pathogenic")),
    map_signif_level = TRUE, y_position = max_y * 1.05, size = 0.4, tip_length = 0.03, textsize = 4, family="Helvetica", color = "black") +
  # chimp vs human within pathogenic
  geom_signif(
    annotations = ifelse(test2$p.value < 0.001, "***", ifelse(test2$p.value < 0.01, "**", ifelse(test2$p.value < 0.05, "*", "ns"))),
    y_position = max_y * 1.15, xmin = 2 - 0.3, xmax = 2 + 0.3, size = 0.4, tip_length = 0.03, textsize = 4, family="Helvetica", color = "black") +
  # chimp vs human within non-pathogenic
  geom_signif(
    annotations = ifelse(test3$p.value < 0.001, "***", ifelse(test3$p.value < 0.01, "**", ifelse(test3$p.value < 0.05, "*", "ns"))),
    y_position = max_y * 1.25, xmin = 1 - 0.3, xmax = 1 + 0.3, size = 0.4, tip_length = 0.03, textsize = 4, family="Helvetica", color = "black") +
  theme(legend.position = "none") +
  guides(fill = "none", colour = "none") +
  scale_x_discrete(labels = c("non_pathogenic" = "Non-pathogenic", "pathogenic" = "Pathogenic")) +
  scale_y_continuous(name = "Genetic diversity") +
  scale_fill_manual(values = mycolors) +
  scale_colour_manual(values = mycolors) +
  theme_bw() +
  theme(
    axis.line = element_line(linewidth = 1, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    text = element_text(family="Helvetica"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(colour = "black", size = 8, family="Helvetica"),
    axis.text.x = element_text(colour = "black", size = 8, family="Helvetica"),
    axis.text.y = element_text(colour = "black", size = 8, family="Helvetica"))

ggsave("het_pathogenic_trs.pdf", plot = p, device = cairo_pdf, width = 3, height = 6)
