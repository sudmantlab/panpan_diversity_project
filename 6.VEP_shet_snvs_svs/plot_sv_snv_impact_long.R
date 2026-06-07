setwd("/Users/joanocha/Google Drive/My Drive/POSTDOC/PANPAN/analysis/Figure3_VEP_snvs_svs")
library(tidyverse)
library(data.table)
library(pals)
library(MetBrewer)
library(cowplot)
library(knitr)
library(vcfR)
library(janitor)
library(ggupset)
library(ggh4x)
library(patchwork)

snv_data <- fread("snvs_vep/long_only/snv_impact_long.tsv.gz") %>%
  mutate(variant_type = "SNVs")

sv_data <- fread("svs_vep/sv_impact_long.tsv") %>%
  mutate(variant_type = "SVs")

combined_data <- bind_rows(snv_data, sv_data)
impact_levels <- c("HIGH", "MODERATE", "LOW", "MODIFIER")
variant_type_levels <- c("SNVs", "SVs")

################ STRACKED PLOT FOR PROPORTIONS OF DATA" ####################
counted_data <- combined_data %>%
  count(species, variant_type, impact, name = "n")
# 1: Proportion within each Variant Type
plot_data1 <- counted_data %>%
  group_by(species, variant_type) %>%
  mutate(
    proportion = n / sum(n),
    impact = factor(impact, levels = impact_levels)
  ) %>%
  ungroup() %>%
  filter(impact != "MODIFIER")

# 2: Proportion of ALL variants per species
plot_data2 <- counted_data %>%
  group_by(species) %>%
  mutate(
    proportion = n / sum(n),
    impact = factor(impact, levels = impact_levels)
  ) %>%
  ungroup() %>%
  filter(impact != "MODIFIER")

create_variant_plot <- function(data, y_lab) {
  ggplot(data, aes(x = variant_type, y = proportion, fill = impact)) +
    geom_col(color = "black", position = "stack") +
    facet_wrap(~ species) +
    scale_fill_manual(values = c("HIGH" = "darkred", "MODERATE" = "darkorange", "LOW" = "lightgrey")) +
    labs(
      x = "",
      y = y_lab,
      fill = "VEP Impact"
    ) +
    theme_bw() + # Using theme_bw() for consistency as it was in your second plot
    theme(
      axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color = "grey", fill = NA),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      strip.background = element_rect(fill = "grey90")
    )
}


variant_plot1 <- create_variant_plot(
  data = plot_data1,
  y_lab = "Proportion SVs and SNVs with functional impact"
)

variant_plot2 <- create_variant_plot(
  data = plot_data2,
  y_lab = "Proportion of all variants with functional impact"
)

combined_plot <- variant_plot2 + variant_plot1 +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'top')
print(combined_plot)
ggsave("variant_impact_combined.pdf", plot = combined_plot, width = 6.2, height = 5, dpi = 300)


#### length distribution of SVs with different impact levels
sv_data %>%
  mutate(impact=fct_relevel(impact, c("MODIFIER", "LOW","MODERATE", "HIGH"))) %>%
  mutate(svtype=fct_relevel(svtype, c("INS", "DEL", "DUP", "INV"))) %>%
  ggplot(aes(x=abs(svlen), fill = svtype)) +
  geom_histogram(color="black", bins=50, linewidth=0.2) +
  scale_x_log10() +
  facet_grid(impact~species, scales = "free_y") +
  labs(x="SV length") +
  theme_classic()
sv_data %>%
  mutate(type_simplified=fct_relevel(type_simplified, c("SINE/Alu", "LINE/L1", "Retroposon/SVA", "Simple/tandem repeat", "Other repeats", "Non-repeat"))) %>%
  filter(svtype %in% c("INS", "DEL")) %>%
  mutate(impact=fct_relevel(impact, c("MODIFIER", "LOW","MODERATE", "HIGH"))) %>%
  ggplot(aes(x=abs(svlen), fill = type_simplified)) +
  geom_histogram(color="black", bins=50, linewidth=0.2) +
  scale_fill_manual(values = pals::cols25(n=25)[c(1:3, 7, 8, 13)], breaks=c("SINE/Alu", "LINE/L1", "Retroposon/SVA", "Simple/tandem repeat", "Other repeats", "Non-repeat"), drop=FALSE) +
  scale_x_log10() +
  facet_grid(impact~species, scales = "free_y") +
  labs(x="SV length") +
  theme_classic()



#### SFS of SVs with different impact levels ########################
sv_impact_sfs_tmp <- crossing(species="bonobo", total_allele_count=1:(10-1), impact_simplified=c("modifier/low", "moderate/high")) %>%
  bind_rows(crossing(species="chimpanzee", total_allele_count=1:(48-1), impact_simplified=c("modifier/low", "moderate/high"))) %>%
  bind_rows(crossing(species="human", total_allele_count=1:94, impact_simplified=c("modifier/low", "moderate/high")))
sv_impact_sfs <- sv_data  %>%
  filter(species!="bonobo" | total_allele_count!=10) %>%
  filter(species!="chimpanzee" | total_allele_count!=48) %>%
  mutate(impact_simplified = case_when(
    impact == "HIGH"     ~ "moderate/high",
    impact == "MODERATE" ~ "moderate/high",
    impact == "LOW"      ~ "modifier/low",
    impact == "MODIFIER" ~ "modifier/low",
    TRUE                 ~ NA_character_
  )) %>%
  count(species, total_allele_count, impact_simplified) %>%
  left_join(sv_impact_sfs_tmp, .) %>%
  mutate(n=ifelse(is.na(n), 0, n)) %>%
  mutate(p=n/sum(n), .by=c(species, impact_simplified))
my_sv_colors <- c(
  "moderate/high" = "darkred",
  "modifier/low" = "lightgrey"
)
custom_x <- list(
  scale_x_continuous(breaks = c(1, 9)),
  scale_x_continuous(breaks = c(1, (1:5)*10)),
  scale_x_continuous(breaks = c(1, (1:4)*25))
)

# Plot 1: Line and Ribbon Plot
sv_line_plot <- sv_impact_sfs %>%
  ggplot(aes(x=total_allele_count, y=p)) +
  geom_line(aes(color=impact_simplified)) +
  geom_point(aes(color=impact_simplified), size=0.2) +
  geom_ribbon(aes(fill=impact_simplified, ymax=p), ymin=0, alpha=0.2) +
  scale_color_manual(values = my_sv_colors) + # Use the new palette
  scale_fill_manual(values = my_sv_colors) +  # Use the new palette
  labs(x="number of non-reference alleles", y="proportion of SVs", color="VEP impact", fill="VEP impact") +
  scale_y_sqrt() +
  facet_grid(.~species, scale="free", space="free_x")+
  ggh4x::facetted_pos_scales(x = custom_x) +
  theme_bw() +
  theme(legend.position = "top", panel.grid = element_blank())
sv_col_plot <- sv_impact_sfs %>%
  ggplot(aes(x=total_allele_count, y=p)) +
  geom_col(aes(fill=impact_simplified), position=position_dodge(), color="black", linewidth=0.1) +
  scale_fill_manual(values = my_sv_colors) + # Use the new palette
  labs(x="number of non-reference alleles", y="proportion of SVs", color="VEP impact", fill="VEP impact") +
  scale_y_sqrt() +
  facet_grid(.~species, scale="free", space="free_x")+
  ggh4x::facetted_pos_scales(x = custom_x) +
  theme_classic() +
  theme(legend.position = "top", panel.grid = element_blank())
#print(sv_line_plot)
print(sv_col_plot)

######### SFS SNVS WITH DIFFERENT IMPACT LEVELS ###########
snv_impact_sfs_tmp <- crossing(species="bonobo", total_allele_count=1:(10-1), impact_simplified=c("modifier/low", "moderate/high")) %>%
  bind_rows(crossing(species="chimpanzee", total_allele_count=1:(48-1), impact_simplified=c("modifier/low", "moderate/high"))) %>%
  bind_rows(crossing(species="human", total_allele_count=1:94, impact_simplified=c("modifier/low", "moderate/high")))
# Process the snv_data to calculate proportions
snv_impact_sfs <- snv_data %>%
  filter(species!="bonobo" | total_allele_count!=10) %>%
  filter(species!="chimpanzee" | total_allele_count!=48) %>%
  mutate(impact_simplified = case_when(
    impact == "HIGH"     ~ "moderate/high",
    impact == "MODERATE" ~ "moderate/high",
    impact == "LOW"      ~ "modifier/low",
    impact == "MODIFIER" ~ "modifier/low",
    TRUE                 ~ NA_character_
  )) %>%
  count(species, total_allele_count, impact_simplified) %>%
  left_join(sv_impact_sfs_tmp, .) %>%
  mutate(n=ifelse(is.na(n), 0, n)) %>%
  mutate(p=n/sum(n), .by=c(species, impact_simplified))
my_colors <- c(
  "moderate/high" = "darkred",
  "modifier/low" = "lightgrey"
)
custom_x <- list(
  scale_x_continuous(breaks = c(1, 9)),
  scale_x_continuous(breaks = c(1, (1:5)*10)),
  scale_x_continuous(breaks = c(1, (1:4)*25))
)
# Plot 1: Line and Ribbon Plot
snv_line_plot <- snv_impact_sfs %>%
  ggplot(aes(x=total_allele_count, y=p)) +
  geom_line(aes(color=impact_simplified)) +
  geom_point(aes(color=impact_simplified), size=0.2) +
  geom_ribbon(aes(fill=impact_simplified, ymax=p), ymin=0, alpha=0.2) +
  scale_color_manual(values = my_sv_colors) + # Use the new palette
  scale_fill_manual(values = my_sv_colors) +  # Use the new palette
  labs(x="number of non-reference alleles", y="proportion of SVs", color="VEP impact", fill="VEP impact") +
  scale_y_sqrt() +
  facet_grid(.~species, scale="free", space="free_x")+
  ggh4x::facetted_pos_scales(x = custom_x) +
  theme_bw() +
  theme(legend.position = "top", panel.grid = element_blank())
snv_col_plot <- snv_impact_sfs %>%
  ggplot(aes(x = total_allele_count, y = p)) +
  geom_col(aes(fill = impact_simplified), position = position_dodge(), color = "black", linewidth = 0.1) +
  scale_fill_manual(values = my_colors) + # Use the new palette
  labs(x = "number of non-reference alleles", y = "proportion of SNVs", color = "VEP impact", fill = "VEP impact") +
  scale_y_sqrt() +
  facet_grid(. ~ species, scale = "free", space = "free_x") +
  ggh4x::facetted_pos_scales(x = custom_x) +
  theme_classic() +
  theme(legend.position = "top", panel.grid = element_blank())
print(snv_col_plot)

combined_plot2 <- snv_col_plot / sv_col_plot +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'right')
print(combined_plot2)
ggsave('SFS_SNVs_SVs.pdf', combined_plot2, width =11, height =4)

combined_plot2 <- snv_line_plot / sv_line_plot +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'top')

