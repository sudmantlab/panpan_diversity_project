library(ggplot2)
library(dplyr)
library(data.table)
library(tidyr)
library(tibble)

# 1. LOAD DATA
afr     <- fread("diversity_stats/combined_pi_stats_AFR.tsv")
non_afr <- fread("diversity_stats/combined_pi_stats_nonAFR.tsv")
pantros <- fread("diversity_stats/pantros_combined_pi_stats.tsv")
ppan    <- fread("diversity_stats/ppan_combined_pi_stats.tsv")
metadata <- fread("PANPAN_PradoMartine_deManuel_METADATA.tsv")

# 2. HUMAN CALIBRATION
# Using the callable genome size: 2,875,101,262 bp
#The human scaling now uses the NCBI value for the GRCh38 non-gap sequence length (2,875,101,262$ bp)
#Method: $(\pi_{raw} \times 130,475,909) / 2,875,101,262.
callable_human_g <- 2875101262

# Calibrated Pi = (raw_pi * n_sites) / callable_human_g
avg_afr     <- mean(afr$mean_pi * afr$n_sites / callable_human_g)
avg_non_afr <- mean(non_afr$mean_pi * non_afr$n_sites / callable_human_g)

# 3. APE DATA PROCESSING
apes_raw <- bind_rows(pantros, ppan)

# 1. Deduplicate the metadata so there is only 1 row per individual
meta_unique <- metadata %>%
  distinct(sample, .keep_all = TRUE)

# 2. Join the raw stats with the unique metadata
apes_processed <- apes_raw %>%
  left_join(meta_unique, by = "sample")
head(apes_processed)
# Define Evolutionary Population Order
pop_order <- c("Bonobo", "Central", "Eastern", "Nigeria-Cameroon", "Western", "Western x Central")
apes_processed$CommonName <- factor(apes_processed$CommonName, levels = pop_order)

# Extract specific color mapping from the metadata "colors" column
color_map <- apes_processed %>%
  select(CommonName, colors) %>%
  distinct() %>%
  deframe()

# central '#4c5d4c'
# western '#9dced9'
# bonobo '#e9d7cb' 
# eastern '#ffb35a'
# hybrid '#bbc671'
# --- PLOT 1: Individual Barplot (Evolutionary Order + Diversity Sorted) ---
# Logic: Order by pop_order first, then by pi value within each pop
plot1_data <- apes_processed %>%
  arrange(CommonName, mean_pi) %>%
  mutate(sample_factor = factor(sample, levels = unique(sample)))


plot1 <- ggplot(plot1_data, aes(x = sample_factor, y = mean_pi, fill = CommonName, alpha = type3)) +
  geom_bar(stat = "identity") +
  # Reference Lines (Dashed)
  geom_hline(yintercept = avg_afr, color = "black", linetype = "dotted", size = 0.7) +
  geom_hline(yintercept = avg_non_afr, color = "black", linetype = "dotted", size = 0.7) +
  # Text annotations for the human lines
  annotate("text", x = 5, y = avg_afr, label = "Human AFR", vjust = -0.5, color = "black", fontface="bold") +
  annotate("text", x = 5, y = avg_non_afr, label = "non-AFR", vjust = 1.5, color = "black", fontface="bold") +
  scale_alpha_manual(values = c("long-read" = 1, "short-read" = 0.5), name = "Read Type") +
  scale_fill_manual(values = color_map) +
  labs(x = "Individuals", y = "SNV heterozygosity", title = "Individual Diversity (Ordered by Population & Pi)") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
plot1
# --- PLOT 2: Population Box Plot ---
plot2_data <- apes_processed %>% filter(CommonName != "Western x Central")
plot2 <- ggplot(plot2_data, aes(x = CommonName, y = mean_pi, fill = CommonName)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1.5) +
  scale_fill_manual(values = color_map) +
  labs(x = "", y = "SNV heterozygosity") +
  theme_classic() +
  coord_flip() + # <--- This flips the plot horizontally
  theme(
    legend.position = "none", # Optional: hides legend since Y-axis already has names
    axis.text.y = element_text(angle = 0, hjust = 1) # Changed to axis.text.y
  )

plot2
# --- PLOT 3: Faceted by CommonName (Short vs Long side-by-side) ---


plot3_data <- apes_processed %>% 
  filter(!CommonName %in% c( "Nigeria-Cameroon", "Western x Central"))

plot3 <- ggplot(plot3_data, aes(x = type3, y = mean_pi)) +
  
  # Base box plots (removed trim = FALSE as it only applies to violins)
  geom_boxplot(fill = "white", color = "gray30", alpha = 0.8) +
  
  # Jittered dots colored by CommonName (population)
  geom_jitter(aes(color = CommonName), width = 0.15, alpha= 0.8, size = 1.5) +
  
  # Keeps the plots horizontal
  coord_flip() + 
  
  # Facet labels on the left, stacked in one column
  facet_wrap(~CommonName, ncol = 1, strip.position = "left") +
  
  # Apply your specific metadata colors to the dots
  scale_color_manual(values = color_map, name = "Population") +
  
  labs(x = "", 
       y = "SNV heterozygosity") +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9),
    legend.position = "right",
    
    # --- FACET LABEL FORMATTING ---
    strip.placement = "outside",          # Keeps them acting like axis labels
    
    # 1. "In boxes": Adds the black border and white background back
    strip.background = element_rect(fill = "white", color = "black"), 
    
    # 2. "Angled 90 degrees" and "less bold": angle = 90, face = "plain"
    strip.text.y.left = element_text(angle = 90, size = 9, face = "plain") 
  )

print(plot3)
# ggsave("Pi_Bars_Evolutionary_Order.pdf", plot1, width = 12, height = 5)
# ggsave("Pi_Population_Boxplots.pdf", plot2, width = 7, height = 5)
# ggsave("Pi_Faceted_ReadType.pdf", plot3, width = 10, height = 5)