setwd("/Users/joanocha/Library/CloudStorage/GoogleDrive-joana.laranjeira.rocha@gmail.com/My Drive/POSTDOC/PANPAN/analysis/Figure1_Diversity_Assembly_stats/PANPAN_assembly_stats/")
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(ggrepel)
library(dplyr)
library(tidyr)
library(gghighlight)
library(ggplot2)
library(dplyr)
library(data.table)
# 1. Load both datasets
df_stats <- fread("genome_stats_summary.tsv")
df_meta  <- fread("PanPan_Metadata.tsv")
meta_unique <- df_meta %>%
  select(genomeName, CommonName, colors) %>%
  distinct()
df_stats <- df_stats %>%
  left_join(meta_unique, by = c("ID" = "genomeName"))
df_stats <- df_stats %>%
  group_by(sex) %>%
  mutate(sex_label = paste0(sex, " (n=", n(), ")")) %>%
  ungroup()
label_levels <- sort(unique(df_stats$sex_label))

color_mapping <- df_stats %>% 
  select(CommonName, colors) %>% 
  distinct() %>% 
  filter(!is.na(CommonName))

named_colors <- setNames(color_mapping$colors, color_mapping$CommonName)

# --- NEW: Calculate the averages for both axes ---
avg_depth <- mean(df_stats$Coverage_depth, na.rm = TRUE)
avg_kmer  <- mean(df_stats$Coverage_Kmer, na.rm = TRUE)


# 4. Create the Plot
plot_coverage <- ggplot(data = df_stats, 
                        aes(x = Coverage_depth, 
                            y = Coverage_Kmer)) +
  
  # Add a 1:1 dashed reference line (y = x) so you can easily see deviations
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "grey70", linewidth = 0.6) +
  
  # --- NEW: Add the vertical and horizontal average lines ---
  # Using a darker color ("grey30") so they stand out against the 1:1 line
  geom_vline(xintercept = avg_depth, linetype = "dashed", color = "grey30", linewidth = 0.7, alpha = 0.8) +
  geom_hline(yintercept = avg_kmer, linetype = "dashed", color = "grey30", linewidth = 0.7, alpha = 0.8) +
  
# Add the points mapping Color to CommonName and Shape to sex
geom_point(aes(color = CommonName, shape = sex_label), size = 3.5, stroke = 1.2, alpha = 0.8) +
  
  scale_shape_manual(values = setNames(c(4, 15), label_levels), name = "Sex") +
  
  # Apply your dataset's exact custom colors
  scale_color_manual(values = named_colors, name = "Population") +
  
  # Labels
  xlab("Read Depth Coverage") +
  ylab("K-mer Based Coverage") +
  
  # Theme matching the QV plot aesthetics
  theme_classic() +
  theme(
    # Note: Because there are now TWO legends (Color and Shape), placing them inside 
    # the plot might block some data points. "right" is the safest option here.
    legend.position = "right", 
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11, color = "black"),
    legend.background = element_rect(fill = "transparent", color = NA),
    
    # Creates the clean, squared-off black border around the plot area
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(), 
    
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 13)
  )

print(plot_coverage)
ggsave('plot_coverage.pdf', plot_coverage, width=5, height = 3, dpi=300)

assemblies<-read.delim("allPan_references_corrected.tsv")
dataset<-fread("PanPan_Metadata.tsv")
meta_statistics <- merge(assemblies,dataset,by=c("genomeName", "type"), all.x=TRUE)
meta_statistics<- meta_statistics[meta_statistics$genomeName != "AG18352_2", ]


meta_statistics <- meta_statistics %>%
  mutate(
    method_rank = case_when(
      grepl("Reference", method, ignore.case = TRUE) ~ 1, # Protects references
      method == "HiFi+HiC+ONT" ~ 2,
      method == "HiFi+HiC" ~ 3,
      method == "HiFi" ~ 4,
      TRUE ~ 5 # Catch-all for any other unexpected methods
    )
  ) %>%
  # 2. Group by genomeName
  group_by(genomeName) %>%
  # 3. Filter to keep only the rows matching the highest priority (lowest number) for that genome
  filter(method_rank == min(method_rank)) %>%
  # 4. Remove the grouping and the temporary rank column
  ungroup() %>%
  select(-method_rank)

NGx<-meta_statistics[meta_statistics $stat == 'NGx',]
NGx_50<-NGx[NGx$v1 == 50,]
NG50_means <- NGx_50 %>%
  # Filter to include only the three methods you care about
  filter(method %in% c("HiFi", "HiFi+HiC", "HiFi+HiC+ONT")) %>%
  
  # Group the data by the method
  group_by(method) %>%
  
  # Calculate the mean of v2 (the NG50 length)
  summarise(
    mean_NG50 = mean(v2, na.rm = TRUE),
    count_of_assemblies = n() # Counts how many samples are in this group
  )

print(NG50_means)

NG50_means_combined <- NGx_50 %>%
  # 1. Filter to include only the three methods
  filter(method %in% c("HiFi", "HiFi+HiC", "HiFi+HiC+ONT")) %>%
  
  # 2. Create a new column that lumps HiFi and HiFi+HiC together
  mutate(
    method_group = if_else(
      method %in% c("HiFi", "HiFi+HiC"), 
      "HiFi & HiFi+HiC (Combined)", 
      "HiFi+HiC+ONT"
    )
  ) %>%
  
  # 3. Group the data by this NEW merged label
  group_by(method_group) %>%
  
  # 4. Calculate the mean and count
  summarise(
    mean_NG50 = mean(v2, na.rm = TRUE),
    count_of_assemblies = n() # This will easily verify your 28+22 = 50 count!
  )

print(NG50_means_combined)

NGx<-meta_statistics[meta_statistics $stat == 'NGx',]
NGx_50<-NGx[NGx$v1 == 50,]
NGx_50_primary<-NGx_50[NGx_50$type == "primary",]
NGx_50_hap1<-NGx_50[NGx_50$type == "hap1",]
NGx_50_hap2<-NGx_50[NGx_50$type == "hap2",]
auN<-meta_statistics [meta_statistics $stat == 'auN',]
auNG<-meta_statistics [meta_statistics $stat == 'auNG',]
# Calculate averages by population, excluding T2T haplotypes
#NGx$type <- factor(NGx$type, levels = c("primary", "hap1", "hap2"))
# 1. Reshape the auNG data to have hap1 and hap2 side-by-side
auNG_wide <- auNG %>%
  # Keep only the haplotypes
  filter(type %in% c("hap1", "hap2")) %>%
  filter(!genomeName %in% c("mPanPan1", "mPanTro3")) %>%
  # (Assuming the auNG value is stored in 'v2' just like NGx)
  select(genomeName, type, v2, CommonName, sex, colors) %>% 
  
  # Pivot so hap1 and hap2 become their own distinct columns containing the v2 values
  pivot_wider(names_from = type, values_from = v2) 


# 2. Dynamically create the legend labels to include counts (e.g., "female (n=14)")
auNG_wide <- auNG_wide %>%
  group_by(sex) %>%
  mutate(sex_label = paste0(sex, " (n=", n(), ")")) %>%
  ungroup()

# Extract unique shapes sorted alphabetically (female first, male second)
label_levels <- sort(unique(auNG_wide$sex_label))

# 3. Automatically extract the exact name-to-color mapping directly from your metadata
color_mapping <- auNG_wide %>% 
  select(CommonName, colors) %>% 
  distinct() %>% 
  filter(!is.na(CommonName))

named_colors <- setNames(color_mapping$colors, color_mapping$CommonName)

# 4. Calculate the averages for the dashed reference lines specific to auNG
avg_hap1 <- mean(auNG_wide$hap1, na.rm = TRUE)
avg_hap2 <- mean(auNG_wide$hap2, na.rm = TRUE)
tmp_avg_hap1 <- median(auNG_wide$hap1, na.rm = TRUE)
tmo_avg_hap2 <- median(auNG_wide$hap2, na.rm = TRUE)
# 5. Create the Plot
library(scales)
plot_auNG <- ggplot(data = auNG_wide, 
                    aes(x = hap1, 
                        y = hap2)) +
  
  # Add a 1:1 reference line (y = x) so you can easily see deviations
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "grey70", linewidth = 0.6) +
  
  # Add the vertical and horizontal average lines
  geom_vline(xintercept = avg_hap1, linetype = "dashed", color = "grey30", linewidth = 0.7, alpha = 0.8) +
  geom_hline(yintercept = avg_hap2, linetype = "dashed", color = "grey30", linewidth = 0.7, alpha = 0.8) +
  
  # Add the points mapping Color to CommonName and Shape to sex
  geom_point(aes(color = CommonName, shape = sex_label), size = 3.5, stroke = 1.2, alpha = 0.8) +
  
  # Assign shapes matching your Figure 1 plot: 4 = 'x' cross (female), 15 = solid square (male)
  scale_shape_manual(values = setNames(c(4, 15), label_levels), name = "Sex") +
  
  # Apply your dataset's exact custom colors
  scale_color_manual(values = named_colors, name = "Population") +
  
  # Formatting X and Y axes to display commas for large sequence lengths (e.g., 100,000,000)
  scale_x_continuous(labels = label_comma()) +
  scale_y_continuous(labels = label_comma()) +
  
  # Labels
  xlab("H1 auNG (bp)") +
  ylab("H2 auNG (bp)") +
  
  theme_classic() +
  theme(
    legend.position = "right", 
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11, color = "black"),
    legend.background = element_rect(fill = "transparent", color = NA),
    
    # Creates the clean, squared-off black border around the plot area
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(), 
    
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 13)
  )

print(plot_auNG)
ggsave("auNG_Comparison_Hap1_vs_Hap2.pdf", plot_auNG, width = 5, height = 3, dpi = 300)

# ---------------------------------------------------------
# Calculate the mean auNG for each specific method
# ---------------------------------------------------------
auNG_means <- auNG %>%
  # Filter to include only the three methods you care about
  filter(method %in% c("HiFi", "HiFi+HiC", "HiFi+HiC+ONT")) %>%
  group_by(method) %>%
  summarise(
    mean_auNG = mean(v2, na.rm = TRUE),
    count_of_assemblies = n()
  )

print(auNG_means)

# ---------------------------------------------------------
# Calculate the combined mean auNG (HiFi & HiFi+HiC grouped together)
# ---------------------------------------------------------
auNG_means_combined <- auNG %>%
  filter(method %in% c("HiFi", "HiFi+HiC", "HiFi+HiC+ONT")) %>%
  mutate(
    method_group = if_else(
      method %in% c("HiFi", "HiFi+HiC"), 
      "HiFi & HiFi+HiC (Combined)", 
      "HiFi+HiC+ONT"
    )
  ) %>%
  group_by(method_group) %>%
  summarise(
    mean_auNG = mean(v2, na.rm = TRUE),
    count_of_assemblies = n() 
  )

print(auNG_means_combined)


library(ggplot2)
library(dplyr)
library(scales)

# 1. Calculate the overall averages (pooling hap1 and hap2)
non_ref_auNG_haps <- auNG %>%
  filter(type %in% c("hap1", "hap2")) %>%
  filter(!grepl("Reference", method, ignore.case = TRUE))

avg_auNG_hifi_hic_all <- mean(non_ref_auNG_haps$v2[non_ref_auNG_haps$method %in% c("HiFi", "HiFi+HiC")], na.rm = TRUE)
avg_auNG_ont_all      <- mean(non_ref_auNG_haps$v2[non_ref_auNG_haps$method == "HiFi+HiC+ONT"], na.rm = TRUE)

# 2. Create the clean text labels for the plot (rounded to 1 decimal place)
label_hifi_hic <- paste0("Avg HiFi & HiFi+HiC: ", round(avg_auNG_hifi_hic_all / 1e6, 1), " Mbp")
label_ont      <- paste0("Avg HiFi+HiC+ONT: ", round(avg_auNG_ont_all / 1e6, 1), " Mbp")

# 3. Create the Plot
plot_auNG <- ggplot(data = auNG_wide, 
                    aes(x = hap1 / 1e6, 
                        y = hap2 / 1e6)) +
  
  # Add a 1:1 reference line (y = x)
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "grey70", linewidth = 0.6) +
  
  # --- NEW: HiFi & HiFi+HiC Average Lines (Dark Blue) ---
  geom_vline(xintercept = avg_auNG_hifi_hic_all / 1e6, linetype = "dotdash", color = "darkblue", linewidth = 0.7, alpha = 0.6) +
  geom_hline(yintercept = avg_auNG_hifi_hic_all / 1e6, linetype = "dotdash", color = "darkblue", linewidth = 0.7, alpha = 0.6) +
  
  # Annotation for HiFi & HiFi+HiC (anchored to the right side of the plot)
  annotate("text", x = Inf, y = avg_auNG_hifi_hic_all / 1e6, label = label_hifi_hic, 
           vjust = -0.6, hjust = 1.05, color = "darkblue", size = 3.5, fontface = "italic") +
  
  # --- NEW: HiFi+HiC+ONT Average Lines (Dark Red) ---
  geom_vline(xintercept = avg_auNG_ont_all / 1e6, linetype = "dotdash", color = "darkred", linewidth = 0.7, alpha = 0.6) +
  geom_hline(yintercept = avg_auNG_ont_all / 1e6, linetype = "dotdash", color = "darkred", linewidth = 0.7, alpha = 0.6) +
  
  # Annotation for ONT (anchored to the right side of the plot)
  annotate("text", x = Inf, y = avg_auNG_ont_all / 1e6, label = label_ont, 
           vjust = -0.6, hjust = 1.05, color = "darkred", size = 3.5, fontface = "italic") +
  # --------------------------------------------------------

# Add the points mapping Color to CommonName and Shape to sex
geom_point(aes(color = CommonName, shape = sex_label), size = 3.5, stroke = 1.2, alpha = 0.8) +
  
  # Assign shapes matching your Figure 1 plot
  scale_shape_manual(values = setNames(c(4, 15), label_levels), name = "Sex") +
  scale_color_manual(values = named_colors, name = "Population") +
  
  # Formatting axes
  scale_x_continuous(labels = label_comma()) +
  scale_y_continuous(labels = label_comma()) +
  
  # Labels
  xlab("H1 auNG (Mbp)") +
  ylab("H2 auNG (Mbp)") +
  
  # Theme
  theme_classic() +
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

print(plot_auNG)
ggsave("auNG_Comparison_Hap1_vs_Hap2_Mbp_extra.pdf", plot_auNG, width = 5, height = 3, dpi = 300)

library(scales)

# 1. First, create the clean text labels for the averages
# (This rounds the numbers to 1 decimal place and adds "Mbp")
label_h1 <- paste0("Avg: ", round(avg_hap1 / 1e6, 1), " Mbp")
label_h2 <- paste0("Avg: ", round(avg_hap2 / 1e6, 1), " Mbp")

# 2. Create the Plot
plot_auNG <- ggplot(data = auNG_wide, 
                    aes(x = hap1 / 1e6, 
                        y = hap2 / 1e6)) +
  
  # Add a 1:1 reference line (y = x)
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "grey70", linewidth = 0.6) +
  
  # Add the vertical and horizontal average lines
  geom_vline(xintercept = avg_hap1 / 1e6, linetype = "dashed", color = "grey30", linewidth = 0.7, alpha = 0.8) +
  geom_hline(yintercept = avg_hap2 / 1e6, linetype = "dashed", color = "grey30", linewidth = 0.7, alpha = 0.8) +
  
  # --- NEW: Write the exact numbers near the axes ---
  # H1 Average (Text near the bottom X-axis, slightly to the right of the vertical line)
  annotate("text", x = avg_hap1 / 1e6, y = -Inf, label = label_h1, 
           vjust = -1.5, hjust = -0.1, color = "grey30", size = 3.5, fontface = "italic") +
  
  # H2 Average (Text near the left Y-axis, slightly above the horizontal line)
  annotate("text", x = -Inf, y = avg_hap2 / 1e6, label = label_h2, 
           vjust = -0.5, hjust = -0.1, color = "grey30", size = 3.5, fontface = "italic") +
  # --------------------------------------------------

# Add the points mapping Color to CommonName and Shape to sex
geom_point(aes(color = CommonName, shape = sex_label), size = 3.5, stroke = 1.2, alpha = 0.8) +
  
  # Assign shapes matching your Figure 1 plot
  scale_shape_manual(values = setNames(c(4, 15), label_levels), name = "Sex") +
  scale_color_manual(values = named_colors, name = "Population") +
  
  scale_x_continuous(labels = label_comma()) +
  scale_y_continuous(labels = label_comma()) +
  
  # Labels
  xlab("H1 auNG (Mbp)") +
  ylab("H2 auNG (Mbp)") +
  
  theme_classic() +
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

print(plot_auNG)

ggsave("auNG_Comparison_Hap1_vs_Hap2_Mbp_simple.pdf", plot_auNG, width = 5, height = 3, dpi = 300)



library(ggplot2)
library(gghighlight)
library(scales)
mypallet <- c(
  "#7d7ea3",
  "#ff8d68",
  "#858058", 
  "#e9d7cb",
  "#F9AB47", 
  "#D96242", 
  "#9dced9",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",  
  "grey",
  "grey")
  
mypallet2 <- c(
  "#7d7ea3",
  "#ff8d68",
  "#858058", 
  "#e9d7cb",
  "#F9AB47", 
  "#D96242", 
  "#9dced9",
  
  "salmon",
  "#683257", 

  "#07bbc2", 
  "#026078", 
  "#3B97B6", 
  "#114155", 
  "#026078", 
  "#3B97B6", 
  "#9dced9",
  "#026078", 
  "#3B97B6",
  "#114155",
  "#9dced9",
  "#026078", 
  
  "salmon",
  "#683257", 
  
  
  "#F9AB47", 
  "#D96242",
  

  "#736f3f",
  "darkgreen",
  
  
  "#a8ca6b", 
  "forestgreen",

  "salmon",
  "#683257", 
  "#80104a",
  "#e9d7cb", 
  "#9dced9",
  "#026078", 
  "#3B97B6",
  "#114155",
  "#07bbc2", 
  "#026078", 
  "#3B97B6", 
  "#114155", 
  "#026078", 
  "#3B97B6", 
  "#9dced9",
  "#026078", 
  "#3B97B6",
  "#114155",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey",
  "grey"
)

library(dplyr)
library(ggplot2)
library(scales)

# 1. Create the new column: genomeName + (method)
# Updated line:
NGx$NewName2 <- paste0(NGx$SpeciesCode, ":", NGx$genomeName, " (", NGx$method, ")")

# 2. Define the sorting logic to order the legend
# Extract unique combinations of the necessary columns
legend_order <- NGx %>%
  select(method, NewName2) %>%
  distinct() %>%
  # Assign priority rank (1 is highest priority/top of legend)
  mutate(
    method_priority = case_when(
      grepl("Reference", method, ignore.case = TRUE) ~ 1,
      method == "HiFi+HiC+ONT" ~ 2,
      method == "HiFi+HiC" ~ 3,
      method == "HiFi" ~ 4,
      TRUE ~ 5 # Catch-all for anything else
    )
  ) %>%
  # Sort by Method Priority, then by Common Name, then by Genome Name (PR/AG ID)
  arrange(method_priority,  NewName2)

# 3. Apply the ordered levels to the NewName2 column
NGx$NewName2 <- factor(NGx$NewName2, levels = legend_order$NewName2)

# 4. Set factor levels for 'type' so linetypes are strictly assigned
NGx$type <- factor(NGx$type, levels = c("primary", "hap1", "hap2"))

# 5. Create the Plot
plot_NGx_all_colors <- ggplot(data = NGx, 
                              aes(x = v1, y = v2, 
                                  color = NewName2,       # Now an ordered factor!
                                  linetype = type,        # Solid, dotted, dashed based on type
                                  group = interaction(NewName2, type))) + # Keeps each line distinct
  
  geom_line(linewidth = 0.8, alpha = 0.9) +
  
  # Linetypes exactly as requested
  scale_linetype_manual(values = c("primary" = "solid", "hap1" = "dotted", "hap2" = "dashed")) +
  
  # Apply custom palette. Because NewName2 is a factor, colors are assigned in the new order.
  scale_color_manual(values = rep(mypallet, length.out = length(unique(NGx$NewName2)))) +
  
  # X-axis: matches the 0, 25, 50, 75, 100 layout
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0.02, 0)) +
  
  # Labels and visual theme
  xlab("NG(x)") + 
  ylab("Sequence length (bp)") +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_blank(), 
    legend.text = element_text(size = 9),
    legend.key.width = unit(1.5, "cm"), 
    
    # Creates the clean, squared-off box around the plot area
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(), 
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 12)
  )

print(plot_NGx_all_colors)
ggsave('plot_NGx_all_colors_non_redundant.pdf', plot_NGx_all_colors, width = 12, height=4, dpi=300)






library(dplyr)
library(ggplot2)
library(scales)

# 1. Calculate the averages from your NGx_50 dataframe
# Strictly EXCLUDE the reference genomes
non_ref_NG50 <- NGx_50[!grepl("Reference", NGx_50$method, ignore.case = TRUE), ]

# Average NG50 for 'HiFi' and 'HiFi+HiC' combined
avg_NG50_hifi_hic <- mean(non_ref_NG50$v2[non_ref_NG50$method %in% c("HiFi", "HiFi+HiC")], na.rm = TRUE)

# Average NG50 for 'HiFi+HiC+ONT'
avg_NG50_ont <- mean(non_ref_NG50$v2[non_ref_NG50$method == "HiFi+HiC+ONT"], na.rm = TRUE)

# Format the labels to include the exact numbers (formatted with commas)
label_hifi_hic <- paste0("Avg HiFi/HiFi+HiC: ", format(round(avg_NG50_hifi_hic), big.mark=","), " bp")
label_ont      <- paste0("Avg HiFi+HiC+ONT: ", format(round(avg_NG50_ont), big.mark=","), " bp")

# 2. Create the Plot
plot_NGx_all_colors <- ggplot(data = NGx, 
                              aes(x = v1, y = v2, 
                                  color = NewName2,       
                                  linetype = type,        
                                  group = interaction(NewName2, type))) + 
  
  geom_line(linewidth = 0.8, alpha = 0.9) +
  
  # --- AVERAGE NG50 LINES & ANNOTATIONS ---
  # Vertical dashed line at NG50 (x = 50)
  geom_vline(xintercept = 50, linetype = "dashed", color = "grey30", linewidth = 0.6) +
  
  # Horizontal line + Annotation for HiFi & HiFi+HiC
  geom_hline(yintercept = avg_NG50_hifi_hic, linetype = "dotdash", color = "darkblue", alpha = 0.7) +
  annotate("text", x = 52, y = avg_NG50_hifi_hic, 
           label = label_hifi_hic, # Uses the pasted string with the number
           color = "darkblue", hjust = 0, vjust = -0.5, size = 3.5) +
  
  # Horizontal line + Annotation for HiFi+HiC+ONT
  geom_hline(yintercept = avg_NG50_ont, linetype = "dotdash", color = "darkred", alpha = 0.7) +
  annotate("text", x = 52, y = avg_NG50_ont, 
           label = label_ont,      # Uses the pasted string with the number
           color = "darkred", hjust = 0, vjust = -0.5, size = 3.5) +
  # ----------------------------------------

# Linetypes exactly as requested
scale_linetype_manual(values = c("primary" = "solid", "hap1" = "dotted", "hap2" = "dashed")) +
  
  # Apply custom palette
  scale_color_manual(values = rep(mypallet, length.out = length(unique(NGx$NewName2)))) +
  
  # X-axis: matches the 0, 25, 50, 75, 100 layout
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0.02, 0)) +
  
  # Y-axis: Linear Scale (Replaced log10 scale)
  # Added comma formatting to the Y-axis so the large numbers are readable
  #scale_y_continuous(
  #  limits = c(0, NA), # Starts at 0 to ground the linear scale naturally
  #  labels = label_comma()
  #) +
  
  # Labels and visual theme
  xlab("NGx (%)") + 
  ylab("Sequence length (bp)") +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_blank(), 
    legend.text = element_text(size = 9),
    legend.key.width = unit(1.5, "cm"), 
    
    # Creates the clean, squared-off box around the plot area
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(), 
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 12)
  )

print(plot_NGx_all_colors)

# Save the plot
ggsave('plot_NGx_all_colors_linear.pdf', plot_NGx_all_colors, width = 12, height = 4, dpi = 300)






