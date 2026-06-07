setwd("/Users/joanocha/Google Drive/My Drive/POSTDOC/PANPAN/analysis/Figure4_SINGER/Yun_regions/")
library(dplyr)
library(readr)
library(tidyverse)



original_data <- read_tsv("old_regions_overlap.txt", col_names = FALSE)

# Give the columns meaningful names for easier use
# Columns 4 and 13-15 seem to be the key ones
colnames(original_data)[c(1:4, 13:15)] <- c("gene_chr", "gene_start", "gene_end", "ensembl_transcript_id_version", "window_chr", "window_start", "window_end")

# Create a clean version of the transcript ID without the version number
original_data <- original_data %>%
  mutate(ensembl_transcript_id = sub("\\..*$", "", ensembl_transcript_id_version))

# Get a unique list of the IDs you need to look up
transcript_ids_to_query <- unique(original_data$ensembl_transcript_id)


gencode47<-read_tsv("gencode.v47.basic.annotation.gtf", comment = "##", col_names = FALSE)
View(gencode47)
gene_map <- gencode47 %>%
  filter(X3 == "transcript") %>%
  mutate(
    ensembl_transcript_id_version = str_extract(X9, 'transcript_id "([^"]+)"', group = 1),
    ensembl_gene_id_version = str_extract(X9, 'gene_id "([^"]+)"', group = 1),
    gene_name = str_extract(X9, 'gene_name "([^"]+)"', group = 1),
    # **NEW**: Also extract the gene_type
    gene_type = str_extract(X9, 'gene_type "([^"]+)"', group = 1)
  ) %>%
  # Remove any rows where the extraction might have failed
  filter(!is.na(ensembl_transcript_id_version), !is.na(gene_name)) %>%
  
  # Select all four columns
  dplyr::select(ensembl_transcript_id_version, ensembl_gene_id_version, gene_name, gene_type) %>%
  
  distinct()

# View the result to confirm it worked
View(gene_map)
                       

final_annotated_data <- original_data %>%
  left_join(gene_map, by = "ensembl_transcript_id_version")  %>% 
  dplyr::select(last_col(offset = 6):last_col()) %>%
  distinct(window_chr, window_start, window_end, ensembl_gene_id_version, .keep_all = TRUE)
View(final_annotated_data)

diversity_data <- read_tsv("filtered_diversity_scan.txt", col_names = FALSE)
diversity_df_prepared <- diversity_data %>%
  rename(
    window_start = X1,
    avg_tmrca_generations = X2,
    chromosome_num = X3
  ) %>%
  mutate(window_start = as.integer(window_start)) %>%
  mutate(window_chr = paste0("chr", chromosome_num)) %>%
  dplyr::select(window_chr, window_start, avg_tmrca_generations)
diversity_df_prepared $avg_tmrca_myr <- diversity_df_prepared $avg_tmrca_generations * 28 / 1e6
View(diversity_df_prepared)


final_merged_data <- left_join(
  final_annotated_data, 
  diversity_df_prepared, 
  by = c("window_chr", "window_start")
) %>%  # <-- Make sure this is the pipe operator %>%
  
  # Arrange the results in descending order
  arrange(desc(avg_tmrca_myr))

# View the final, sorted table
View(final_merged_data)


gene_summary_wide <- final_merged_data %>%
  filter(!is.na(gene_name)) %>%
  group_by(gene_name, gene_type) %>%
  summarise(
    max_tmrca_generations = max(avg_tmrca_generations, na.rm = TRUE),
    max_tmrca_mya = max(avg_tmrca_myr, na.rm = TRUE)
  ) %>%
  arrange(desc(max_tmrca_mya))

# View the new wide-format table
View(gene_summary_wide)


unique_gene_list_filtered <- final_merged_data %>%
  # First, get the unique, non-NA gene names
  filter(!is.na(gene_name)) %>%
  distinct(gene_name) %>%
  
  # **NEW**: Filter out genes containing "LIN", "LOC", or "ENSG"
  filter(!str_detect(gene_name, "LIN|LOC|ENSG")) %>%
  
  # Pull the final list into a vector
  pull(gene_name)

# Print the first few gene names to the console to check
print(head(unique_gene_list_filtered))
write_lines(unique_gene_list_filtered, "unique_gene_names_filtered.txt")

write_csv(final_merged_data,"Africans_1kGP_high_pairwiseTMRCA_NatGenetics.csv") 
write_csv(gene_summary_wide,"Africans_1kGP_high_pairwiseTMRCA_NatGenetics_WIDE.csv") 

