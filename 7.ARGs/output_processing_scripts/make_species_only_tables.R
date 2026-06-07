library(tidyverse)

setwd("/Users/joanocha/Desktop/singer_tmp")
outdir <- "new_output"

# ── 1. Load data ──────────────────────────────────────────────────────────────

chimp <- read_csv("SuppTable_avg_tmrca_Singer_chimp.csv",
                  show_col_types = FALSE, name_repair = "universal") %>%
  mutate(display = trimws(genes_renamed)) %>%
  filter(!is.na(display), display != "", display != ".") %>%
  group_by(display) %>%
  summarise(
    descriptions = first(descriptions),
    Chimp_Pt     = max(avg_pairwise_myr, na.rm = TRUE),
    .groups = "drop"
  )

bonobo <- read_csv("SuppTable_avg_tmrca_Singer_bonobo.csv",
                   show_col_types = FALSE, name_repair = "universal") %>%
  mutate(display = trimws(genes_renamed)) %>%
  filter(!is.na(display), display != "", display != ".") %>%
  group_by(display) %>%
  summarise(
    descriptions  = first(descriptions),
    gene_biotypes = first(gene_biotypes),
    Bonobo_Ppa    = max(avg_pairwise_myr, na.rm = TRUE),
    .groups = "drop"
  )

human <- read_csv("human_outliers_gMYA_NatGenetics_WIDE.csv",
                  show_col_types = FALSE) %>%
  rename(display = gene_name, AFR_hsa = max_tmrca_mya)

# ── 2. Gene sets & species sharing ───────────────────────────────────────────

chimp_genes  <- unique(chimp$display)
bonobo_genes <- unique(bonobo$display)
human_genes  <- unique(human$display)

all_genes <- union(union(chimp_genes, bonobo_genes), human_genes)

summary_tbl <- tibble(display = all_genes) %>%
  mutate(
    in_chimp  = display %in% chimp_genes,
    in_bonobo = display %in% bonobo_genes,
    in_human  = display %in% human_genes,
    n_species = in_chimp + in_bonobo + in_human,
    species_sharing = case_when(
      in_chimp & in_bonobo & in_human ~ "all_three",
      in_chimp & in_bonobo            ~ "bonobo_chimp",
      in_bonobo & in_human            ~ "bonobo_human",
      in_chimp & in_human             ~ "chimp_human",
      in_bonobo                       ~ "bonobo_only",
      in_chimp                        ~ "chimp_only",
      TRUE                            ~ "human_only"
    )
  ) %>%
  left_join(chimp  %>% select(display, descriptions, Chimp_Pt), by = "display") %>%
  left_join(bonobo %>% select(display, descriptions_b = descriptions,
                               gene_biotypes, Bonobo_Ppa),       by = "display") %>%
  left_join(human  %>% select(display, gene_type, AFR_hsa),      by = "display") %>%
  mutate(
    descriptions  = coalesce(descriptions, descriptions_b),
    gene_biotypes = coalesce(gene_biotypes, gene_type)
  ) %>%
  select(
    genes_renamed   = display,
    descriptions,
    gene_biotypes,
    n_species,
    species_sharing,
    Bonobo_Ppa,
    Chimp_Pt,
    AFR_hsa
  ) %>%
  arrange(desc(n_species), genes_renamed)

write_csv(summary_tbl,
          file.path(outdir, "SuppTable_species_overlap_summary_6MYA_NatGenetics.csv"))
message("Summary: ", nrow(summary_tbl), " genes")

# ── 3. Counts table ───────────────────────────────────────────────────────────

category_order <- c(
  "All three species"   = "all_three",
  "Bonobo + Chimp only" = "bonobo_chimp",
  "Bonobo + Human only" = "bonobo_human",
  "Chimp + Human only"  = "chimp_human",
  "Bonobo only"         = "bonobo_only",
  "Chimp only"          = "chimp_only",
  "Human only (AFR)"    = "human_only"
)

counts_rows <- list()

for (cat_label in names(category_order)) {
  cat_key <- category_order[[cat_label]]
  genes_in_cat <- summary_tbl %>%
    filter(species_sharing == cat_key) %>%
    arrange(genes_renamed)
  n <- nrow(genes_in_cat)

  counts_rows[[length(counts_rows) + 1]] <- tibble(
    category = cat_label, n_genes = n, genes = NA_character_
  )
  for (i in seq_len(n)) {
    counts_rows[[length(counts_rows) + 1]] <- tibble(
      category = NA_character_, n_genes = NA_integer_,
      genes = genes_in_cat$genes_renamed[i]
    )
  }
}

counts_tbl <- bind_rows(counts_rows)
write_csv(counts_tbl,
          file.path(outdir, "SuppTable_species_overlap_counts_6MYA_NatGenetics.csv"),
          na = "")
message("Counts: ", nrow(counts_tbl), " rows")

# ── 4. Summary ────────────────────────────────────────────────────────────────
cat("\nGenes per category:\n")
summary_tbl %>% count(species_sharing) %>% arrange(desc(n)) %>% print()

cat("\nShared genes:\n")
summary_tbl %>%
  filter(n_species > 1) %>%
  select(genes_renamed, species_sharing, Bonobo_Ppa, Chimp_Pt, AFR_hsa) %>%
  print(n = 50)
