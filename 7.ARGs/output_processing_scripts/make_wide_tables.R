library(tidyverse)

setwd("/Users/joanocha/Desktop/singer_tmp")

# ── 1. Load data ──────────────────────────────────────────────────────────────

chimp <- read_csv("SuppTable_chimp_wide_genes_Twithin_6MYA.csv", show_col_types = FALSE) %>%
  mutate(display = trimws(genes_renamed)) %>%
  # When multiple original IDs rename to the same symbol, keep the one with most pops
  # then highest Central value as tiebreak
  group_by(display) %>%
  slice_max(num_populations, n = 1, with_ties = FALSE) %>%
  ungroup()

human <- read_csv("SuppTable_humans_wide_6myt_genes_Twithin.csv", show_col_types = FALSE) %>%
  mutate(display = trimws(genes))

# Bonobo is window-level: one gene per row, collapse to per-gene max
bonobo_raw <- read_csv("SuppTable_avg_tmrca_Singer_bonobo.csv",
                       show_col_types = FALSE, name_repair = "universal") %>%
  rename_with(~ str_remove(.x, "^\\."), starts_with("."))  # strip BOM dot

bonobo <- bonobo_raw %>%
  mutate(display = trimws(genes_renamed)) %>%
  group_by(display) %>%
  summarise(
    descriptions  = first(descriptions),
    gene_biotypes = first(gene_biotypes),
    Bonobo_Ppa    = max(avg_pairwise_myr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(display), display != "", display != ".")

# ── 2. Build per-species presence sets ───────────────────────────────────────

chimp_genes  <- unique(chimp$display)
human_genes  <- unique(human$display)
bonobo_genes <- unique(bonobo$display)

all_genes <- union(union(chimp_genes, human_genes), bonobo_genes)

# ── 3. Population-sharing helpers ────────────────────────────────────────────

chimp_pop_label <- function(gene) {
  row <- chimp %>% filter(display == gene) %>% slice(1)
  if (nrow(row) == 0) return(NA_character_)
  pops <- c(
    if (!is.na(row$Central[1])) "Central",
    if (!is.na(row$Eastern[1])) "Eastern",
    if (!is.na(row$Western[1])) "Western"
  )
  if (length(pops) == 3) "all_chimp_pops" else paste(pops, collapse = "+")
}

human_pop_label <- function(gene) {
  row <- human %>% filter(display == gene) %>% slice(1)
  if (nrow(row) == 0) return(NA_character_)
  pops <- c(
    if (!is.na(row$AFR[1])) "AFR",
    if (!is.na(row$EAS[1])) "EAS",
    if (!is.na(row$EUR[1])) "EUR",
    if (!is.na(row$SAS[1])) "SAS",
    if (!is.na(row$AMR[1])) "AMR"
  )
  if (length(pops) == 5) "all_human_pops" else paste(pops, collapse = "+")
}

# ── 4. Build summary table ───────────────────────────────────────────────────

summary_tbl <- tibble(display = all_genes) %>%
  mutate(
    in_chimp  = display %in% chimp_genes,
    in_human  = display %in% human_genes,
    in_bonobo = display %in% bonobo_genes,
    n_species = in_chimp + in_human + in_bonobo,
    species_sharing = case_when(
      in_chimp & in_human & in_bonobo ~ "all_three",
      in_chimp & in_bonobo            ~ "bonobo_chimp",
      in_human & in_bonobo            ~ "bonobo_human",
      in_chimp & in_human             ~ "chimp_human",
      in_bonobo                       ~ "bonobo_only",
      in_chimp                        ~ "chimp_only",
      TRUE                            ~ "human_only"
    )
  ) %>%
  # Join values from each species
  left_join(
    chimp %>% select(display, descriptions, gene_biotypes,
                     chimp_n_pops = num_populations,
                     Central_Pt = Central, Eastern_Pt = Eastern, Western_Pt = Western),
    by = "display"
  ) %>%
  left_join(
    human %>% select(display, descriptions_h = descriptions,
                     human_n_pops = num_populations,
                     AFR_hsa = AFR, EAS_hsa = EAS, EUR_hsa = EUR,
                     SAS_hsa = SAS, AMR_hsa = AMR),
    by = "display"
  ) %>%
  left_join(
    bonobo %>% select(display, descriptions_b = descriptions,
                      gene_biotypes_b = gene_biotypes, Bonobo_Ppa),
    by = "display"
  ) %>%
  mutate(
    descriptions  = coalesce(descriptions, descriptions_h, descriptions_b),
    gene_biotypes = coalesce(gene_biotypes, gene_biotypes_b),
    chimp_n_pops  = replace_na(as.integer(chimp_n_pops), 0L),
    human_n_pops  = replace_na(as.integer(human_n_pops), 0L),
    chimp_pop_sharing = map_chr(display, chimp_pop_label),
    human_pop_sharing = map_chr(display, human_pop_label)
  ) %>%
  select(
    genes_renamed    = display,
    descriptions,
    gene_biotypes,
    n_species,
    species_sharing,
    chimp_n_pops,
    chimp_pop_sharing,
    human_n_pops,
    human_pop_sharing,
    Bonobo_Ppa,
    Eastern_Pt,
    Western_Pt,
    Central_Pt,
    AFR_hsa,
    EAS_hsa,
    EUR_hsa,
    SAS_hsa,
    AMR_hsa
  ) %>%
  arrange(desc(n_species), genes_renamed)

write_csv(summary_tbl, "SuppTable_species+pops_overlap_summary_6MYA.csv")
message("Summary table: ", nrow(summary_tbl), " genes")

# ── 5. Build counts table ─────────────────────────────────────────────────────

category_order <- c(
  "All three species"   = "all_three",
  "Bonobo + Chimp only" = "bonobo_chimp",
  "Bonobo + Human only" = "bonobo_human",
  "Chimp + Human only"  = "chimp_human",
  "Bonobo only"         = "bonobo_only",
  "Chimp only"          = "chimp_only",
  "Human only"          = "human_only"
)

counts_rows <- list()

for (cat_label in names(category_order)) {
  cat_key <- category_order[[cat_label]]
  genes_in_cat <- summary_tbl %>%
    filter(species_sharing == cat_key) %>%
    arrange(genes_renamed)

  n <- nrow(genes_in_cat)

  # Header row for this category
  counts_rows[[length(counts_rows) + 1]] <- tibble(
    category      = cat_label,
    n_genes       = n,
    genes         = NA_character_,
    chimp_pop_sharing = NA_character_,
    human_pop_sharing = NA_character_
  )

  # One row per gene
  for (i in seq_len(n)) {
    g <- genes_in_cat[i, ]
    counts_rows[[length(counts_rows) + 1]] <- tibble(
      category      = NA_character_,
      n_genes       = NA_integer_,
      genes         = g$genes_renamed,
      chimp_pop_sharing = if_else(g$chimp_n_pops > 0, g$chimp_pop_sharing, NA_character_),
      human_pop_sharing = if_else(g$human_n_pops > 0, g$human_pop_sharing, NA_character_)
    )
  }
}

counts_tbl <- bind_rows(counts_rows)
write_csv(counts_tbl, "SuppTable_species+pops_overlap_counts_6MYA.csv", na = "")
message("Counts table: ", nrow(counts_tbl), " rows")

# ── 6. Quick check ────────────────────────────────────────────────────────────
cat("\nGenes per category:\n")
summary_tbl %>% count(species_sharing) %>% arrange(desc(n)) %>% print()

cat("\nSample summary rows (shared genes):\n")
summary_tbl %>%
  filter(n_species > 1) %>%
  select(genes_renamed, species_sharing, chimp_pop_sharing, human_pop_sharing,
         Bonobo_Ppa, Eastern_Pt, Central_Pt) %>%
  head(10) %>%
  print()
