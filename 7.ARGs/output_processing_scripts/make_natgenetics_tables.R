library(tidyverse)

setwd("/Users/joanocha/Desktop/singer_tmp")
outdir <- "new_output"

# ── 1. Load data ──────────────────────────────────────────────────────────────

# Chimp wide: per-pop (Central, Eastern, Western)
chimp <- read_csv("SuppTable_chimp_wide_genes_Twithin_6MYA.csv",
                  show_col_types = FALSE) %>%
  mutate(display = trimws(genes_renamed)) %>%
  group_by(display) %>%
  slice_max(num_populations, n = 1, with_ties = FALSE) %>%
  ungroup()

# Bonobo long: one row per window → collapse to per-gene max
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

# Human NatGenetics WIDE: AFR single population, use as authoritative gene list
human <- read_csv("human_outliers_gMYA_NatGenetics_WIDE.csv",
                  show_col_types = FALSE) %>%
  rename(display = gene_name, AFR_hsa = max_tmrca_mya)

# ── 2. Gene sets ──────────────────────────────────────────────────────────────

chimp_genes  <- unique(chimp$display)
bonobo_genes <- unique(bonobo$display)
human_genes  <- unique(human$display)

all_genes <- union(union(chimp_genes, bonobo_genes), human_genes)

# ── 3. Population-sharing helpers ─────────────────────────────────────────────

chimp_pop_label <- function(gene) {
  row <- chimp %>% filter(display == gene) %>% slice(1)
  if (nrow(row) == 0) return(NA_character_)
  pops <- c(
    if (!is.na(row$Central[1])) "Central",
    if (!is.na(row$Eastern[1])) "Eastern",
    if (!is.na(row$Western[1])) "Western"
  )
  if (length(pops) == 3) "all_chimp_pops"
  else if (length(pops) == 0) NA_character_
  else paste(pops, collapse = "+")
}

# ── 4. Build summary table ────────────────────────────────────────────────────

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
  left_join(
    chimp %>% select(display, descriptions,
                     chimp_n_pops = num_populations,
                     Central_Pt = Central, Eastern_Pt = Eastern, Western_Pt = Western),
    by = "display"
  ) %>%
  left_join(
    bonobo %>% select(display, descriptions_b = descriptions,
                      gene_biotypes_b = gene_biotypes, Bonobo_Ppa),
    by = "display"
  ) %>%
  left_join(
    human %>% select(display, gene_type, AFR_hsa),
    by = "display"
  ) %>%
  mutate(
    descriptions  = coalesce(descriptions, descriptions_b),
    gene_biotypes = coalesce(gene_biotypes_b, gene_type),
    chimp_n_pops  = replace_na(as.integer(chimp_n_pops), 0L),
    chimp_pop_sharing = map_chr(display, chimp_pop_label)
  ) %>%
  select(
    genes_renamed     = display,
    descriptions,
    gene_biotypes,
    n_species,
    species_sharing,
    chimp_n_pops,
    chimp_pop_sharing,
    Bonobo_Ppa,
    Central_Pt,
    Eastern_Pt,
    Western_Pt,
    AFR_hsa
  ) %>%
  arrange(desc(n_species), genes_renamed)

write_csv(summary_tbl,
          file.path(outdir, "SuppTable_species+pops_overlap_summary_6MYA_NatGenetics.csv"))
message("Summary table: ", nrow(summary_tbl), " genes")

# ── 5. Build counts table ─────────────────────────────────────────────────────

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
    category          = cat_label,
    n_genes           = n,
    genes             = NA_character_,
    chimp_pop_sharing = NA_character_,
    human_pop_sharing = NA_character_
  )

  for (i in seq_len(n)) {
    g <- genes_in_cat[i, ]
    counts_rows[[length(counts_rows) + 1]] <- tibble(
      category          = NA_character_,
      n_genes           = NA_integer_,
      genes             = g$genes_renamed,
      chimp_pop_sharing = if_else(g$chimp_n_pops > 0, g$chimp_pop_sharing, NA_character_),
      human_pop_sharing = if_else(!is.na(g$AFR_hsa), "AFR", NA_character_)
    )
  }
}

counts_tbl <- bind_rows(counts_rows)
write_csv(counts_tbl,
          file.path(outdir, "SuppTable_species+pops_overlap_counts_6MYA_NatGenetics.csv"),
          na = "")
message("Counts table: ", nrow(counts_tbl), " rows")

# ── 6. Summary ────────────────────────────────────────────────────────────────
cat("\nGenes per category:\n")
summary_tbl %>% count(species_sharing) %>% arrange(desc(n)) %>% print()

cat("\nShared genes (n_species > 1):\n")
summary_tbl %>%
  filter(n_species > 1) %>%
  select(genes_renamed, species_sharing, chimp_pop_sharing, Bonobo_Ppa, AFR_hsa) %>%
  print(n = 50)
