suppressPackageStartupMessages(library(tidyverse))

gi <- read_tsv("GO_analysis/backgrounds/Homo_sapiens.gene_info.gz",
               show_col_types = FALSE, na = c("", "-")) %>%
  rename(tax_id = `#tax_id`) %>%
  mutate(best_symbol = coalesce(Symbol_from_nomenclature_authority, Symbol)) %>%
  filter(!str_detect(best_symbol, "^LOC"))

norm <- function(x) str_to_lower(str_trim(x))
pipe_split <- fixed("|")

map_desc <- gi %>%
  mutate(k = norm(description)) %>% filter(!is.na(k)) %>%
  distinct(k, .keep_all = TRUE) %>% select(k, sym = best_symbol)

map_fullname <- gi %>%
  mutate(k = norm(Full_name_from_nomenclature_authority)) %>% filter(!is.na(k)) %>%
  distinct(k, .keep_all = TRUE) %>% select(k, sym = best_symbol)

map_other <- gi %>%
  filter(!is.na(Other_designations)) %>%
  separate_rows(Other_designations, sep = pipe_split) %>%
  mutate(k = norm(Other_designations)) %>% filter(k != "") %>%
  distinct(k, .keep_all = TRUE) %>% select(k, sym = best_symbol)

cat("Lookup map sizes:  desc=", nrow(map_desc),
    "  full_name=", nrow(map_fullname),
    "  other_designations=", nrow(map_other), "\n\n")

wide <- read_tsv("gene_impact_wide.tsv", show_col_types = FALSE)
loc_wide <- wide %>% filter(str_detect(name, "^LOC"))

resolve_one <- function(desc) {
  if (is.na(desc) || desc == "") return(NA_character_)
  k0 <- norm(desc)
  k1 <- str_trim(str_remove(k0, "-like( protein.*)?$"))
  k2 <- str_trim(str_remove(k0, " protein.*$"))
  keys <- unique(c(k0, k1, k2))
  for (k in keys) {
    m <- map_desc$sym[match(k, map_desc$k)]
    if (!is.na(m)) return(m)
    m <- map_fullname$sym[match(k, map_fullname$k)]
    if (!is.na(m)) return(m)
    m <- map_other$sym[match(k, map_other$k)]
    if (!is.na(m)) return(m)
  }
  NA_character_
}

loc_wide$new_symbol <- vapply(loc_wide$description, resolve_one, character(1))

cat("Total LOC* rows:", nrow(loc_wide), "\n")
cat("Successfully mapped:", sum(!is.na(loc_wide$new_symbol)),
    " (", round(100 * mean(!is.na(loc_wide$new_symbol)), 1), "%)\n")
cat("Unmapped:", sum(is.na(loc_wide$new_symbol)), "\n\n")

cat("Sample successful renames:\n")
loc_wide %>% filter(!is.na(new_symbol)) %>%
  select(LOC = name, NewSymbol = new_symbol, description) %>%
  slice_sample(n = 15) %>% print()

write_tsv(loc_wide %>% select(loc_name = name, new_symbol, description),
          "GO_analysis/backgrounds/loc_to_hgnc_map.tsv")
cat("\nWrote: GO_analysis/backgrounds/loc_to_hgnc_map.tsv\n")
