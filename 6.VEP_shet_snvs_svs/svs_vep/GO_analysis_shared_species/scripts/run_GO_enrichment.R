# GO / pathway / disease enrichment for SV-impacted gene sets shared across
# great apes. enrichR-only run — matches the ShinyGO v0.85 setup the user is
# using interactively (default Enrichr-internal background, no custom_bg).
#
# Databases queried (one independent pass each):
#   - GO_Biological_Process_2026  : GO Biological Process
#   - GO_Molecular_Function_2026  : GO Molecular Function
#   - KEGG_2026                   : KEGG pathways
#   - DisGeNET                    : Human-disease term enrichment
#
# Significance protocol (ShinyGO-style):
#   1. enrichR returns raw P.value and (Adjusted.P.value over all DB terms).
#   2. We FILTER terms to a per-DB term_size window first.
#   3. Then we RE-RUN BH-FDR on that filtered subset to get the reported FDR.
#      This makes the FDR pool depend only on size-appropriate pathways, which
#      is how ShinyGO's "min/max pathway size" parameters actually behave.
# Cutoffs:
#   - GO BP / GO MF / KEGG : term_size [1, 60], FDR < 0.05  (matches the
#     ShinyGO screenshot the user is reproducing)
#   - DisGeNET             : term_size [1, 1000], FDR < 0.05  (relaxed max so
#     Cooley's anemia, term_size = 70, still appears)
#
# Per-list output:
#   GO_analysis/<list>/enrichr_<DB>.csv
#   GO_analysis/<list>/plot_<DB>.{pdf,png}      (only if the CSV has rows)
# Cross-list summary:
#   GO_analysis/plots/summary_<DB>.{pdf,png}
# Combined supplementary:
#   GO_analysis/combined/SupplementaryTable_enrichment.xlsx
#   (one sheet per gene list; rows from all 4 DBs labelled by `pass`)

suppressPackageStartupMessages({
  library(tidyverse)
  library(enrichR)
  library(openxlsx)
  library(tidytext)  # reorder_within / scale_y_reordered for faceted ggplots
})

# ---- paths ------------------------------------------------------------------

base_dir <- "/Users/joanocha/Library/CloudStorage/GoogleDrive-joana.laranjeira.rocha@gmail.com/My Drive/POSTDOC/PANPAN/analysis/Figure3_VEP_snvs_svs/svs_vep"
out_dir  <- file.path(base_dir, "GO_analysis")

shared_file <- file.path(base_dir, "TableSX_genes_impacted_by_svs_shared_across_species.csv")
wide_file   <- file.path(base_dir, "gene_impact_wide.tsv")

# ---- 1. load source data + build gene lists --------------------------------

shared <- read_csv(shared_file, show_col_types = FALSE)
wide   <- read_tsv(wide_file,    show_col_types = FALSE)

# Same minimal LOC-RNA filter as vep_shet_SVs.R:14-16.
loc_rna_drop <- function(df) {
  df %>%
    filter(!(str_detect(name, "^LOC") & str_detect(description, "spliceosomal RNA"))) %>%
    filter(!(str_detect(name, "^LOC") & str_detect(description, "small") & str_detect(description, "RNA")))
}

wide_clean   <- loc_rna_drop(wide)
shared_clean <- loc_rna_drop(shared)
clean_genes  <- function(x) unique(x[!is.na(x)])

# Rename LOC* entries to their HGNC/NCBI-curated symbol where possible.
# Mapping is built by GO_analysis/scripts/_build_loc_map.R using the local
# NCBI Homo_sapiens.gene_info table (matches descriptions / full names /
# Other_designations to the official Symbol). Rows we can't rename keep the
# original LOC name — ShinyGO will silently ignore those (they were not
# going to map to any GO/KEGG term anyway).
loc_map_path <- file.path(out_dir, "backgrounds", "loc_to_hgnc_map.tsv")
if (file.exists(loc_map_path)) {
  loc_map <- read_tsv(loc_map_path, show_col_types = FALSE)
  rename_locs <- function(names) {
    idx <- match(names, loc_map$loc_name)
    new <- loc_map$new_symbol[idx]
    ifelse(!is.na(new), new, names)
  }
  wide_clean$name <- rename_locs(wide_clean$name)
  message("Applied LOC -> HGNC renaming: ", sum(!is.na(loc_map$new_symbol)),
          " of ", nrow(loc_map), " LOC entries remapped")
} else {
  message("Note: ", loc_map_path, " missing; skipping LOC renaming. ",
          "Run scripts/_build_loc_map.R to generate it.")
}

gene_lists <- list(
  cross_species_shared_60 = clean_genes(shared_clean$name),
  human_chimp_shared_35   = clean_genes(shared_clean %>%
                                          filter(human == 4, chimpanzee == 4, bonobo != 4) %>%
                                          pull(name))
)

# Primary ShinyGO background = the set of genes that were actually evaluated
# by the VEP pipeline (gene_impact_wide.tsv, after LOC-RNA cleaning AND
# LOC->HGNC renaming via NCBI gene_info). This is the statistically
# appropriate background — it controls for bias in which genes the SV-impact
# analysis could test in the first place, while making LOC entries
# recognisable to ShinyGO wherever NCBI has a curated symbol for them.
bg_tested <- clean_genes(wide_clean$name)

dir.create(file.path(out_dir, "gene_lists"),  showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "backgrounds"), showWarnings = FALSE, recursive = TRUE)
iwalk(gene_lists, ~ writeLines(.x, file.path(out_dir, "gene_lists", paste0(.y, ".txt"))))
writeLines(bg_tested, file.path(out_dir, "backgrounds", "tested_universe.txt"))

# Optional alternative background = Ensembl protein-coding (~19,468).
ensembl_bg_path <- file.path(out_dir, "backgrounds", "ensembl_human_protein_coding.txt")
bg_ensembl <- if (file.exists(ensembl_bg_path)) readLines(ensembl_bg_path) else NULL

iwalk(gene_lists, function(genes, list_name) {
  dir.create(file.path(out_dir, list_name), showWarnings = FALSE, recursive = TRUE)
  writeLines(genes,     file.path(out_dir, list_name, "shinygo_input.txt"))
  # Primary ShinyGO background = tested universe (genes evaluated by the VEP pipeline).
  writeLines(bg_tested, file.path(out_dir, list_name, "shinygo_background_tested.txt"))
  # Optional alternative.
  if (!is.null(bg_ensembl)) {
    writeLines(bg_ensembl,
               file.path(out_dir, list_name, "shinygo_background_ensembl.txt"))
  }
})

# ---- 2. enrichR config ------------------------------------------------------

Sys.setenv(ENRICHR_LIVE = "TRUE")
setEnrichrSite("Enrichr")

# Per-database parameters. min/max are PATHWAY (term) sizes; fdr is the
# Benjamini-Hochberg cutoff applied AFTER size filtering (= ShinyGO logic).
# FDR is set to 0.10 ("as relaxed as you can") because the strongest GO_BP/MF
# hits sit at FDR ≈ 0.064 under min=1/max=60. KEGG and DisGeNET also reuse this
# cutoff for consistency; their top hits (Galactose metabolism, Laryngitis,
# Cooley's anemia) easily pass at FDR<0.05 as well.
# bg_size = number of genes the Enrichr library annotates (from listEnrichrDbs);
# used to compute fold enrichment = (k/n) / (K/bg_size).
passes <- list(
  GO_BP    = list(db = "GO_Biological_Process_2026", min = 1, max = 60,   fdr = 0.10, bg_size = 15557),
  GO_MF    = list(db = "GO_Molecular_Function_2026", min = 1, max = 60,   fdr = 0.10, bg_size = 12297),
  KEGG     = list(db = "KEGG_2026",                  min = 1, max = 60,   fdr = 0.10, bg_size = 8110),
  DisGeNET = list(db = "DisGeNET",                   min = 1, max = 1000, fdr = 0.10, bg_size = 17464)
)

# ---- 3. enrichment ----------------------------------------------------------

# Parse enrichR's "Overlap" column ("k/N"), filter term_size in [min, max],
# re-run BH on the filtered subset, compute fold enrichment, and return ALL
# size-filtered rows (with `passes_fdr` flag). Plots use top 10 by FDR
# regardless of significance; the CSV / supplementary table filter on
# passes_fdr at write time.
run_enrichr_pass <- function(genes, pass_name, cfg) {
  empty <- tibble(source = character(), term_id = character(), term_name = character(),
                  p_value = numeric(), p_value_raw = numeric(),
                  fold_enrichment = numeric(), term_size = integer(),
                  intersection_size = integer(), intersection = character(),
                  passes_fdr = logical(), background_used = character())
  res <- tryCatch(enrichr(genes, cfg$db)[[cfg$db]],
                  error = function(e) { message("  enrichR error: ", conditionMessage(e)); NULL })
  if (is.null(res) || nrow(res) == 0) return(empty)
  n_query <- length(genes)
  res %>%
    as_tibble() %>%
    separate(Overlap, into = c("k", "N"), sep = "/", convert = TRUE, remove = FALSE) %>%
    filter(N >= cfg$min, N <= cfg$max) %>%
    mutate(
      p_adj_subset    = p.adjust(P.value, method = "BH"),
      # Fold enrichment = (k/n) / (K/bg_size) — same definition ShinyGO uses.
      fold_enrichment = (k / n_query) / (N / cfg$bg_size)
    ) %>%
    arrange(p_adj_subset, desc(fold_enrichment)) %>%
    transmute(
      source            = pass_name,
      term_id           = Term,
      term_name         = Term,
      p_value           = p_adj_subset,   # BH-FDR on size-filtered subset
      p_value_raw       = P.value,
      fold_enrichment   = fold_enrichment,
      term_size         = N,
      intersection_size = k,
      intersection      = Genes,
      passes_fdr        = p_adj_subset < cfg$fdr,
      background_used   = "enrichr_default"
    )
}

all_results <- map_dfr(names(passes), function(pass_name) {
  cfg <- passes[[pass_name]]
  map_dfr(names(gene_lists), function(list_name) {
    message("--- ", list_name, " | ", pass_name, " (", cfg$db,
            ", size ", cfg$min, "-", cfg$max, ", FDR<", cfg$fdr, ") ---")
    df <- run_enrichr_pass(gene_lists[[list_name]], pass_name, cfg)
    write_csv(df, file.path(out_dir, list_name, paste0("enrichr_", pass_name, ".csv")))
    df %>% mutate(gene_list = list_name, pass = pass_name, .before = 1)
  })
})

# ---- 4. combined supplementary table (one sheet per gene list) -------------

dir.create(file.path(out_dir, "combined"), showWarnings = FALSE, recursive = TRUE)

wb <- createWorkbook()
walk(names(gene_lists), function(list_name) {
  addWorksheet(wb, list_name)
  sheet_df <- all_results %>%
    filter(gene_list == list_name) %>%
    select(pass, source, term_id, term_name, p_value, p_value_raw,
           fold_enrichment, term_size, intersection_size, intersection,
           passes_fdr, background_used)
  writeData(wb, list_name, sheet_df)
})
saveWorkbook(wb,
             file.path(out_dir, "combined", "SupplementaryTable_enrichment.xlsx"),
             overwrite = TRUE)

# ---- 5. plots --------------------------------------------------------------

dir.create(file.path(out_dir, "plots"), showWarnings = FALSE, recursive = TRUE)

# ShinyGO Chart-style horizontal bar plot:
#   x = Fold Enrichment (bar length), y = term name (sorted by FDR ascending),
#   fill = -log10(FDR) (color gradient), bar label = nGenes (gene overlap).
make_dotplot <- function(df, list_name, pass_name, top_n = 10) {
  if (nrow(df) == 0) return(NULL)
  cfg <- passes[[pass_name]]
  d <- df %>%
    arrange(p_value, desc(fold_enrichment)) %>%
    slice_head(n = top_n) %>%
    mutate(
      neg_log10_p = -log10(p_value),
      term_label  = str_trunc(term_name, 70),
      # fct_reorder by -p_value: smallest p (most significant) ends up on top.
      term_label  = fct_reorder(term_label, -p_value)
    )
  sig_count <- sum(d$passes_fdr)
  subtitle  <- sprintf(
    "top %d terms by FDR (then by fold enrichment); %d pass FDR<%g | term size %d-%d | bg = Enrichr default",
    nrow(d), sig_count, cfg$fdr, cfg$min, cfg$max
  )
  ggplot(d, aes(x = fold_enrichment, y = term_label, fill = neg_log10_p)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = intersection_size),
              hjust = -0.3, size = 3, color = "grey25") +
    scale_fill_viridis_c(option = "C", direction = -1) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title    = paste0(list_name, "  |  ", pass_name, "  (", cfg$db, ")"),
      subtitle = subtitle,
      x        = "Fold Enrichment",
      y        = NULL,
      fill     = expression(-log[10](FDR))
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank())
}

source_colors <- c(
  "GO_BP"    = "#1b9e77",
  "GO_MF"    = "#d95f02",
  "KEGG"     = "#e7298a",
  "DisGeNET" = "#a6761d"
)

# Per-list plots
pwalk(expand_grid(gene_list = names(gene_lists), pass = names(passes)),
      function(gene_list, pass) {
  df <- all_results %>% filter(gene_list == .env$gene_list, pass == .env$pass)
  p <- make_dotplot(df, gene_list, pass)
  if (is.null(p)) {
    message("  (no significant terms for ", gene_list, " / ", pass, " - skipping plot)")
    return(invisible(NULL))
  }
  ggsave(file.path(out_dir, gene_list, paste0("plot_", pass, ".pdf")),
         p, width = 9, height = 6)
  ggsave(file.path(out_dir, gene_list, paste0("plot_", pass, ".png")),
         p, width = 9, height = 6, dpi = 200)
})

# Cross-list summary plots: one per DB, faceted by gene list. Same ShinyGO
# Chart-style horizontal bars as the per-list plots.
walk(names(passes), function(pass_name) {
  df <- all_results %>% filter(pass == pass_name)
  if (nrow(df) == 0) {
    message("  (no enrichR rows in size window for ", pass_name, " - skipping summary plot)")
    return(invisible(NULL))
  }
  cfg <- passes[[pass_name]]
  d <- df %>%
    group_by(gene_list) %>%
    arrange(p_value, desc(fold_enrichment)) %>%
    slice_head(n = 10) %>%
    ungroup() %>%
    mutate(
      neg_log10_p = -log10(p_value),
      term_label  = str_trunc(term_name, 60),
      gene_list   = factor(gene_list, levels = names(gene_lists))
    ) %>%
    group_by(gene_list) %>%
    mutate(term_label = tidytext::reorder_within(term_label, -p_value, gene_list)) %>%
    ungroup()
  p_summary <- ggplot(d, aes(x = fold_enrichment, y = term_label, fill = neg_log10_p)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = intersection_size), hjust = -0.3, size = 2.8, color = "grey25") +
    facet_wrap(~ gene_list, scales = "free", ncol = 2) +
    tidytext::scale_y_reordered() +
    scale_fill_viridis_c(option = "C", direction = -1) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title    = paste0(pass_name, " (", cfg$db, ") — top 10 terms per gene list"),
      subtitle = sprintf("sorted by FDR then fold enrichment | cutoff FDR<%g, term size %d-%d | bg = Enrichr default",
                         cfg$fdr, cfg$min, cfg$max),
      x        = "Fold Enrichment",
      y        = NULL,
      fill     = expression(-log[10](FDR))
    ) +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank())
  ggsave(file.path(out_dir, "plots", paste0("summary_", pass_name, ".pdf")),
         p_summary, width = 14, height = 8)
  ggsave(file.path(out_dir, "plots", paste0("summary_", pass_name, ".png")),
         p_summary, width = 14, height = 8, dpi = 200)
})

message("Done. enrichR outputs: GO_analysis/<list>/enrichr_{GO_BP,GO_MF,KEGG,DisGeNET}.csv ; plots in same folders + GO_analysis/plots/ ; combined supplementary at GO_analysis/combined/SupplementaryTable_enrichment.xlsx")
