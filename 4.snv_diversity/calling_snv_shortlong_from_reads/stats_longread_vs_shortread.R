
suppressPackageStartupMessages({
  library(dplyr); library(data.table); library(tidyr); library(tibble)
})

pantros  <- fread("diversity_stats/pantros_combined_pi_stats.tsv")
ppan     <- fread("diversity_stats/ppan_combined_pi_stats.tsv")
metadata <- fread("PANPAN_PradoMartine_deManuel_METADATA.tsv")

apes <- bind_rows(pantros, ppan) %>%
  left_join(metadata %>% distinct(sample, .keep_all = TRUE), by = "sample") %>%
  mutate(type3 = trimws(type3)) %>%
  filter(!is.na(CommonName), !CommonName %in% c("Nigeria-Cameroon", "Western x Central"))

# Per-population means and % difference
summary_tbl <- apes %>%
  group_by(CommonName, type3) %>%
  summarise(n = n(), mean_pi = mean(mean_pi), .groups = "drop") %>%
  pivot_wider(names_from = type3, values_from = c(n, mean_pi)) %>%
  mutate(pct_higher_long = 100 * (`mean_pi_long-read` - `mean_pi_short-read`) /
                                 `mean_pi_short-read`)

cat("=== Per-population means (mean of individual pi) ===\n")
print(as.data.frame(summary_tbl), digits = 4)

# Per-pop significance tests + multiple-testing correction
tests_tbl <- apes %>%
  group_by(CommonName) %>%
  summarise(
    wilcox_p = tryCatch(wilcox.test(mean_pi ~ type3)$p.value, error = function(e) NA),
    ttest_p  = tryCatch(t.test(mean_pi ~ type3)$p.value,      error = function(e) NA),
    .groups  = "drop"
  ) %>%
  mutate(
    wilcox_p_bonferroni = p.adjust(wilcox_p, method = "bonferroni"),
    wilcox_p_BH         = p.adjust(wilcox_p, method = "BH"),
    ttest_p_bonferroni  = p.adjust(ttest_p,  method = "bonferroni"),
    ttest_p_BH          = p.adjust(ttest_p,  method = "BH")
  )

out_tbl <- summary_tbl %>% left_join(tests_tbl, by = "CommonName")
dir.create("figures", showWarnings = FALSE)
fwrite(out_tbl, "figures/EDFig1H_longread_vs_shortread_stats.csv")
cat("\nWrote figures/EDFig1H_longread_vs_shortread_stats.csv\n")
cat("\n=== Per-pop tests with multiple-testing correction ===\n")
print(as.data.frame(tests_tbl), digits = 3)

# Overall (pooling all populations)
overall <- apes %>%
  group_by(type3) %>%
  summarise(n = n(), mean_pi = mean(mean_pi), .groups = "drop")
cat("\n=== Overall pooled ===\n")
print(as.data.frame(overall), digits = 4)
pct_overall <- 100 * (overall$mean_pi[overall$type3 == "long-read"] -
                      overall$mean_pi[overall$type3 == "short-read"]) /
                     overall$mean_pi[overall$type3 == "short-read"]
cat(sprintf("\nOverall long-read is %.2f%% higher than short-read (pooled mean-of-means)\n",
            pct_overall))

# Mean of per-population % differences
cat(sprintf("Mean of per-pop %% differences: %.2f%% (range %.2f%% - %.2f%%)\n",
            mean(summary_tbl$pct_higher_long, na.rm = TRUE),
            min(summary_tbl$pct_higher_long,  na.rm = TRUE),
            max(summary_tbl$pct_higher_long,  na.rm = TRUE)))

# Significance tests per population
cat("\n=== Wilcoxon & t-tests per population (long-read vs short-read) ===\n")
tests <- apes %>%
  group_by(CommonName) %>%
  summarise(
    n_long  = sum(type3 == "long-read"),
    n_short = sum(type3 == "short-read"),
    wilcox_p = tryCatch(wilcox.test(mean_pi ~ type3)$p.value, error = function(e) NA),
    ttest_p  = tryCatch(t.test(mean_pi ~ type3)$p.value,      error = function(e) NA),
    .groups  = "drop"
  )
print(as.data.frame(tests), digits = 4)

# Overall test, stratified by population (paired-like): linear model with pop as covariate
cat("\n=== Linear model: mean_pi ~ type3 + CommonName ===\n")
m <- lm(mean_pi ~ type3 + CommonName, data = apes)
print(summary(m)$coefficients)

# Overall Wilcoxon ignoring population
cat("\n=== Overall Wilcoxon (ignoring pop) ===\n")
print(wilcox.test(mean_pi ~ type3, data = apes))
cat("\n=== Overall t-test (ignoring pop) ===\n")
print(t.test(mean_pi ~ type3, data = apes))
