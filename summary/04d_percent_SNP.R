#!/usr/bin/env Rscript
# Percent of twin heritability explained by SNP (MASH/AdjHE-RE) per phenotype set
# For each phenotype: pct = 100 * h2_MASH / h2_Twin
# Summarize per Set (gordon, probaConns, SA) with mean, SD, SE, median, etc.
# Input: results/summary/mash_twin_wide.csv (Twin_h2, h2_gordon_AdjHE_RE, h2_proba_AdjHE_RE, h2_SA_AdjHE_RE)
# Output: results/summary/percent_SNP_explained.csv + per-set details
# Also used in quarto 06_results_catalog.qmd

library(tidyverse)

ROOT <- "/users/4/coffm049/papers/functionalBrainHerit"
WIDE <- file.path(ROOT, "results/summary/mash_twin_wide.csv")
OUT  <- file.path(ROOT, "results/summary/percent_SNP_explained.csv")
# Portable fallback for local Windows
if (!file.exists(WIDE)) {
  WIDE <- file.path("results", "summary", "mash_twin_wide.csv")
  OUT  <- file.path("results", "summary", "percent_SNP_explained.csv")
  ROOT <- normalizePath(".", winslash = "/")
}
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)

wide <- tryCatch(read_csv(WIDE, show_col_types = FALSE), error = function(e) tibble())
if (nrow(wide)==0) stop("mash_twin_wide.csv not found or empty: ", WIDE)

# Map Set -> MASH column
# gordon -> h2_gordon_AdjHE_RE, probaConns -> h2_proba_AdjHE_RE, SA -> h2_SA_AdjHE_RE
# Note: column names use underscore (AdjHE_RE) for files, display as AdjHE-RE (hyphen)
mash_col_for_set <- c(
  "gordon" = "h2_gordon_AdjHE_RE",
  "probaConns" = "h2_proba_AdjHE_RE",
  "SA" = "h2_SA_AdjHE_RE"
)

# For SA, we want 14 networks excl. 4/6/17 per 01-07 labeling (as in 03_plot)
sa_keep <- c(1,2,3,5,7,8,9,10,11,12,13,14,15,16)

results <- list()
details <- list()

for (set_name in names(mash_col_for_set)) {
  mash_col <- mash_col_for_set[[set_name]]
  if (!mash_col %in% names(wide)) {
    # Try alternative naming (proba vs probaConns)
    alt <- gsub("probaConns", "proba", mash_col)
    if (alt %in% names(wide)) mash_col <- alt else next
  }
  sub <- wide %>% filter(Set == set_name) %>% select(Pheno, Twin_h2, all_of(mash_col)) %>% rename(MASH_h2 = all_of(mash_col))
  # For SA, filter to 14 networks (excl. 4/6/17) to avoid shift, per 05-topoViz.R
  if (set_name == "SA") {
    sub <- sub %>% mutate(pheno_num = parse_number(Pheno)) %>% filter(pheno_num %in% sa_keep) %>% select(-pheno_num)
  }
  # Filter to finite, Twin_h2 > 0.01 to avoid near-zero denominator blow-up
  # Keep all for main, but also report filtered
  sub <- sub %>% mutate(
    Twin_h2 = as.numeric(Twin_h2),
    MASH_h2 = as.numeric(MASH_h2),
    pct = 100 * MASH_h2 / Twin_h2
  )
  # Define valid pct for summary: finite Twin and MASH, 0 < Twin_h2 <= 1, 0 <= MASH_h2 <= 1.2 (cap at 120% to avoid extreme outliers)
  # Keep uncapped for mean, but also report median which is robust
  valid <- sub %>% filter(is.finite(Twin_h2), is.finite(MASH_h2), Twin_h2 > 0.001, Twin_h2 <= 1, MASH_h2 >= 0, MASH_h2 <= 1.2, is.finite(pct))
  # Also keep all finite pct for n
  all_valid_n <- nrow(valid)
  if (all_valid_n == 0) {
    message("No valid pct for set ", set_name)
    next
  }
  # Summary stats
  mean_pct <- mean(valid$pct, na.rm = TRUE)
  sd_pct <- sd(valid$pct, na.rm = TRUE)
  se_pct <- sd_pct / sqrt(nrow(valid))
  median_pct <- median(valid$pct, na.rm = TRUE)
  q1 <- quantile(valid$pct, 0.25, na.rm = TRUE)
  q3 <- quantile(valid$pct, 0.75, na.rm = TRUE)
  min_pct <- min(valid$pct, na.rm = TRUE)
  max_pct <- max(valid$pct, na.rm = TRUE)
  # Also compute overall mean h2 for context
  mean_twin <- mean(valid$Twin_h2, na.rm = TRUE)
  mean_mash <- mean(valid$MASH_h2, na.rm = TRUE)
  # Ratio of means vs mean of ratios
  ratio_of_means <- 100 * mean_mash / mean_twin

  results[[set_name]] <- tibble(
    Set = set_name,
    n_pheno = nrow(sub),
    n_valid = nrow(valid),
    mean_Twin_h2 = mean_twin,
    mean_MASH_h2 = mean_mash,
    ratio_of_means_pct = ratio_of_means,
    mean_pct = mean_pct,
    median_pct = median_pct,
    sd_pct = sd_pct,
    se_pct = se_pct,
    q1_pct = as.numeric(q1),
    q3_pct = as.numeric(q3),
    min_pct = min_pct,
    max_pct = max_pct
  )
  # Save per-pheno details for quarto inspection
  details[[set_name]] <- valid %>% mutate(Set = set_name) %>% select(Set, Pheno, Twin_h2, MASH_h2, pct)
}

if (length(results)==0) stop("No results produced")

out <- bind_rows(results) %>% mutate(across(where(is.numeric), ~round(.x, 2)))
write_csv(out, OUT)
message("Wrote ", OUT)
print(out)

# Also write per-pheno details
details_out <- file.path(dirname(OUT), "percent_SNP_explained_per_pheno.csv")
bind_rows(details) %>% mutate(across(where(is.numeric), ~round(.x, 4))) %>% write_csv(details_out)
message("Wrote ", details_out)

# Print for log
cat("\n=== Percent SNP of twin (MASH/Twin*100) by Set ===\n")
print(out %>% select(Set, n_valid, mean_pct, se_pct, median_pct, sd_pct, ratio_of_means_pct) %>% arrange(Set))
