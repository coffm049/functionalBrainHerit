#!/usr/bin/env Rscript
# Compare R2 between methods for each phenotype set, and SA w/ vs wo total
# Input: results/summary/mash_twin_wide.csv  (and SA w_total files if present)
# Outputs:
#   results/summary/r2_by_set.csv
#   results/summary/r2_sa_w_vs_wo_total.csv
#   results/summary/plots/r2_scatter_by_set.png
#   results/summary/plots/r2_sa_w_vs_wo.png
# Run: Rscript summary/compare_r2.R

library(tidyverse)

ROOT <- "/users/4/coffm049/papers/functionalBrainHerit"
WIDE <- file.path(ROOT, "results/summary/mash_twin_wide.csv")
OUT_DIR <- file.path(ROOT, "results/summary")
PLOT_DIR <- file.path(OUT_DIR, "plots")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

wide <- read_csv(WIDE, show_col_types = FALSE)

# Helper to compute R2 (Pearson r squared) and Spearman
r2_stats <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 3) return(tibble(n = length(x), Pearson_r = NA_real_, Spearman_rho = NA_real_, R2_pearson = NA_real_, R2_spearman = NA_real_))
  r <- suppressWarnings(cor(x, y, method = "pearson"))
  rho <- suppressWarnings(cor(x, y, method = "spearman"))
  tibble(n = length(x), Pearson_r = r, Spearman_rho = rho, R2_pearson = r^2, R2_spearman = rho^2)
}

# -------------------------------------------------------------------------
# 1. R2 between Twin and each MASH method, per Set (gordon, probaConns, SA)
# -------------------------------------------------------------------------
melt <- wide %>%
  select(Set, Pheno, Twin_h2, starts_with("h2_")) %>%
  pivot_longer(starts_with("h2_"), names_to = "MashMethod", values_to = "Mash_h2")

r2_by_set <- melt %>%
  group_by(Set, MashMethod) %>%
  summarise(r2_stats(Twin_h2, Mash_h2), .groups = "drop")

# Overall (collapsed sets, for reference)
r2_overall <- melt %>%
  group_by(MashMethod) %>%
  summarise(r2_stats(Twin_h2, Mash_h2), .groups = "drop") %>%
  mutate(Set = "ALL")

r2_all <- bind_rows(r2_by_set, r2_overall) %>%
  arrange(Set, MashMethod)

write_csv(r2_all, file.path(OUT_DIR, "r2_by_set.csv"))
message("Wrote ", file.path(OUT_DIR, "r2_by_set.csv"))
print(r2_all)

# Scatter with R2 annotation per Set
p1 <- melt %>%
  filter(!is.na(Twin_h2), !is.na(Mash_h2)) %>%
  ggplot(aes(x = Twin_h2, y = Mash_h2)) +
  geom_point(alpha = 0.3, size = 0.7) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 0.5) +
  facet_wrap(~Set + MashMethod, scales = "free", ncol = 3) +
  labs(x = "Twin h2", y = "MASH h2", title = "Twin vs MASH h2 (R2 = Pearson r^2)") +
  theme_minimal() +
  geom_text(data = r2_all %>% filter(Set != "ALL"),
            aes(label = sprintf("R2=%.2f\nn=%d", R2_pearson, n)),
            x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, size = 2.5, inherit.aes = FALSE)

ggsave(file.path(PLOT_DIR, "r2_scatter_by_set.png"), p1, width = 9, height = 6, dpi = 300)
message("Wrote ", file.path(PLOT_DIR, "r2_scatter_by_set.png"))

# -------------------------------------------------------------------------
# 2. SA: w_total vs wo_total (total surface controlled vs not)
# -------------------------------------------------------------------------
# Try to read w_total vs wo_total from wide if compare_mash_twin.R has been updated
# to include both streams. Fallback: read SA results directly.

sa_w_vs_wo <- NULL

# Attempt 1: wide has distinct columns for w_total vs wo_total
if (any(grepl("w_total", names(wide)))) {
  sa_w_vs_wo <- wide %>%
    filter(Set == "SA") %>%
    select(Pheno, Twin_h2, matches("h2_SA.*w_total"), matches("h2_SA.*wo_total")) %>%
    pivot_longer(cols = -c(Pheno, Twin_h2), names_to = "SA_Method", values_to = "h2") %>%
    # Expect columns like h2_SA_AdjHE_RE and h2_SA_AdjHE_RE_w_total
    mutate(wo_total = if_else(grepl("wo_total", SA_Method), h2, NA_real_),
           w_total  = if_else(grepl("w_total") & !grepl("wo_total", SA_Method), h2, NA_real_)) %>%
    group_by(Pheno) %>%
    summarise(wo_total = first(na.omit(wo_total)), w_total = first(na.omit(w_total)), Twin_h2 = first(Twin_h2), .groups = "drop")
} else {
  # Attempt 2: read SA files directly
  # Twin wo_total vs w_total
  twin_wo_path <- file.path(ROOT, "results/SA/twinEsts/herit_wo_total.Rds")
  twin_w_path  <- file.path(ROOT, "results/SA/twinEsts/herit_w_total.Rds")
  # MASH wo_total vs w_total
  mash_wo_path <- file.path(ROOT, "results/SA/AdjHE_RE_wo_total.csv")
  mash_w_path  <- file.path(ROOT, "results/SA/AdjHE_RE.csv")

  # Helper to extract h2 per pheno from RDS/CSV
  extract_twin_h2 <- function(path) {
    if (!file.exists(path)) return(NULL)
    df <- readRDS(path)
    # df has columns phenotype (network_surfarea*) and herit list
    # Reuse compare_mash_twin.R::extract_twin logic simplified
    df %>%
      mutate(
        Pheno = df[[setdiff(names(df), c("herit","data"))[1]]],
        h2 = map_dbl(herit, function(x) {
          if (is.null(x) || inherits(x, "twinlm_error")) return(NA_real_)
          cf <- tryCatch(x$coef, error = function(e) NULL)
          if (is.null(cf) || nrow(cf) < 3) return(NA_real_)
          A <- as.numeric(cf[1,1]); C <- as.numeric(cf[2,1]); E <- as.numeric(cf[3,1])
          Ttot <- A + C + E
          if (Ttot > 0) A / Ttot else NA_real_
        })
      ) %>%
      select(Pheno, h2)
  }

  extract_mash_h2 <- function(path) {
    if (!file.exists(path)) return(NULL)
    d <- read_csv(path, show_col_types = FALSE)
    if ("pheno" %in% names(d)) d <- rename(d, Pheno = pheno)
    # Keep max PCs if present
    if ("PCs" %in% names(d)) {
      d <- d %>% group_by(Pheno) %>% filter(PCs == max(PCs, na.rm = TRUE)) %>% ungroup()
    }
    select(d, Pheno, h2)
  }

  twin_wo <- extract_twin_h2(twin_wo_path) %>% rename(wo_total = h2)
  twin_w  <- extract_twin_h2(twin_w_path)  %>% rename(w_total = h2)
  mash_wo <- extract_mash_h2(mash_wo_path) %>% rename(wo_total = h2)
  mash_w  <- extract_mash_h2(mash_w_path)  %>% rename(w_total = h2)

  if (!is.null(twin_wo) && !is.null(twin_w)) {
    sa_w_vs_wo_twin <- full_join(twin_wo, twin_w, by = "Pheno") %>%
      mutate(Set = "SA", Method = "Twin")
  } else sa_w_vs_wo_twin <- NULL

  if (!is.null(mash_wo) && !is.null(mash_w)) {
    sa_w_vs_wo_mash <- full_join(mash_wo, mash_w, by = "Pheno") %>%
      mutate(Set = "SA", Method = "AdjHE_RE")
  } else sa_w_vs_wo_mash <- NULL

  sa_w_vs_wo <- bind_rows(sa_w_vs_wo_twin, sa_w_vs_wo_mash)
}

if (!is.null(sa_w_vs_wo) && nrow(sa_w_vs_wo) > 0) {
  r2_sa <- sa_w_vs_wo %>%
    group_by(Set, Method) %>%
    summarise(r2_stats(wo_total, w_total), .groups = "drop")

  write_csv(r2_sa, file.path(OUT_DIR, "r2_sa_w_vs_wo_total.csv"))
  message("Wrote ", file.path(OUT_DIR, "r2_sa_w_vs_wo_total.csv"))
  print(r2_sa)

  p2 <- sa_w_vs_wo %>%
    ggplot(aes(x = wo_total, y = w_total)) +
    geom_point(alpha = 0.8, size = 2) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 0.5) +
    facet_wrap(~Method, scales = "free") +
    labs(x = "SA h2 wo_total (no total surface)", y = "SA h2 w_total (with total surface)",
         title = "SA: with vs without total surface control") +
    theme_minimal() +
    geom_text(data = r2_sa,
              aes(label = sprintf("R2=%.2f\nn=%d", R2_pearson, n)),
              x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, size = 3, inherit.aes = FALSE)

  ggsave(file.path(PLOT_DIR, "r2_sa_w_vs_wo.png"), p2, width = 7, height = 4, dpi = 300)
  message("Wrote ", file.path(PLOT_DIR, "r2_sa_w_vs_wo.png"))
} else {
  message("SA w_total vs wo_total files not found — run the w_total SLURM jobs first (SA/SLURM_AdjHE_RE_with_total.sh and SA/twinEsts/SLURM_twin_with_total.sh) and rebuild mash_twin_wide.csv")
}
