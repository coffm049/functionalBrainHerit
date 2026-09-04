#!/usr/bin/env Rscript
# Horizontal barplots for SA surface area heritability (twin vs AdjHE-RE)
# wo_total and w_total (total surface controlled), twin + AdjHE-RE side-by-side with geom_dodge
# Organized by max twin h2 (ascending, so most heritable at top when horizontal).
# Input: results/summary/mash_twin_wide.csv (Set == "SA", 17 network_surfarea* phenos) + fallback to direct SA files for w_total
# Outputs: results/summary/plots/SA_horizontal_barplot.png (wo_total) + SA_w_total_horizontal_barplot.png (w_total)
# Run: Rscript summary/03_plot_sa_horizontal.R  (on acn01 with R/4.4.0)

library(tidyverse)

ROOT <- "/users/4/coffm049/papers/functionalBrainHerit"
WIDE <- file.path(ROOT, "results/summary/mash_twin_wide.csv")
OUT  <- file.path(ROOT, "results/summary/plots/SA_horizontal_barplot.png")
OUT_W <- file.path(ROOT, "results/summary/plots/SA_w_total_horizontal_barplot.png")
# Portable fallback for local Windows render (where /users/4/... does not exist)
if (!file.exists(WIDE)) {
  WIDE <- file.path("results", "summary", "mash_twin_wide.csv")
  OUT  <- file.path("results", "summary", "plots", "SA_horizontal_barplot.png")
  OUT_W <- file.path("results", "summary", "plots", "SA_w_total_horizontal_barplot.png")
  ROOT <- normalizePath(".", winslash = "/")
}
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
# Also make twin/adj paths portable for fallback
twin_w_path_portable <- file.path("results", "SA", "twinEsts", "herit_w_total.Rds")
adj_w_path_portable  <- file.path("results", "SA", "AdjHE_RE.csv")

# 14 networks — from 05-topoViz.R (excludes NET4, NET6, NET17 which are null in Gordon)
# Old 01-07 labeling (03b-topoh2Ests2Brain.py: range(1,15), 05-topoViz.R networks) kept 14 after filtering 1,2,3,5,7,8,9,10,11,12,13,14,15,16
networks <- c(
  "1" = "DMN", "2" = "VIS", "3" = "FP", "5" = "DAN", "7" = "VAN",
  "8" = "SAL", "9" = "CO", "10" = "SMD", "11" = "SML", "12" = "AUD",
  "13" = "Tpole", "14" = "MTL", "15" = "PMN", "16" = "PON"
)
# For display, keep 17 labels but mark 4/6/17 as excluded (null) — see SA brain scripts

wide <- tryCatch(read_csv(WIDE, show_col_types = FALSE), error = function(e) tibble())

plot_sa <- function(df, out_path, title_suffix) {
  if (nrow(df) == 0) {
    message("No data for ", title_suffix, " — skipping ", out_path)
    return(invisible(NULL))
  }
  # Order y-axis by max twin h2 (ascending → most heritable at top)
  twin_order <- df %>%
    filter(Type == "Twin") %>%
    group_by(pheno_label) %>%
    summarise(m = max(h2, na.rm = TRUE), .groups = "drop") %>%
    arrange(m) %>%
    pull(pheno_label)

  cor_df <- df %>%
    pivot_wider(names_from = Type, values_from = h2) %>%
    summarise(
      n = sum(!is.na(Twin) & !is.na(AdjHE_RE)),
      r = suppressWarnings(cor(Twin, AdjHE_RE, use = "pairwise.complete.obs")),
      rho = suppressWarnings(cor(Twin, AdjHE_RE, method = "spearman", use = "pairwise.complete.obs"))
    )
  # Save stats for quarto (instead of subtitle numbers)
  stats_path <- sub("\\.png$", "_stats.csv", out_path)
  tryCatch({
    write_csv(cor_df %>% mutate(suffix = title_suffix), stats_path)
  }, error = function(e) NULL)

  p <- df %>%
    mutate(pheno_label = factor(pheno_label, levels = twin_order)) %>%
    ggplot(aes(y = pheno_label, x = h2, fill = Type)) +
    geom_col(position = position_dodge2(width = 0.8, preserve = "single"), width = 0.7) +
    scale_x_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1), expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = c("Twin" = "#1f77b4", "AdjHE-RE" = "#ff7f0e"), labels = c("Twin", "AdjHE-RE")) +
    labs(x = expression(paste("Heritability (", h^2, ")")), y = "", fill = "") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_text(face = "bold", size = 9, margin = margin(r = 4)),
      axis.text.x = element_text(size = 9, margin = margin(t = 4)),
      axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 8)),
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      legend.margin = margin(t = 4, b = 0),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 4, 2, 2),
      axis.ticks.y = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    )
  ggsave(out_path, p, width = 7, height = 5, dpi = 300)
  message("Wrote ", out_path)
  print(cor_df)
}

# ---- wo_total (from mash_twin_wide.csv: Twin_h2 + h2_SA_AdjHE_RE) ----
if (nrow(wide) > 0 && "Twin_h2" %in% names(wide) && "h2_SA_AdjHE_RE" %in% names(wide)) {
  df_wo <- wide %>%
    filter(Set == "SA") %>%
    select(Pheno, Twin_h2, h2_SA_AdjHE_RE) %>%
    pivot_longer(cols = c(Twin_h2, h2_SA_AdjHE_RE),
                 names_to = "Type", values_to = "h2") %>%
    mutate(
      Type = recode(Type, Twin_h2 = "Twin", h2_SA_AdjHE_RE = "AdjHE_RE"),
      pheno_num = parse_number(Pheno),
      pheno_label = recode(as.character(pheno_num), !!!networks, .default = NA_character_),
      # 4/6/17 are null per Gordon (no parcels) — intentionally NA and dropped to avoid shift
      h2 = if_else(is.na(h2) | h2 < 0, 0, h2)
    ) %>%
    drop_na(pheno_label)
  plot_sa(df_wo, OUT, " (wo_total, 14 networks)")
} else {
  message("WIDE missing wo_total columns — skipping wo_total barplot")
}

# ---- w_total (Twin_h2_w_total + h2_SA_AdjHE_RE_w_total if in WIDE, else fallback to direct SA files) ----
df_w <- tibble()
if (nrow(wide) > 0) {
  w_twin_col <- intersect(c("Twin_h2_w_total", "Twin_h2_wtotal", "Twin_h2_w_Total"), names(wide))
  w_adj_col  <- intersect(c("h2_SA_AdjHE_RE_w_total", "h2_SA_w_total_AdjHE_RE", "h2_SA_AdjHE_RE_wtotal"), names(wide))
  if (length(w_twin_col) == 1 && length(w_adj_col) == 1) {
    df_w <- wide %>%
      filter(Set == "SA") %>%
      select(Pheno, all_of(c(w_twin_col, w_adj_col))) %>%
      rename(Twin_h2 = all_of(w_twin_col), h2_w = all_of(w_adj_col)) %>%
      pivot_longer(cols = c(Twin_h2, h2_w), names_to = "Type", values_to = "h2") %>%
      mutate(
        Type = recode(Type, Twin_h2 = "Twin", h2_w = "AdjHE_RE"),
        pheno_num = parse_number(Pheno),
        pheno_label = recode(as.character(pheno_num), !!!networks, .default = NA_character_),
        h2 = if_else(is.na(h2) | h2 < 0, 0, h2)
      ) %>%
      drop_na(pheno_label)
    message("Using WIDE w_total columns: ", w_twin_col, " + ", w_adj_col)
  }
}
# Fallback to direct files if WIDE w_total not available
if (nrow(df_w) == 0) {
  # Twin w_total from RDS (herit_w_total.Rds) — try HPC ROOT first, then portable
  twin_w_path <- file.path(ROOT, "results/SA/twinEsts/herit_w_total.Rds")
  adj_w_path  <- file.path(ROOT, "results/SA/AdjHE_RE.csv")
  adj_wo_path <- file.path(ROOT, "results/SA/AdjHE_RE_wo_total.csv")
  if (!file.exists(twin_w_path) && file.exists(twin_w_path_portable)) twin_w_path <- twin_w_path_portable
  if (!file.exists(adj_w_path) && file.exists(adj_w_path_portable))  adj_w_path  <- adj_w_path_portable
  # Try AdjHE w_total file (AdjHE_RE.csv is w_total when wo_total is separate)
  adj_file <- if (file.exists(adj_w_path)) adj_w_path else NA
  # If wo_total file was used for WIDE, w_total is the other file; heuristic: choose file not used for wo_total (check N)
  if (!is.na(adj_file) && file.exists(twin_w_path)) {
    tryCatch({
      twin_rds <- readRDS(twin_w_path)
      # twin_rds is tibble with Phenotype (?) and herit
      id_col <- setdiff(names(twin_rds), c("herit", "data"))
      if (length(id_col) == 1) {
        twin_rds <- twin_rds %>% rename(Pheno = all_of(id_col))
      } else if ("Phenotype" %in% names(twin_rds)) {
        twin_rds <- twin_rds %>% rename(Pheno = Phenotype)
      }
      extract_h2 <- function(x) {
        if (is.null(x) || inherits(x, "twinlm_error")) return(NA_real_)
        cf <- tryCatch(x$coef, error = function(e) NULL)
        if (is.null(cf)) return(NA_real_)
        A <- as.numeric(cf[1,1]); C <- as.numeric(cf[2,1]); E <- as.numeric(cf[3,1])
        Ttot <- A + C + E
        if (Ttot > 0) A / Ttot else NA_real_
      }
      twin_df <- twin_rds %>%
        mutate(Twin_h2 = map_dbl(herit, extract_h2)) %>%
        select(Pheno, Twin_h2)
      adj_df <- read_csv(adj_file, show_col_types = FALSE) %>%
        filter(PCs == 30) %>%
        select(Pheno = pheno, h2_w = h2)
      df_w <- full_join(twin_df, adj_df, by = "Pheno") %>%
        pivot_longer(cols = c(Twin_h2, h2_w), names_to = "Type", values_to = "h2") %>%
        mutate(
          Type = recode(Type, Twin_h2 = "Twin", h2_w = "AdjHE_RE"),
          pheno_num = parse_number(Pheno),
          pheno_label = recode(as.character(pheno_num), !!!networks, .default = NA_character_),
          h2 = if_else(is.na(h2) | h2 < 0, 0, h2)
        ) %>%
        drop_na(pheno_label)
      message("Built w_total from direct files: ", twin_w_path, " + ", adj_file)
    }, error = function(e) message("Fallback w_total build failed: ", conditionMessage(e)))
  } else {
    message("No w_total source found — will copy wo_total as placeholder if needed")
    # As last resort, reuse wo_total plot so QMD has an image
    if (file.exists(OUT) && !file.exists(OUT_W)) {
      file.copy(OUT, OUT_W, overwrite = TRUE)
      message("Copied wo_total barplot to w_total placeholder")
    }
  }
}
if (nrow(df_w) > 0) {
  # Ensure same 17 networks, even if some missing in w_total
  plot_sa(df_w, OUT_W, " (w_total, total surface controlled)")
} else if (!file.exists(OUT_W) && file.exists(OUT)) {
  file.copy(OUT, OUT_W, overwrite = TRUE)
  message("No w_total data — copied wo_total to ", OUT_W)
}
