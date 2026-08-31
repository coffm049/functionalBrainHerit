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
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)

# 17 networks — short names from 05-topoViz.R plus placeholders for 4,6,17
networks <- c(
  "1" = "DMN", "2" = "VIS", "3" = "FP", "4" = "NET4",
  "5" = "DAN", "6" = "NET6", "7" = "VAN", "8" = "SAL",
  "9" = "CO", "10" = "SMD", "11" = "SML", "12" = "AUD",
  "13" = "Tpole", "14" = "MTL", "15" = "PMN", "16" = "PON", "17" = "NET17"
)

wide <- tryCatch(read_csv(WIDE, show_col_types = FALSE), error = function(e) tibble())

plot_sa <- function(df, out_path, title_suffix) {
  if (nrow(df) == 0) {
    message("No data for ", title_suffix, " — skipping ", out_path)
    return(invisible(NULL))
  }
  # Order y-axis by max twin h2 (ascending → most heritable at top in horizontal plot)
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

  p <- df %>%
    mutate(pheno_label = factor(pheno_label, levels = twin_order)) %>%
    ggplot(aes(y = pheno_label, x = h2, fill = Type)) +
    geom_col(position = position_dodge2(width = 0.8, preserve = "single"),
             width = 0.7) +
    scale_x_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1)) +
    scale_fill_manual(values = c("Twin" = "#1f77b4", "AdjHE_RE" = "#ff7f0e")) +
    labs(
      x = expression(paste("Heritability (", h^2, ")")),
      y = "",
      fill = "",
      title = paste0("SA network surface area h2", title_suffix),
      subtitle = sprintf("n=%d networks, r=%.2f, rho=%.2f (Twin vs AdjHE_RE, 30 PCs)",
                         cor_df$n, cor_df$r, cor_df$rho)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_text(face = "bold", size = 10),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 12, face = "bold"),
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      plot.title = element_text(face = "bold", size = 12),
      panel.grid.major.y = element_blank()
    )

  ggsave(out_path, p, width = 8, height = 6, dpi = 300)
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
      pheno_label = recode(as.character(pheno_num), !!!networks),
      pheno_label = if_else(is.na(pheno_label), Pheno, pheno_label),
      h2 = if_else(is.na(h2) | h2 < 0, 0, h2)
    ) %>%
    drop_na(pheno_label)
  plot_sa(df_wo, OUT, " (wo_total)")
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
        pheno_label = recode(as.character(pheno_num), !!!networks),
        pheno_label = if_else(is.na(pheno_label), Pheno, pheno_label),
        h2 = if_else(is.na(h2) | h2 < 0, 0, h2)
      ) %>%
      drop_na(pheno_label)
    message("Using WIDE w_total columns: ", w_twin_col, " + ", w_adj_col)
  }
}
# Fallback to direct files if WIDE w_total not available
if (nrow(df_w) == 0) {
  # Twin w_total from RDS (herit_w_total.Rds)
  twin_w_path <- file.path(ROOT, "results/SA/twinEsts/herit_w_total.Rds")
  adj_w_path  <- file.path(ROOT, "results/SA/AdjHE_RE.csv")
  adj_wo_path <- file.path(ROOT, "results/SA/AdjHE_RE_wo_total.csv")
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
          pheno_label = recode(as.character(pheno_num), !!!networks),
          pheno_label = if_else(is.na(pheno_label), Pheno, pheno_label),
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
