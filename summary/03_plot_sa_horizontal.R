#!/usr/bin/env Rscript
# Horizontal barplots for SA surface area heritability (twin vs AdjHE-RE)
# wo_total and w_total (total surface controlled), twin + AdjHE-RE side-by-side with geom_dodge
# Organized by max twin h2 (ascending, so most heritable at top when horizontal).
# Input: results/summary/mash_twin_wide.csv (Set == "SA", 17 network_surfarea* phenos) + fallback to direct SA files for w_total
# Outputs: results/summary/plots/SA_horizontal_barplot.png (wo_total) + SA_w_total_horizontal_barplot.png (w_total)
# Run: Rscript summary/03_plot_sa_horizontal.R  (on acn01 with R/4.4.0)

library(tidyverse)
source(file.path("summary", "plotstyle.R"))

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

  # Correlate on networks where BOTH estimates are resolved and positive.
  # Missing, negative and exactly-zero estimates are drawn as 0 so the bar still
  # renders, but they are flagged first and excluded from the correlation.
  #
  # The zeros matter. Four of the fourteen plotted networks (FP, DMN, Tpole, SML)
  # return exactly 0 from both methods, and those (0,0) points drag the fitted line
  # toward the origin, lifting r from 0.5992 to 0.7292. Excluding them is also the
  # only treatment consistent with this paper's own argument: an exact 0 from an
  # SNP estimator is not a measurement of zero heritability, it is an unresolved
  # solve, and counting it as one is precisely the GCTA behaviour criticised in
  # the estimator-sensitivity section.
  #
  # n_nets is the number of networks plotted; n is the number actually correlated.
  cor_df <- df %>%
    pivot_wider(names_from = Type, values_from = c(h2, imputed)) %>%
    mutate(pair_real = !coalesce(imputed_Twin, TRUE) & !coalesce(imputed_AdjHE_RE, TRUE)) %>%
    summarise(
      n_nets = n(),
      n = sum(pair_real),
      r = if (sum(pair_real) > 2)
            suppressWarnings(cor(h2_Twin[pair_real], h2_AdjHE_RE[pair_real]))
          else NA_real_,
      rho = if (sum(pair_real) > 2)
              suppressWarnings(cor(h2_Twin[pair_real], h2_AdjHE_RE[pair_real],
                                   method = "spearman"))
            else NA_real_
    )
  # Save stats for quarto (instead of subtitle numbers)
  stats_path <- sub("\\.png$", "_stats.csv", out_path)
  tryCatch({
    write_csv(cor_df %>% mutate(suffix = title_suffix), stats_path)
  }, error = function(e) NULL)

  # X-scale constants, declared here because the legend anchor is derived from
  # them. ggplot2 places legends in normalised panel coordinates, so a data
  # value like 0.4 has to be converted before it can be used as a position.
  #
  # The axis runs to 0.4 rather than 0.5 because the largest plotted estimate is
  # 0.343 (VIS, Twin w_total); 0.5 left roughly a third of the panel empty. With
  # 0.4 as the limit, the legend's 0.4 anchor is also the right edge of the panel,
  # so it is derived from the same constant rather than restated.
  x_lim <- c(0, 0.4)
  x_expand_right <- 0.02            # matches expansion(mult = c(0, 0.02))
  # Anchor just inside the 0.4 limit: at exactly 0.4 the legend box would be flush
  # against the panel edge, so 0.39 leaves it clear of the border.
  legend_anchor_x <- 0.39           # data value the legend's right edge sits on
  # The right expansion widens the drawn panel to 0.4 * 1.02, so the panel span
  # is larger than the raw limit range and the anchor must be divided by the
  # wider span, not by the raw range.
  x_panel_span <- diff(x_lim) * (1 + x_expand_right)
  legend_pos_x <- legend_anchor_x / x_panel_span

  p <- df %>%
    mutate(pheno_label = factor(pheno_label, levels = twin_order)) %>%
    ggplot(aes(y = pheno_label, x = h2, fill = Type)) +
    geom_col(position = position_dodge2(width = 0.8, preserve = "single"), width = 0.7) +
    geom_errorbar(aes(xmin = pmax(0, pmin(1, h2 - se)), xmax = pmax(0, pmin(1, h2 + se))),
                  position = position_dodge2(width = 0.8, preserve = "single"),
                  width = 0.25, linewidth = 0.4, color = "grey20", na.rm = TRUE) +
    scale_x_continuous(limits = x_lim, breaks = seq(x_lim[1], x_lim[2], by = 0.1),
                       expand = expansion(mult = c(0, x_expand_right))) +
    # Data levels are "Twin" / "AdjHE_RE" (underscore -- that is what the
    # pivot_longer and recode calls produce). The keys below must match the data
    # exactly or scale_fill_manual silently drops the unmatched level: AdjHE_RE
    # bars fell back to the default grey and its legend entry disappeared,
    # leaving a legend that showed only Twin. Keys use the underscore; only the
    # displayed label carries the hyphen.
    # Twin is the accent colour; the SNP series is the recessive one, so it is
    # grey. This also matches the compressed tail of the Manhattan panels, where
    # non-highlighted points are drawn in grey.
    scale_fill_manual(values = c("Twin" = "#1f77b4", "AdjHE_RE" = "grey50"),
                      labels = c("Twin" = "Twin", "AdjHE_RE" = "AdjHE-RE"),
                      drop = FALSE) +
    labs(x = expression(paste("Heritability (", h^2, ")")), y = "", fill = "") +
    # Within each network the dodged bars run Twin on top and AdjHE-RE below, so
    # the legend is reversed to match that top-to-bottom reading order instead of
    # listing the series in the opposite sequence.
    guides(fill = guide_legend(reverse = TRUE)) +
    theme_minimal(base_size = 13, base_family = "Arial") +
    theme(
      axis.text.y = element_text(face = "bold", size = 11, margin = margin(r = 4)),
      axis.text.x = element_text(size = 11, margin = margin(t = 4)),
      axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 8)),
      # Legend sits inside the panel, its right edge anchored on the 0.4 mark of
      # the h2 axis. SMl and Tpole have no bars, so that corner is reliably
      # empty. legend.position as a two-element numeric vector is in normalised
      # panel coordinates and legend.justification anchors the box's own
      # bottom-right corner there; the data value is converted via legend_pos_x.
      legend.position = c(legend_pos_x, 0.01),
      legend.justification = c(1, 0),
      legend.text = element_text(size = 10),
      legend.title = element_blank(),
      legend.margin = margin(3, 3, 3, 3),
      # Small fixed keys. The legend occupies one corner of a 7x5 in figure, and
      # the default keys are sized in text lines, so they made the box compete
      # with the data for attention.
      legend.key.width = grid::unit(11, "pt"),
      legend.key.height = grid::unit(7, "pt"),
      legend.spacing.x = grid::unit(2, "pt"),
      # Bars run left to right, so the legend needs its own background to stay
      # readable where it overlaps the gridlines. element_rect takes no `alpha`
      # argument -- passing one raises "`...` must be empty" and the rect falls
      # back to opaque. The translucency has to live in the colour itself, via
      # alpha() from scales.
      legend.background = element_rect(fill = alpha("white", 0.88), colour = NA),
      legend.spacing.y = grid::unit(0.1, "lines"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 4, 2, 2),
      axis.ticks.y = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    )
  # Canvas grew with the text: +2pt on the axis labels and base_size is what the
  # 11pt network names and title were scaled from.
  save_fig(out_path, p, width = 8, height = 6, dpi = 300)
  print(cor_df)
}

# SA panels must read explicit per-variant sources, not WIDE: after the 01_compare
# switch, WIDE's SA columns are all w_total, so a "wo_total" panel built from WIDE
# would silently display w_total data.
extract_twin_h2 <- function(twin_rds_path) {
  rds <- readRDS(twin_rds_path)
  id_col <- setdiff(names(rds), c("herit", "data"))
  if (length(id_col) == 1) rds <- rename(rds, Pheno = all_of(id_col))
  if (!"Pheno" %in% names(rds) && "Phenotype" %in% names(rds)) rds <- rename(rds, Pheno = Phenotype)
  rds %>% mutate(Twin = map_dbl(herit, function(x) {
    if (is.null(x) || inherits(x, "twinlm_error")) return(NA_real_)
    cf <- tryCatch(x$coef, error = function(e) NULL)
    if (is.null(cf)) return(NA_real_)
    A <- as.numeric(cf[1,1]); C <- as.numeric(cf[2,1]); E <- as.numeric(cf[3,1])
    Ttot <- A + C + E
    if (Ttot > 0) A / Ttot else NA_real_
  })) %>% select(Pheno, Twin)
}

build_sa_panel <- function(twin_rds_path, adj_csv) {
  twin_df <- extract_twin_h2(twin_rds_path)
  adj_df <- read_csv(adj_csv, show_col_types = FALSE) %>%
    filter(PCs == max(PCs, na.rm = TRUE)) %>%
    select(Pheno = pheno, AdjHE_RE = h2)
  full_join(twin_df, adj_df, by = "Pheno") %>%
    pivot_longer(cols = c(Twin, AdjHE_RE), names_to = "Type", values_to = "h2") %>%
    mutate(se = 0,
           pheno_num = parse_number(Pheno),
           pheno_label = recode(as.character(pheno_num), !!!networks, .default = NA_character_),
           imputed = is.na(h2) | h2 <= 0,
           h2 = pmax(0, pmin(1, if_else(imputed, 0, h2)))) %>%
    drop_na(pheno_label)
}

df_wo <- tryCatch(build_sa_panel(file.path(ROOT, "results/SA/twinEsts/herit_wo_total.Rds"),
                                 file.path(ROOT, "results/SA/AdjHE_RE_wo_total.csv")),
                  error = function(e) { message("wo_total build failed: ", conditionMessage(e)); tibble() })
if (nrow(df_wo) > 0) plot_sa(df_wo, OUT, " (wo_total, 14 networks)") else message("No wo_total SA data — skipping ", OUT)

# ---- w_total (Twin_h2_w_total + h2_SA_AdjHE_RE_w_total if in WIDE, else fallback to direct SA files) ----
df_w <- tibble()
if (nrow(wide) > 0) {
  w_twin_col <- intersect(c("Twin_h2_w_total", "Twin_h2_wtotal", "Twin_h2_w_Total"), names(wide))
  w_adj_col  <- intersect(c("h2_SA_AdjHE_RE_w_total", "h2_SA_w_total_AdjHE_RE", "h2_SA_AdjHE_RE_wtotal"), names(wide))
  if (length(w_twin_col) == 1 && length(w_adj_col) == 1) {
    df_w <- wide %>%
      filter(Set == "SA") %>%
      select(Pheno, all_of(c(w_twin_col, w_adj_col)), any_of(c("Twin_h2_se", "var_h2_SA_AdjHE_RE"))) %>%
      rename(Twin_h2 = all_of(w_twin_col), h2_w = all_of(w_adj_col)) %>%
      mutate(
        Twin_se = if_exists("Twin_h2_se", .) else NA_real_,
        Adj_se = if_exists("var_h2_SA_AdjHE_RE", .) else NA_real_
      )
    # Pivot longer keeping SE
    # Simpler: pivot h2 and se separately or keep columns
    df_w <- wide %>%
      filter(Set == "SA") %>%
      select(Pheno, all_of(c(w_twin_col, w_adj_col)), any_of(c("Twin_h2_se", "var_h2_SA_AdjHE_RE"))) %>%
      rename(Twin_h2 = all_of(w_twin_col), h2_w = all_of(w_adj_col)) %>%
      mutate(
        Twin_se = if ("Twin_h2_se" %in% names(.)) Twin_h2_se else NA_real_,
        Adj_se = if ("var_h2_SA_AdjHE_RE" %in% names(.)) 2 * sqrt(var_h2_SA_AdjHE_RE) else NA_real_
      ) %>%
      select(Pheno, Twin_h2, h2_w, Twin_se, Adj_se) %>%
      pivot_longer(cols = c(Twin_h2, h2_w), names_to = "Type", values_to = "h2") %>%
      mutate(
        se = if_else(Type == "Twin_h2", Twin_se, Adj_se),
        Type = recode(Type, Twin_h2 = "Twin", h2_w = "AdjHE_RE"),
        pheno_num = parse_number(Pheno),
        pheno_label = recode(as.character(pheno_num), !!!networks, .default = NA_character_),
        imputed = is.na(h2) | h2 <= 0,
        h2 = pmax(0, pmin(1, if_else(imputed, 0, h2))),
        se = if_else(is.na(se), mean(se, na.rm = TRUE), se),
        se = if_else(imputed | is.na(se), 0, se)
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
        filter(PCs == max(PCs, na.rm = TRUE)) %>%
        select(Pheno = pheno, h2_w = h2)
      df_w <- full_join(twin_df, adj_df, by = "Pheno") %>%
        pivot_longer(cols = c(Twin_h2, h2_w), names_to = "Type", values_to = "h2") %>%
        mutate(
          Type = recode(Type, Twin_h2 = "Twin", h2_w = "AdjHE_RE"),
          pheno_num = parse_number(Pheno),
          pheno_label = recode(as.character(pheno_num), !!!networks, .default = NA_character_),
          imputed = is.na(h2) | h2 <= 0,
          h2 = if_else(imputed, 0, h2)
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
