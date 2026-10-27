#!/usr/bin/env Rscript
# Simplified horizontal barplots for SA surface area heritability (Twin vs AdjHE-RE) with SE error bars
# Input: results/summary/mash_twin_wide.csv
# Outputs: results/summary/plots/SA_w_total_horizontal_barplot.png + .pdf

library(tidyverse)
library(data.table)
source(file.path("summary", "plotstyle.R"))

WIDE <- "results/summary/mash_twin_wide.csv"
OUT_W <- "results/summary/plots/SA_w_total_horizontal_barplot.png"

dir.create(dirname(OUT_W), recursive = TRUE, showWarnings = FALSE)

networks <- c(
  "1" = "DMN", "2" = "VIS", "3" = "FP", "5" = "DAN", "7" = "VAN",
  "8" = "SAL", "9" = "CO", "10" = "SMD", "11" = "SML", "12" = "AUD",
  "13" = "Tpole", "14" = "MTL", "15" = "PMN", "16" = "PON"
)

# Read data strictly (will fail if file missing)
df_raw <- read_csv(WIDE, show_col_types = FALSE)
stopifnot(nrow(df_raw) > 0)

df_sa <- df_raw %>%
  filter(Set == "SA") %>%
  select(Pheno, Twin_h2, Twin_h2_se, h2_SA_AdjHE_RE, var_h2_SA_AdjHE_RE)

# Compute SEs and impute missing Twin SE with mean of available Twin SEs
mean_twin_se <- mean(df_sa$Twin_h2_se, na.rm = TRUE)

df_long <- df_sa %>%
  mutate(
    Twin_se = if_else(is.na(Twin_h2_se), mean_twin_se, Twin_h2_se),
    Adj_se = 2 * sqrt(var_h2_SA_AdjHE_RE)
  ) %>%
  select(Pheno, Twin_h2, h2_SA_AdjHE_RE, Twin_se, Adj_se) %>%
  pivot_longer(
    cols = c(Twin_h2, h2_SA_AdjHE_RE),
    names_to = "Type",
    values_to = "h2"
  ) %>%
  mutate(
    se = if_else(Type == "Twin_h2", Twin_se, Adj_se),
    Type = if_else(Type == "Twin_h2", "Twin", "AdjHE_RE"),
    pheno_num = parse_number(Pheno),
    pheno_label = recode(as.character(pheno_num), !!!networks, .default = NA_character_),
    imputed = is.na(h2) | h2 <= 0,
    h2 = pmax(0, pmin(1, if_else(imputed, 0, h2))),
    se = if_else(imputed | is.na(se), 0, se)
  ) %>%
  drop_na(pheno_label)

# Order y-axis by max twin h2 (ascending -> most heritable at top)
twin_order <- df_long %>%
  filter(Type == "Twin") %>%
  group_by(pheno_label) %>%
  summarise(m = max(h2, na.rm = TRUE), .groups = "drop") %>%
  arrange(m) %>%
  pull(pheno_label)

x_lim <- c(0, 0.4)
x_expand_right <- 0.02
legend_anchor_x <- 0.39
x_panel_span <- diff(x_lim) * (1 + x_expand_right)
legend_pos_x <- legend_anchor_x / x_panel_span

p <- df_long %>%
  mutate(pheno_label = factor(pheno_label, levels = twin_order)) %>%
  ggplot(aes(y = pheno_label, x = h2, fill = Type)) +
  geom_col(position = position_dodge2(width = 0.8, preserve = "single"), width = 0.7) +
  geom_errorbar(aes(xmin = pmax(0, pmin(1, h2 - se)), xmax = pmax(0, pmin(1, h2 + se))),
                position = position_dodge2(width = 0.8, preserve = "single"),
                width = 0.25, linewidth = 0.4, color = "grey20", na.rm = TRUE) +
  scale_x_continuous(limits = x_lim, breaks = seq(x_lim[1], x_lim[2], by = 0.1),
                     expand = expansion(mult = c(0, x_expand_right))) +
  scale_fill_manual(values = c("Twin" = "#1f77b4", "AdjHE_RE" = "grey50"),
                    labels = c("Twin" = "Twin", "AdjHE_RE" = "AdjHE-RE"),
                    drop = FALSE) +
  labs(x = expression(paste("Heritability (", h^2, ")")), y = "", fill = "") +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_text(face = "bold", size = 11, margin = margin(r = 4)),
    axis.text.x = element_text(size = 11, margin = margin(t = 4)),
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 8)),
    legend.position = c(legend_pos_x, 0.01),
    legend.justification = c(1, 0),
    legend.text = element_text(size = 10),
    legend.title = element_blank(),
    legend.margin = margin(3, 3, 3, 3),
    legend.key.width = grid::unit(11, "pt"),
    legend.key.height = grid::unit(7, "pt"),
    legend.spacing.x = grid::unit(2, "pt"),
    legend.background = element_rect(fill = alpha("white", 0.88), colour = NA),
    legend.spacing.y = grid::unit(0.1, "lines"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(4, 4, 2, 2),
    axis.ticks.y = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

save_fig(OUT_W, p, width = 8, height = 6, dpi = 300)
message("Successfully generated simplified SA horizontal barplot with SE bars at: ", OUT_W)
