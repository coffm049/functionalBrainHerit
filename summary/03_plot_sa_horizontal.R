#!/usr/bin/env Rscript
# Horizontal barplot for SA surface area heritability (twin vs AdjHE_RE, w/ total surface controlled)
# Organized by max twin h2 (ascending, so most heritable at top when horizontal).
# Input: results/summary/mash_twin_wide.csv  (Set == "SA", 17 network_surfarea* phenos)
# Output: results/summary/plots/SA_horizontal_barplot.png
# Run: Rscript summary/plot_sa_horizontal.R  (on acn01 with R/4.4.0)

library(tidyverse)

ROOT <- "/users/4/coffm049/papers/functionalBrainHerit"
WIDE <- file.path(ROOT, "results/summary/mash_twin_wide.csv")
OUT  <- file.path(ROOT, "results/summary/plots/SA_horizontal_barplot.png")
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)

# 17 networks — short names from 05-topoViz.R plus placeholders for 4,6,17
networks <- c(
  "1" = "DMN", "2" = "VIS", "3" = "FP", "4" = "NET4",
  "5" = "DAN", "6" = "NET6", "7" = "VAN", "8" = "SAL",
  "9" = "CO", "10" = "SMD", "11" = "SML", "12" = "AUD",
  "13" = "Tpole", "14" = "MTL", "15" = "PMN", "16" = "PON", "17" = "NET17"
)

wide <- read_csv(WIDE, show_col_types = FALSE)

# SA rows: Pheno like "network_surfarea7" (see SA/AdjHE_RE.json mpheno 1..17)
df <- wide %>%
  filter(Set == "SA") %>%
  select(Pheno, Twin_h2, h2_SA_AdjHE_RE) %>%
  pivot_longer(cols = c(Twin_h2, h2_SA_AdjHE_RE),
               names_to = "Type", values_to = "h2") %>%
  mutate(
    Type = recode(Type, Twin_h2 = "Twin", h2_SA_AdjHE_RE = "AdjHE_RE"),
    # extract numeric suffix for ordering/labeling
    pheno_num = parse_number(Pheno),
    pheno_label = recode(as.character(pheno_num), !!!networks),
    # keep original Pheno for any missing mapping
    pheno_label = if_else(is.na(pheno_label), Pheno, pheno_label),
    h2 = if_else(is.na(h2) | h2 < 0, 0, h2)
  ) %>%
  drop_na(pheno_label)

# Order y-axis by max twin h2 (ascending → most heritable at top in horizontal plot)
twin_order <- df %>%
  filter(Type == "Twin") %>%
  group_by(pheno_label) %>%
  summarise(m = max(h2, na.rm = TRUE), .groups = "drop") %>%
  arrange(m) %>%
  pull(pheno_label)

# Also compute correlation for subtitle
cor_df <- df %>%
  pivot_wider(names_from = Type, values_from = h2) %>%
  summarise(
    n = sum(!is.na(Twin) & !is.na(AdjHE_RE)),
    r = cor(Twin, AdjHE_RE, use = "pairwise.complete.obs"),
    rho = cor(Twin, AdjHE_RE, method = "spearman", use = "pairwise.complete.obs")
  )

p <- df %>%
  mutate(pheno_label = factor(pheno_label, levels = twin_order)) %>%
  ggplot(aes(y = pheno_label, x = h2, fill = Type)) +
  geom_col(position = position_dodge2(width = 0.8, preserve = "single"),
           width = 0.7) +
  scale_x_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1)) +
  labs(
    x = expression(paste("Heritability (", h^2, ")")),
    y = "",
    fill = "",
    title = "SA network surface area h2 (total surface controlled)",
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

ggsave(OUT, p, width = 8, height = 6, dpi = 300)
message("Wrote ", OUT)
print(cor_df)
