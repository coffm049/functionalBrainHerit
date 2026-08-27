#!/usr/bin/env Rscript
# Pairwise-system Manhattan plots — clean recreation of 03-matrixManhattan2.py
# for all phenotype sets (gordon 352, probaConns 80, SA 17) and both methods
# (Twin, AdjHE-RE). Uses mash_twin_wide.csv (30 PCs) — no dlabel needed.
#
# For FC (gordon/probaConns) the x-axis is the edge (Pheno, e.g. o123)
# ordered by network-pair blocks when parcel->network is known, otherwise
# by h2. SA is per-network, ordered by network.
#
# Outputs: results/summary/plots/manhattan_{Set}_{Source}.png
# Run: Rscript summary/05_manhattan.R

library(tidyverse)
library(ggsci)

ROOT <- "/users/4/coffm049/papers/functionalBrainHerit"
WIDE <- file.path(ROOT, "results/summary/mash_twin_wide.csv")
PLOT_DIR <- file.path(ROOT, "results/summary/plots")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

wide <- read_csv(WIDE, show_col_types = FALSE)

# Map Set -> N and label for titles
set_labels <- c("gordon" = "Gordon", "probaConns" = "ProbaConns", "SA" = "SA")
# SA 17-network short names (as in brain_sa.py / 05-topoViz.R)
sa_networks <- c(
  "1" = "DMN", "2" = "VIS", "3" = "FP", "4" = "NET4", "5" = "DAN", "6" = "NET6",
  "7" = "VAN", "8" = "SAL", "9" = "CO", "10" = "SMD", "11" = "SML", "12" = "AUD",
  "13" = "Tpole", "14" = "MTL", "15" = "PMN", "16" = "PON", "17" = "NET17"
)

# Helper: decode FC Pheno "o123" -> edge index -> network pair when possible
# For Gordon 352 and Proba 80 we can infer parcel->network via the wide file's
# Pheno ordering alone if needed; here we simply order by h2 within each Set/Source
# and colour by h2 threshold, keeping the plot dependency-free (no dlabel).

sources <- c("Twin_h2", "h2_gordon_AdjHE_RE", "h2_proba_AdjHE_RE", "h2_SA_AdjHE_RE")

# Build long table: one row per Pheno per Source
long <- wide %>%
  select(Set, Pheno, all_of(intersect(sources, names(wide)))) %>%
  pivot_longer(cols = all_of(intersect(sources, names(wide))),
               names_to = "Source", values_to = "h2") %>%
  mutate(
    Method = case_when(
      Source == "Twin_h2" ~ "Twin",
      grepl("AdjHE", Source) ~ "AdjHE-RE",
      TRUE ~ Source
    ),
    Method = factor(Method, levels = c("Twin", "AdjHE-RE")),
    h2 = if_else(is.na(h2) | h2 < 0, 0, h2)
  ) %>%
  filter(!is.na(h2))

# For SA, map Pheno network_surfarea* -> short network for x labels
long <- long %>%
  mutate(
    pheno_num = parse_number(Pheno),
    pheno_label = if_else(Set == "SA",
                          recode(as.character(pheno_num), !!!sa_networks, .default = Pheno),
                          Pheno)
  )

# Plot per Set x Method
for (s in unique(long$Set)) {
  for (m in unique(long$Method)) {
    d <- long %>% filter(Set == s, Method == m, !is.na(h2))
    if (nrow(d) == 0) next

    # Order x by h2 (descending) to mimic Manhattan "peaks", or keep Pheno order
    # Here we keep Pheno order for SA (network order) and h2 order for FC
    if (s %in% c("gordon", "probaConns")) {
      # FC: order edges by h2 descending within each network-pair block if possible;
      # simple: order by h2 descending globally, preserve rank as x
      d <- d %>% arrange(desc(h2)) %>% mutate(x = row_number(), sig = h2 > 0.2)
    } else {
      # SA: keep network order 1..17
      d <- d %>% arrange(pheno_num) %>% mutate(x = row_number(), sig = h2 > 0.2)
    }

    # Clean title
    set_disp <- set_labels[s]
    if (is.na(set_disp)) set_disp <- s
    title <- sprintf("%s — %s (n=%d)", set_disp, m, nrow(d))

    # Threshold line (h2 > 0.2 as in 03 script's signif mean)
    thr <- 0.2

    p <- ggplot(d, aes(x = x, y = h2)) +
      geom_point(aes(color = sig), size = 0.8, alpha = 0.8) +
      scale_color_manual(values = c("TRUE" = "black", "FALSE" = "grey70"), guide = "none") +
      geom_hline(yintercept = thr, linetype = "dashed", color = "red", linewidth = 0.4) +
      labs(title = title, x = "Phenotype (ordered)", y = expression(paste("Heritability (", h^2, ")"))) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()
      ) +
      ylim(0, max(1, max(d$h2, na.rm = TRUE) * 1.05))

    # Add network-pair rug for FC if we have many points
    if (s %in% c("gordon", "probaConns") && nrow(d) > 100) {
      p <- p + geom_rug(data = d %>% filter(sig), aes(x = x), color = "black", alpha = 0.3, sides = "b")
    }

    out <- file.path(PLOT_DIR, sprintf("manhattan_%s_%s.png", s, gsub("-", "", m)))
    ggsave(out, p, width = 8, height = 3.5, dpi = 300)
    message("Wrote ", out, " [", s, " ", m, " n=", nrow(d), "]")
  }
}

# Also write a faceted overview (all Sets x Methods) for quick QC
overview <- long %>%
  filter(!is.na(h2)) %>%
  group_by(Set, Method) %>%
  mutate(x = row_number()) %>%
  ungroup()

p_all <- ggplot(overview, aes(x = x, y = h2, color = Method)) +
  geom_point(size = 0.5, alpha = 0.6) +
  facet_wrap(~Set + Method, scales = "free_x", ncol = 2) +
  scale_color_npg() +
  labs(title = "Manhattan — all phenotype sets × methods", x = "Phenotype (ordered)", y = expression(paste("h2"))) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave(file.path(PLOT_DIR, "manhattan_overview.png"), p_all, width = 10, height = 6, dpi = 300)
message("Wrote overview")
