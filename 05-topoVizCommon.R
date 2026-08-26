library(tidyverse)
library(cowplot)
library(ggsci)
scaleFUN <- function(x) sprintf("%.2f", x)
theme_set(theme_bw())
options(
    ggplot2.discrete.colour = ggsci::scale_colour_d3,
    ggplot2.discrete.fill = ggsci::scale_fill_d3,
    axis.text = element_text(face = "bold")
)


# Order of networks from Sanju
networks <- c(
    "1" = "DMN", "2" = "VIS", "3" = "FP", "5" = "DAN", "7" = "VAN",
    "8" = "SAL", "9" = "CO", "10" = "SMD", "11" = "SML", "12" = "AUD", "13" = "Tpole",
    "14" = "MTL", "15" = "PMN", "16" = "PON"
)
df <- readRDS("../results/FCs100824/Topography/SA/herit_w_total.Rds") %>%
    filter(!grepl("anthro", phenotype)) %>%
    mutate(
        pheno = parse_number(phenotype),
        pheno = factor(networks[as.character(pheno)]),
        h2 = map(herit, function(x) x$heritability[1]),
        C = map(herit, function(x) x$coef[2, 1])
    ) %>%
    unnest(C, h2)

newOrder <- df %>%
    reframe(m = max(h2, na.rm = T), .by = pheno) %>%
    arrange(m) %>%
    pull(pheno)




df %>%
    arrange(h2) %>%
    mutate(pheno = factor(pheno, levels = newOrder)) %>%
    ggplot(aes(y = pheno, x = C)) +
    geom_bar(stat = "identity", position = "dodge2") +
    # xlim(0, 0.6) +
    theme_minimal() +
    scale_color_npg() +
    scale_fill_npg() +
    labs(x = "Common Effect", y = "") +
    theme(
        axis.title.x = element_text(size = 24, face = "bold"),
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 24, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.line = element_line(linewidth = 0.5),
        legend.text = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 24, face = "bold")
    )




# USE npg theme
ggsave2("figures/topo/SAsCommon.png", dpi = 300)

df %>%
    arrange(h2) %>%
    select(phenotype, pheno, h2, C) %>%
    write_csv("formattedData/topoCommon.csv")
