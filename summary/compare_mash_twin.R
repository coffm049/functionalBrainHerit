library(tidyverse)

ROOT <- "/users/4/coffm049/papers/functionalBrainHerit"
NPC  <- 30
OUT  <- file.path(ROOT, "results", "summary")
dir.create(OUT, showWarnings = FALSE)

empty_mash <- tibble(Set = character(), Pheno = character(), h2 = numeric(),
                     var_h2 = numeric(), stream = character())
empty_row  <- tibble(Pheno = character(), h2 = numeric(), var_h2 = numeric())

read_mash_stream <- function(pattern, label) {
  files <- Sys.glob(file.path(ROOT, pattern))
  if (length(files) == 0) { message("No MASH files: ", label); return(empty_mash) }
  dfs <- map_dfr(files, function(f) {
    d <- suppressWarnings(read_csv(f, show_col_types = FALSE))
    if ("pheno" %in% names(d)) d <- rename(d, Pheno = pheno)
    if (!"Pheno" %in% names(d) || !"h2" %in% names(d)) return(empty_row)
    d <- d %>% mutate(var_h2 = as.numeric(`var(h2)`)) %>%
      select(Pheno, h2, var_h2, any_of("PCs"))
    d
  })
  if (nrow(dfs) == 0) return(empty_mash)
  if (!"PCs" %in% names(dfs) || all(is.na(dfs$PCs))) {
    message("stream ", label, ": no/NA PCs column; assuming alternating npc 0/30")
    dfs <- dfs %>% mutate(PCs = rep(c(0, 30), length.out = nrow(dfs)))
  }
  avail <- unique(dfs$PCs)
  use_npc <- if (NPC %in% avail) NPC else max(avail, na.rm = TRUE)
  if (!(NPC %in% avail)) message("stream ", label, ": npc ", NPC, " absent; using npc=", use_npc)
  dfs %>% filter(PCs == use_npc) %>% distinct(Pheno, .keep_all = TRUE) %>% mutate(stream = label)
}

streams <- list(
  SA_GCTA         = c(pattern = "results/SA/GCTA_wo_total.csv",
                      label = "SA_GCTA"),
  SA_AdjHE_RE     = c(pattern = "results/SA/AdjHE_RE_wo_total.csv",
                      label = "SA_AdjHE_RE"),
  gordon_AdjHE_FE = c(pattern = "results/FCs/gordon/pconns.AdjHE.FE.*.csv",
                      label = "gordon_AdjHE_FE"),
  gordon_AdjHE_RE = c(pattern = "results/FCs/gordon/pconns.AdjHE.RE.*.csv",
                      label = "gordon_AdjHE_RE"),
  proba_AdjHE_FE  = c(pattern = "results/FCs/probaConns/probaConns.AdjHE.FE.*.csv",
                      label = "proba_AdjHE_FE"),
  proba_AdjHE_RE  = c(pattern = "results/FCs/probaConns/probaConns.AdjHE.RE.*.csv",
                      label = "proba_AdjHE_RE"),
  gordon_GCTA     = c(pattern = "results/FCs/gordon/pconns.GCTA.FE.*.csv",
                      label = "gordon_GCTA"),
  proba_GCTA      = c(pattern = "results/FCs/probaConns/probaConns.GCTA.FE.*.csv",
                      label = "proba_GCTA"),
  gordon_HEREG    = c(pattern = "results/FCs/gordon/pconns.HEreg.FE.*.csv",
                      label = "gordon_HEREG"),
  proba_HEREG     = c(pattern = "results/FCs/probaConns/probaConns.HEreg.FE.*.csv",
                      label = "proba_HEREG")
)

mash <- map_dfr(streams, function(s) read_mash_stream(s[["pattern"]], s[["label"]]),
                .id = "stream_key") %>%
  mutate(Set = case_when(grepl("^SA", stream) ~ "SA",
                         grepl("gordon", stream) ~ "gordon",
                         grepl("proba", stream) ~ "probaConns"))

extract_twin <- function(x) {
  if (is.null(x) || inherits(x, "twinlm_error"))
    return(tibble(A = NA_real_, C = NA_real_, E = NA_real_,
                  Twin_h2 = NA_real_, Twin_h2_se = NA_real_))
  cf <- tryCatch(x$coef, error = function(e) NULL)
  if (is.null(cf) || !is.matrix(cf) || nrow(cf) < 3)
    return(tibble(A = NA_real_, C = NA_real_, E = NA_real_,
                  Twin_h2 = NA_real_, Twin_h2_se = NA_real_))
  A <- as.numeric(cf[1, 1]); C <- as.numeric(cf[2, 1]); E <- as.numeric(cf[3, 1])
  vA <- as.numeric(cf[1, 2])^2; vC <- as.numeric(cf[2, 2])^2; vE <- as.numeric(cf[3, 2])^2
  Ttot <- A + C + E
  Twin_h2 <- if_else(Ttot > 0, A / Ttot, NA_real_)
  var_h2 <- ((C + E)^2 * vA + A^2 * vC + A^2 * vE) / (Ttot^4)
  tibble(A = A, C = C, E = E, Twin_h2 = Twin_h2, Twin_h2_se = sqrt(var_h2))
}

read_twin_set <- function(patterns, set) {
  files <- unlist(lapply(patterns, function(p) Sys.glob(file.path(ROOT, p))))
  if (length(files) == 0) { message("No twin RDS for ", set); return(tibble()) }
  bind_rows(lapply(files, readRDS)) %>%
    { .tmp <- .
      id_col <- setdiff(names(.tmp), c("herit", "data"))
      if (length(id_col) != 1) stop(paste("twin id col ambiguous:", paste(names(.tmp), collapse = ",")))
      rename(.tmp, Phenotype = all_of(id_col)) } %>%
    mutate(est = map(herit, extract_twin)) %>%
    select(Phenotype, est) %>% unnest(est) %>% mutate(Set = set)
}

twin <- bind_rows(
  read_twin_set("results/SA/twinEsts/herit_wo_total.Rds", "SA"),
  read_twin_set(c("results/FCs/gordon/herit_*.Rds",
                  "results/FCs/gordon/twinEstResults/herit_*.Rds"), "gordon"),
  read_twin_set(c("results/FCs/probaConns/herit_*.Rds",
                  "results/FCs/probaConns/twinEstResults/herit_*.Rds"), "probaConns")
)

mash_cols <- mash %>% split(.$stream) %>% imap(function(d, nm) {
  d <- d %>% select(Set, Pheno, h2, var_h2) %>%
    rename(!!sym(paste0("h2_", nm)) := h2,
           !!sym(paste0("var_h2_", nm)) := var_h2)
  d
})
wide <- reduce(mash_cols, full_join, by = c("Set", "Pheno"))

final <- twin %>% rename(Pheno = Phenotype) %>%
  full_join(wide, by = c("Set", "Pheno")) %>%
  arrange(Set, Pheno)
write_csv(final, file.path(OUT, "mash_twin_wide.csv"))

long <- final %>% pivot_longer(cols = starts_with("h2_") | all_of("Twin_h2"),
                               names_to = "Source", values_to = "h2") %>%
  mutate(Type = if_else(Source == "Twin_h2", "Twin", "MASH")) %>%
  arrange(Set, Pheno, Type, Source)
write_csv(long, file.path(OUT, "unified_long.csv"))

melt <- final %>% select(Set, Pheno, starts_with("h2_"), Twin_h2) %>%
  pivot_longer(starts_with("h2_"), names_to = "MashMethod", values_to = "Mash_h2")

corr <- melt %>% group_by(Set, MashMethod) %>%
  summarise(n = sum(!is.na(Mash_h2) & !is.na(Twin_h2)),
            Pearson_r = cor(Mash_h2, Twin_h2, method = "pearson", use = "pairwise.complete.obs"),
            Spearman_rho = cor(Mash_h2, Twin_h2, method = "spearman", use = "pairwise.complete.obs"),
            .groups = "drop")
write_csv(corr, file.path(OUT, "correlations_by_set.csv"))

corr_all <- melt %>% group_by(MashMethod) %>%
  summarise(n = sum(!is.na(Mash_h2) & !is.na(Twin_h2)),
            Pearson_r = cor(Mash_h2, Twin_h2, method = "pearson", use = "pairwise.complete.obs"),
            Spearman_rho = cor(Mash_h2, Twin_h2, method = "spearman", use = "pairwise.complete.obs"),
            .groups = "drop")
write_csv(corr_all, file.path(OUT, "correlations_overall.csv"))

plot_dir <- file.path(OUT, "plots"); dir.create(plot_dir, showWarnings = FALSE)
for (sm in unique(melt$MashMethod)) {
  d <- melt %>% filter(MashMethod == sm, !is.na(Mash_h2), !is.na(Twin_h2))
  p <- ggplot(d, aes(x = Twin_h2, y = Mash_h2)) + geom_point(alpha = 0.3, size = 0.5) +
    geom_abline(intercept = 0, slope = 1, color = "red") + facet_wrap(~Set) +
    labs(title = sm, x = "Twin h2 = A/(A+C+E)", y = "MASH h2") + theme_minimal()
  ggsave(file.path(plot_dir, paste0("scatter_", sm, ".png")), p, width = 6, height = 4)
  d2 <- d %>% mutate(diff = Mash_h2 - Twin_h2, avg = (Mash_h2 + Twin_h2) / 2)
  p2 <- ggplot(d2, aes(x = avg, y = diff)) + geom_point(alpha = 0.3, size = 0.5) +
    geom_hline(yintercept = mean(d2$diff, na.rm = TRUE), color = "red") + facet_wrap(~Set) +
    labs(title = paste("Bland-Altman", sm), x = "Mean h2", y = "MASH - Twin") + theme_minimal()
  ggsave(file.path(plot_dir, paste0("blandaltman_", sm, ".png")), p2, width = 6, height = 4)
}

for (s in c("SA", "gordon", "probaConns")) {
  d <- twin %>% filter(Set == s) %>% mutate(t = 1:n()) %>%
    pivot_longer(c(A, C, E), names_to = "Comp", values_to = "val")
  p <- ggplot(d, aes(x = t, y = val, fill = Comp)) + geom_area() +
    labs(title = paste("Twin ACE composition", s)) + theme_minimal()
  ggsave(file.path(plot_dir, paste0("twin_composition_", s, ".png")), p, width = 6, height = 4)
}

cat("\n=== MASH vs Twin correlations (by set) ===\n")
print(corr)
cat("\n=== Row counts ===\n")
print(final %>% count(Set, name = "n_pheno"))
cat("\nDone. Outputs in", OUT, "\n")
