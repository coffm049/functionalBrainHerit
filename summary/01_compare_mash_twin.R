library(tidyverse)

ROOT <- "/users/4/coffm049/papers/functionalBrainHerit"
NPC  <- 20
OUT  <- file.path(ROOT, "results", "summary")
dir.create(OUT, showWarnings = FALSE)

source(file.path(ROOT, "summary", "mash_flags.R"))

empty_mash <- tibble(Set = character(), Pheno = character(), h2 = numeric(),
                     var_h2 = numeric(), stream = character())
empty_row  <- tibble(Pheno = character(), h2 = numeric(), var_h2 = numeric())

read_mash_stream <- function(pattern, label) {
  files <- Sys.glob(file.path(ROOT, pattern))
  if (length(files) == 0) { message("No MASH files: ", label); return(empty_mash) }
  dropped_flagged <- 0
  dropped_nonpositive <- 0
  dfs <- map_dfr(files, function(f) {
    # Read every column as text and coerce explicitly. read_csv's type guessing
    # turns an h2 column that is mostly "0" but occasionally "nan"/"" into
    # character, which then fails every numeric comparison below without raising
    # an error.
    d <- suppressWarnings(read_csv(f, show_col_types = FALSE,
                                   col_types = cols(.default = col_character())))
    if ("pheno" %in% names(d)) d <- rename(d, Pheno = pheno)
    if (!"Pheno" %in% names(d) || !"h2" %in% names(d)) return(empty_row)
    d <- d %>% mutate(h2 = suppressWarnings(as.numeric(h2)))

    # AdjHE emits a flag column; GCTA and HEreg do not, so a negative sigma_g
    # clamped by the estimator is indistinguishable there from a genuine 0 and
    # would survive a flag-based filter. Where there is no flag to justify the
    # row, drop h2 <= 0 instead, so every method is summarised over the same
    # "positive genetic estimate" subset. Same spirit as the h2 > Twin_h2 censor
    # applied further down.
    if ("flag" %in% names(d)) {
      keep <- mash_row_kept(d$flag)
      dropped_flagged <<- dropped_flagged + sum(!keep)
      d <- d %>% filter(keep)
    } else {
      nonpos <- !is.na(d$h2) & d$h2 <= 0
      dropped_nonpositive <<- dropped_nonpositive + sum(nonpos)
      d <- d %>% filter(is.na(h2) | h2 > 0)
    }

    # var(h2) is absent from some streams. This used to be an unguarded
    # backtick reference, so one stream missing the column aborted the whole
    # script and left the previous mash_twin_wide.csv in place -- a stale table
    # with a fresh-looking mtime is exactly what we could not diagnose.
    if ("var(h2)" %in% names(d))
      d <- d %>% mutate(var_h2 = suppressWarnings(as.numeric(`var(h2)`)))
    else
      d <- d %>% mutate(var_h2 = NA_real_)
    if ("PCs" %in% names(d))
      d <- d %>% mutate(PCs = suppressWarnings(as.numeric(PCs)))

    d %>% select(Pheno, h2, var_h2, any_of("PCs"))
  })
  if (dropped_flagged > 0)
    message(sprintf("stream %s: dropped %d rows with a bad flag", label, dropped_flagged))
  if (dropped_nonpositive > 0) {
    message(sprintf("stream %s: no flag column; dropped %d rows with h2 <= 0", label, dropped_nonpositive))
  }
  if (nrow(dfs) == 0) { message("stream ", label, ": EMPTY after filtering"); return(empty_mash) }
  if (!"PCs" %in% names(dfs) || all(is.na(dfs$PCs))) {
    message("stream ", label, ": no/NA PCs column; assuming npc 20")
    dfs <- dfs %>% mutate(PCs = 20)
  }
  avail <- unique(dfs$PCs)
  use_npc <- if (NPC %in% avail) NPC else max(avail, na.rm = TRUE)
  if (!(NPC %in% avail)) message("stream ", label, ": npc ", NPC, " absent; using npc=", use_npc)
  dfs %>% filter(PCs == use_npc) %>% distinct(Pheno, .keep_all = TRUE) %>% mutate(stream = label)
}

streams <- list(
  SA_AdjHE_FE        = c(pattern = "results/SA/AdjHE_FE_wo_total.csv",
                         label = "SA_AdjHE_FE"),
  SA_AdjHE_RE        = c(pattern = "results/SA/AdjHE_RE_wo_total.csv",
                         label = "SA_AdjHE_RE"),
  SA_GCTA            = c(pattern = "results/SA/GCTA.csv",
                         label = "SA_GCTA"),
  SA_HEreg           = c(pattern = "results/SA/HE.csv",
                         label = "SA_HEreg"),
  gordon_AdjHE_FE    = c(pattern = "results/FCs/gordon/pconns.AdjHE.FE.*.csv",
                         label = "gordon_AdjHE_FE"),
  gordon_AdjHE_RE    = c(pattern = "results/FCs/gordon/pconns.AdjHE.RE.*.csv",
                         label = "gordon_AdjHE_RE"),
  gordon_GCTA        = c(pattern = "results/FCs/gordon/pconns.GCTA.GCTA.*.csv",
                         label = "gordon_GCTA"),
  gordon_HEreg       = c(pattern = "results/FCs/gordon/pconns.HEreg.HEreg.*.csv",
                         label = "gordon_HEreg"),
  proba_AdjHE_FE     = c(pattern = "results/FCs/probaConns/probaConns.AdjHE.FE.*.csv",
                         label = "proba_AdjHE_FE"),
  proba_AdjHE_RE     = c(pattern = "results/FCs/probaConns/probaConns.AdjHE.RE.*.csv",
                         label = "proba_AdjHE_RE"),
  proba_GCTA         = c(pattern = "results/FCs/probaConns/probaConns.GCTA.GCTA.*.csv",
                         label = "proba_GCTA"),
  proba_HEreg        = c(pattern = "results/FCs/probaConns/probaConns.HEreg.HEreg.*.csv",
                         label = "proba_HEreg")
)

mash <- map_dfr(streams, function(s) {
  # One bad stream must not abort the run. Before this guard a single unguarded
  # column reference killed the script mid-way, leaving the previous wide table
  # in place with no non-zero exit to notice.
  tryCatch(read_mash_stream(s[["pattern"]], s[["label"]]),
           error = function(e) {
             warning(sprintf("stream %s FAILED: %s", s[["label"]], conditionMessage(e)))
             empty_mash
           })
}, .id = "stream_key") %>%
  mutate(Set = case_when(grepl("^SA", stream) ~ "SA",
                         grepl("gordon", stream) ~ "gordon",
                         grepl("proba", stream) ~ "probaConns"))

# Provenance, written before the wide table: how many phenotypes each stream
# actually contributed. A later reader can then tell a current table from a
# leftover instead of having to guess from an mtime.
stream_counts <- map_dfr(names(streams), function(k)
  tibble(stream = k,
         n_files = length(Sys.glob(file.path(ROOT, streams[[k]][["pattern"]]))),
         n_rows = sum(mash$stream_key == k),
         n_pheno = n_distinct(mash$Pheno[mash$stream_key == k])))
write_csv(stream_counts, file.path(OUT, "stream_row_counts.csv"))
if (any(stream_counts$n_pheno == 0))
  warning("streams contributing no phenotypes: ",
          paste(stream_counts$stream[stream_counts$n_pheno == 0], collapse = ", "))

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
  rds <- bind_rows(lapply(files, readRDS))
  # Gate: warn loudly if mostly all-error
  if ("herit" %in% names(rds)) {
    err_count <- sum(sapply(rds$herit, function(h) inherits(h, "twinlm_error")))
    tot <- length(rds$herit)
    if (tot > 0 && err_count / tot >= 0.8)
      warning(sprintf("%s: %d/%d twin estimates are errors (%.0f%%); likely failed run", set, err_count, tot, 100*err_count/tot))
  }
  rds %>%
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

final <- final %>% mutate(across(starts_with("h2_"),
                                 ~ if_else(!is.na(Twin_h2) & .x > Twin_h2, NA_real_, .x)))

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

twin_labels <- c("SA" = "SA", "gordon" = "Gordon", "probaConns" = "ProbaConns")
for (s in c("SA", "gordon", "probaConns")) {
  d <- twin %>%
    filter(Set == s) %>%
    arrange(Phenotype) %>%
    mutate(t = row_number()) %>%
    pivot_longer(c(A, C, E), names_to = "Comp", values_to = "val") %>%
    mutate(Comp = factor(Comp, levels = c("A", "C", "E"),
                         labels = c("Additive (A)", "Common (C)", "Unique (E)")))
  # drop fully-NA phenotypes for cleaner x-axis
  keep <- d %>% group_by(t) %>% summarise(m = max(val, na.rm = TRUE), .groups = "drop") %>% filter(is.finite(m))
  d <- d %>% filter(t %in% keep$t)
  p <- ggplot(d, aes(x = t, y = val, fill = Comp)) +
    geom_area(color = "white", linewidth = 0.15, alpha = 0.9) +
    scale_fill_manual(values = c("Additive (A)" = "#1f77b4", "Common (C)" = "#ff7f0e", "Unique (E)" = "#2ca02c")) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    labs(title = paste0(twin_labels[s], " — Twin ACE composition"),
         x = "Phenotype (ordered)", y = "Variance component") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5),
          panel.grid = element_blank())
  ggsave(file.path(plot_dir, paste0("twin_composition_", s, ".png")), p, width = 7, height = 3.5, dpi = 300)
}

cat("\n=== MASH vs Twin correlations (by set) ===\n")
print(corr)
cat("\n=== Row counts ===\n")
print(final %>% count(Set, name = "n_pheno"))
cat("\nDone. Outputs in", OUT, "\n")
