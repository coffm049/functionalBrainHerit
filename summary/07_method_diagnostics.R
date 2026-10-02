#!/usr/bin/env Rscript
# Cross-method diagnostics for every phenotype set: AdjHE-FE, AdjHE-RE, GCTA,
# HEreg and Twin.
#
# This reads the RAW MASH CSVs rather than mash_twin_wide.csv on purpose. The
# wide table is produced by 01_compare_mash_twin.R, and its GCTA/HEreg columns
# never appeared because that script aborted partway through and left the
# previous table in place. Going straight to the source files means a stream that
# 01_compare drops still shows up here, and the per-stream schema table below is
# what explains why.
#
# Run from the repo root:  Rscript summary/07_method_diagnostics.R
#
# The row filter is sourced from summary/mash_flags.R, shared with
# 01_compare_mash_twin.R, so the counts below are the counts the pipeline keeps
# rather than a second, drifting copy of the rules.
#
# Run from the repo root:  Rscript summary/07_method_diagnostics.R
#
# Outputs (results/summary/):
#   method_diagnostics_streams.csv    per stream: files, rows, columns, flag/npc/N
#   method_diagnostics_summary.csv     per set x method: h2 distribution
#   method_diagnostics_pairwise.csv    per set: Pearson R between every method pair
#   method_diagnostics_divergence.csv  per set: FE vs RE agreement and flag rate
#   method_diagnostics_flag_vs_twin.csv flag rate vs Twin h2 (selection check)

library(tidyverse)

ROOT <- "/users/4/coffm049/papers/functionalBrainHerit"
NPC  <- 20
OUT  <- file.path(ROOT, "results", "summary")
dir.create(OUT, showWarnings = FALSE)

source(file.path(ROOT, "summary", "mash_flags.R"))

streams <- tribble(
  ~set,         ~method,         ~pattern,
  "gordon",     "AdjHE-FE",      "results/FCs/gordon/pconns.AdjHE.FE.*.csv",
  "gordon",     "AdjHE-RE",      "results/FCs/gordon/pconns.AdjHE.RE.*.csv",
  "gordon",     "GCTA",          "results/FCs/gordon/pconns.GCTA.GCTA.*.csv",
  "gordon",     "HEreg",         "results/FCs/gordon/pconns.HEreg.HEreg.*.csv",
  "probaConns", "AdjHE-FE",      "results/FCs/probaConns/probaConns.AdjHE.FE.*.csv",
  "probaConns", "AdjHE-RE",      "results/FCs/probaConns/probaConns.AdjHE.RE.*.csv",
  "probaConns", "GCTA",          "results/FCs/probaConns/probaConns.GCTA.GCTA.*.csv",
  "probaConns", "HEreg",         "results/FCs/probaConns/probaConns.HEreg.HEreg.*.csv",
  # SA is one CSV per run. The wo_total files are the ones the report uses;
  # the w_total files are the same spec with totalNetworkSurface as a covariate.
  "SA", "AdjHE-FE",      "results/SA/AdjHE_FE_wo_total.csv",
  "SA", "AdjHE-FE-wtot", "results/SA/AdjHE_FE.csv",
  "SA", "AdjHE-RE",      "results/SA/AdjHE_RE_wo_total.csv",
  "SA", "AdjHE-RE-wtot", "results/SA/AdjHE_RE.csv",
  "SA", "GCTA",          "results/SA/GCTA_wo_total.csv",
  "SA", "GCTA-wtot",     "results/SA/GCTA.csv",
  "SA", "HEreg",         "results/SA/HE.csv"
)

# Mirrors read_mash_stream() in 01_compare_mash_twin.R so the "analysis" counts
# here line up with what the pipeline actually keeps: drop rows carrying a
# bad flag, and for streams with no flag column at all drop h2 <= 0 instead.
read_stream <- function(pattern) {
  files <- Sys.glob(file.path(ROOT, pattern))
  if (length(files) == 0) {
    return(list(status = "MISSING", data = NULL, schema = NA_character_,
                n_flagged = NA_integer_, n_nonpos = NA_integer_))
  }
  first <- suppressWarnings(read_csv(files[1], show_col_types = FALSE, n_max = 5))
  schema <- paste(names(first), collapse = ",")

  dfs <- suppressWarnings(map_dfr(files, function(f)
    read_csv(f, show_col_types = FALSE, col_types = cols(.default = col_character()))))

  if (!"h2" %in% names(dfs))
    return(list(status = "NO_H2", data = NULL, schema = schema,
                n_flagged = NA_integer_, n_nonpos = NA_integer_))

  dfs <- dfs %>% mutate(h2 = suppressWarnings(as.numeric(h2)))
  # Guard every optional column: a stream with an unexpected schema is exactly
  # what this script exists to surface, so it must report rather than crash.
  if ("PCs" %in% names(dfs))
    dfs <- dfs %>% mutate(PCs = suppressWarnings(as.numeric(PCs))) %>%
      filter(is.na(PCs) | PCs == NPC)
  if ("N" %in% names(dfs)) dfs <- dfs %>% mutate(N = suppressWarnings(as.numeric(N)))
  if ("pheno" %in% names(dfs)) dfs <- rename(dfs, Pheno = pheno)

  npc_avail <- if ("PCs" %in% names(dfs)) sort(unique(dfs$PCs[!is.na(dfs$PCs)])) else numeric()

  if ("flag" %in% names(dfs)) {
    keep <- mash_row_kept(dfs$flag)
    n_flagged <- sum(!keep)
  } else {
    n_flagged <- 0L
    keep <- is.na(dfs$h2) | dfs$h2 > 0
  }
  n_nonpos <- sum(!is.na(dfs$h2) & dfs$h2 <= 0)

  # Without a phenotype column the rows cannot be aligned to Twin or to another
  # method, so report the counts and stop rather than failing downstream.
  if (!"Pheno" %in% names(dfs))
    return(list(status = "NO_PHENO", data = NULL, schema = schema,
                npc_avail = npc_avail, n_flagged = n_flagged, n_nonpos = n_nonpos))

  kept <- dfs %>% filter(keep) %>% filter(!is.na(Pheno))

  list(status = "OK",
       data = kept %>% select(Pheno, h2, N),
       schema = schema,
       npc_avail = npc_avail,
       n_flagged = n_flagged,
       n_nonpos = n_nonpos)
}

# Re-read a stream keeping the flag column, so we can ask whether the flag fires
# preferentially on phenotypes with low Twin h2. read_stream() has already thrown
# the flagged rows away, so this cannot be recovered from the data it returns.
read_flagged_twin <- function(files, tw) {
  d <- suppressWarnings(map_dfr(files, function(f)
    read_csv(f, show_col_types = FALSE, col_types = cols(.default = col_character()))))
  if (!all(c("pheno", "h2") %in% names(d))) return(NULL)
  d <- d %>% rename(Pheno = pheno) %>%
    mutate(h2 = suppressWarnings(as.numeric(h2))) %>%
    filter(!is.na(h2), !is.na(Pheno))
  if (!nrow(d)) return(NULL)
  flagged <- if ("flag" %in% names(d)) !mash_row_kept(d$flag) else rep(FALSE, nrow(d))
  j <- inner_join(tibble(Pheno = d$Pheno, flagged = flagged), tw, by = "Pheno")
  j %>% filter(!is.na(Twin_h2)) %>%
    mutate(bin = cut(Twin_h2, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE))
}

cat("=== 1. STREAM INVENTORY AND SCHEMA ===\n\n")
inv <- list(); long <- list(); dist <- list()

for (i in seq_len(nrow(streams))) {
  s <- streams[i, ]
  r <- read_stream(s$pattern)
  label <- paste(s$set, s$method, sep = " / ")

  inv[[i]] <- tibble(
    set = s$set, method = s$method, status = r$status,
    n_files = length(Sys.glob(file.path(ROOT, s$pattern))),
    columns = if (is.na(r$schema)) NA_character_ else r$schema,
    has_pheno = if (is.na(r$schema)) NA else grepl("pheno", r$schema),
    has_h2    = if (is.na(r$schema)) NA else grepl("(^|,)h2(,|$)", r$schema),
    has_varh2 = if (is.na(r$schema)) NA else grepl("var\\(h2\\)", r$schema),
    has_flag  = if (is.na(r$schema)) NA else grepl("flag", r$schema),
    npc_avail = if (is.null(r$npc_avail)) NA else paste(r$npc_avail, collapse = "/"),
    n_flagged = r$n_flagged, n_h2_nonpos = r$n_nonpos
  )

  if (r$status != "OK") {
    cat(sprintf("%-28s %-8s  files=%d\n", label, r$status, inv[[i]]$n_files))
    next
  }

  d <- r$data
  long[[length(long) + 1]] <- d %>% mutate(set = s$set, method = s$method)

  n_min <- if ("N" %in% names(d)) min(d$N, na.rm = TRUE) else NA_real_
  n_max <- if ("N" %in% names(d)) max(d$N, na.rm = TRUE) else NA_real_

  dist[[length(dist) + 1]] <- tibble(
    set = s$set, method = s$method, n = sum(!is.na(d$h2)),
    frac_pos = mean(d$h2 > 0, na.rm = TRUE),
    median = median(d$h2, na.rm = TRUE),
    q25 = quantile(d$h2, 0.25, na.rm = TRUE),
    q75 = quantile(d$h2, 0.75, na.rm = TRUE),
    mean = mean(d$h2, na.rm = TRUE), max = max(d$h2, na.rm = TRUE),
    N_min = n_min, N_max = n_max
  )

  cat(sprintf("%-28s %-8s  n=%5d  med=%.4f  frac>0=%.3f  N=[%s,%s]\n",
              label, "OK", sum(!is.na(d$h2)), median(d$h2, na.rm = TRUE),
              mean(d$h2 > 0, na.rm = TRUE), n_min, n_max))
}

inv  <- bind_rows(inv)
long <- bind_rows(long)
dist <- bind_rows(dist)

write_csv(inv,  file.path(OUT, "method_diagnostics_streams.csv"))
write_csv(dist, file.path(OUT, "method_diagnostics_summary.csv"))

cat("\n=== 2. h2 DISTRIBUTION BY SET AND METHOD ===\n\n")
print(dist %>% mutate(across(where(is.numeric), ~round(.x, 4))), n = Inf)

# Twin h2 comes from the wide table, which is trustworthy for that column: the
# missing streams are all SNP methods. Printed with its mtime so a stale file is
# obvious rather than silent.
cat("\n=== 3. TWIN REFERENCE ===\n\n")
wide_path <- file.path(OUT, "mash_twin_wide.csv")
if (file.exists(wide_path)) {
  cat("mash_twin_wide.csv mtime:", format(file.info(wide_path)$mtime), "\n")
  twin <- suppressWarnings(read_csv(wide_path, show_col_types = FALSE)) %>%
    select(set = Set, Pheno, Twin_h2)
  print(twin %>% filter(!is.na(Twin_h2)) %>% group_by(set) %>%
          summarise(n = n(), median = round(median(Twin_h2), 4),
                    q25 = round(quantile(Twin_h2, .25), 4),
                    q75 = round(quantile(Twin_h2, .75), 4), .groups = "drop"))
  long <- bind_rows(long, twin %>% transmute(set, Pheno, h2 = Twin_h2, method = "Twin"))
} else {
  twin <- NULL
  cat("mash_twin_wide.csv not found - run 01_compare_mash_twin.R first\n")
}

cat("\n=== 4. PAIRWISE PEARSON R WITHIN SET ===\n\n")
pw <- list()
for (s in unique(long$set)) {
  sub <- long %>% filter(set == s) %>% filter(!is.na(h2)) %>%
    select(method, Pheno, h2) %>% distinct(method, Pheno, .keep_all = TRUE)
  ms <- sort(unique(sub$method))
  if (length(ms) < 2) next
  w <- sub %>% pivot_wider(names_from = method, values_from = h2)
  for (i in seq_along(ms)) for (j in seq_along(ms)) if (i < j) {
    a <- ms[i]; b <- ms[j]
    if (!all(c(a, b) %in% names(w))) next
    cc <- complete.cases(w[[a]], w[[b]])
    pw[[length(pw) + 1]] <- tibble(
      set = s, method_a = a, method_b = b, n_complete = sum(cc),
      pearson_r = if (sum(cc) >= 3) cor(w[[a]][cc], w[[b]][cc]) else NA_real_)
  }
}
pw <- bind_rows(pw)
if (!is.null(pw)) {
  write_csv(pw, file.path(OUT, "method_diagnostics_pairwise.csv"))
  print(pw %>% mutate(pearson_r = round(pearson_r, 3)), n = Inf)
}

# The AdjHE-FE and AdjHE-RE panels are shown side by side, but they retain very
# different fractions of phenotypes (flag rates differ by ~2x), so they are not
# estimating the same thing over the same rows. Two things need separating before
# the pair is read as a method contrast:
#   1. how much of the disagreement is a pure scale factor (RE ~ k x FE), versus
#      a genuine re-ranking of which connections look heritable;
#   2. whether the flag filter is selecting on the outcome.
# (2) matters most: if flagged rows are systematically the low-Twin-h2 ones, then
# FE and RE are not "AdjHE with different site handling" but different subsets
# selected on h2 itself, and no amount of side-by-side plotting fixes that.
cat("\n=== 5. AdjHE-FE vs AdjHE-RE: SCALE vs RE-RANKING ===\n\n")

inv2 <- inv %>% select(set, method, n_files, npc_avail, n_flagged, n_h2_nonpos)
print(inv2, n = Inf)

div <- list()
for (s in sort(unique(long$set))) {
  sub <- long %>% filter(set == s, !is.na(h2)) %>%
    select(method, Pheno, h2) %>% distinct(method, Pheno, .keep_all = TRUE)
  fe <- sub %>% filter(method == "AdjHE-FE") %>% transmute(Pheno, fe = h2)
  re <- sub %>% filter(method == "AdjHE-RE") %>% transmute(Pheno, re = h2)
  if (nrow(fe) < 3 || nrow(re) < 3) next
  j <- inner_join(fe, re, by = "Pheno")
  if (nrow(j) < 3) next
  cc <- complete.cases(j$fe, j$re)
  j <- j[cc, ]
  if (nrow(j) < 3) next
  # OLS of RE on FE: slope ~1 means "same ranking, different units".
  fit <- tryCatch(lm(re ~ fe, data = j), error = function(e) NULL)
  div[[length(div) + 1]] <- tibble(
    set = s, n_complete = nrow(j),
    pearson_r  = cor(j$fe, j$re),
    spearman_r = cor(j$fe, j$re, method = "spearman"),
    ols_slope_re_on_fe = if (is.null(fit)) NA_real_ else unname(coef(fit)[2]),
    median_ratio_re_fe = median(j$re / j$fe[j$fe > 0]),
    # Share of the FE variance a pure rescale would explain: r^2. Near 1 means
    # the panels differ only in scale; near 0 means they rank connections
    # differently and a shared x-axis invites a false read.
    r2_after_rescale = cor(j$fe, j$re)^2)
}
div <- bind_rows(div)
if (!is.null(div)) {
  write_csv(div, file.path(OUT, "method_diagnostics_divergence.csv"))
  print(div %>% mutate(across(where(is.numeric), ~round(.x, 4))), n = Inf)
}

# Flag rate as a function of Twin h2. If the AdjHE flag filter is not independent
# of the phenotype, the retained subset is outcome-selected and the FE/RE
# difference above cannot be read as a method contrast.
cat("\n=== 6. IS THE FLAG FILTER SELECTING ON h2? ===\n\n")
if (!is.null(twin) && nrow(twin)) {
  flagsel <- list()
  for (s in sort(unique(long$set))) {
    for (m in c("AdjHE-FE", "AdjHE-RE")) {
      pat <- streams %>% filter(set == s, method == m) %>% pull(pattern)
      if (!length(pat)) next
      files <- Sys.glob(file.path(ROOT, pat))
      if (!length(files)) next
      r <- read_flagged_twin(files, twin %>% filter(set == s) %>% select(Pheno, Twin_h2))
      if (!is.null(r)) {
        r$set <- s; r$method <- m
        flagsel[[length(flagsel) + 1]] <- r
      }
    }
  }
  fs <- bind_rows(flagsel)
  if (!is.null(fs) && nrow(fs)) {
    agg <- fs %>% group_by(set, method, bin) %>%
      summarise(n = n(), frac_flagged = round(mean(flagged), 4),
                median_twin_h2 = round(median(Twin_h2), 4), .groups = "drop")
    write_csv(agg, file.path(OUT, "method_diagnostics_flag_vs_twin.csv"))
    print(agg, n = Inf)
    # A flag that fires more often on low-Twin-h2 phenotypes means retention is
    # outcome-dependent; report the slope so it is not left to the eye.
    cat("\nflag rate vs Twin h2 (logistic slope; negative = flags low-h2 rows):\n")
    slopes <- fs %>% group_by(set, method) %>%
      summarise(n = n(),
                slope = {
                  dd <- data.frame(flagged = flagged, Twin_h2 = Twin_h2)
                  if (nrow(dd) >= 20 && length(unique(dd$Twin_h2)) > 2)
                    tryCatch(unname(coef(glm(flagged ~ Twin_h2, family = binomial,
                                             data = dd))[2]),
                             error = function(e) NA_real_) else NA_real_
                },
                .groups = "drop")
    print(slopes, n = Inf)
  }
} else {
  cat("no Twin reference; skipping flag-selection check\n")
}

cat("\nDone. Outputs in", OUT, "\n")
