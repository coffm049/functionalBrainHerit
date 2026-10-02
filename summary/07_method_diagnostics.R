#!/usr/bin/env Rscript
# Cross-method diagnostics for every phenotype set: AdjHE-FE, AdjHE-RE, GCTA,
# HEreg and Twin.
#
# This reads the RAW MASH CSVs rather than mash_twin_wide.csv on purpose. The
# wide table is produced by 01_compare_mash_twin.R, and its GCTA/HEreg columns
# have never appeared, so anything built on it silently inherits that gap. Going
# straight to the source files means a stream that 01_compare drops still shows up
# here, and the per-stream schema table below is what explains why.
#
# Run from the repo root:  Rscript summary/07_method_diagnostics.R
#
# Outputs (results/summary/):
#   method_diagnostics_streams.csv    per stream: files, rows, columns, flag/npc/N
#   method_diagnostics_summary.csv     per set x method: h2 distribution
#   method_diagnostics_pairwise.csv    per set: Pearson R between every method pair

library(tidyverse)

ROOT <- "/users/4/coffm049/papers/functionalBrainHerit"
NPC  <- 20
OUT  <- file.path(ROOT, "results", "summary")
dir.create(OUT, showWarnings = FALSE)

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

  has_flag <- "flag" %in% names(dfs) && any(!is.na(dfs$flag) & nzchar(dfs$flag))
  if (has_flag) {
    bad <- paste(c("neg_sigma_g", "h2_neg_clamped", "nonpositive_h2"), collapse = "|")
    n_flagged <- sum(!is.na(dfs$flag) & grepl(bad, dfs$flag))
    keep <- is.na(dfs$flag) | !nzchar(dfs$flag) | !grepl(bad, dfs$flag)
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
