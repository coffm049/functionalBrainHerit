# Shared MASH row-filter rules.
#
# Sourced by 01_compare_mash_twin.R and 07_method_diagnostics.R. These two
# scripts must agree on which rows they keep: when the lists were duplicated
# inline they drifted, so 07 reported retention rates that no longer matched what
# the pipeline actually keeps and its "analysis" counts became untrustworthy.
# Requires dplyr to be attached (coalesce).

# Flag values meaning "this estimate is not usable". Matched as substrings via
# grepl, so a stream that spells a flag slightly differently is still caught --
# but a genuinely new flag has to be added here deliberately rather than
# silently changing one script's numbers and not the other's.
MASH_BAD_FLAGS <- c(
  "neg_sigma_g",
  "h2_neg_clamped",
  "nonpositive_h2",
  "nonpos_det",
  "ill_conditioned",
  "singular",
  "nan_solve",
  "h2_gt_1_invalid"
)

MASH_BAD_FLAG_RE <- paste(MASH_BAD_FLAGS, collapse = "|")

# Keep rows with an empty/NA flag, and rows whose flag is not a known-bad value.
# An unrecognised non-empty flag is kept on the assumption that it is
# informational; anything fatal belongs in MASH_BAD_FLAGS above.
mash_row_kept <- function(flag) {
  f <- coalesce(as.character(flag), "")
  f == "" | !grepl(MASH_BAD_FLAG_RE, f, perl = TRUE)
}
