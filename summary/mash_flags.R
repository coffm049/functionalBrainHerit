# Shared MASH row-filter rules.
#
# Sourced by 01_compare_mash_twin.R and 07_method_diagnostics.R. These two
# scripts must agree on which rows they keep: when the lists were duplicated
# inline they drifted, so 07 reported retention rates that no longer matched what
# the pipeline actually keeps and its "analysis" counts became untrustworthy.
# Requires dplyr to be attached (coalesce).

# Flag values meaning "this estimate is not usable".
#
# These are the complete, literal set of values MASH can emit, taken from
# MASH/src/Estimate/estimators/AdjHE.py (commit ffb888d). Only AdjHE emits a
# flag column at all -- GCTA and HEreg have no equivalent (see the methods
# section of 06_results_catalog.qmd).
#
#   2-component solve (_adjhe_2comp; random_groups null -> AdjHE-FE):
#     nonpos_det            determinant of the projected normal equations <= 0
#     neg_sigma_g           sigma_g < 0
#   3-component solve (_adjhe_3comp; random_groups = abcd_site -> AdjHE-RE):
#     ill_conditioned       cond2(XtX) > 1e10, or non-finite
#     singular              numpy.linalg.solve raised LinAlgError
#     nan_solve             NaN in the solved variance components
#     neg_sigma_g           sigma_g < 0
#   packaging (_package, both paths):
#     h2_neg_clamped        h2 < 0, clamped to 0
#     h2_gt_1_invalid       h2 > 1, returned as NaN (deliberately not clamped to 1)
#
# ill_conditioned / singular / nan_solve can ONLY come from the 3-component path,
# so their absence from an FE stream is expected rather than reassuring: the FE
# solve has no conditioning screen at all.
MASH_BAD_FLAGS <- c(
  # 2-component
  "nonpos_det",
  # 3-component
  "ill_conditioned",
  "singular",
  "nan_solve",
  # both
  "neg_sigma_g",
  "h2_neg_clamped",
  "h2_gt_1_invalid"
)

MASH_BAD_FLAG_RE <- paste(MASH_BAD_FLAGS, collapse = "|")

# Keep rows with an empty/NA flag, and rows whose flag is not a known-bad value.
# Matched as substrings because _package COMPOUNDS flags -- a row can read
# "neg_sigma_g;h2_neg_clamped" -- so exact equality would miss half of them.
# An unrecognised non-empty flag is kept on the assumption that it is
# informational; anything fatal belongs in MASH_BAD_FLAGS above.
mash_row_kept <- function(flag) {
  f <- coalesce(as.character(flag), "")
  f == "" | !grepl(MASH_BAD_FLAG_RE, f, perl = TRUE)
}
