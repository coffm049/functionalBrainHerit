# Shared ggplot styling for the results catalog.
#
# Imaging Neuroscience requires Arial (or Helvetica) on every figure, so the font
# is declared explicitly rather than left to the device default. This matters:
# the default ggsave() device silently ignores theme_minimal(base_family=), and
# on this Windows host the device's default sans happens to be Arial, so the
# setting looks like it works. On the HPC the default is whatever fontconfig
# picks, which is generally not Arial. ragg resolves the family through
# systemfonts, so base_family and theme(text=) are actually honoured there.
#
# Source this before building any plot, and pass device = agg_dev to ggsave().

BASE_FAMILY <- "Arial"

# ragg is used when available because it honours the requested family. On a host
# without ragg we fall back to the default device and say so, since the font then
# depends on the device default and must be checked by eye.
.has_ragg <- requireNamespace("ragg", quietly = TRUE)

if (.has_ragg) {
  agg_dev <- function(filename, width, height, units = "in", res = 300, ...) {
    ragg::agg_png(filename, width = width, height = height,
                  units = units, res = res, background = "white", ...)
  }
} else {
  agg_dev <- function(filename, width, height, units = "in", res = 300, ...) {
    grDevices::png(filename, width = width, height = height,
                   units = units, res = res, type = "cairo", bg = "white", ...)
  }
}

# Applied globally so bare theme_minimal() calls inherit the family too.
suppressPackageStartupMessages(library(ggplot2))
theme_set(theme_minimal(base_family = BASE_FAMILY))
theme_update(base_family = BASE_FAMILY,
             text = element_text(family = BASE_FAMILY))

# Report what the family will actually resolve to, so the log records the real
# font rather than the requested one.
#
# Arial is proprietary and is generally NOT installed on a Linux cluster, and
# systemfonts::match_fonts() cannot be used to detect that: given a family that
# does not exist at all it still returns a fallback font rather than failing. So
# presence is checked directly against the enumerated family list instead.
#
# What matters for the journal is that the rendered glyphs match Arial.
# Liberation Sans / Arimo are metric-compatible with Arial and are what
# fontconfig substitutes for it by default, so those are reported as acceptable.
# Only a total absence of both Arial and a substitute is a real failure, because
# then the figures fall back to DejaVu Sans, which does not meet the rule.
.report_font <- function() {
  if (!requireNamespace("systemfonts", quietly = TRUE)) {
    message("plotstyle: systemfonts not installed, cannot verify the font.")
    return(invisible(NULL))
  }
  fams <- tryCatch(unique(systemfonts::system_fonts()$family),
                   error = function(e) character())
  has_arial <- any(tolower(fams) == "arial")

  # Metric-compatible substitutes, best first.
  subs <- c("liberation sans", "arimo", "nimbus sans", "helvetica")
  # which() over subs (not over fams) so the indices address the candidate list,
  # and so a single match cannot recycle and pick the wrong entry.
  hit <- subs[which(subs %in% tolower(fams))]

  if (has_arial) {
    message("plotstyle: Arial present; figures will render in Arial.")
  } else if (length(hit)) {
    message("plotstyle: Arial itself is not installed. fontconfig will render ",
            hit[1], ", which is metric-compatible with Arial and is the normal ",
            "Linux substitute. The font embedded in the output will be named ",
            hit[1], ", not Arial.")
  } else {
    message("plotstyle: WARNING neither Arial nor a metric-compatible ",
            "substitute (", paste(subs, collapse = "/"), ") is installed. ",
            "Figures will fall back to the default sans (usually DejaVu Sans), ",
            "which does NOT meet the Arial/Helvetica requirement.")
  }
  invisible(NULL)
}
.report_font()