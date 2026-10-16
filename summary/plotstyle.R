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

# Report what the family actually resolves to, so the log records the real font
# rather than the requested one. systemfonts will substitute if Arial is absent,
# which is worth knowing before submission.
.report_font <- function() {
  if (!requireNamespace("systemfonts", quietly = TRUE)) return(invisible(NULL))
  fams <- tryCatch(systemfonts::system_fonts()$family, error = function(e) character())
  hit <- grep("^arial$", fams, ignore.case = TRUE, value = TRUE)
  if (length(hit)) {
    message("plotstyle: Arial present; figures will render in Arial.")
  } else {
    sub <- grep("liberation ?sans|helvetica|nimbus ?sans|dejavu ?sans",
                fams, ignore.case = TRUE, value = TRUE)
    message("plotstyle: WARNING Arial not installed; closest match(es): ",
            paste(utils::head(unique(sub), 3), collapse = ", "))
  }
  invisible(NULL)
}
.report_font()