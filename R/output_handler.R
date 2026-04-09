# =============================================================================
# output_handler.R — Internal helper shared by all plot_*() functions
# =============================================================================

#' Handle session vs png output for plot lists
#' @keywords internal
.handle_output <- function(plots, output, outdir, prefix,
                            width, height, dpi) {
  if (output == "session") {
    return(plots)
  }

  # output == "png"
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  paths <- character(length(plots))
  nms   <- names(plots)

  for (i in seq_along(plots)) {
    p     <- plots[[i]]
    nm    <- nms[i]
    fname <- file.path(outdir, paste0(prefix, "_", nm, ".png"))
    ggplot2::ggsave(fname, plot = p, width = width, height = height, dpi = dpi)
    message(sprintf("  saved: %s", fname))
    paths[i] <- fname
  }

  invisible(stats::setNames(paths, nms))
}

# ── NULL coalescing operator (if not already defined) ─────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b
