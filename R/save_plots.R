# =============================================================================
# save_plots.R — Pipe-friendly plot saver
# =============================================================================

#' Save a named list of ggplots to disk
#'
#' Designed to be used at the end of a pipe after any \code{plot_*()} call
#' that returned with \code{output = "session"}.
#'
#' @param plots Named list of ggplot objects (as returned by \code{plot_*()}).
#' @param plotname character. Base name for files. Each plot is saved as
#'   \code{<plotname>_<name>.<extension>}. If \code{""}, uses the list name only.
#' @param extension character. File extension without dot.
#'   Passed to \code{ggplot2::ggsave()}. Default \code{"png"}.
#' @param outdir character. Output directory. Created if it does not exist.
#'   Default \code{"."}.
#' @param width,height,dpi Dimensions passed to \code{ggplot2::ggsave()}.
#'   Default 8 × 6 at 300 dpi.
#'
#' @return The input \code{plots} list, invisibly (pipe-friendly).
#' @export
save_plots <- function(plots,
                       plotname  = "",
                       extension = "png",
                       outdir    = ".",
                       width     = 8, height = 6, dpi = 300) {

  if (!is.list(plots))
    stop("`plots` must be a named list of ggplot objects.")

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  purrr::iwalk(plots, function(p, nm) {
    if (!inherits(p, "ggplot")) {
      warning(sprintf("Element '%s' is not a ggplot, skipping.", nm))
      return(invisible(NULL))
    }
    stem  <- if (nchar(plotname) > 0) paste0(plotname, "_", nm) else nm
    fname <- file.path(outdir, paste0(stem, ".", extension))
    ggplot2::ggsave(fname, plot = p, width = width, height = height, dpi = dpi)
    message(sprintf("  saved: %s", fname))
  })

  invisible(plots)
}
