# =============================================================================
# plot_envelope.R — Raw envelope of all turnover curves
# =============================================================================

#' Plot raw envelope of turnover curves faceted by variable
#'
#' @param gf_result A \code{gf_result} object.
#' @param alpha numeric. Transparency of individual lines. Default 0.02.
#' @param output character. \code{"session"} or \code{"png"}.
#' @param outdir character. Output directory.
#' @param width,height,dpi plot dimensions.
#'
#' @return Named list with one ggplot (\code{"envelope"}).
#' @export
plot_envelope <- function(gf_result,
                          alpha  = 0.02,
                          output = c("session", "png"),
                          outdir = ".",
                          width  = 10, height = 8, dpi = 300) {

  output <- match.arg(output)
  if (!inherits(gf_result, "gf_result"))
    stop("`gf_result` must be a gf_result object.")

  tf <- gf_result$turnover_full

  p <- ggplot2::ggplot(tf, ggplot2::aes(x = x, y = F_x, group = SNP_ID)) +
    ggplot2::geom_line(alpha = alpha, linewidth = 0.2) +
    ggplot2::facet_wrap(~ variable, scales = "free_x") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x     = "Environmental gradient",
      y     = "Turnover",
      title = "Turnover curve envelope by variable"
    )

  plots <- list(envelope = p)
  .handle_output(plots, output, outdir, prefix = "envelope", width, height, dpi)
}
