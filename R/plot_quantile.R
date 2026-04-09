# =============================================================================
# plot_quantile.R — Quantile envelope of turnover curves
# =============================================================================

#' Plot quantile envelope of turnover curves
#'
#' Draws ribbon bands at Q5/Q95 (light) and Q25/Q75 (dark) with median line.
#'
#' @param gf_result A \code{gf_result} object.
#' @param probs numeric vector of length 5: lower outer, lower inner, median,
#'   upper inner, upper outer. Default \code{c(0.05, 0.25, 0.50, 0.75, 0.95)}.
#' @param output character. \code{"session"} or \code{"png"}.
#' @param outdir character. Output directory.
#' @param width,height,dpi plot dimensions.
#'
#' @return Named list with one ggplot (\code{"quantile"}).
#' @export
plot_quantile <- function(gf_result,
                          probs  = c(0.05, 0.25, 0.50, 0.75, 0.95),
                          output = c("session", "png"),
                          outdir = ".",
                          width  = 10, height = 8, dpi = 300) {

  output <- match.arg(output)
  if (!inherits(gf_result, "gf_result"))
    stop("`gf_result` must be a gf_result object.")
  if (length(probs) != 5)
    stop("`probs` must have exactly 5 values.")

  tf <- gf_result$turnover_full

  qc <- tf |>
    dplyr::group_by(variable, x) |>
    dplyr::summarise(
      q_lo_out = stats::quantile(F_x, probs[1], na.rm = TRUE),
      q_lo_in  = stats::quantile(F_x, probs[2], na.rm = TRUE),
      q_mid    = stats::quantile(F_x, probs[3], na.rm = TRUE),
      q_hi_in  = stats::quantile(F_x, probs[4], na.rm = TRUE),
      q_hi_out = stats::quantile(F_x, probs[5], na.rm = TRUE),
      .groups  = "drop"
    )

  p <- ggplot2::ggplot(qc, ggplot2::aes(x = x)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = q_lo_out, ymax = q_hi_out),
                         fill = "grey80") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = q_lo_in,  ymax = q_hi_in),
                         fill = "grey60") +
    ggplot2::geom_line(ggplot2::aes(y = q_mid), colour = "black") +
    ggplot2::facet_wrap(~ variable, scales = "free") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x        = "Environmental gradient",
      y        = "Turnover",
      title    = "Quantile envelope of turnover",
      subtitle = sprintf("Bands: %s / %s  |  Line: median",
                         paste0(probs[c(1,5)] * 100, "%", collapse = "–"),
                         paste0(probs[c(2,4)] * 100, "%", collapse = "–"))
    )

  plots <- list(quantile = p)
  .handle_output(plots, output, outdir, prefix = "quantile", width, height, dpi)
}
