# =============================================================================
# plot_diagnostic.R — Overfitting diagnostic plots
# =============================================================================

#' Plot overfitting diagnostics for gradient forest models
#'
#' Produces diagnostic plots based on turnover concentration (frac_peak)
#' and curve stability (correlation to mean curve).
#'
#' @param gf_summary A \code{gf_summary} object.
#' @param type character vector. Any combination of:
#'   \code{"scatter"} (frac_peak vs cor_to_mean with zone colours),
#'   \code{"hist_concentration"} (histogram of frac_peak),
#'   \code{"boxplot_concentration"} (frac_peak boxplot by variable),
#'   \code{"boxplot_stability"} (cor_to_mean boxplot by variable),
#'   \code{"dominance"} (dominant SNP fraction per variable),
#'   \code{"r2_vs_turnover"} (R² vs total turnover scatter).
#'   Default: all six.
#' @param thr_frac_good,thr_frac_bad,thr_cor_good,thr_cor_bad numeric.
#'   Zone thresholds used for the scatter plot. Should match those used in
#'   \code{\link{gf_summarise}}.
#' @param output character. \code{"session"} or \code{"png"}.
#' @param outdir character. Output directory.
#' @param width,height,dpi plot dimensions.
#'
#' @return Named list of ggplots.
#' @export
plot_diagnostic <- function(gf_summary,
                            type = c("scatter", "hist_concentration",
                                     "boxplot_concentration",
                                     "boxplot_stability",
                                     "dominance", "r2_vs_turnover"),
                            thr_frac_good = 0.2,
                            thr_frac_bad  = 0.4,
                            thr_cor_good  = 0.6,
                            thr_cor_bad   = 0.4,
                            output = c("session", "png"),
                            outdir = ".",
                            width  = 8, height = 6, dpi = 300) {

  output <- match.arg(output)
  type   <- match.arg(type, several.ok = TRUE)
  if (!inherits(gf_summary, "gf_summary"))
    stop("`gf_summary` must be a gf_summary object.")

  tc   <- gf_summary$turnover_concentration
  cc   <- gf_summary$curve_cor
  cxy  <- gf_summary$curve_xy
  snp  <- gf_summary$snp_summary
  vi   <- gf_summary$var_snp_matrix
  plots <- list()

  # ── scatter: frac_peak vs cor_to_mean ────────────────────────────────────────
  if ("scatter" %in% type) {
    plots[["scatter"]] <- ggplot2::ggplot(
      cxy, ggplot2::aes(x = frac_peak, y = cor_to_mean)
    ) +
      ggplot2::geom_rect(
        xmin = 0, xmax = thr_frac_good, ymin = thr_cor_good, ymax = 1,
        fill = "darkseagreen2", alpha = 0.3, inherit.aes = FALSE
      ) +
      ggplot2::geom_rect(
        xmin = thr_frac_good, xmax = thr_frac_bad,
        ymin = thr_cor_bad,   ymax = thr_cor_good,
        fill = "khaki1", alpha = 0.3, inherit.aes = FALSE
      ) +
      ggplot2::geom_rect(
        xmin = thr_frac_bad, xmax = 1, ymin = 0, ymax = thr_cor_bad,
        fill = "indianred1", alpha = 0.3, inherit.aes = FALSE
      ) +
      ggplot2::geom_point(ggplot2::aes(colour = zone), alpha = 0.8, size = 0.8) +
      ggplot2::scale_colour_manual(
        values = c(good = "darkgreen", borderline = "goldenrod3",
                   overfit = "firebrick")
      ) +
      ggplot2::geom_vline(
        xintercept = c(thr_frac_good, thr_frac_bad), linetype = "dashed"
      ) +
      ggplot2::geom_hline(
        yintercept = c(thr_cor_good, thr_cor_bad), linetype = "dashed"
      ) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        x        = "Turnover concentration (frac_peak)",
        y        = "Correlation to mean curve",
        title    = "Overfitting diagnostic — Gradient Forest",
        subtitle = "Zones indicate per-SNP fit stability",
        colour   = "Zone"
      )
  }

  # ── histogram: frac_peak ─────────────────────────────────────────────────────
  if ("hist_concentration" %in% type) {
    plots[["hist_concentration"]] <- ggplot2::ggplot(
      tc, ggplot2::aes(x = frac_peak)
    ) +
      ggplot2::geom_histogram(bins = 40, fill = "grey30", colour = "white") +
      ggplot2::theme_bw() +
      ggplot2::labs(
        x     = "Fraction of turnover in peak bin",
        y     = "Number of SNPs",
        title = "Turnover concentration (overfitting diagnostic)"
      )
  }

  # ── boxplot: frac_peak by variable ───────────────────────────────────────────
  if ("boxplot_concentration" %in% type) {
    plots[["boxplot_concentration"]] <- ggplot2::ggplot(
      tc, ggplot2::aes(x = variable, y = frac_peak)
    ) +
      ggplot2::geom_boxplot(outlier.colour = "firebrick") +
      ggplot2::coord_flip() +
      ggplot2::theme_bw() +
      ggplot2::labs(
        x     = "Environmental variable",
        y     = "Frac. of turnover at peak",
        title = "Turnover concentration by variable"
      )
  }

  # ── boxplot: cor_to_mean by variable ─────────────────────────────────────────
  if ("boxplot_stability" %in% type) {
    plots[["boxplot_stability"]] <- ggplot2::ggplot(
      cc, ggplot2::aes(x = variable, y = cor_to_mean)
    ) +
      ggplot2::geom_boxplot(outlier.colour = "firebrick") +
      ggplot2::coord_flip() +
      ggplot2::theme_bw() +
      ggplot2::labs(
        x     = "Environmental variable",
        y     = "Correlation to mean curve",
        title = "Curve stability (GF diagnostic)"
      )
  }

  # ── dominance: fraction explained by dominant SNP per variable ───────────────
  if ("dominance" %in% type) {
    dom_df <- vi |>
      dplyr::group_by(variable) |>
      dplyr::mutate(frac = total_turnover / sum(total_turnover)) |>
      dplyr::summarise(top_frac = max(frac), .groups = "drop")

    plots[["dominance"]] <- ggplot2::ggplot(
      dom_df, ggplot2::aes(
        x = stats::reorder(variable, top_frac),
        y = top_frac
      )
    ) +
      ggplot2::geom_col(fill = "firebrick") +
      ggplot2::coord_flip() +
      ggplot2::theme_bw() +
      ggplot2::labs(
        x     = "Environmental variable",
        y     = "Fraction explained by dominant SNP",
        title = "SNP dominance per variable"
      )
  }

  # ── r² vs total turnover ─────────────────────────────────────────────────────
  if ("r2_vs_turnover" %in% type) {
    plots[["r2_vs_turnover"]] <- ggplot2::ggplot(
      snp, ggplot2::aes(x = mean_r2, y = total_turnover)
    ) +
      ggplot2::geom_point(alpha = 0.4) +
      ggplot2::geom_smooth(method = "lm", se = FALSE, colour = "#0072B2") +
      ggplot2::theme_bw() +
      ggplot2::labs(
        x     = "Mean R²",
        y     = "Total turnover",
        title = "R² vs total turnover per SNP"
      )
  }

  .handle_output(plots, output, outdir, prefix = "diagnostic", width, height, dpi)
}
