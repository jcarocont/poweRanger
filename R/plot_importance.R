# =============================================================================
# plot_importance.R — Variable importance plots
# =============================================================================

#' Plot variable importance summaries
#'
#' Produces up to four bar plots from a \code{gf_summary}:
#' total importance, percent of SNPs, and optionally top-N frequency counts.
#'
#' @param gf_summary A \code{gf_summary} object.
#' @param type character vector. Any combination of:
#'   \code{"total"} (total importance by variable),
#'   \code{"percent"} (% of SNPs per variable),
#'   \code{"snp_turnover"} (SNPs ranked by total turnover).
#'   Default: all three.
#' @param output character. \code{"session"} or \code{"png"}.
#' @param outdir character. Output directory.
#' @param width,height,dpi plot dimensions.
#'
#' @return Named list of ggplots.
#' @export
plot_importance <- function(gf_summary,
                            type   = c("total", "percent", "snp_turnover"),
                            output = c("session", "png"),
                            outdir = ".",
                            width  = 8, height = 6, dpi = 300) {

  output <- match.arg(output)
  type   <- match.arg(type, several.ok = TRUE)
  if (!inherits(gf_summary, "gf_summary"))
    stop("`gf_summary` must be a gf_summary object.")

  vi  <- gf_summary$var_importance
  snp <- gf_summary$snp_summary
  plots <- list()

  if ("total" %in% type) {
    plots[["total"]] <- ggplot2::ggplot(
      vi, ggplot2::aes(
        x = stats::reorder(variable, total_importance),
        y = total_importance
      )
    ) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::theme_bw() +
      ggplot2::labs(
        x     = "Environmental variable",
        y     = "Total importance",
        title = "Global variable importance"
      )
  }

  if ("percent" %in% type) {
    plots[["percent"]] <- ggplot2::ggplot(
      vi, ggplot2::aes(
        x = stats::reorder(variable, percent_snps),
        y = percent_snps
      )
    ) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::theme_classic() +
      ggplot2::labs(
        x     = "Variable",
        y     = "% of SNPs",
        title = "Proportion of SNPs associated per variable"
      )
  }

  if ("snp_turnover" %in% type) {
    plots[["snp_turnover"]] <- ggplot2::ggplot(
      snp, ggplot2::aes(
        x = stats::reorder(SNP_ID, total_turnover),
        y = total_turnover
      )
    ) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.y = ggplot2::element_blank()) +
      ggplot2::labs(
        x     = "SNP",
        y     = "Total turnover",
        title = "SNPs ranked by environmental turnover"
      )
  }

  .handle_output(plots, output, outdir, prefix = "importance", width, height, dpi)
}
