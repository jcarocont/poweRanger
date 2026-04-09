# =============================================================================
# plot_heatmap.R — SNP × variable turnover heatmap
# =============================================================================

#' Plot a heatmap of total turnover per SNP × variable
#'
#' @param gf_result A \code{gf_result} object.
#' @param gf_summary Optional \code{gf_summary} object. If provided, SNPs are
#'   ordered by total_turnover from \code{snp_summary}.
#' @param output character. \code{"session"} or \code{"png"}.
#' @param outdir character. Output directory.
#' @param width,height,dpi plot dimensions.
#'
#' @return Named list with one ggplot (\code{"heatmap"}).
#' @export
plot_heatmap <- function(gf_result,
                         gf_summary = NULL,
                         output     = c("session", "png"),
                         outdir     = ".",
                         width      = 8, height = 6, dpi = 300) {

  output <- match.arg(output)
  if (!inherits(gf_result, "gf_result"))
    stop("`gf_result` must be a gf_result object.")

  tf <- gf_result$turnover_full

  mat <- tf |>
    dplyr::group_by(SNP_ID, variable) |>
    dplyr::summarise(total_turnover = sum(F_x, na.rm = TRUE), .groups = "drop")

  if (!is.null(gf_summary) && inherits(gf_summary, "gf_summary")) {
    snp_order <- gf_summary$snp_summary |>
      dplyr::arrange(dplyr::desc(total_turnover)) |>
      dplyr::pull(SNP_ID)
    mat$SNP_ID <- factor(mat$SNP_ID, levels = snp_order)
  }

  p <- ggplot2::ggplot(mat,
    ggplot2::aes(x = variable, y = SNP_ID, fill = total_turnover)
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(option = "viridis", name = "Turnover") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)
    ) +
    ggplot2::labs(
      x     = "Variable",
      y     = "SNP",
      title = "Turnover — SNP × Variable"
    )

  plots <- list(heatmap = p)
  .handle_output(plots, output, outdir, prefix = "heatmap", width, height, dpi)
}
