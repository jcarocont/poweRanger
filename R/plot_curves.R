# =============================================================================
# plot_curves.R — Per-variable turnover curves (individual + composite)
# =============================================================================

#' Plot turnover curves per environmental variable
#'
#' For each variable, draws all individual SNP curves (grey) plus the
#' composite (sum) curve in blue.
#'
#' @param gf_result A \code{gf_result} object.
#' @param vars character vector of variables to plot. Default: all.
#' @param output character. \code{"session"} returns a named list of ggplots;
#'   \code{"png"} saves files and returns paths invisibly.
#' @param outdir character. Output directory (used when \code{output="png"}).
#' @param width,height,dpi plot dimensions.
#'
#' @return Named list of ggplot objects (invisibly if \code{output="png"}).
#' @export
plot_curves <- function(gf_result,
                        vars   = NULL,
                        output = c("session", "png"),
                        outdir = ".",
                        width  = 8, height = 6, dpi = 300) {

  output <- match.arg(output)
  if (!inherits(gf_result, "gf_result"))
    stop("`gf_result` must be a gf_result object.")

  tf   <- gf_result$turnover_full
  vars <- vars %||% unique(tf$variable)

  plots <- lapply(vars, function(v) {
    df <- tf[tf$variable == v & !is.na(tf$x), ]
    if (nrow(df) == 0) {
      warning(sprintf("No data for variable '%s', skipping.", v))
      return(NULL)
    }
    composite <- dplyr::summarise(
      dplyr::group_by(df, x),
      F_x_sum = sum(F_x, na.rm = TRUE),
      .groups = "drop"
    )
    ggplot2::ggplot() +
      ggplot2::geom_line(
        data  = df,
        ggplot2::aes(x = x, y = F_x, group = SNP_ID),
        alpha = 0.12, colour = "grey50"
      ) +
      ggplot2::geom_line(
        data      = composite,
        ggplot2::aes(x = x, y = F_x_sum),
        linewidth = 1.2, colour = "#0072B2"
      ) +
      ggplot2::labs(
        title    = paste("Turnover curves —", v),
        subtitle = sprintf("N SNPs = %d", dplyr::n_distinct(df$SNP_ID)),
        x        = v,
        y        = expression(Sigma * R^2)
      ) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5),
        plot.subtitle    = ggplot2::element_text(hjust = 0.5),
        panel.grid.minor = ggplot2::element_blank()
      )
  })
  names(plots) <- vars
  plots <- Filter(Negate(is.null), plots)

  .handle_output(plots, output, outdir, prefix = "curves", width, height, dpi)
}
