# =============================================================================
# classes.R — S3 constructors for powerRanger
# =============================================================================

#' Constructor for gf_result
#'
#' @param turnover_full long tibble: SNP_ID x variable x x x F_x x importance x r_squared
#' @param vars character vector of variables used
#' @param n_models number of models processed
#' @keywords internal
new_gf_result <- function(turnover_full, vars, n_models) {
  structure(
    list(
      turnover_full = turnover_full,
      vars          = vars,
      n_models      = n_models
    ),
    class = "gf_result"
  )
}

#' Constructor for gf_summary
#'
#' @keywords internal
new_gf_summary <- function(snp_summary, var_importance, var_snp_matrix,
                            turnover_concentration, curve_cor, curve_xy,
                            snps_keep, zone_counts, mean_curve) {
  structure(
    list(
      snp_summary            = snp_summary,
      var_importance         = var_importance,
      var_snp_matrix         = var_snp_matrix,
      turnover_concentration = turnover_concentration,
      curve_cor              = curve_cor,
      curve_xy               = curve_xy,
      snps_keep              = snps_keep,
      zone_counts            = zone_counts,
      mean_curve             = mean_curve
    ),
    class = "gf_summary"
  )
}

#' @export
print.gf_result <- function(x, ...) {
  cat("── gf_result ──────────────────────────────\n")
  cat(sprintf("  Models    : %d\n", x$n_models))
  cat(sprintf("  Variables : %s\n", paste(x$vars, collapse = ", ")))
  cat(sprintf("  SNPs      : %d\n", dplyr::n_distinct(x$turnover_full$SNP_ID)))
  cat(sprintf("  Rows      : %d\n", nrow(x$turnover_full)))
  invisible(x)
}

#' @export
print.gf_summary <- function(x, ...) {
  cat("── gf_summary ─────────────────────────────\n")
  cat(sprintf("  SNPs total   : %d\n", nrow(x$snp_summary)))
  cat(sprintf("  SNPs passing : %d\n", nrow(x$snps_keep)))
  cat("  Zone counts:\n")
  print(x$zone_counts)
  invisible(x)
}
# =============================================================================
# filter_models.R — Filter and prepare a ranger model list for GF analysis
# =============================================================================

#' Filter and prepare a list of ranger models for gradient forest analysis
#'
#' Validates, filters by out-of-bag R², optionally caches \code{treeInfo},
#' and returns a cleaned named list with a summary attribute.
#'
#' This function should be run before \code{\link{turnover_curves}}.
#' It makes the filtering step explicit and auditable.
#'
#' @param models Named list. Each element is either a \code{ranger} object
#'   or a list with a \code{$model} ranger element (as produced by
#'   \code{random_forest_model.R}).
#' @param min_r2 numeric. Minimum OOB R² to retain a model. Models with
#'   \code{r.squared <= min_r2} are dropped. Default \code{0} (drops only
#'   non-positive R²). Use e.g. \code{0.01} to be stricter.
#' @param cache_treeinfo logical. If \code{TRUE} (default), pre-computes and
#'   stores \code{ranger::treeInfo()} for each model inside the list entry
#'   under \code{.__treeInfo_cache__}. This makes \code{turnover_curves()}
#'   faster at the cost of memory.
#' @param verbose logical. Print per-model diagnostics. Default \code{FALSE}.
#'
#' @return A named list of the same structure as \code{models}, containing
#'   only models that passed the filter. The list has two attributes:
#'   \describe{
#'     \item{\code{filter_summary}}{a \code{data.frame} with columns
#'       \code{SNP_ID}, \code{r_squared}, \code{passed}, \code{reason}.}
#'     \item{\code{n_passed}}{integer: number of retained models.}
#'   }
#'
#' @examples
#' \dontrun{
#' models_clean <- filter_models(models, min_r2 = 0.01, cache_treeinfo = TRUE)
#' attr(models_clean, "filter_summary")
#'
#' gf <- turnover_curves(models_clean, env, vars = c("bio1", "bio12"))
#' }
#'
#' @export
filter_models <- function(models,
                          min_r2          = 0,
                          cache_treeinfo  = TRUE,
                          verbose         = FALSE) {

  # ── basic validation ─────────────────────────────────────────────────────────
  if (!is.list(models) || is.null(names(models)))
    stop("`models` must be a named list.")
  if (!is.numeric(min_r2) || length(min_r2) != 1 || min_r2 < 0 || min_r2 >= 1)
    stop("`min_r2` must be a single numeric in [0, 1).")

  snp_ids <- names(models)
  n_total <- length(snp_ids)
  message(sprintf("[filter_models] Input: %d models", n_total))

  # ── per-model evaluation ─────────────────────────────────────────────────────
  summary_rows <- vector("list", n_total)
  keep         <- logical(n_total)

  for (i in seq_along(snp_ids)) {
    sid   <- snp_ids[i]
    entry <- models[[i]]
    m     <- .get_ranger(entry)

    # not a ranger object
    if (is.null(m)) {
      reason <- "not a ranger object"
      summary_rows[[i]] <- data.frame(
        SNP_ID    = sid,
        r_squared = NA_real_,
        passed    = FALSE,
        reason    = reason,
        stringsAsFactors = FALSE
      )
      if (verbose) message(sprintf("  SKIP [%s]: %s", sid, reason))
      next
    }

    r2 <- m$r.squared

    # missing or non-finite r²
    if (is.null(r2) || !is.finite(r2)) {
      reason <- "r.squared is NULL or non-finite"
      summary_rows[[i]] <- data.frame(
        SNP_ID    = sid,
        r_squared = NA_real_,
        passed    = FALSE,
        reason    = reason,
        stringsAsFactors = FALSE
      )
      if (verbose) message(sprintf("  SKIP [%s]: %s", sid, reason))
      next
    }

    # below threshold
    if (r2 <= min_r2) {
      reason <- sprintf("r.squared = %.4f <= min_r2 = %.4f", r2, min_r2)
      summary_rows[[i]] <- data.frame(
        SNP_ID    = sid,
        r_squared = r2,
        passed    = FALSE,
        reason    = reason,
        stringsAsFactors = FALSE
      )
      if (verbose) message(sprintf("  SKIP [%s]: %s", sid, reason))
      next
    }

    # passed — optionally cache treeInfo
    if (cache_treeinfo) {
      ti <- tryCatch(
        ranger::treeInfo(m),
        error = function(e) {
          warning(sprintf("treeInfo failed for '%s': %s", sid, e$message))
          NULL
        }
      )
      if (!is.null(ti)) {
        if (inherits(entry, "ranger")) {
          models[[i]]$.__treeInfo_cache__ <- ti
        } else {
          models[[i]]$model$.__treeInfo_cache__ <- ti
        }
      }
    }

    keep[i] <- TRUE
    summary_rows[[i]] <- data.frame(
      SNP_ID    = sid,
      r_squared = r2,
      passed    = TRUE,
      reason    = "ok",
      stringsAsFactors = FALSE
    )
  }

  # ── assemble output ──────────────────────────────────────────────────────────
  filter_summary <- do.call(rbind, summary_rows)
  n_passed  <- sum(keep)
  n_skipped <- n_total - n_passed

  message(sprintf("[filter_models] Passed: %d / %d  |  Skipped: %d",
                  n_passed, n_total, n_skipped))

  if (n_passed == 0)
    stop("No models passed the filter. Lower `min_r2` or check your model list.")

  # r² distribution of retained models
  r2_vals <- filter_summary$r_squared[filter_summary$passed]
  message(sprintf("[filter_models] R² range (passed): %.4f – %.4f  |  median: %.4f",
                  min(r2_vals), max(r2_vals), stats::median(r2_vals)))

  out <- models[keep]
  attr(out, "filter_summary") <- filter_summary
  attr(out, "n_passed")       <- n_passed

  out
}
# =============================================================================
# gf_summarise.R — Diagnostic summary tables for a gf_result
# =============================================================================

#' Compute diagnostic summary tables from a gf_result
#'
#' Produces per-SNP and per-variable summaries, turnover concentration
#' metrics, curve stability (correlation to mean), and zone classification
#' (good / borderline / overfit).
#'
#' @param gf_result An object of class \code{gf_result} produced by
#'   \code{\link{turnover_curves}}.
#' @param thr_frac_good numeric. Concentration threshold below which a SNP
#'   is considered well-distributed. Default 0.2.
#' @param thr_frac_bad numeric. Concentration threshold above which a SNP
#'   is considered potentially overfit. Default 0.4.
#' @param thr_cor_good numeric. Correlation-to-mean-curve threshold above
#'   which a SNP is considered stable. Default 0.6.
#' @param thr_cor_bad numeric. Correlation threshold below which a SNP
#'   is considered unstable. Default 0.4.
#' @param filter_frac_peak numeric. Hard filter: SNPs with frac_peak above
#'   this are excluded from \code{snps_keep}. Default 0.3.
#' @param filter_cor numeric. Hard filter: SNPs with cor_to_mean below this
#'   are excluded from \code{snps_keep}. Default 0.5.
#'
#' @return An object of class \code{gf_summary}.
#'
#' @export
gf_summarise <- function(gf_result,
                         thr_frac_good    = 0.2,
                         thr_frac_bad     = 0.4,
                         thr_cor_good     = 0.6,
                         thr_cor_bad      = 0.4,
                         filter_frac_peak = 0.3,
                         filter_cor       = 0.5) {

  if (!inherits(gf_result, "gf_result"))
    stop("`gf_result` must be a gf_result object.")

  tf <- gf_result$turnover_full

  # ── snp_summary ─────────────────────────────────────────────────────────────
  snp_summary <- tf |>
    dplyr::group_by(SNP_ID) |>
    dplyr::summarise(
      total_turnover = sum(F_x,       na.rm = TRUE),
      mean_r2        = mean(r_squared, na.rm = TRUE),
      .groups        = "drop"
    )

  # ── variable importance ──────────────────────────────────────────────────────
  var_importance <- tf |>
    dplyr::distinct(SNP_ID, variable, importance) |>
    dplyr::group_by(variable) |>
    dplyr::summarise(
      total_importance = sum(importance, na.rm = TRUE),
      n_snps           = dplyr::n(),
      .groups          = "drop"
    ) |>
    dplyr::mutate(
      percent_snps = 100 * n_snps / dplyr::n_distinct(tf$SNP_ID)
    ) |>
    dplyr::arrange(dplyr::desc(total_importance))

  # ── var_snp_matrix ───────────────────────────────────────────────────────────
  var_snp_matrix <- tf |>
    dplyr::group_by(SNP_ID, variable) |>
    dplyr::summarise(total_turnover = sum(F_x, na.rm = TRUE), .groups = "drop")

  # ── turnover_concentration ───────────────────────────────────────────────────
  turnover_concentration <- tf |>
    dplyr::group_by(SNP_ID, variable) |>
    dplyr::summarise(
      total_fx  = sum(F_x,  na.rm = TRUE),
      peak_fx   = max(F_x,  na.rm = TRUE),
      frac_peak = peak_fx / total_fx,
      .groups   = "drop"
    )

  # ── mean_curve & curve_cor ───────────────────────────────────────────────────
  mean_curve <- tf |>
    dplyr::group_by(variable, x) |>
    dplyr::summarise(mean_fx = mean(F_x, na.rm = TRUE), .groups = "drop")

  curve_cor <- tf |>
    dplyr::left_join(mean_curve, by = c("variable", "x")) |>
    dplyr::group_by(variable, SNP_ID) |>
    dplyr::summarise(
      cor_to_mean = stats::cor(F_x, mean_fx, use = "complete.obs"),
      .groups     = "drop"
    )

  # ── zone classification ──────────────────────────────────────────────────────
  curve_xy <- curve_cor |>
    dplyr::left_join(turnover_concentration, by = c("SNP_ID", "variable")) |>
    dplyr::group_by(SNP_ID) |>
    dplyr::summarise(
      frac_peak   = max(frac_peak,    na.rm = TRUE),
      cor_to_mean = mean(cor_to_mean, na.rm = TRUE),
      .groups     = "drop"
    ) |>
    dplyr::mutate(
      zone = dplyr::case_when(
        frac_peak < thr_frac_good & cor_to_mean > thr_cor_good ~ "good",
        frac_peak > thr_frac_bad  & cor_to_mean < thr_cor_bad  ~ "overfit",
        TRUE ~ "borderline"
      )
    )

  zone_counts <- curve_xy |>
    dplyr::distinct(SNP_ID, zone) |>
    dplyr::count(zone) |>
    dplyr::mutate(pct = 100 * n / sum(n))

  # ── snps_keep ────────────────────────────────────────────────────────────────
  snps_keep <- curve_cor |>
    dplyr::left_join(turnover_concentration, by = c("SNP_ID", "variable")) |>
    dplyr::filter(frac_peak < filter_frac_peak, cor_to_mean > filter_cor) |>
    dplyr::distinct(SNP_ID)

  message(sprintf("[gf_summarise] SNPs: %d total | %d passing filter | zones: %s",
                  nrow(snp_summary),
                  nrow(snps_keep),
                  paste(sprintf("%s=%d", zone_counts$zone, zone_counts$n),
                        collapse = ", ")))

  new_gf_summary(
    snp_summary            = snp_summary,
    var_importance         = var_importance,
    var_snp_matrix         = var_snp_matrix,
    turnover_concentration = turnover_concentration,
    curve_cor              = curve_cor,
    curve_xy               = curve_xy,
    snps_keep              = snps_keep,
    zone_counts            = zone_counts,
    mean_curve             = mean_curve
  )
}
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
# =============================================================================
# turnover_curves.R — Main entry point for powerRanger
# =============================================================================

#' Compute gradient forest turnover curves from a list of ranger models
#'
#' For each model (one per locus/SNP) and each environmental variable,
#' computes the normalised cumulative turnover curve F(x) weighted by
#' the per-predictor R² contribution.
#'
#' @param models Named list of ranger models (or named list of lists with a
#'   `$model` ranger element, as produced by `random_forest_model.R`).
#'   Names are used as SNP/locus identifiers.
#' @param env data.frame or matrix of environmental variables.
#'   Rows = sites/individuals, columns = predictors.
#'   Column names must match predictor names used during training.
#' @param vars character vector of predictor names to compute turnover for.
#'   Defaults to all predictors found in the first model.
#' @param num_bins integer. Number of grid points along each environmental
#'   gradient. Default 101.
#' @param min_r2 numeric. Models with overall r.squared <= this value are
#'   skipped. Default 0.
#' @param .progress logical. Show a progress bar (requires cli). Default TRUE.
#'
#' @return An object of class \code{gf_result}.
#'
#' @export
turnover_curves <- function(models, env, vars = NULL,
                            num_bins = 101L, min_r2 = 0,
                            .progress = TRUE) {

  # ── input checks ────────────────────────────────────────────────────────────
  ok <- .validate_models(models)
  models <- models[ok]

  if (!is.data.frame(env) && !is.matrix(env))
    stop("`env` must be a data.frame or matrix.")
  env <- as.data.frame(env)

  # ── filter by r² ────────────────────────────────────────────────────────────
  r2_vec  <- vapply(models, function(e) {
    m <- .get_ranger(e); if (is.null(m)) NA_real_ else m$r.squared
  }, numeric(1))
  models  <- models[!is.na(r2_vec) & r2_vec > min_r2]

  if (length(models) == 0)
    stop("No models passed the r² filter.")

  message(sprintf("[turnover_curves] %d models after r² > %g filter",
                  length(models), min_r2))

  # ── default vars ────────────────────────────────────────────────────────────
  first_model <- .get_ranger(models[[1]])
  all_preds   <- names(first_model$variable.importance)

  if (is.null(vars)) {
    vars <- all_preds
  } else {
    missing_v <- setdiff(vars, colnames(env))
    if (length(missing_v) > 0)
      stop(sprintf("Variables not found in `env`: %s",
                   paste(missing_v, collapse = ", ")))
  }

  # ── pre-cache treeInfo ───────────────────────────────────────────────────────
  message("[turnover_curves] Caching treeInfo for all models...")
  ti_cache <- lapply(models, function(e) {
    m <- .get_ranger(e)
    ranger::treeInfo(m)
  })

  # ── main loop ───────────────────────────────────────────────────────────────
  snp_ids <- names(models)
  n       <- length(snp_ids)

  rows <- vector("list", n * length(vars))
  idx  <- 1L

  for (i in seq_along(snp_ids)) {
    sid   <- snp_ids[i]
    model <- .get_ranger(models[[i]])
    ti    <- ti_cache[[i]]
    vi    <- model$variable.importance
    r2    <- model$r.squared

    if (.progress && i %% max(1L, n %/% 20L) == 0)
      message(sprintf("  [%d / %d] %s", i, n, sid))

    for (v in vars) {
      if (!(v %in% colnames(env))) next
      x_vec <- env[[v]]

      curve <- .turnover_curve_one(model, ti, x_vec, v, num_bins)
      if (is.null(curve)) next

      imp <- if (v %in% names(vi)) vi[[v]] else NA_real_

      rows[[idx]] <- tibble::tibble(
        SNP_ID    = sid,
        variable  = v,
        x         = curve$x,
        F_x       = curve$F_x,
        importance = imp,
        r_squared  = r2
      )
      idx <- idx + 1L
    }
  }

  turnover_full <- dplyr::bind_rows(rows[seq_len(idx - 1L)])

  message(sprintf("[turnover_curves] Done. %d rows | %d SNPs | %d variables",
                  nrow(turnover_full),
                  dplyr::n_distinct(turnover_full$SNP_ID),
                  dplyr::n_distinct(turnover_full$variable)))

  new_gf_result(turnover_full, vars = vars, n_models = length(models))
}
# =============================================================================
# utils.R — Internal helpers for powerRanger
# =============================================================================

#' Extract R² contribution of a single predictor from a ranger model
#' @keywords internal
.extract_r2_contribution <- function(model, predictor_name) {
  r2 <- model$r.squared
  if (is.null(r2) || r2 <= 0) return(0)
  vi <- model$variable.importance
  imp <- vi[predictor_name]
  if (is.null(imp) || is.na(imp) || imp <= 0) return(0)
  r2 * (imp / sum(vi, na.rm = TRUE))
}

#' Extract split density for one predictor using a pre-cached treeInfo tibble
#'
#' @param treeinfo_cache tibble returned by ranger::treeInfo(), pre-computed
#' @param x_vec numeric vector of env values for this predictor
#' @param predictor_name character
#' @param num_bins integer
#' @keywords internal
.extract_split_density_cached <- function(treeinfo_cache, x_vec,
                                          predictor_name, num_bins) {
  splitvals <- treeinfo_cache[treeinfo_cache$splitvarName == predictor_name,
                              "splitval", drop = TRUE]
  if (length(splitvals) < 2) return(NULL)

  rng    <- range(x_vec, na.rm = TRUE)
  grid_x <- seq(rng[1], rng[2], length.out = num_bins)

  dens_splits <- stats::density(splitvals, from = rng[1],
                                to = rng[2], n = num_bins)$y
  dens_data   <- stats::density(x_vec, from = rng[1],
                                to = rng[2], n = num_bins)$y

  list(grid_x = grid_x, dens_splits = dens_splits, dens_data = dens_data)
}

#' Compute normalised cumulative turnover curve for one SNP × predictor
#' @keywords internal
.turnover_curve_one <- function(model, treeinfo_cache, x_vec,
                                predictor_name, num_bins) {
  r2c  <- .extract_r2_contribution(model, predictor_name)
  if (r2c == 0) return(NULL)

  dens <- .extract_split_density_cached(treeinfo_cache, x_vec,
                                        predictor_name, num_bins)
  if (is.null(dens)) return(NULL)

  epsilon  <- 1e-6
  raw_rate <- dens$dens_splits / (dens$dens_data + epsilon)
  cumul    <- cumsum(raw_rate)
  mx       <- max(cumul, na.rm = TRUE)
  F_x      <- if (mx == 0) rep(0, num_bins) else cumul * (r2c / mx)

  tibble::tibble(x = dens$grid_x, F_x = F_x)
}

#' Validate that models is a named list of ranger objects (or wrapped)
#' @keywords internal
.validate_models <- function(models) {
  if (!is.list(models) || is.null(names(models)))
    stop("`models` must be a named list.")

  is_ranger <- function(entry) {
    if (inherits(entry, "ranger")) return(TRUE)
    if (is.list(entry) && inherits(entry$model, "ranger")) return(TRUE)
    FALSE
  }

  ok <- vapply(models, is_ranger, logical(1))
  if (!any(ok)) stop("No ranger models found in `models`.")
  if (!all(ok)) warning(sprintf("%d entries in `models` are not ranger objects and will be skipped.",
                                sum(!ok)))
  invisible(ok)
}

#' Safely extract the ranger model from a list entry (or return as-is)
#' @keywords internal
.get_ranger <- function(entry) {
  if (inherits(entry, "ranger")) return(entry)
  if (is.list(entry) && inherits(entry$model, "ranger")) return(entry$model)
  NULL
}
