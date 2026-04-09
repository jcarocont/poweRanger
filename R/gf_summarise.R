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
