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
