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
