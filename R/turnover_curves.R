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
