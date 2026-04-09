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
