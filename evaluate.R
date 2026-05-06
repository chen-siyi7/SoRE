## =============================================================================
## Credible-set evaluation utilities.
##
## evaluate_cs() reports coverage and mean credible-set size against a known
## causal index (used in simulations). It returns NA for cs_size when the
## method produces no credible sets at all, which avoids biasing replicate
## means downward by treating empty-CS replicates as size-0.
## =============================================================================

#' Evaluate credible sets against the true causal index.
#'
#' For each known causal variant, coverage is the fraction of causal variants
#' contained in at least one returned credible set. Mean credible-set size is
#' the mean of |CS_l| across the returned credible sets.
#'
#' @param causal_idx Integer vector of true causal variant indices.
#' @param cs_list    List of credible sets, each an integer vector of variant
#'                   indices. Empty CS entries are dropped before evaluation.
#' @return Named numeric vector with components:
#'         \describe{
#'           \item{coverage}{Fraction of causal variants covered.}
#'           \item{cs_size}{Mean CS size, or NA if no CS were returned.}
#'           \item{n_cs}{Number of non-empty credible sets returned.}
#'         }
#' @export
evaluate_cs <- function(causal_idx, cs_list) {
  cs_list <- cs_list[lengths(cs_list) > 0]
  n_cs <- length(cs_list)
  if (n_cs == 0L)
    return(c(coverage = 0.0, cs_size = NA_real_, n_cs = 0))
  cov <- mean(vapply(causal_idx,
                     function(j) any(vapply(cs_list,
                                            function(cs) j %in% cs,
                                            logical(1))),
                     logical(1)))
  sz <- mean(lengths(cs_list))
  c(coverage = as.numeric(cov),
    cs_size  = as.numeric(sz),
    n_cs     = as.integer(n_cs))
}


#' Summarise a V-SoRE result as a per-CS data frame.
#'
#' Useful for the Alzheimer's-disease application table or for any per-locus
#' reporting. Each row corresponds to one returned credible set.
#'
#' @param res            A list returned by vsore().
#' @param variant_names  Optional length-p character vector of variant names.
#' @return data.frame with columns: cs_index, cs_size, lead_pip, lead_variant,
#'         hat_tau2.
#' @export
summarise_vsore <- function(res, variant_names = NULL) {
  cs <- res$cs
  pip <- res$pip
  p <- length(pip)
  nms <- if (!is.null(variant_names)) variant_names else as.character(seq_len(p))
  if (length(nms) != p)
    stop("length(variant_names) must equal length(res$pip).")

  if (length(cs) == 0L)
    return(data.frame(cs_index     = integer(0),
                      cs_size      = integer(0),
                      lead_pip     = numeric(0),
                      lead_variant = character(0),
                      hat_tau2     = numeric(0),
                      stringsAsFactors = FALSE))

  rows <- lapply(seq_along(cs), function(i) {
    idx <- cs[[i]]
    lead_j <- idx[which.max(pip[idx])]
    data.frame(cs_index     = i,
               cs_size      = length(idx),
               lead_pip     = pip[lead_j],
               lead_variant = as.character(nms[lead_j]),
               hat_tau2     = res$hat_tau2,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}


#' Aggregate evaluate_cs() across simulation replicates with Monte Carlo SEs.
#'
#' Useful for reproducing the simulation tables of the manuscript with
#' Monte Carlo standard errors as suggested by the simulation-reporting
#' guidance of Morris, White & Crowther (Stat Med 2019).
#'
#' @param replicate_results A list of evaluate_cs() outputs, length R.
#' @return Named numeric vector with mean and Monte Carlo SE for coverage and
#'         mean CS size, plus the number of valid replicates.
#' @export
aggregate_replicates <- function(replicate_results) {
  if (!length(replicate_results))
    return(c(coverage_mean = NA_real_, coverage_mcse = NA_real_,
             cs_size_mean = NA_real_, cs_size_mcse = NA_real_,
             n_valid = 0))
  cov_v <- vapply(replicate_results, function(r) as.numeric(r["coverage"]), numeric(1))
  sz_v  <- vapply(replicate_results, function(r) as.numeric(r["cs_size"]),  numeric(1))
  ## Drop replicates where cs_size is NA (no CS returned) when computing the
  ## CS-size mean and its SE; coverage is well-defined even then (it is 0).
  sz_valid <- !is.na(sz_v)
  R_cov  <- length(cov_v)
  R_sz   <- sum(sz_valid)
  cov_mean  <- mean(cov_v)
  cov_mcse  <- sqrt(cov_mean * (1 - cov_mean) / R_cov)
  if (R_sz >= 2L) {
    sz_mean <- mean(sz_v[sz_valid])
    sz_mcse <- stats::sd(sz_v[sz_valid]) / sqrt(R_sz)
  } else {
    sz_mean <- if (R_sz == 1L) sz_v[sz_valid] else NA_real_
    sz_mcse <- NA_real_
  }
  c(coverage_mean = cov_mean,
    coverage_mcse = cov_mcse,
    cs_size_mean  = sz_mean,
    cs_size_mcse  = sz_mcse,
    n_valid       = as.integer(R_cov))
}
