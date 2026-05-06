## =============================================================================
## V-SoRE: empirical-Bayes annotation-ranked fine-mapping
##
## Implements Algorithm 2 of:
##   Chen S. (2026) Annotation-Ranked Bayesian Fine-Mapping of GWAS Signals
##   with Data-Adaptive Reversion to Exchangeability.
##
## Dependencies: susieR (>= 0.12.0).
## =============================================================================

#' Run V-SoRE on summary statistics.
#'
#' Wraps SuSiE-RSS as an inner solver inside an outer empirical-Bayes
#' reweighting loop driven by the agreement between annotation ranks and
#' current PIPs.
#'
#' Algorithm 2 of the manuscript begins with a uniform-weight bypass: when the
#' anchor weights are uniform, no reweighting is possible and V-SoRE returns
#' the SuSiE-RSS solution directly with hat_tau2 = NA. This implementation
#' enforces the bypass.
#'
#' @param z       Numeric vector of marginal z-scores, length p.
#' @param R       p x p ridge-regularized LD matrix.
#' @param omega0  Numeric vector of strictly positive anchor weights, length p.
#' @param n_eff   Effective sample size (use compute_neff() for binary GWAS).
#' @param L       SuSiE-RSS number of single-effect components. Default 10.
#' @param K0      Top-K0 window size used for the Spearman correlation between
#'                anchor ranks and current PIPs. Default 20.
#' @param T_max   Maximum number of outer reweighting iterations. Default 15.
#' @param delta   Convergence tolerance on max change in log-weights. Default 1e-3.
#' @param verbose If TRUE, print per-iteration progress. Default FALSE.
#'
#' @return List with components:
#'   \describe{
#'     \item{fit}{The final SuSiE-RSS fit object.}
#'     \item{pip}{Posterior inclusion probabilities, length p.}
#'     \item{hat_tau2}{Estimated annotation-data agreement parameter, or NA
#'                     under the uniform-weight bypass.}
#'     \item{cs}{List of 95\% credible sets (each a vector of variant indices).}
#'     \item{n_iter}{Number of outer iterations actually run.}
#'     \item{converged}{Logical, TRUE if outer loop converged within T_max.}
#'     \item{bypassed}{Logical, TRUE if the uniform-weight bypass was triggered.}
#'   }
#'
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   p <- 100
#'   R <- build_ld(p, n_blocks = 5, rho = 0.85)
#'   b <- numeric(p); b[c(15, 60)] <- 1.5
#'   z <- as.numeric(R %*% b + t(chol(R)) %*% rnorm(p))
#'   omega <- runif(p, 0.5, 2); omega[c(15, 60)] <- 5
#'   res <- vsore(z, R, omega, n_eff = 50000)
#'   res$pip; res$hat_tau2; res$cs
#' }
#' @export
vsore <- function(z, R, omega0, n_eff,
                  L = 10, K0 = 20, T_max = 15, delta = 1e-3,
                  verbose = FALSE) {

  ## --- Input validation -----------------------------------------------------
  if (!requireNamespace("susieR", quietly = TRUE))
    stop("Package 'susieR' is required. Install with install.packages('susieR').")

  z <- as.numeric(z)
  R <- as.matrix(R)
  omega0 <- as.numeric(omega0)
  p <- length(z)

  if (nrow(R) != p || ncol(R) != p)
    stop("R must be a p x p matrix matching length(z) = ", p, ".")
  if (length(omega0) != p)
    stop("omega0 must have length p = ", p, ".")
  if (any(!is.finite(z)))
    stop("z contains non-finite values.")
  if (any(!is.finite(omega0)) || any(omega0 <= 0))
    stop("omega0 must contain strictly positive finite values.")
  if (!is.numeric(n_eff) || length(n_eff) != 1L || n_eff <= 0)
    stop("n_eff must be a single positive number.")
  if (L < 1L || K0 < 1L || T_max < 1L || delta <= 0)
    stop("L, K0, T_max must be >= 1 and delta must be > 0.")

  K0 <- min(K0, p)

  ## --- Uniform-anchor bypass (Algorithm 2 step 1, Proposition 3) -----------
  if (.is_uniform(omega0)) {
    if (verbose) message("Uniform anchor weights detected: bypassing to SuSiE-RSS.")
    fit <- susieR::susie_rss(z = z, R = R, n = n_eff, L = L)
    return(.format_result(fit, hat_tau2 = NA_real_, n_iter = 0L,
                          converged = TRUE, bypassed = TRUE))
  }

  ## --- Initial SuSiE-RSS pass with anchor priors ---------------------------
  log_om0 <- log(omega0)
  I_K0 <- order(omega0, decreasing = TRUE)[seq_len(K0)]

  fit <- susieR::susie_rss(z = z, R = R, n = n_eff, L = L,
                           prior_weights = omega0 / sum(omega0))
  gamma <- susieR::susie_get_pip(fit)

  log_om_prev <- log_om0
  hat_tau2 <- NA_real_
  converged <- FALSE
  n_iter <- 0L

  ## --- Outer empirical-Bayes reweighting loop ------------------------------
  for (t in seq_len(T_max)) {
    n_iter <- t

    ## Step a: Spearman correlation on the top-K0 window.
    rs <- .spearman_r(omega0[I_K0], gamma[I_K0])

    ## Step b: tau-hat update on [0.05, 1.5].
    hat_tau2 <- max(0.05, 1.5 * (1 - rs))

    ## Step c: log-weight convex combination of anchors and current PIPs.
    gamma_clip <- pmin(pmax(gamma, 1e-6), 1 - 1e-6)
    log_om_t <- (log_om0 + hat_tau2 * log(gamma_clip)) / (1 + hat_tau2)
    omega_t <- exp(log_om_t)

    ## Step d: re-run SuSiE-RSS with the updated weights.
    fit <- susieR::susie_rss(z = z, R = R, n = n_eff, L = L,
                             prior_weights = omega_t / sum(omega_t))
    gamma <- susieR::susie_get_pip(fit)

    ## Step e: convergence check (after a minimum of 3 iterations).
    dw <- max(abs(log_om_t - log_om_prev))
    log_om_prev <- log_om_t

    if (verbose)
      message(sprintf("  iter %2d  rs=%6.3f  hat_tau2=%5.3f  dlog_w=%.2e",
                      t, rs, hat_tau2, dw))

    if (t >= 3L && dw < delta) {
      converged <- TRUE
      break
    }
  }

  .format_result(fit, hat_tau2 = hat_tau2, n_iter = n_iter,
                 converged = converged, bypassed = FALSE)
}


#' V-SoRE with square-root-flattened anchor weights.
#'
#' Used in the IIA/PL-sensitivity analysis (Supplementary Section S13).
#' Despite the historical name, this is NOT the Mallows ranking model;
#' it is V-SoRE with Plackett-Luce weights replaced by sqrt(omega0), which
#' compresses the anchor-weight scale and reduces the influence of
#' independence-of-irrelevant-alternatives in regions of high LD.
#'
#' Forwards all other arguments unchanged to vsore().
#'
#' @inheritParams vsore
#' @return Same structure as vsore().
#' @export
vsore_flattened <- function(z, R, omega0, n_eff,
                            L = 10, K0 = 20, T_max = 15, delta = 1e-3,
                            verbose = FALSE) {
  if (any(omega0 <= 0))
    stop("omega0 must be strictly positive for sqrt flattening.")
  vsore(z = z, R = R, omega0 = sqrt(as.numeric(omega0)),
        n_eff = n_eff, L = L, K0 = K0, T_max = T_max,
        delta = delta, verbose = verbose)
}


## =============================================================================
## Internal helpers (not exported)
## =============================================================================

#' Test whether a numeric vector is constant up to numerical tolerance.
#' @keywords internal
.is_uniform <- function(x, tol = .Machine$double.eps^0.5) {
  if (length(x) <= 1L) return(TRUE)
  diff_range <- max(x) - min(x)
  scale <- max(abs(x), 1)
  diff_range / scale < tol
}

#' Ties-aware Spearman correlation that returns 0 when undefined.
#' Uses base R cor() so that midranks handle ties correctly.
#' @keywords internal
.spearman_r <- function(x, y) {
  if (length(x) < 3L) return(0)
  if (.is_uniform(x) || .is_uniform(y)) return(0)
  r <- suppressWarnings(stats::cor(x, y, method = "spearman"))
  if (is.na(r)) 0 else r
}

#' Format the V-SoRE return value uniformly.
#' @keywords internal
.format_result <- function(fit, hat_tau2, n_iter, converged, bypassed) {
  pip <- susieR::susie_get_pip(fit)
  cs_obj <- susieR::susie_get_cs(fit)
  cs <- cs_obj$cs
  cs <- cs[lengths(cs) > 0]
  list(fit = fit,
       pip = pip,
       hat_tau2 = hat_tau2,
       cs = cs,
       n_iter = n_iter,
       converged = converged,
       bypassed = bypassed)
}
