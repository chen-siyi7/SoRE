# =============================================================================
# vsore.R — V-SoRE outer reweighting loop around SuSiE-RSS
# =============================================================================

#' V-SoRE: variational empirical-Bayes fine-mapping with annotation ranking
#'
#' Runs the V-SoRE outer reweighting loop around SuSiE-RSS. The outer loop
#' alternates between estimating the annotation-data agreement parameter
#' tau^2 from the Spearman rank correlation between the anchor weights and
#' the current posterior inclusion probabilities (PIPs), and re-running
#' SuSiE-RSS with updated weights.
#'
#' @param z Numeric vector of marginal z-scores, length p.
#' @param R p-by-p LD matrix (ridge-regularized reference panel estimate).
#' @param anchor_weights Numeric vector of length p; external annotation
#'   anchor weights tilde-omega_j. Set all entries to 1 to recover SuSiE-RSS.
#' @param n Effective sample size for the SuSiE-RSS call.
#' @param L Maximum number of single-effect components (default 10).
#' @param K0 Top-K0 variants by anchor weight used for the Spearman
#'   diagnostic (default 20).
#' @param T_max Maximum number of outer iterations (default 15).
#' @param tol Outer-loop convergence tolerance on log-weight change
#'   (default 1e-3).
#' @param eps PIP clipping value to avoid log(0) (default 1e-6).
#' @param n_inner Inner SuSiE-RSS iterations per outer pass (default 100).
#' @param verbose Logical; print progress at each outer iteration.
#'
#' @return A list with elements
#'   \describe{
#'     \item{\code{pip}}{p-vector of final posterior inclusion probabilities}
#'     \item{\code{tau2}}{Final annotation-data agreement parameter}
#'     \item{\code{weights}}{p-vector of final weights tilde-omega^{(t)}}
#'     \item{\code{n_iter}}{Number of outer iterations completed}
#'     \item{\code{converged}}{Logical, whether tol was reached}
#'     \item{\code{susie_fit}}{Final susieR object from the inner loop}
#'   }
#'
#' @details
#' If all entries of \code{anchor_weights} equal 1, the function returns
#' the SuSiE-RSS solution exactly with \code{tau2 = NA}.
#'
#' tau^2 is mapped from the Spearman correlation rho_s by
#' tau^2 = max(0.05, 1.5 * (1 - rho_s)).
#' The weight update is the convex combination on the log scale
#' log(tilde-omega^(t)) = (log(tilde-omega^(0)) + tau^2 * log(gamma^(t-1))) /
#'                       (1 + tau^2).
#'
#' Posterior contraction guarantees and the relationship to SuSiE-RSS are
#' described in Sections 3 and 4 of the manuscript.
#'
#' @references
#' Chen, S. (2026). Bayesian Fine-Mapping with an Annotation Ranking Prior
#' Using GWAS Summary Statistics. Biometrics.
#'
#' Zou, Y., Carbonetto, P., Wang, G., Stephens, M. (2022). Fine-mapping from
#' summary data with the Sum of Single Effects model. PLOS Genetics.
#'
#' @examples
#' \dontrun{
#'   # Simulated example
#'   set.seed(1)
#'   p <- 200; n <- 50000; s0 <- 2
#'   R <- 0.9 ^ abs(outer(1:p, 1:p, "-"))
#'   causal <- c(50, 150)
#'   b <- numeric(p); b[causal] <- 5
#'   z <- as.numeric(R %*% b + crossprod(chol(R), rnorm(p)))
#'   omega <- rep(1, p)
#'   omega[causal] <- 100  # strong annotation
#'   fit <- vsore(z, R, omega, n = n)
#'   sum(fit$pip > 0.5)
#'   fit$tau2
#' }
#'
#' @export
vsore <- function(z, R, anchor_weights, n,
                  L = 10, K0 = 20, T_max = 15,
                  tol = 1e-3, eps = 1e-6, n_inner = 100,
                  verbose = FALSE) {
  if (!requireNamespace("susieR", quietly = TRUE)) {
    stop("The 'susieR' package is required. Install with: install.packages('susieR').")
  }

  z <- as.numeric(z)
  p <- length(z)
  if (!is.matrix(R) || nrow(R) != p || ncol(R) != p) {
    stop("R must be a p x p matrix matching length(z).")
  }
  if (length(anchor_weights) != p || any(anchor_weights <= 0)) {
    stop("anchor_weights must be a positive numeric vector of length p.")
  }

  uniform_weights <- isTRUE(all.equal(diff(range(anchor_weights)), 0,
                                      tolerance = 1e-10))

  fit <- susieR::susie_rss(z = z, R = R, n = n, L = L, max_iter = n_inner,
                           prior_weights = anchor_weights /
                             sum(anchor_weights))
  pips <- susie_pips(fit)

  if (uniform_weights) {
    return(list(pip = pips, tau2 = NA_real_,
                weights = anchor_weights, n_iter = 0,
                converged = TRUE, susie_fit = fit))
  }

  IK0 <- order(anchor_weights, decreasing = TRUE)[seq_len(min(K0, p))]
  log_w0 <- log(anchor_weights)
  log_w <- log_w0

  converged <- FALSE
  tau2 <- NA_real_

  for (t in seq_len(T_max)) {
    rho_s <- spearman_top_k(anchor_weights[IK0], pips[IK0])
    tau2  <- max(0.05, 1.5 * (1 - rho_s))

    pips_clip <- pmin(pmax(pips, eps), 1 - eps)
    log_w_new <- (log_w0 + tau2 * log(pips_clip)) / (1 + tau2)

    if (verbose) {
      message(sprintf("[V-SoRE] iter %02d  rho_s=%.3f  tau2=%.3f", t, rho_s, tau2))
    }

    fit <- susieR::susie_rss(z = z, R = R, n = n, L = L,
                             max_iter = n_inner,
                             prior_weights = exp(log_w_new) /
                               sum(exp(log_w_new)))
    pips <- susie_pips(fit)

    if (max(abs(log_w_new - log_w)) < tol && t >= 3) {
      converged <- TRUE
      log_w <- log_w_new
      break
    }
    log_w <- log_w_new
  }

  list(pip = pips, tau2 = tau2,
       weights = exp(log_w), n_iter = t,
       converged = converged, susie_fit = fit)
}

#' Extract PIPs from a susieR fit object
#'
#' Wrapper that handles different susieR versions.
#' @keywords internal
#' @noRd
susie_pips <- function(fit) {
  if (!is.null(fit$pip)) {
    as.numeric(fit$pip)
  } else if (!is.null(fit$alpha)) {
    as.numeric(1 - apply(1 - fit$alpha, 2, prod))
  } else {
    stop("Unable to extract PIPs from susieR fit object.")
  }
}

#' Spearman rank correlation on the top-K subset
#'
#' Compute Spearman's rank correlation between two vectors restricted to
#' their top-K positions by the first vector. Used internally by
#' \code{\link{vsore}} to estimate tau^2.
#'
#' @param x Numeric vector of weights or scores.
#' @param y Numeric vector of PIPs.
#'
#' @return Scalar Spearman correlation in [-1, 1].
#' @export
spearman_top_k <- function(x, y) {
  if (length(x) != length(y)) stop("x and y must have the same length.")
  cor(rank(x, ties.method = "average"),
      rank(y, ties.method = "average"),
      method = "pearson")
}

#' Alias for \code{\link{vsore}}
#'
#' Provided for backward compatibility with earlier scripts.
#' @inheritParams vsore
#' @export
vsore_fit <- function(z, R, anchor_weights, n, ...) {
  vsore(z, R, anchor_weights, n, ...)
}
