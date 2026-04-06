#' Run V-SoRE: variational SoRE fine-mapping from GWAS summary statistics
#'
#' V-SoRE wraps SuSiE-RSS with an outer reweighting loop driven by a
#' Plackett-Luce annotation prior. The annotation-data agreement parameter
#' \eqn{\hat\tau^2} is estimated per locus via the Spearman correlation between
#' anchor weights and current PIPs within the top-\eqn{K_0} annotation window.
#'
#' @param z Numeric vector of length \eqn{p}: marginal \eqn{z}-scores.
#' @param R Numeric \eqn{p \times p} LD (correlation) matrix, ridge-regularised.
#' @param omega0 Numeric vector of length \eqn{p}: anchor annotation weights
#'   \eqn{\tilde\omega_j > 0}. Need not sum to 1.
#' @param n_eff Positive integer: effective GWAS sample size passed to
#'   \code{susie_rss} as \code{n}. For binary traits use \eqn{N_\mathrm{eff}}
#'   from \code{\link{compute_neff}}.
#' @param L Integer: number of SuSiE single-effect components (default 10).
#' @param K0 Integer: annotation window size for \eqn{\hat\tau^2} estimation
#'   (default 20). Should satisfy \eqn{K_0 \geq 2s_0 + 5}.
#' @param T_max Integer: maximum outer iterations (default 15).
#' @param delta Numeric: convergence threshold on max log-weight change
#'   (default \eqn{10^{-3}}).
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{fit}}{SuSiE fit object from the final outer iteration.}
#'     \item{\code{pip}}{Numeric vector of length \eqn{p}: posterior inclusion
#'       probabilities.}
#'     \item{\code{hat_tau2}}{Scalar: estimated \eqn{\hat\tau^2} at convergence.
#'       Values near 0 indicate strong annotation-data agreement; values above 1
#'       indicate the data substantially revised the annotation ordering.}
#'     \item{\code{cs}}{List of 95\% credible sets from \code{susie_get_cs}.}
#'     \item{\code{n_iter}}{Integer: number of outer iterations completed.}
#'   }
#'
#' @references
#' Chen, S. (2025). Bayesian Fine-Mapping with an Annotation Ranking Prior
#' Using GWAS Summary Statistics. \emph{Biometrics}.
#'
#' Wang, G., Sarkar, A., Carbonetto, P., and Stephens, M. (2020). A simple new
#' approach to variable selection in regression, with application to genetic
#' fine-mapping. \emph{JRSS-B}, 82, 1273-1300.
#'
#' @seealso \code{\link{make_weights}}, \code{\link{compute_neff}},
#'   \code{\link{vsore_mallows}}
#'
#' @examples
#' \dontrun{
#' R    <- build_ld(p = 200)
#' z    <- rnorm(200)
#' w    <- make_weights(causal_idx = c(50, 51), p = 200,
#'                      pve = 0.01, s0 = 2, N = 50000, R = R)
#' res  <- vsore(z, R, w, n_eff = 50000)
#' res$hat_tau2
#' susieR::susie_plot(res$fit, y = "PIP")
#' }
#'
#' @importFrom susieR susie_rss susie_get_pip susie_get_cs
#' @export
vsore <- function(z, R, omega0, n_eff, L = 10, K0 = 20,
                  T_max = 15, delta = 1e-3) {

  p           <- length(z)
  omega0      <- as.numeric(omega0)
  log_om0     <- log(omega0 + 1e-300)
  I_K0        <- order(-omega0)[seq_len(min(K0, p))]

  fit         <- susieR::susie_rss(z, R, n = n_eff, L = L,
                                    prior_weights = omega0 / sum(omega0))
  gamma       <- susieR::susie_get_pip(fit)
  hat_tau2    <- 0.5
  log_om_prev <- log_om0
  n_iter      <- 0L

  for (t in seq_len(T_max)) {
    n_iter   <- t
    rs       <- .spearman_r(omega0[I_K0], gamma[I_K0])
    ht       <- max(0.05, 1.5 * (1 - rs))
    gc       <- pmin(pmax(gamma, 1e-6), 1 - 1e-6)
    log_ot   <- (log_om0 + ht * log(gc)) / (1 + ht)
    ot       <- exp(log_ot)

    fit_new  <- susieR::susie_rss(z, R, n = n_eff, L = L,
                                   prior_weights = ot / sum(ot))
    gamma    <- susieR::susie_get_pip(fit_new)
    dw       <- max(abs(log_ot - log_om_prev))
    log_om_prev <- log_ot
    fit      <- fit_new
    hat_tau2 <- ht

    if (t >= 3L && dw < delta) break
  }

  cs <- susieR::susie_get_cs(fit)$cs
  cs <- cs[lengths(cs) > 0]

  list(fit      = fit,
       pip      = gamma,
       hat_tau2 = hat_tau2,
       cs       = cs,
       n_iter   = n_iter)
}


#' Run M-SoRE: Mallows-approximated V-SoRE via square-root-flattened weights
#'
#' M-SoRE runs V-SoRE with \eqn{\sqrt{\tilde\omega_j}} in place of
#' \eqn{\tilde\omega_j}, reducing annotation contrast by half on the log scale.
#' This approximates the Mallows prior's symmetry among near-equal-weight
#' variants and may be preferable in regions of near-perfect LD such as the MHC.
#'
#' @inheritParams vsore
#' @inherit vsore return
#' @seealso \code{\link{vsore}}
#' @export
vsore_mallows <- function(z, R, omega0, n_eff, L = 10, K0 = 20,
                           T_max = 15, delta = 1e-3) {
  vsore(z, R, sqrt(as.numeric(omega0) + 1e-300),
        n_eff = n_eff, L = L, K0 = K0, T_max = T_max, delta = delta)
}


## internal Spearman helper (not exported)
.spearman_r <- function(x, y) {
  n <- length(x)
  if (n < 3L) return(0)
  1 - 6 * sum((rank(x) - rank(y))^2) / (n * (n^2 - 1))
}
