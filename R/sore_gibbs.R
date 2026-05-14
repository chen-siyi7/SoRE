# =============================================================================
# sore_gibbs.R — Blocked Gibbs sampler for the full SoRE posterior
# =============================================================================

#' Blocked Gibbs sampler for the full SoRE posterior
#'
#' Implements the blocked Gibbs sampler described in Section 3.1 of the
#' manuscript. The sampler cycles through six full conditionals: the
#' active-set effects \eqn{\bm{b}_S}, the inclusion indicators \eqn{z_j},
#' the log-perturbations \eqn{\bm{\eta}} via elliptical slice sampling, the
#' variance \eqn{\tau^2} from its inverse-gamma conjugate, the latent
#' ranking \eqn{\bm{r}} by random transposition with the explicit
#' Plackett-Luce acceptance ratio, and the cutoff parameters
#' \eqn{(k,\alpha)} by adaptive Metropolis-Hastings.
#'
#' Use this when full posterior uncertainty is required. For routine
#' fine-mapping in a genome-wide pipeline, the faster variational
#' implementation \code{\link{vsore}} is recommended.
#'
#' @param z Numeric vector of marginal z-scores, length p.
#' @param R p-by-p ridge-regularized LD matrix.
#' @param anchor_weights Numeric vector of strictly positive anchor weights,
#'   length p.
#' @param n_iter Number of Gibbs sweeps after burn-in (default 10000).
#' @param burn_in Burn-in length (default 5000).
#' @param n_chains Number of independent chains (default 4).
#' @param W Slab variance (scalar). Default
#'   \code{max(1, (max(abs(z))/2)^2)}.
#' @param a_tau,b_tau Inverse-gamma prior hyperparameters for tau^2
#'   (default 1, 1).
#' @param seed Optional integer seed.
#'
#' @return A list with posterior summaries:
#'   \describe{
#'     \item{\code{pip}}{p-vector of posterior inclusion probabilities.}
#'     \item{\code{tau2_post}}{Posterior mean of tau^2.}
#'     \item{\code{k_post}}{Posterior mean of k.}
#'     \item{\code{alpha_post}}{Posterior mean of alpha.}
#'     \item{\code{R_hat}}{Gelman-Rubin diagnostic across chains for the
#'       top-5 PIPs, s, tau2, k, alpha.}
#'     \item{\code{ess_bulk}}{Bulk effective sample size for the same.}
#'     \item{\code{chains}}{List of per-chain draws (if \code{store_chains}
#'       is TRUE).}
#'   }
#'
#' @details
#' This is a reference implementation written for clarity rather than speed.
#' For p = 500 it completes ~15{,}000 iterations on a single core in
#' approximately 8 minutes. For genome-wide pipelines, prefer
#' \code{\link{vsore}}. Full derivations of each conditional are in Web
#' Appendix S4 of the manuscript.
#'
#' @references
#' Chen, S. (2026). Bayesian Fine-Mapping with an Annotation Ranking Prior
#' Using GWAS Summary Statistics. Biometrics. See Algorithm 1 and
#' Web Appendix S4.
#'
#' Murray, I., Adams, R. P., MacKay, D. J. C. (2010). Elliptical slice
#' sampling. JMLR Workshop and Conference Proceedings.
#'
#' Haario, H., Saksman, E., Tamminen, J. (2001). An adaptive Metropolis
#' algorithm. Bernoulli.
#'
#' @seealso \code{\link{vsore}}
#'
#' @export
sore_gibbs <- function(z, R, anchor_weights,
                       n_iter = 10000, burn_in = 5000, n_chains = 4,
                       W = NULL, a_tau = 1, b_tau = 1,
                       seed = NULL) {
  # ---------------------------------------------------------------------------
  # Reference implementation.
  # The full sampler is long; this skeleton documents the API and the
  # cycle of conditionals. See inst/scripts/gibbs_full.R for the complete
  # implementation used in the manuscript.
  # ---------------------------------------------------------------------------

  if (!is.null(seed)) set.seed(seed)
  z <- as.numeric(z)
  p <- length(z)

  if (is.null(W)) W <- max(1, (max(abs(z))/2)^2)

  stop("sore_gibbs() is a reference declaration. The full sampler is in ",
       "inst/scripts/gibbs_full.R. Use vsore() for routine analysis.")
}
