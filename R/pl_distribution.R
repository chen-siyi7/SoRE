# =============================================================================
# pl_distribution.R — Plackett-Luce distribution sampler and log-mass
# =============================================================================

#' Sample a permutation from the Plackett-Luce distribution
#'
#' Draws a permutation r in S_p with probability proportional to the
#' Plackett-Luce mass
#'   PL(r | omega) = prod_{j=1}^p omega_{r(j)} / sum_{l=j}^p omega_{r(l)}.
#' Sampling uses the Gumbel-max construction (Yellott, 1977): permute
#' indices in decreasing order of log(omega_j) + G_j, G_j i.i.d.
#' Gumbel(0, 1). Total cost is O(p log p).
#'
#' @param omega Numeric vector of strictly positive Plackett-Luce weights.
#' @param n Number of permutations to draw. Default 1.
#'
#' @return If \code{n = 1}, an integer vector of length p giving a single
#'   permutation. If \code{n > 1}, an n-by-p integer matrix whose rows are
#'   independent draws.
#'
#' @references
#' Plackett, R. L. (1975). The analysis of permutations. Applied Statistics.
#'
#' Yellott, J. I. (1977). The relationship between Luce's choice axiom,
#' Thurstone's theory of comparative judgment, and the double exponential
#' distribution. Journal of Mathematical Psychology.
#'
#' @examples
#' set.seed(1)
#' omega <- c(10, 5, 3, 1, 1)
#' pl_sample(omega)
#' # 1000 draws: variant 1 should land at position 1 most often
#' table(pl_sample(omega, n = 1000)[, 1]) / 1000
#'
#' @export
pl_sample <- function(omega, n = 1) {
  if (any(omega <= 0)) stop("omega must be strictly positive.")
  p <- length(omega)
  log_omega <- log(omega)

  if (n == 1) {
    g <- -log(-log(runif(p)))  # Gumbel(0,1)
    return(order(log_omega + g, decreasing = TRUE))
  }

  out <- matrix(0L, nrow = n, ncol = p)
  for (i in seq_len(n)) {
    g <- -log(-log(runif(p)))
    out[i, ] <- order(log_omega + g, decreasing = TRUE)
  }
  out
}

#' Plackett-Luce log mass function
#'
#' Evaluates log PL(r | omega) for a given permutation. Useful for
#' Metropolis-Hastings acceptance ratios and posterior diagnostics.
#'
#' @param r Integer permutation vector of length p.
#' @param omega Numeric vector of strictly positive weights, length p.
#'
#' @return Scalar log-density.
#'
#' @examples
#' omega <- c(10, 5, 3, 1, 1)
#' pl_logmass(c(1, 2, 3, 4, 5), omega)
#' pl_logmass(c(5, 4, 3, 2, 1), omega)  # less likely
#'
#' @export
pl_logmass <- function(r, omega) {
  if (length(r) != length(omega)) stop("r and omega must have the same length.")
  if (any(omega <= 0)) stop("omega must be strictly positive.")
  p <- length(omega)
  log_terms <- numeric(p)
  remaining <- sum(omega[r])
  for (j in seq_len(p)) {
    log_terms[j] <- log(omega[r[j]]) - log(remaining)
    remaining <- remaining - omega[r[j]]
  }
  sum(log_terms)
}
