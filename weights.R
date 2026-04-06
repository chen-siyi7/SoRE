#' Construct anchor annotation weights for V-SoRE
#'
#' Builds the anchor weight vector \eqn{\tilde\omega_j = \exp(0.8 a_j)
#' (1 + |\hat{z}_j^{(2)}|)} from a categorical annotation score \eqn{a_j}
#' and an auxiliary \eqn{z}-score from a correlated phenotype or prior GWAS.
#'
#' @param p Integer: number of variants in the locus window.
#' @param annot_score Numeric vector of length \eqn{p}: annotation scores
#'   \eqn{a_j \geq 0}. Higher values indicate greater prior probability of
#'   causality. Typical sources: binary coding/chromatin indicators (set to 0/1
#'   and rescaled), continuous conservation scores, or combinations thereof.
#' @param aux_z Numeric vector of length \eqn{p}: auxiliary \eqn{z}-scores
#'   from a correlated phenotype or a prior GWAS (e.g. Bellenguez et al. 2022
#'   for an AD analysis). Set to \code{rep(0, p)} to disable the auxiliary term.
#' @param log_additive Logical: if \code{TRUE} (default), combine annotation
#'   sources in \code{annot_score} on the log scale (i.e. \eqn{a_j} is the
#'   sum of log source weights). If \code{FALSE}, treat \code{annot_score}
#'   directly as \eqn{a_j}.
#'
#' @return Numeric vector of length \eqn{p}: anchor weights \eqn{\tilde\omega_j > 0}.
#'
#' @details
#' The weight formula \eqn{\tilde\omega_j = \exp(0.8 a_j)(1 + |\hat{z}_j^{(2)}|)}
#' combines a categorical annotation component (functional categories, chromatin
#' accessibility, conservation) with a continuous auxiliary association signal.
#' The factor 0.8 was chosen so that the expected log-weight ratio between a
#' high-annotation and null variant is approximately 2.4, corresponding to an
#' 11-fold weight ratio. The auxiliary term \eqn{(1 + |\hat{z}_j^{(2)}|)} ensures
#' all weights are strictly positive even when \eqn{a_j = 0}.
#'
#' For real GWAS analyses, \code{annot_score} can be constructed by
#' \code{\link{build_annot_score}} from multiple annotation sources.
#'
#' @seealso \code{\link{build_annot_score}}, \code{\link{vsore}}
#'
#' @examples
#' p   <- 500
#' a   <- runif(p)             # placeholder annotation scores
#' z2  <- rnorm(p)             # placeholder auxiliary z-scores
#' w   <- make_weights(p, a, z2)
#' summary(w)
#'
#' @export
make_weights <- function(p, annot_score, aux_z = rep(0, p),
                          log_additive = TRUE) {
  stopifnot(length(annot_score) == p, length(aux_z) == p)
  a <- as.numeric(annot_score)
  exp(0.8 * a) * (1 + abs(as.numeric(aux_z)))
}


#' Build a composite annotation score from multiple sources
#'
#' Combines binary and continuous annotation sources into a single log-additive
#' score \eqn{a_j = \sum_k c_k \cdot \mathbf{1}[\text{source}_k]} for use in
#' \code{\link{make_weights}}.
#'
#' @param annot_mat Logical or numeric matrix of dimensions \eqn{p \times K}:
#'   each column is one annotation source. Binary sources (0/1) are most common;
#'   continuous sources are accepted but treated as-is.
#' @param weights Numeric vector of length \eqn{K}: log-scale weights for each
#'   annotation source. Default values reflect approximate enrichment factors
#'   used in the paper: coding variants 3.91 (\eqn{\log 50}), open chromatin
#'   1.61 (\eqn{\log 5}), conservation 1.10 (\eqn{\log 3}), prior GWAS 1.0.
#'   Any positive values may be supplied.
#' @param source_names Optional character vector of length \eqn{K}: names for
#'   the annotation sources, used only for messages.
#'
#' @return Numeric vector of length \eqn{p}: composite annotation score
#'   \eqn{a_j = \sum_k \mathrm{weights}_k \cdot \mathrm{annot\_mat}_{jk}}.
#'
#' @seealso \code{\link{make_weights}}
#'
#' @examples
#' p   <- 500
#' K   <- 3
#' mat <- matrix(rbinom(p * K, 1, 0.1), p, K)
#' a   <- build_annot_score(mat, weights = c(3.91, 1.61, 1.10))
#' range(a)
#'
#' @export
build_annot_score <- function(annot_mat, weights = NULL,
                               source_names = NULL) {
  annot_mat <- as.matrix(annot_mat)
  p <- nrow(annot_mat)
  K <- ncol(annot_mat)

  if (is.null(weights)) {
    ## defaults from the paper: coding, chromatin, conservation, prior GWAS
    weights <- c(log(50), log(5), log(3), 1.0)[seq_len(K)]
    if (K > 4) weights <- c(weights, rep(1.0, K - 4))
  }
  stopifnot(length(weights) == K)

  if (!is.null(source_names))
    colnames(annot_mat) <- source_names

  as.numeric(annot_mat %*% weights)
}


#' Compute effective sample size for binary GWAS traits
#'
#' Maps case-control GWAS \eqn{z}-scores to the continuous-trait RSS scale
#' using the liability threshold model (Lee et al. 2012).
#'
#' @param n_cases Integer: number of cases.
#' @param n_controls Integer: number of controls.
#' @param prevalence Numeric in (0, 1): population disease prevalence \eqn{K}.
#'
#' @return Scalar: effective sample size \eqn{N_\mathrm{eff}}.
#'
#' @references
#' Lee, S. H., Wray, N. R., Goddard, M. E., and Visscher, P. M. (2012).
#' Estimating missing heritability for disease from GWAS.
#' \emph{AJHG}, 88, 294-305.
#'
#' @examples
#' compute_neff(n_cases = 21982, n_controls = 41944, prevalence = 0.05)
#'
#' @export
compute_neff <- function(n_cases, n_controls, prevalence) {
  N    <- n_cases + n_controls
  ybar <- n_cases / N
  K    <- prevalence
  phi  <- dnorm(qnorm(K))
  4 * n_cases * n_controls / N *
    K^2 * (1 - K)^2 / (phi^2 * ybar * (1 - ybar))
}
