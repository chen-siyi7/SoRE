## =============================================================================
## Annotation weights and effective sample size for binary GWAS traits.
## =============================================================================

#' Build a composite annotation weight on the log-additive scale.
#'
#' Combines per-source annotation indicators or scores into a single positive
#' weight per variant. Source-specific log-weights default to the values used
#' in the AD application of the manuscript Section 6.1: log(50) for coding
#' consequence, log(5) for chromatin accessibility, log(3) for conservation,
#' and 1.0 for the (additive) prior-GWAS term.
#'
#' Note: the prior-GWAS term in the AD application is added on the log scale as
#' log(1 + |z_prior|), which is constructed by the caller (see weights.R::
#' build_prior_gwas_term). This function is the linear combiner.
#'
#' @param annot_mat   p x K matrix of per-variant annotation values.
#' @param weights     Length-K vector of source-specific log-weights. If NULL,
#'                    defaults are used for K up to 4.
#' @param source_names Optional character vector naming the K columns.
#' @return Numeric vector of length p, on the log-weight scale (caller can
#'         exponentiate to get omega0 for vsore()).
#' @export
build_annot_score <- function(annot_mat, weights = NULL, source_names = NULL) {
  annot_mat <- as.matrix(annot_mat)
  K <- ncol(annot_mat)
  if (is.null(weights)) {
    defaults <- c(log(50), log(5), log(3), 1.0)
    weights  <- if (K <= 4L) defaults[seq_len(K)]
                else c(defaults, rep(1.0, K - 4L))
  }
  if (length(weights) != K)
    stop("length(weights) must equal ncol(annot_mat) = ", K)
  if (!is.null(source_names)) {
    if (length(source_names) != K)
      stop("length(source_names) must equal ncol(annot_mat) = ", K)
    colnames(annot_mat) <- source_names
  }
  as.numeric(annot_mat %*% weights)
}


#' Construct the additive prior-GWAS log-term used in the AD application.
#'
#' Returns log(1 + |z_prior|), which is non-negative, monotone in |z|, and
#' bounded above by log(1 + max|z|), preventing extreme prior-GWAS signals
#' from dominating the composite log-weight.
#'
#' @param z_prior Numeric vector of prior-GWAS z-scores (length p).
#' @return Numeric vector of length p.
#' @export
build_prior_gwas_term <- function(z_prior) {
  log1p(abs(as.numeric(z_prior)))
}


#' Convert an additive log-weight into a positive anchor weight for vsore().
#'
#' @param log_omega Numeric vector of log-anchor-weights, length p.
#' @return Numeric vector of strictly positive anchor weights, length p.
#' @export
make_anchor_weights <- function(log_omega) {
  exp(as.numeric(log_omega))
}


#' Effective sample size for binary GWAS traits under the liability-threshold model.
#'
#' Implements eq. (7) of the manuscript:
#'
#'   N_eff = 4 n1 n0 / N  *  K^2 (1-K)^2 / [phi(qnorm(K))^2 ybar (1-ybar)],
#'
#' where K is the population prevalence and ybar = n1/N is the sample case
#' fraction. With n1 = 21,982, n0 = 41,944, K = 0.05 the result is 54,238;
#' with K = 0.025 (an alternative early-AD prevalence) it is 41,856.
#'
#' This function returns the value computed from the formula. Users should
#' verify the prevalence value used matches the trait being analysed and
#' should report the resulting N_eff in their analysis.
#'
#' @param n_cases    Number of cases.
#' @param n_controls Number of controls.
#' @param prevalence Population prevalence K of the binary trait, in (0, 1).
#' @return Numeric scalar, the effective sample size.
#' @export
compute_neff <- function(n_cases, n_controls, prevalence) {
  if (n_cases <= 0 || n_controls <= 0)
    stop("n_cases and n_controls must be positive.")
  if (prevalence <= 0 || prevalence >= 1)
    stop("prevalence must be in (0, 1).")
  N <- n_cases + n_controls
  ybar <- n_cases / N
  K <- prevalence
  phi <- stats::dnorm(stats::qnorm(K))
  4 * n_cases * n_controls / N *
    K^2 * (1 - K)^2 / (phi^2 * ybar * (1 - ybar))
}


#' Effective sample size: simpler 4 n1 n0 / N formula (no liability-threshold correction).
#'
#' Some annotation pipelines and reference papers use the unadjusted form.
#' Provided for comparison; not used in the manuscript.
#'
#' @param n_cases    Number of cases.
#' @param n_controls Number of controls.
#' @return Numeric scalar.
#' @export
compute_neff_unadjusted <- function(n_cases, n_controls) {
  if (n_cases <= 0 || n_controls <= 0)
    stop("n_cases and n_controls must be positive.")
  4 * n_cases * n_controls / (n_cases + n_controls)
}
