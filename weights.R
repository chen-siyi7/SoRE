make_weights <- function(p, annot_score, aux_z = rep(0, p),
                          log_additive = TRUE) {
  stopifnot(length(annot_score) == p, length(aux_z) == p)
  a <- as.numeric(annot_score)
  exp(0.8 * a) * (1 + abs(as.numeric(aux_z)))
}


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

compute_neff <- function(n_cases, n_controls, prevalence) {
  N    <- n_cases + n_controls
  ybar <- n_cases / N
  K    <- prevalence
  phi  <- dnorm(qnorm(K))
  4 * n_cases * n_controls / N *
    K^2 * (1 - K)^2 / (phi^2 * ybar * (1 - ybar))
}
