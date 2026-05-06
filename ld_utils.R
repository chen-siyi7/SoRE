## =============================================================================
## LD matrix utilities for the simulation studies of Sections 5 and S12.
##
## build_ld() constructs a block AR(1) LD matrix used as the in-sample LD.
## make_mismatch_ld() applies sub-Gaussian element-wise noise of order
## 1/sqrt(n_ref) to simulate reference-panel mismatch.
## =============================================================================

#' Build a block AR(1) LD matrix.
#'
#' Each block of the LD matrix is an AR(1) correlation matrix with parameter
#' rho. The full p x p LD matrix is block-diagonal with n_blocks blocks of
#' size p / n_blocks each. A small ridge is then added to ensure positive
#' definiteness.
#'
#' @param p        Total number of variants. Must be divisible by n_blocks.
#' @param n_blocks Number of AR(1) blocks.
#' @param rho      AR(1) correlation parameter, in (-1, 1).
#' @param ridge    Diagonal ridge added to the LD matrix. Default 0.01.
#' @return p x p positive-definite LD matrix.
#' @export
build_ld <- function(p = 500, n_blocks = 5, rho = 0.9, ridge = 0.01) {
  if (p %% n_blocks != 0L)
    stop("p must be divisible by n_blocks; got p = ", p,
         ", n_blocks = ", n_blocks)
  if (abs(rho) >= 1)
    stop("rho must be strictly between -1 and 1.")
  if (ridge < 0 || ridge > 1)
    stop("ridge must be in [0, 1].")

  block <- p %/% n_blocks
  d <- seq_len(block)
  ar1 <- rho ^ abs(outer(d, d, "-"))

  if (requireNamespace("Matrix", quietly = TRUE)) {
    R <- as.matrix(Matrix::bdiag(replicate(n_blocks, ar1, simplify = FALSE)))
  } else {
    ## Fallback when Matrix is not available.
    R <- matrix(0, p, p)
    for (b in seq_len(n_blocks)) {
      idx <- ((b - 1L) * block + 1L):(b * block)
      R[idx, idx] <- ar1
    }
  }
  (1 - ridge) * R + ridge * diag(p)
}


#' Build a reference-panel-mismatched LD matrix from a true LD matrix.
#'
#' Adds element-wise sub-Gaussian noise of order 1/sqrt(n_ref) to the off-diagonal
#' entries of R_true, symmetrizes the result, ensures positive definiteness by
#' a small eigen-shift, re-scales to a correlation matrix, and applies a final
#' ridge regularization. Used for the reference-panel mismatch experiment in
#' Supplementary Section S12.
#'
#' Does NOT call set.seed() inside the function; pass a generator-state
#' argument or call set.seed() in the caller for reproducibility.
#'
#' @param R_true  True p x p LD matrix.
#' @param n_ref   Reference-panel sample size; controls the noise magnitude.
#' @param ridge   Diagonal ridge applied after the noise-and-symmetrize steps.
#' @param eps_min Minimum eigenvalue floor before ridging. Default 0.01.
#' @return p x p reference-panel-mismatched LD matrix (correlation form).
#' @export
make_mismatch_ld <- function(R_true, n_ref, ridge = 0.01, eps_min = 0.01) {
  R_true <- as.matrix(R_true)
  p <- nrow(R_true)
  if (ncol(R_true) != p)
    stop("R_true must be square.")
  if (n_ref <= 0)
    stop("n_ref must be positive.")

  noise <- matrix(stats::rnorm(p * p, sd = 1 / sqrt(n_ref)), p, p)
  noise <- (noise + t(noise)) / 2

  R_m <- R_true + noise
  R_m <- (R_m + t(R_m)) / 2  # numerical symmetrization

  ev_min <- min(eigen(R_m, symmetric = TRUE, only.values = TRUE)$values)
  if (ev_min < eps_min)
    R_m <- R_m + (eps_min - ev_min) * diag(p)

  d <- sqrt(diag(R_m))
  R_m <- R_m / outer(d, d)
  (1 - ridge) * R_m + ridge * diag(p)
}


#' Sample z-scores from the working RSS likelihood.
#'
#' Generates z = R %*% b + xi with xi ~ N(0, R), the RSS likelihood used in
#' the manuscript Section 2.3.
#'
#' @param R  Positive-definite p x p LD matrix.
#' @param b  Length-p vector of true scaled effects.
#' @return Numeric vector of length p.
#' @export
sample_z <- function(R, b) {
  R <- as.matrix(R)
  b <- as.numeric(b)
  p <- length(b)
  if (nrow(R) != p || ncol(R) != p)
    stop("R must be a p x p matrix matching length(b) = ", p, ".")
  L <- t(chol(R))
  as.numeric(R %*% b + L %*% stats::rnorm(p))
}
