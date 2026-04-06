build_ld <- function(p = 500, n_blocks = 5, rho = 0.9, ridge = 0.01) {
  block <- p / n_blocks
  d     <- seq_len(block)
  ar1   <- rho ^ abs(outer(d, d, "-"))
  R     <- as.matrix(Matrix::bdiag(replicate(n_blocks, ar1, simplify = FALSE)))
  (1 - ridge) * R + ridge * diag(p)
}



make_mismatch_ld <- function(R_true, nref, ridge = 0.01, seed = 99) {
  p <- nrow(R_true)
  set.seed(seed)
  noise <- matrix(rnorm(p * p, sd = 1 / sqrt(nref)), p, p)
  noise <- (noise + t(noise)) / 2
  R_m   <- R_true + noise
  R_m   <- (R_m + t(R_m)) / 2
  ev    <- eigen(R_m, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < 0.01)
    R_m <- R_m + (0.01 - min(ev)) * diag(p)
  d     <- sqrt(diag(R_m))
  R_m   <- R_m / outer(d, d)
  (1 - ridge) * R_m + ridge * diag(p)
}
