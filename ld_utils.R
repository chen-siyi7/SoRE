#' Build a block-diagonal AR(1) LD matrix (simulation utility)
#'
#' Constructs the block-diagonal LD matrix used in the SoRE simulation study:
#' \eqn{\bm{R} = (1-\lambda)\bigoplus_{b=1}^B \bm{A}_b + \lambda \bm{I}_p}
#' where \eqn{(\bm{A}_b)_{ij} = \rho^{|i-j|}}.
#'
#' @param p Integer: total number of variants (default 500).
#' @param n_blocks Integer: number of equal-size LD blocks (default 5).
#' @param rho Numeric in (0, 1): AR(1) correlation parameter within each block
#'   (default 0.9).
#' @param ridge Numeric: ridge regularisation constant added to the diagonal
#'   (default 0.01).
#'
#' @return Numeric \eqn{p \times p} positive-definite correlation matrix.
#'
#' @seealso \code{\link{make_mismatch_ld}}
#'
#' @examples
#' R <- build_ld(p = 200, n_blocks = 4, rho = 0.85)
#' range(eigen(R, only.values = TRUE)$values)
#'
#' @importFrom Matrix bdiag
#' @export
build_ld <- function(p = 500, n_blocks = 5, rho = 0.9, ridge = 0.01) {
  block <- p / n_blocks
  d     <- seq_len(block)
  ar1   <- rho ^ abs(outer(d, d, "-"))
  R     <- as.matrix(Matrix::bdiag(replicate(n_blocks, ar1, simplify = FALSE)))
  (1 - ridge) * R + ridge * diag(p)
}


#' Generate a mismatched reference LD matrix (simulation utility)
#'
#' Adds symmetric Gaussian noise to a true LD matrix to simulate the effect
#' of estimating LD from a finite reference panel of size \code{nref}, then
#' re-normalises to a valid correlation matrix and applies ridge regularisation.
#'
#' @param R_true Numeric \eqn{p \times p}: true LD matrix.
#' @param nref Integer: reference panel size. The noise standard deviation is
#'   \eqn{1/\sqrt{N_\mathrm{ref}}}, giving a relative Frobenius error of
#'   approximately \eqn{p/\sqrt{N_\mathrm{ref}}} before regularisation.
#' @param ridge Numeric: ridge constant applied after normalisation (default 0.01).
#' @param seed Integer: random seed (default 99).
#'
#' @return Numeric \eqn{p \times p} positive-definite correlation matrix.
#'
#' @seealso \code{\link{build_ld}}
#'
#' @examples
#' R      <- build_ld(p = 200)
#' R_mis  <- make_mismatch_ld(R, nref = 500)
#' norm(R_mis - R, "F") / norm(R, "F")   # relative Frobenius error
#'
#' @export
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
