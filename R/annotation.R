# =============================================================================
# annotation.R — anchor weight construction from external annotations
# =============================================================================

#' Build SoRE anchor weights from functional annotation sources
#'
#' Constructs the multiplicative anchor weight vector tilde-omega_j used as
#' input to V-SoRE or SoRE Gibbs. The default construction matches equation
#' (8) in the manuscript and combines four annotation sources: coding
#' consequence, cell-type-specific chromatin accessibility, evolutionary
#' conservation, and prior-GWAS z-statistics.
#'
#' @param coding Optional integer vector of length p giving coding category:
#'   \code{3} = protein-truncating (LoF), \code{2} = missense,
#'   \code{1} = synonymous or other. Default NULL.
#' @param chrom_open Optional logical vector of length p indicating overlap
#'   with a cell-type-specific open-chromatin peak. Default NULL.
#' @param conserved Optional logical vector of length p indicating evolutionary
#'   conservation (e.g., GERP > 4). Default NULL.
#' @param z_prior Optional numeric vector of length p giving marginal
#'   z-scores from a related prior GWAS. Default NULL.
#' @param p Optional explicit length. Required if all inputs are NULL.
#'
#' @return Numeric vector of strictly positive anchor weights, length p.
#'
#' @details
#' The default coding multiplier is \code{c(50, 10, 1)} corresponding to
#' protein-truncating, missense, synonymous. The chromatin multiplier is 5
#' for in-peak variants and 1 otherwise. The conservation multiplier is 3
#' for conserved variants and 1 otherwise. The prior-GWAS contribution is
#' \code{1 + |z_prior|}. Any subset of inputs may be omitted; the missing
#' factors default to 1. Setting all inputs to NULL returns a uniform
#' weight vector, which recovers SuSiE-RSS.
#'
#' @examples
#' p <- 200
#' set.seed(1)
#' w <- build_anchor_weights(
#'   coding     = sample(1:3, p, replace = TRUE, prob = c(0.95, 0.04, 0.01)),
#'   chrom_open = runif(p) < 0.05,
#'   conserved  = runif(p) < 0.10,
#'   z_prior    = rnorm(p, sd = 0.5)
#' )
#' summary(w)
#' length(w)
#'
#' @export
build_anchor_weights <- function(coding = NULL, chrom_open = NULL,
                                  conserved = NULL, z_prior = NULL,
                                  p = NULL) {
  lengths <- c(length(coding), length(chrom_open),
               length(conserved), length(z_prior))
  lengths <- lengths[lengths > 0]
  if (length(lengths) == 0) {
    if (is.null(p)) stop("Provide p when all annotation arguments are NULL.")
  } else {
    if (length(unique(lengths)) > 1) {
      stop("All annotation inputs must have the same length.")
    }
    p_inferred <- unique(lengths)
    if (!is.null(p) && p != p_inferred) {
      stop("p does not match the length of the annotation inputs.")
    }
    p <- p_inferred
  }

  w <- rep(1, p)

  if (!is.null(coding)) {
    cod_map <- c(`1` = 1, `2` = 10, `3` = 50)
    w <- w * cod_map[as.character(coding)]
  }
  if (!is.null(chrom_open)) {
    w <- w * ifelse(as.logical(chrom_open), 5, 1)
  }
  if (!is.null(conserved)) {
    w <- w * ifelse(as.logical(conserved), 3, 1)
  }
  if (!is.null(z_prior)) {
    w <- w * (1 + abs(as.numeric(z_prior)))
  }

  as.numeric(w)
}
