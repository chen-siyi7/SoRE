# =============================================================================
# credible_sets.R — greedy credible-set construction from PIPs
# =============================================================================

#' Construct credible sets from posterior inclusion probabilities
#'
#' Builds credible sets greedily by sorting variants in decreasing PIP order
#' and adding them until the cumulative sum first exceeds the nominal level.
#' This is the standard construction used by SuSiE-RSS.
#'
#' @param pip Numeric vector of posterior inclusion probabilities, length p.
#' @param coverage Nominal credible-set coverage (default 0.95).
#' @param min_purity Optional minimum CS purity (median minimum absolute
#'   pairwise LD within the CS). If specified, an LD matrix must be supplied.
#' @param R Optional p-by-p LD matrix, required if \code{min_purity} is set.
#'
#' @return A list with elements
#'   \describe{
#'     \item{\code{cs}}{Integer vector of variant indices in the credible set.}
#'     \item{\code{size}}{Length of the credible set.}
#'     \item{\code{coverage}}{Achieved cumulative PIP.}
#'     \item{\code{purity}}{Median minimum |LD| within the CS, if R supplied.}
#'     \item{\code{lead_pip}}{Maximum PIP in the credible set.}
#'   }
#'
#' @examples
#' set.seed(1)
#' pips <- c(0.6, 0.3, 0.05, 0.04, 0.01)
#' credible_sets(pips, coverage = 0.95)
#'
#' @export
credible_sets <- function(pip, coverage = 0.95,
                          min_purity = NULL, R = NULL) {
  pip <- as.numeric(pip)
  if (any(pip < 0 | pip > 1)) stop("pip must be in [0, 1].")
  if (coverage <= 0 || coverage >= 1) stop("coverage must be in (0, 1).")

  ord <- order(pip, decreasing = TRUE)
  cum_pip <- cumsum(pip[ord])
  k <- which(cum_pip >= coverage)[1]
  if (is.na(k)) k <- length(pip)

  cs_idx <- sort(ord[seq_len(k)])
  purity <- NA_real_
  if (!is.null(min_purity)) {
    if (is.null(R)) stop("R must be supplied if min_purity is set.")
    if (length(cs_idx) == 1) {
      purity <- 1
    } else {
      sub <- abs(R[cs_idx, cs_idx])
      diag(sub) <- NA
      purity <- median(apply(sub, 1, min, na.rm = TRUE), na.rm = TRUE)
    }
  }

  list(cs = cs_idx,
       size = length(cs_idx),
       coverage = cum_pip[k],
       purity = purity,
       lead_pip = max(pip[cs_idx]))
}
