# =============================================================================
# diagnostics.R — tau-squared interpretation utility
# =============================================================================

#' Interpret the V-SoRE tau-squared diagnostic
#'
#' Translates the estimated tau-squared from a V-SoRE fit into a categorical
#' interpretation following the thresholds calibrated in the manuscript:
#' tau-squared below 0.2 indicates strong annotation-data concordance,
#' values in [0.2, 1] indicate moderate concordance, and values above 1
#' indicate that the data have substantially overridden the annotation
#' ordering.
#'
#' @param tau2 Numeric scalar or vector of tau-squared values.
#'
#' @return A character vector of length \code{length(tau2)} with values
#'   \code{"strong"}, \code{"moderate"}, \code{"uninformative"}, or
#'   \code{"no annotation"} (when tau2 is NA).
#'
#' @details
#' These thresholds were calibrated using the simulation sweep reported in
#' Web Appendix S11 of the manuscript: a tau-squared near 1.0 corresponds to
#' an annotation log-elevation U around 1.0, where the rank-correlation
#' diagnostic transitions from informative to uninformative.
#'
#' @examples
#' tau2_diagnostic(c(0.1, 0.5, 1.3, NA))
#'
#' @export
tau2_diagnostic <- function(tau2) {
  out <- character(length(tau2))
  out[is.na(tau2)] <- "no annotation"
  out[!is.na(tau2) & tau2 < 0.2] <- "strong"
  out[!is.na(tau2) & tau2 >= 0.2 & tau2 <= 1.0] <- "moderate"
  out[!is.na(tau2) & tau2 > 1.0] <- "uninformative"
  out
}
