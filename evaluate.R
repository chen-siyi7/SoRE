#' Evaluate credible-set coverage and size
#'
#' Given a set of known causal variants and a list of credible sets, computes
#' per-signal coverage (fraction of causal variants captured by at least one
#' credible set) and mean credible-set size.
#'
#' @param causal_idx Integer vector: indices of true causal variants (1-based).
#' @param cs_list List of integer vectors: credible sets, as returned by
#'   \code{susieR::susie_get_cs(fit)$cs}. Empty sets are ignored.
#'
#' @return Named numeric vector with elements:
#'   \describe{
#'     \item{\code{coverage}}{Proportion of causal variants in \code{causal_idx}
#'       contained in at least one credible set. Values in \eqn{[0, 1]}.}
#'     \item{\code{cs_size}}{Mean number of variants across all non-empty
#'       credible sets.}
#'   }
#'   Returns \code{c(coverage = 0, cs_size = 0)} when \code{cs_list} is empty
#'   or contains only empty sets.
#'
#' @examples
#' cs <- list(c(48, 49, 50, 51), c(120, 121))
#' evaluate_cs(causal_idx = c(50, 120), cs_list = cs)
#'
#' @export
evaluate_cs <- function(causal_idx, cs_list) {
  cs_list <- cs_list[lengths(cs_list) > 0]
  if (length(cs_list) == 0)
    return(c(coverage = 0.0, cs_size = 0.0))
  cov <- mean(sapply(causal_idx,
                     function(j) any(sapply(cs_list, function(cs) j %in% cs))))
  sz  <- mean(sapply(cs_list, length))
  c(coverage = as.numeric(cov), cs_size = as.numeric(sz))
}


#' Summarise V-SoRE results for a single locus
#'
#' Extracts the most useful summaries from a \code{vsore} output list: the
#' lead variant (highest PIP), credible-set sizes, purity, and \eqn{\hat\tau^2}.
#'
#' @param res List: output from \code{\link{vsore}}.
#' @param variant_names Optional character vector of length \eqn{p}: variant
#'   identifiers (rsIDs or chr:pos). If \code{NULL}, variants are numbered 1 to p.
#'
#' @return A data frame with one row per credible set and columns:
#'   \describe{
#'     \item{\code{cs_index}}{Credible set index.}
#'     \item{\code{cs_size}}{Number of variants in the credible set.}
#'     \item{\code{lead_pip}}{PIP of the highest-PIP variant in the set.}
#'     \item{\code{lead_variant}}{Name or index of the lead variant.}
#'     \item{\code{hat_tau2}}{Locus-level \eqn{\hat\tau^2}.}
#'   }
#'   Returns a zero-row data frame if no credible sets were detected.
#'
#' @seealso \code{\link{vsore}}, \code{\link{evaluate_cs}}
#'
#' @examples
#' \dontrun{
#' res <- vsore(z, R, omega0, n_eff = 50000)
#' summarise_vsore(res)
#' }
#'
#' @export
summarise_vsore <- function(res, variant_names = NULL) {
  cs   <- res$cs
  pip  <- res$pip
  p    <- length(pip)
  nms  <- if (!is.null(variant_names)) variant_names else seq_len(p)

  if (length(cs) == 0)
    return(data.frame(cs_index    = integer(0),
                      cs_size     = integer(0),
                      lead_pip    = numeric(0),
                      lead_variant = character(0),
                      hat_tau2    = numeric(0)))

  rows <- lapply(seq_along(cs), function(i) {
    idx      <- cs[[i]]
    lead_j   <- idx[which.max(pip[idx])]
    data.frame(cs_index     = i,
               cs_size      = length(idx),
               lead_pip     = pip[lead_j],
               lead_variant = as.character(nms[lead_j]),
               hat_tau2     = res$hat_tau2,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
