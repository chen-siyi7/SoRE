
evaluate_cs <- function(causal_idx, cs_list) {
  cs_list <- cs_list[lengths(cs_list) > 0]
  if (length(cs_list) == 0)
    return(c(coverage = 0.0, cs_size = 0.0))
  cov <- mean(sapply(causal_idx,
                     function(j) any(sapply(cs_list, function(cs) j %in% cs))))
  sz  <- mean(sapply(cs_list, length))
  c(coverage = as.numeric(cov), cs_size = as.numeric(sz))
}



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
