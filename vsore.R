vsore <- function(z, R, omega0, n_eff, L = 10, K0 = 20,
                  T_max = 15, delta = 1e-3) {

  p           <- length(z)
  omega0      <- as.numeric(omega0)
  log_om0     <- log(omega0 + 1e-300)
  I_K0        <- order(-omega0)[seq_len(min(K0, p))]

  fit         <- susieR::susie_rss(z, R, n = n_eff, L = L,
                                    prior_weights = omega0 / sum(omega0))
  gamma       <- susieR::susie_get_pip(fit)
  hat_tau2    <- 0.5
  log_om_prev <- log_om0
  n_iter      <- 0L

  for (t in seq_len(T_max)) {
    n_iter   <- t
    rs       <- .spearman_r(omega0[I_K0], gamma[I_K0])
    ht       <- max(0.05, 1.5 * (1 - rs))
    gc       <- pmin(pmax(gamma, 1e-6), 1 - 1e-6)
    log_ot   <- (log_om0 + ht * log(gc)) / (1 + ht)
    ot       <- exp(log_ot)

    fit_new  <- susieR::susie_rss(z, R, n = n_eff, L = L,
                                   prior_weights = ot / sum(ot))
    gamma    <- susieR::susie_get_pip(fit_new)
    dw       <- max(abs(log_ot - log_om_prev))
    log_om_prev <- log_ot
    fit      <- fit_new
    hat_tau2 <- ht

    if (t >= 3L && dw < delta) break
  }

  cs <- susieR::susie_get_cs(fit)$cs
  cs <- cs[lengths(cs) > 0]

  list(fit      = fit,
       pip      = gamma,
       hat_tau2 = hat_tau2,
       cs       = cs,
       n_iter   = n_iter)
}



vsore_mallows <- function(z, R, omega0, n_eff, L = 10, K0 = 20,
                           T_max = 15, delta = 1e-3) {
  vsore(z, R, sqrt(as.numeric(omega0) + 1e-300),
        n_eff = n_eff, L = L, K0 = K0, T_max = T_max, delta = delta)
}


## internal Spearman helper (not exported)
.spearman_r <- function(x, y) {
  n <- length(x)
  if (n < 3L) return(0)
  1 - 6 * sum((rank(x) - rank(y))^2) / (n * (n^2 - 1))
}
