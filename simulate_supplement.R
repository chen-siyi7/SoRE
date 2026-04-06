set.seed(42)
suppressPackageStartupMessages({
  library(susieR)   # >= 0.12.35
  library(Matrix)
  library(dplyr)
})

## ── parameters ───────────────────────────────────────────────────────────────
P      <- 500
N_GWAS <- 50000

## ── shared helpers ────────────────────────────────────────────────────────────
build_ld <- function(p = 500, n_blocks = 5, rho = 0.9, ridge = 0.01) {
  block <- p / n_blocks
  d     <- seq_len(block)
  ar1   <- rho ^ abs(outer(d, d, "-"))
  R     <- as.matrix(bdiag(replicate(n_blocks, ar1, simplify = FALSE)))
  (1 - ridge) * R + ridge * diag(p)
}

make_weights <- function(causal_idx, p, pve, s0,
                         N = 50000, rg = 0.6, R = NULL,
                         U_range = c(2.5, 3.5), degraded = FALSE) {
  a <- runif(p)
  if (!degraded)
    for (j in causal_idx)
      a[j] <- min(a[j] + runif(1, U_range[1], U_range[2]), 4.0)
  b2 <- numeric(p)
  if (length(causal_idx) > 0 && !degraded) {
    bp           <- numeric(p)
    bp[causal_idx] <- rnorm(s0, 0, sqrt(pve * N / s0))
    b2             <- rg * bp
  }
  z2 <- if (!is.null(R))
    as.numeric(R %*% b2 + t(chol(R)) %*% rnorm(p))
  else numeric(p)
  exp(0.8 * a) * (1 + abs(z2))
}

spearman_r <- function(x, y) {
  n <- length(x)
  1 - 6 * sum((rank(x) - rank(y))^2) / (n * (n^2 - 1))
}

vsore <- function(z, R, omega0, n_eff, L = 10, K0 = 20,
                  T_max = 15, delta = 1e-3) {
  log_om0     <- log(omega0 + 1e-300)
  I_K0        <- order(-omega0)[seq_len(K0)]
  fit         <- susie_rss(z, R, n = n_eff, L = L,
                           prior_weights = omega0 / sum(omega0))
  gamma       <- susie_get_pip(fit)
  hat_tau2    <- 0.5
  log_om_prev <- log_om0
  
  for (t in seq_len(T_max)) {
    rs       <- spearman_r(omega0[I_K0], gamma[I_K0])
    ht       <- max(0.05, 1.5 * (1 - rs))
    gc       <- pmin(pmax(gamma, 1e-6), 1 - 1e-6)
    log_ot   <- (log_om0 + ht * log(gc)) / (1 + ht)
    ot       <- exp(log_ot)
    fit_new  <- susie_rss(z, R, n = n_eff, L = L,
                          prior_weights = ot / sum(ot))
    gamma    <- susie_get_pip(fit_new)
    dw       <- max(abs(log_ot - log_om_prev))
    log_om_prev <- log_ot
    fit      <- fit_new
    hat_tau2 <- ht
    if (t >= 3 && dw < delta) break
  }
  list(pip = gamma, hat_tau2 = hat_tau2, fit = fit)
}

## M-SoRE: V-SoRE with sqrt-flattened weights
vsore_mallows <- function(z, R, omega0, n_eff, L = 10, K0 = 20,
                          T_max = 15, delta = 1e-3) {
  vsore(z, R, sqrt(omega0 + 1e-300), n_eff = n_eff, L = L,
        K0 = K0, T_max = T_max, delta = delta)
}

evaluate <- function(causal_idx, cs_list) {
  cs_list <- cs_list[lengths(cs_list) > 0]
  if (length(cs_list) == 0)
    return(c(coverage = 0.0, cs_size = 0.0))
  cov <- mean(sapply(causal_idx,
                     function(c) any(sapply(cs_list, function(cs) c %in% cs))))
  sz  <- mean(sapply(cs_list, length))
  c(coverage = as.numeric(cov), cs_size = as.numeric(sz))
}

## ── build LD ─────────────────────────────────────────────────────────────────
cat("Building LD matrix...\n")
R_true <- build_ld(p = P)

## ═══════════════════════════════════════════════════════════════════════════
## S11: tau^2 calibration
## 200 Monte Carlo loci per U_j level (increase N_CAL to 625 for publication)
## ═══════════════════════════════════════════════════════════════════════════
cat("\nS11: tau2 calibration...\n")

U_LEVELS <- c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0)
N_CAL    <- 200    # loci per level; 625 for publication quality
K0_CAL   <- 7      # = 2 * s0 + 5 with s0 = 1

cal_rows <- lapply(seq_along(U_LEVELS), function(k) {
  U <- U_LEVELS[k]
  set.seed(200 + round(U * 10))
  tau2_v <- numeric(N_CAL)
  surr_v <- numeric(N_CAL)
  
  for (rep in seq_len(N_CAL)) {
    causal    <- sort(sample(seq_len(P / 5), 1))
    b         <- numeric(P)
    b[causal] <- rnorm(1, 0, sqrt(0.01 * N_GWAS))
    z         <- as.numeric(R_true %*% b + t(chol(R_true)) %*% rnorm(P))
    
    degraded  <- (U == 0)
    u_range   <- if (U > 0) c(U, min(U + 0.5, 4.0)) else c(0, 0)
    omega0    <- make_weights(causal, P, 0.01, 1, N = N_GWAS, R = R_true,
                              U_range = u_range, degraded = degraded)
    
    res           <- vsore(z, R_true, omega0, n_eff = N_GWAS, L = 10, K0 = K0_CAL)
    tau2_v[rep]   <- res$hat_tau2
    
    I_K0        <- order(-omega0)[seq_len(K0_CAL)]
    rs          <- spearman_r(omega0[I_K0], res$pip[I_K0])
    surr_v[rep] <- max(0.05, 1.5 * (1 - rs))
  }
  
  mae <- mean(abs(tau2_v - surr_v))
  cat(sprintf("  U=%.1f  Gibbs=%.2f  Surr=%.2f  MAE=%.3f\n",
              U, mean(tau2_v), mean(surr_v), mae))
  
  data.frame(U_j            = U,
             mean_Gibbs     = mean(tau2_v),
             mean_surrogate = mean(surr_v),
             MAE            = mae,
             n              = N_CAL)
})

df_cal <- bind_rows(cal_rows)
write.csv(df_cal, "calibration_results.csv", row.names = FALSE)
cat("Saved calibration_results.csv\n")

## LaTeX table S11
lines_s11 <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Calibration of the $\\hat\\tau^2$ surrogate $\\max(0.05,\\,1.5(1-\\rho_s))$",
  "  against the V-SoRE variational posterior mean, stratified by annotation",
  "  informativeness $U_j$. $s_0=1$, PVE$=1\\%$; 200 loci per level.",
  "  MAE~=~mean absolute error.}",
  "\\label{t:tau_cal}",
  "\\small",
  "\\begin{tabular}{ccccc}",
  "\\toprule",
  "$U_j$ & Mean $\\hat\\tau^2_{\\mathrm{Gibbs}}$ & Mean $\\hat\\tau^2_{\\mathrm{surrogate}}$ & MAE & $n$ \\\\",
  "\\midrule",
  sapply(seq_len(nrow(df_cal)), function(i) {
    r <- df_cal[i, ]
    sprintf("%.1f & %.2f & %.2f & %.3f & %d \\\\",
            r$U_j, r$mean_Gibbs, r$mean_surrogate, r$MAE, r$n)
  }),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)
writeLines(lines_s11, "latex_table_S11.txt")
cat("Saved latex_table_S11.txt\n")

## ═══════════════════════════════════════════════════════════════════════════
## S13: IIA / Mallows comparison
## 30 replicates per configuration, C1-C6
## ═══════════════════════════════════════════════════════════════════════════
cat("\nS13: IIA/Mallows comparison...\n")

CONFIGS_IIA <- data.frame(
  label = c("C1","C2","C3","C4","C5","C6"),
  s0    = c(1, 2, 3, 1, 2, 3),
  pve   = c(.005, .005, .005, .02, .02, .02)
)
N_IIA <- 500

iia_rows <- lapply(seq_len(nrow(CONFIGS_IIA)), function(i) {
  cfg <- CONFIGS_IIA[i, ]
  set.seed(300 + i)
  
  corrs <- sz_m <- cov_m <- numeric(N_IIA)
  
  for (rep in seq_len(N_IIA)) {
    causal    <- sort(sample(seq_len(P / 5), cfg$s0))
    b         <- numeric(P)
    b[causal] <- rnorm(cfg$s0, 0, sqrt(cfg$pve * N_GWAS / cfg$s0))
    z         <- as.numeric(R_true %*% b + t(chol(R_true)) %*% rnorm(P))
    omega0    <- make_weights(causal, P, cfg$pve, cfg$s0, N = N_GWAS, R = R_true)
    K0        <- min(20, 2 * cfg$s0 + 5)
    
    res_v     <- vsore(z, R_true, omega0, n_eff = N_GWAS, L = 10, K0 = K0)
    res_m     <- vsore_mallows(z, R_true, omega0, n_eff = N_GWAS, L = 10, K0 = K0)
    
    corrs[rep] <- cor(res_v$pip, res_m$pip)
    cs_m_rep   <- susie_get_cs(res_m$fit)$cs
    ev_m       <- evaluate(causal, cs_m_rep)
    cov_m[rep] <- ev_m[["coverage"]]
    sz_m[rep]  <- ev_m[["cs_size"]]
  }
  
  cat(sprintf("  %s  corr=%.3f  M-SoRE cov=%3.0f%%  sz=%.1f\n",
              cfg$label, mean(corrs, na.rm = TRUE),
              100 * mean(cov_m), mean(sz_m)))
  
  data.frame(
    Config    = cfg$label,
    s0        = cfg$s0,
    PVE       = sprintf("%.1f\\%%", cfg$pve * 100),
    pip_corr  = mean(corrs, na.rm = TRUE),
    msore_cov = 100 * mean(cov_m),
    msore_sz  = mean(sz_m)
  )
})

df_iia <- bind_rows(iia_rows)
write.csv(df_iia, "iia_results.csv", row.names = FALSE)
cat("Saved iia_results.csv\n")
