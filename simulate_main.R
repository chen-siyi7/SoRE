## ============================================================================
## simulate_main.R
## Simulation study for main paper (Section 5, Table 1, Figure 2)
##
## Dependencies:
##   install.packages(c("susieR", "Matrix", "ggplot2", "dplyr", "MASS"))
##
## Outputs:
##   simulation_results.csv
##   fig_coverage.pdf
##   fig_cssize.pdf
##   fig_tau.pdf
##   latex_table.txt
## ============================================================================

set.seed(42)
suppressPackageStartupMessages({
  library(susieR)   # >= 0.12.35
  library(Matrix)
  library(MASS)
  library(ggplot2)
  library(dplyr)
})

## ── parameters ───────────────────────────────────────────────────────────────
P      <- 500
N_GWAS <- 50000
N_REPS <- 500

## ── 1. LD matrix ─────────────────────────────────────────────────────────────
build_ld <- function(p = 500, n_blocks = 5, rho = 0.9, ridge = 0.01) {
  block <- p / n_blocks
  d     <- seq_len(block)
  ar1   <- rho ^ abs(outer(d, d, "-"))
  R     <- as.matrix(bdiag(replicate(n_blocks, ar1, simplify = FALSE)))
  (1 - ridge) * R + ridge * diag(p)
}

## ── 2. Mismatched reference LD ───────────────────────────────────────────────
make_mismatch <- function(R_true, nref, seed = 99) {
  p <- nrow(R_true)
  set.seed(seed)
  noise <- matrix(rnorm(p * p, sd = 1 / sqrt(nref)), p, p)
  noise <- (noise + t(noise)) / 2
  R_m   <- R_true + noise
  R_m   <- (R_m + t(R_m)) / 2
  ev    <- eigen(R_m, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < 0.01)
    R_m <- R_m + (0.01 - min(ev)) * diag(p)
  d   <- sqrt(diag(R_m))
  R_m <- R_m / outer(d, d)
  0.99 * R_m + 0.01 * diag(p)
}

## ── 3. Anchor weights ────────────────────────────────────────────────────────
make_weights <- function(causal_idx, p, pve, s0,
                         N = 50000, rg = 0.6, R = NULL,
                         degraded = FALSE) {
  a <- runif(p)
  if (!degraded)
    for (j in causal_idx)
      a[j] <- min(a[j] + runif(1, 2.5, 3.5), 4.0)
  b2 <- numeric(p)
  if (length(causal_idx) > 0 && !degraded) {
    bp           <- numeric(p)
    bp[causal_idx] <- rnorm(s0, 0, sqrt(pve * N / s0))
    b2             <- rg * bp
  }
  z2 <- if (!is.null(R)) as.numeric(R %*% b2 + t(chol(R)) %*% rnorm(p)) else numeric(p)
  exp(0.8 * a) * (1 + abs(z2))
}

## ── 4. Spearman helper ───────────────────────────────────────────────────────
spearman_r <- function(x, y) {
  n <- length(x)
  1 - 6 * sum((rank(x) - rank(y))^2) / (n * (n^2 - 1))
}

## ── 5. V-SoRE ────────────────────────────────────────────────────────────────
## n_eff: effective sample size passed to susie_rss as `n`
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
  list(fit = fit, pip = gamma, hat_tau2 = hat_tau2)
}

## ── 6. Evaluation ────────────────────────────────────────────────────────────
evaluate <- function(causal_idx, cs_list) {
  cs_list <- cs_list[lengths(cs_list) > 0]   # drop empty credible sets
  if (length(cs_list) == 0)
    return(c(coverage = 0.0, cs_size = 0.0))
  cov <- mean(sapply(causal_idx,
                     function(c) any(sapply(cs_list, function(cs) c %in% cs))))
  sz  <- mean(sapply(cs_list, length))
  c(coverage = as.numeric(cov), cs_size = as.numeric(sz))
}

## ── 7. FINEMAP (univariate Wakefield ABF) ────────────────────────────────────
finemap_simple <- function(z, p, W) {
  log_bf <- -0.5 * log(1 + W) + 0.5 * W * z^2 / (1 + W)
  pip    <- exp(log_bf - log(p))
  pip    <- pip / sum(pip)
  ord    <- order(-pip)
  cum    <- cumsum(pip[ord])
  idx    <- which(cum >= 0.95)
  n_cs   <- if (length(idx) > 0) idx[1] else p   # take all if threshold never reached
  list(pip = pip, cs = list(ord[seq_len(n_cs)]))
}

## ── 8. Single replicate ──────────────────────────────────────────────────────
run_rep <- function(R_true, R_inf, s0, pve, N = 50000, p = 500,
                    degraded = FALSE) {
  causal    <- sort(sample(seq_len(p / 5), s0))
  b         <- numeric(p)
  b[causal] <- rnorm(s0, 0, sqrt(pve * N / s0))
  z         <- as.numeric(R_true %*% b + t(chol(R_true)) %*% rnorm(p))
  W         <- max(1.0, (max(abs(z)) / 2)^2)
  
  omega0 <- make_weights(causal, p, pve, s0, N = N, R = R_true, degraded = degraded)
  
  ## SuSiE-RSS
  fit_s  <- susie_rss(z, R_inf, n = N, L = 10)
  cs_s   <- susie_get_cs(fit_s)$cs
  ev_s   <- evaluate(causal, cs_s)
  
  ## V-SoRE
  res_v  <- vsore(z, R_inf, omega0, n_eff = N, L = 10, K0 = min(20, 2 * s0 + 5))
  cs_v   <- susie_get_cs(res_v$fit)$cs
  ev_v   <- evaluate(causal, cs_v)
  
  ## FINEMAP
  fm     <- finemap_simple(z, p, W)
  ev_f   <- evaluate(causal, fm$cs)
  
  c(cov_s = ev_s[["coverage"]],  sz_s = ev_s[["cs_size"]],
    cov_v = ev_v[["coverage"]],  sz_v = ev_v[["cs_size"]],
    tau2  = res_v$hat_tau2,
    cov_f = ev_f[["coverage"]],  sz_f = ev_f[["cs_size"]])
}

## ── 9. Configurations ────────────────────────────────────────────────────────
CONFIGS <- data.frame(
  label    = c("C1","C2","C3","C4","C5","C6","C7","C8"),
  s0       = c(1, 2, 3, 1, 2, 3, 3, 3),
  pve      = c(.005, .005, .005, .02, .02, .02, .005, .005),
  degraded = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE),
  ld_mis   = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
  stringsAsFactors = FALSE
)

## ── 10. Build LD matrices ─────────────────────────────────────────────────────
cat("Building LD matrices...\n")
R_true <- build_ld(p = P)
R_10k  <- make_mismatch(R_true, nref = 10000, seed = 99)
R_500  <- make_mismatch(R_true, nref = 500,   seed = 100)

## ── 11. Main loop (two reference-panel sizes) ─────────────────────────────────
run_panel <- function(panel_label, R_mis) {
  cat(sprintf("\nPanel N_ref = %s\n", panel_label))
  rows <- vector("list", nrow(CONFIGS))
  
  for (i in seq_len(nrow(CONFIGS))) {
    cfg   <- CONFIGS[i, ]
    R_inf <- if (cfg$ld_mis) R_mis else R_true
    set.seed(42 + i)
    
    mat <- replicate(N_REPS, {
      run_rep(R_true, R_inf, cfg$s0, cfg$pve,
              N = N_GWAS, p = P, degraded = cfg$degraded)
    })
    
    mn <- function(k) mean(mat[k, ])
    se <- function(k) sd(mat[k, ]) / sqrt(N_REPS)
    
    rows[[i]] <- data.frame(
      Panel     = panel_label,
      Config    = cfg$label,
      s0        = cfg$s0,
      PVE       = sprintf("%.1f%%", cfg$pve * 100),
      cov_fm    = 100 * mn("cov_f"),
      sz_fm     = mn("sz_f"),
      cov_susie = 100 * mn("cov_s"),
      sz_susie  = mn("sz_s"),
      cov_vsore = 100 * mn("cov_v"),
      sz_vsore  = mn("sz_v"),
      hat_tau2  = mn("tau2"),
      se_cov_s  = 100 * se("cov_s"),
      se_cov_v  = 100 * se("cov_v"),
      se_sz_s   = se("sz_s"),
      se_sz_v   = se("sz_v"),
      se_tau2   = se("tau2")
    )
    
    cat(sprintf("  %s  SuSiE %3.0f%%/%6.1f  VSoRE %3.0f%%/%6.1f  tau2=%.2f\n",
                cfg$label,
                rows[[i]]$cov_susie, rows[[i]]$sz_susie,
                rows[[i]]$cov_vsore, rows[[i]]$sz_vsore,
                rows[[i]]$hat_tau2))
  }
  bind_rows(rows)
}

df10  <- run_panel("10k",  R_10k)
df500 <- run_panel("500",  R_500)
df    <- bind_rows(df10, df500)

write.csv(df, "simulation_results.csv", row.names = FALSE)
cat("\nSaved simulation_results.csv\n")

## ── 12. Figures ───────────────────────────────────────────────────────────────
cfgs <- CONFIGS$label
cols <- c(FINEMAP = "#4e79a7", `SuSiE-RSS` = "#f28e2b", `V-SoRE` = "#59a14f")

panel_labels <- c("10k" = "N[ref]==10000", "500" = "N[ref]==500")

long_cov <- bind_rows(
  df |> transmute(Panel, Config, Method = "FINEMAP",   val = cov_fm,    se = 0),
  df |> transmute(Panel, Config, Method = "SuSiE-RSS", val = cov_susie, se = se_cov_s),
  df |> transmute(Panel, Config, Method = "V-SoRE",    val = cov_vsore, se = se_cov_v)
) |>
  mutate(Method = factor(Method, levels = c("FINEMAP","SuSiE-RSS","V-SoRE")),
         Config = factor(Config, levels = cfgs),
         Panel  = factor(Panel,  levels = c("10k","500"), labels = panel_labels))

long_sz <- bind_rows(
  df |> transmute(Panel, Config, Method = "FINEMAP",   val = sz_fm,    se = 0),
  df |> transmute(Panel, Config, Method = "SuSiE-RSS", val = sz_susie, se = se_sz_s),
  df |> transmute(Panel, Config, Method = "V-SoRE",    val = sz_vsore, se = se_sz_v)
) |>
  mutate(Method = factor(Method, levels = c("FINEMAP","SuSiE-RSS","V-SoRE")),
         Config = factor(Config, levels = cfgs),
         Panel  = factor(Panel,  levels = c("10k","500"), labels = panel_labels))

long_tau <- df |>
  mutate(Config = factor(Config, levels = cfgs),
         Panel  = factor(Panel,  levels = c("10k","500"), labels = panel_labels))

theme_sore <- theme_classic(base_size = 11) +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size = 11))

## Panel A
p_cov <- ggplot(long_cov, aes(Config, val, fill = Method)) +
  geom_col(position = position_dodge(0.8), width = 0.75, alpha = 0.85) +
  geom_errorbar(aes(ymin = val - se, ymax = val + se),
                position = position_dodge(0.8), width = 0.3, linewidth = 0.7) +
  geom_hline(yintercept = 95, linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~Panel, labeller = label_parsed) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(limits = c(0, 115), breaks = seq(0, 100, 25)) +
  labs(x = NULL, y = "Per-signal coverage (%)", fill = NULL,
       title = "Panel A: Per-signal 95% credible-set coverage") +
  theme_sore
ggsave("fig_coverage.pdf", p_cov, width = 10, height = 4.5)
cat("Saved fig_coverage.pdf\n")

## Panel B
p_sz <- ggplot(long_sz, aes(Config, val, fill = Method)) +
  geom_col(position = position_dodge(0.8), width = 0.75, alpha = 0.85) +
  geom_errorbar(aes(ymin = pmax(val - se, 0.1), ymax = val + se),
                position = position_dodge(0.8), width = 0.3, linewidth = 0.7) +
  facet_wrap(~Panel, labeller = label_parsed) +
  scale_fill_manual(values = cols) +
  scale_y_log10() +
  labs(x = NULL, y = "Mean CS size (log scale)", fill = NULL,
       title = "Panel B: Mean credible-set size") +
  theme_sore
ggsave("fig_cssize.pdf", p_sz, width = 10, height = 4.5)
cat("Saved fig_cssize.pdf\n")

## Panel C
p_tau <- ggplot(long_tau, aes(Config, hat_tau2)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0,   ymax=0.2,
           fill="steelblue", alpha=0.10) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0.2, ymax=1.0,
           fill="goldenrod", alpha=0.10) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=1.0, ymax=2.5,
           fill="firebrick", alpha=0.10) +
  geom_hline(yintercept=0.2, linetype="dashed", colour="steelblue", linewidth=0.9) +
  geom_hline(yintercept=1.0, linetype="dashed", colour="firebrick", linewidth=0.9) +
  geom_errorbar(aes(ymin=hat_tau2-se_tau2, ymax=hat_tau2+se_tau2),
                width=0.25, linewidth=0.9, colour="#59a14f") +
  geom_point(size=3.5, colour="#59a14f") +
  facet_wrap(~Panel, labeller=label_parsed) +
  scale_y_continuous(limits=c(0, 2.0)) +
  labs(x=NULL, y=expression(hat(tau)^2),
       title=expression("Panel C: V-SoRE annotation diagnostic " * hat(tau)^2)) +
  theme_classic(base_size=11) +
  theme(strip.background=element_blank(), strip.text=element_text(size=11))
ggsave("fig_tau.pdf", p_tau, width=9, height=4)
cat("Saved fig_tau.pdf\n")

## ── 13. LaTeX table ───────────────────────────────────────────────────────────
fmt_row <- function(r)
  sprintf("%s & %d & %s & %3.0f & %5.1f & %3.0f & %6.1f & %3.0f & %6.1f & %.2f \\\\",
          r$Config, r$s0, r$PVE,
          r$cov_fm, r$sz_fm,
          r$cov_susie, r$sz_susie,
          r$cov_vsore, r$sz_vsore,
          r$hat_tau2)

lines <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Per-signal 95\\% credible-set coverage (\\%) and mean CS size",
  "  across eight configurations (30 replicates).",
  "  Panel~A: $N_{\\mathrm{ref}}=10{,}000$.  Panel~B: $N_{\\mathrm{ref}}=500$.}",
  "\\label{tab:sim}",
  "\\small",
  "\\begin{tabular}{lcc rr rr rr r}",
  "\\toprule",
  " & & &\\multicolumn{2}{c}{FINEMAP}",
  "  &\\multicolumn{2}{c}{SuSiE-RSS}",
  "  &\\multicolumn{2}{c}{V-SoRE} & \\\\",
  "\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}",
  "Config & $s_0$ & PVE & Cov & Size & Cov & Size & Cov & Size & $\\hat\\tau^2$\\\\",
  "\\midrule",
  "\\multicolumn{10}{l}{\\textit{Panel A: $N_{\\mathrm{ref}} = 10{,}000$}}\\\\[2pt]",
  sapply(seq_len(nrow(df10)), function(i) fmt_row(df10[i, ])),
  "\\midrule",
  "\\multicolumn{10}{l}{\\textit{Panel B: $N_{\\mathrm{ref}} = 500$}}\\\\[2pt]",
  sapply(seq_len(nrow(df500)), function(i) fmt_row(df500[i, ])),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)
writeLines(lines, "latex_table.txt")
cat("Saved latex_table.txt\n\nAll done.\n")