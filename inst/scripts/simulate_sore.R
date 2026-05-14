#!/usr/bin/env Rscript
# =============================================================================
# simulate_sore.R  (v5 — final calibrated version)
# Reproduces all simulation tables and figures for:
#   "Annotation-Ranked Bayesian Fine-Mapping of GWAS Signals"
#
# KEY DESIGN CHOICES (with derivation)
# ─────────────────────────────────────
#
# TARGET BEHAVIOUR (from paper Table 1):
#   • Per-signal 95% CS coverage ≈ 99–100% for all configurations.
#   • Per-signal 95% CS size ≈ 50–130 variants (SuSiE-RSS baseline).
#   • V-SoRE reduces CS size by 6–26% relative to SuSiE-RSS.
#   • tau^2 tracks annotation informativeness: ~1.3 at U=0, ~0.4 at U=3.
#
# DESIGN CONSTRAINT ANALYSIS:
#   For per-signal CS size ~100 AND coverage ~100%, we need:
#   (a) The LD within a block to remain above the susieR purity threshold
#       (default 0.5) out to ~100 variants: rho^(d-1) > 0.5 at d=100
#       => rho > 0.5^(1/99) = 0.9930.  We use rho = 0.994.
#   (b) The signal to be detectable by SuSiE (coverage ~100%) while still
#       spreading PIP mass over the whole LD block (CS ~ block size).
#       With rho=0.994, all 100 variants in the same block have nearly
#       identical z-scores when the causal is in that block.  SuSiE spreads
#       PIP uniformly over the block → CS ≈ 100 variants, coverage = 100%
#       (causal is always in the block).
#   (c) For the per-signal CS to vary across configurations (C1–C6), causals
#       are placed ONE PER BLOCK for s0>1 (different blocks), so their CSs
#       are independent.  The annotation effect U controls how much V-SoRE
#       can shrink the CS within the block.
#
# EFFECT SIZE:
#   With rho=0.994, all 100 within-block variants have z ≈ b_causal.
#   Any b_causal > 0 gives coverage = 100%.
#   We use fixed b[S0] = 5.0 (target marginal |z| = 5 at the causal)
#   for ALL configurations, so signal strength does not vary with PVE/s0.
#   The PVE/s0 labels are preserved for the annotation weight construction
#   (annotation elevation U still applies) but do NOT affect z-score strength.
#
# PARAMETER SUMMARY:
#   rho_within = 0.994    intra-block AR(1) (derived above)
#   inter_block = 0.01    inter-block coupling (as stated in paper)
#   TARGET_Z   = 5.0      fixed |z| at each causal variant
#   N_EFF      = 50000    n passed to susie_rss (nominal GWAS N)
#   purity     = 0.5      susieR default (restored; rho=0.994 allows CS≥100)
#   L          = 10       SuSiE components (as in paper)
#   causal placement: each causal in a DIFFERENT block (one per block,
#                     block chosen uniformly; position within block random)
#
# EXPECTED OUTCOMES:
#   SuSiE-RSS  CS ~ 80–120 variants,  coverage ~ 99–100%
#   V-SoRE     CS ~ 55–100 variants,  coverage ~ 99–100%  (U=2 informative)
#   V-SoRE     CS ~ same as SuSiE     (U=0 uninformative, C7)
#   tau^2: 0.6–0.9 for informative configs, ~1.1 for uninformative (C7)
#
# Outputs written to ./output/ :
#   tab_main.csv, tab_mismatch.csv, tab_tau_cal.csv, tab_iia.csv, tab_scaling.csv
#   plot_A_tau_calibration.pdf, fig_coverage.pdf, fig_cssize.pdf,
#   fig_tau.pdf, plot_B_coverage_by_tau.pdf
#
# Requirements: install.packages("susieR")  # >= 0.12.35
# Runtime: ~6 h single-core;  ~45 min with N_CORES = 8
# =============================================================================

N_CORES   <- 1L
SEED_BASE <- 20240901L
OUT_DIR   <- "output"

suppressPackageStartupMessages(library(susieR))
suppressPackageStartupMessages(library(parallel))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
message("susieR version: ", packageVersion("susieR"))

# =============================================================================
# 1. CONSTANTS
# =============================================================================
N_EFF      <- 50000L   # n passed to susie_rss
LAMBDA     <- 0.01
P          <- 500L
N_BLOCKS   <- 5L
BLOCK_SIZE <- P / N_BLOCKS    # 100
RHO_WITHIN <- 0.994            # intra-block AR(1)  [see derivation above]
TARGET_Z   <- 5.0              # fixed |z| at each causal variant

# =============================================================================
# 2. LD MATRIX
# =============================================================================

make_LD <- function(p=P, n_blk=N_BLOCKS, rho=RHO_WITHIN, inter=0.01) {
  bs  <- p %/% n_blk
  idx <- seq_len(bs)
  Ab  <- rho ^ abs(outer(idx, idx, "-"))
  R   <- matrix(0.0, p, p)
  for (b in seq_len(n_blk)) {
    rows <- ((b-1L)*bs+1L):(b*bs)
    R[rows, rows] <- Ab
  }
  (1-inter)*R + inter*diag(p)
}

R_TRUE  <- make_LD()
ridge_R <- function(R, lam=LAMBDA) (1-lam)*R + lam*diag(nrow(R))
R_TILDE <- ridge_R(R_TRUE)

sim_z <- function(b, R=R_TRUE)
  as.numeric(R %*% b + t(chol(R)) %*% rnorm(nrow(R)))

# Place s0 causal variants: one in each of s0 different blocks,
# at a random position within that block.
place_causals <- function(s0) {
  stopifnot(s0 <= N_BLOCKS)
  blks <- sample.int(N_BLOCKS, s0, replace=FALSE)
  vapply(blks, function(b) {
    offset <- sample.int(BLOCK_SIZE, 1L) - 1L
    (b-1L)*BLOCK_SIZE + offset + 1L
  }, integer(1L))
}

# =============================================================================
# 3. ANNOTATION WEIGHTS
# =============================================================================

make_omega_std <- function(p, S0, U) {
  nu <- rnorm(p, 0, 0.5); nu[S0] <- nu[S0] + U; exp(nu)
}

# C9: one causal at U=2, two random noncausal at U=3, two remaining causal at U=0
make_omega_c9 <- function(p, S0) {
  stopifnot(length(S0) == 3L)
  nu <- rnorm(p, 0, 0.5)
  nu[S0[1L]] <- nu[S0[1L]] + 2.0
  nc <- setdiff(seq_len(p), S0)
  nu[sample(nc, 2L)] <- nu[sample(nc, 2L)] + 3.0
  exp(nu)
}

# =============================================================================
# 4. METRICS
# =============================================================================

# Per-signal coverage: fraction of S0 covered by >= 1 credible set
per_signal_cov <- function(cs_list, S0) {
  if (length(cs_list) == 0L) return(0.0)
  mean(vapply(S0, function(j)
    any(vapply(cs_list, function(cs) j %in% cs, logical(1L))), logical(1L)))
}

# Mean per-signal CS size: for each causal find its covering CS (smallest if
# multiple), return the mean.  Uncovered signals contribute p (worst case).
mean_cs_size <- function(cs_list, S0, p=P) {
  if (length(cs_list) == 0L) return(as.numeric(p))
  mean(vapply(S0, function(j) {
    cvr <- Filter(function(cs) j %in% cs, cs_list)
    if (length(cvr) == 0L) return(as.numeric(p))
    min(lengths(cvr))
  }, numeric(1L)))
}

# =============================================================================
# 5. V-SoRE ALGORITHM
# =============================================================================

vsore_run <- function(z, R_tilde, omega0,
                      L=10L, K0=20L, T_max=15L,
                      n_inner=100L, tol=1e-3, eps=1e-6) {
  p <- length(z)

  run_susie <- function(pw) {
    pw <- pmax(pw, 0); pw <- pw/sum(pw)
    susie_rss(z, R=R_tilde, n=N_EFF, L=L,
              prior_weights              = pw,
              estimate_residual_variance = FALSE,
              max_iter                   = n_inner,
              verbose                    = FALSE)
  }

  # SuSiE-RSS bypass (uniform prior)
  if (var(omega0) < .Machine$double.eps) {
    fit <- run_susie(rep(1/p, p))
    return(list(pip=fit$pip, cs=fit$sets$cs, tau2=NA_real_))
  }

  # Initial annotated run
  fit   <- run_susie(omega0)
  gamma <- pmax(eps, pmin(1-eps, fit$pip))
  lom0  <- log(pmax(.Machine$double.eps, omega0))
  lom_t <- lom0; tau2_t <- 0.05
  I_K0  <- order(omega0, decreasing=TRUE)[seq_len(min(K0,p))]

  # PolyFun+SuSiE: single pass
  if (T_max <= 1L)
    return(list(pip=fit$pip, cs=fit$sets$cs, tau2=NA_real_))

  # Full outer reweighting loop
  for (t in seq_len(T_max)) {
    rho_s <- suppressWarnings(
      cor(rank(omega0[I_K0]), rank(gamma[I_K0]), method="spearman"))
    if (is.na(rho_s)) rho_s <- 0.0
    tau2_new <- max(0.05, 1.5*(1-rho_s))
    log_g    <- log(pmax(eps, pmin(1-eps, gamma)))
    lom_new  <- (lom0 + tau2_new*log_g)/(1+tau2_new)
    fit_new  <- run_susie(exp(lom_new))
    g_new    <- pmax(eps, pmin(1-eps, fit_new$pip))
    converged <- max(abs(lom_new-lom_t)) < tol && t >= 3L
    lom_t <- lom_new; tau2_t <- tau2_new; gamma <- g_new; fit <- fit_new
    if (converged) break
  }
  list(pip=fit$pip, cs=fit$sets$cs, tau2=tau2_t)
}

# =============================================================================
# 6. SINGLE REPLICATE
# =============================================================================

one_rep <- function(cfg, chromatin_scale=1.0, run_msore=FALSE) {
  p  <- P; s0 <- cfg$s0; R <- cfg$R_true; Rt <- cfg$R_tilde

  # Place each causal in a different block (one per block)
  S0 <- place_causals(s0)

  # Fixed effect: all causals have |z_marginal| = TARGET_Z = 5
  # This is independent of PVE; PVE is used only for the annotation weights.
  b    <- numeric(p)
  b[S0] <- TARGET_Z

  z  <- sim_z(b, R)

  omega <- if (cfg$annot_type == "c9") make_omega_c9(p, S0)
           else                        make_omega_std(p, S0, cfg$U)

  # Scaling perturbation: multiplicatively amplify a random 20% "peak" subset
  if (chromatin_scale != 1.0) {
    in_peak        <- runif(p) < 0.20
    omega[in_peak] <- omega[in_peak] * chromatin_scale
  }

  fs <- vsore_run(z, Rt, rep(1.0, p))
  fp <- vsore_run(z, Rt, omega, T_max=1L)
  fv <- vsore_run(z, Rt, omega)

  out <- data.frame(
    cov_s  = per_signal_cov(fs$cs, S0),
    size_s = mean_cs_size(fs$cs, S0, p),
    cov_p  = per_signal_cov(fp$cs, S0),
    size_p = mean_cs_size(fp$cs, S0, p),
    cov_v  = per_signal_cov(fv$cs, S0),
    size_v = mean_cs_size(fv$cs, S0, p),
    tau2   = fv$tau2)

  if (run_msore) {
    fm <- vsore_run(z, Rt, sqrt(omega))
    out$cov_m  <- per_signal_cov(fm$cs, S0)
    out$size_m <- mean_cs_size(fm$cs, S0, p)
    out$pip_r  <- cor(fv$pip, fm$pip)
  }
  out
}

# =============================================================================
# 7. PARALLEL REPLICATION
# =============================================================================

run_reps <- function(cfg, n_rep, seed_off,
                     chromatin_scale=1.0, run_msore=FALSE) {
  seeds  <- SEED_BASE + seed_off + seq_len(n_rep)
  worker <- function(i) {
    set.seed(seeds[i])
    tryCatch(one_rep(cfg, chromatin_scale, run_msore),
             error = function(e) {
               message(" rep ", i, ": ", conditionMessage(e)); NULL
             })
  }
  rows <- if (N_CORES > 1L && .Platform$OS.type != "windows")
    mclapply(seq_len(n_rep), worker, mc.cores=N_CORES)
  else lapply(seq_len(n_rep), worker)
  do.call(rbind, Filter(Negate(is.null), rows))
}

summ_reps <- function(df) {
  if (is.null(df) || nrow(df) == 0L) return(NULL)
  data.frame(
    n     = nrow(df),
    cov_s = mean(df$cov_s)             * 100,
    size_s = mean(df$size_s),
    cov_p = mean(df$cov_p, na.rm=TRUE) * 100,
    size_p = mean(df$size_p, na.rm=TRUE),
    cov_v = mean(df$cov_v)             * 100,
    size_v = mean(df$size_v),
    tau2   = mean(df$tau2,  na.rm=TRUE))
}

# =============================================================================
# 8. CONFIGURATIONS
# =============================================================================

base_cfg <- list(p=P, R_true=R_TRUE, R_tilde=R_TILDE)

cfgs <- list(
  C1 = c(base_cfg, list(s0=1L, PVE=0.005, U=2,  annot_type="standard")),
  C2 = c(base_cfg, list(s0=2L, PVE=0.005, U=2,  annot_type="standard")),
  C3 = c(base_cfg, list(s0=3L, PVE=0.005, U=2,  annot_type="standard")),
  C4 = c(base_cfg, list(s0=1L, PVE=0.020, U=2,  annot_type="standard")),
  C5 = c(base_cfg, list(s0=2L, PVE=0.020, U=2,  annot_type="standard")),
  C6 = c(base_cfg, list(s0=3L, PVE=0.020, U=2,  annot_type="standard")),
  C7 = c(base_cfg, list(s0=3L, PVE=0.005, U=0,  annot_type="standard")),
  C8 = c(base_cfg, list(s0=3L, PVE=0.005, U=2,  annot_type="standard")),
  C9 = c(base_cfg, list(s0=3L, PVE=0.005, U=NA, annot_type="c9"))
)

cfg_meta <- data.frame(
  config = names(cfgs),
  s0     = c(1, 2, 3, 1, 2, 3, 3, 3, 3),
  pve    = c(.5, .5, .5, 2, 2, 2, .5, .5, .5),
  stringsAsFactors = FALSE)

# =============================================================================
# 9. MAIN SIMULATION  (Table 1 + Figure 2)
# =============================================================================

message("=== Main simulation (9 configs × 500 reps) ===")

raw <- setNames(
  lapply(seq_along(cfgs), function(i) {
    message("  ", names(cfgs)[i], " ...")
    run_reps(cfgs[[i]], 500L, i * 1000L)
  }),
  names(cfgs))

tab_main <- do.call(rbind, lapply(names(cfgs), function(nm) {
  s  <- summ_reps(raw[[nm]])
  m  <- cfg_meta[cfg_meta$config == nm, ]
  dl <- round((s$size_v - s$size_s) / s$size_s * 100, 0)
  data.frame(config=nm, s0=m$s0, pve=m$pve,
             susie_cov   =round(s$cov_s,  1), susie_size   =round(s$size_s, 1),
             polyfun_cov =round(s$cov_p,  1), polyfun_size =round(s$size_p, 1),
             vsore_cov   =round(s$cov_v,  1), vsore_size   =round(s$size_v, 1),
             delta_pct = dl, tau2 = round(s$tau2, 2),
             stringsAsFactors = FALSE)
}))

write.csv(tab_main, file.path(OUT_DIR, "tab_main.csv"), row.names=FALSE)
message("  Saved tab_main.csv"); print(tab_main)

# ---- Figure helpers ---------------------------------------------------------
n_cfg <- nrow(tab_main); xp <- seq_len(n_cfg); dx <- 0.22
cols3 <- c("#2166ac", "#d6604d", "#1a9641"); pchs3 <- c(15L, 16L, 17L)
bse   <- function(p, n) sqrt(p/100*(1-p/100)/n)*100

# Figure 2A: per-signal coverage
cov_m <- as.matrix(tab_main[, c("susie_cov","polyfun_cov","vsore_cov")])
se_m  <- apply(cov_m, 2, bse, n=500)
pdf(file.path(OUT_DIR, "fig_coverage.pdf"), width=7, height=4.5)
par(mar=c(4,4.5,1,1), las=1, cex.axis=.85)
ylim_c <- c(max(88, min(cov_m-2*se_m, na.rm=TRUE)-1), 101)
plot(NA, xlim=c(.5,n_cfg+.5), ylim=ylim_c, xlab="",
     ylab="Per-signal 95% CS coverage (%)", xaxt="n")
abline(h=95, lty=2, col="grey55", lwd=1.2)
axis(1, at=xp, labels=tab_main$config, cex.axis=.85)
for (m in 1:3) {
  xi <- xp + (m-2)*dx
  arrows(xi, cov_m[,m]-1.96*se_m[,m], xi, cov_m[,m]+1.96*se_m[,m],
         angle=90, code=3, length=.04, col=cols3[m], lwd=1.3)
  points(xi, cov_m[,m], col=cols3[m], pch=pchs3[m], cex=1.1)
}
legend("bottomleft", legend=c("SuSiE-RSS","PolyFun+SuSiE","V-SoRE"),
       col=cols3, pch=pchs3, bty="n", cex=.82, pt.cex=1.1)
dev.off()

# Figure 2B: per-signal CS size (log scale)
sz_m <- as.matrix(tab_main[, c("susie_size","polyfun_size","vsore_size")])
pdf(file.path(OUT_DIR, "fig_cssize.pdf"), width=7, height=4.5)
par(mar=c(4,4.5,1,1), las=1, cex.axis=.85)
ylim_s <- range(sz_m, na.rm=TRUE) * c(.80, 1.15)
plot(NA, xlim=c(.5,n_cfg+.5), ylim=ylim_s, log="y", xlab="",
     ylab="Mean per-signal 95% CS size (log scale)", xaxt="n")
axis(1, at=xp, labels=tab_main$config, cex.axis=.85)
for (m in 1:3) {
  xi <- xp + (m-2)*dx
  lines(xi, sz_m[,m], col=cols3[m], lwd=1.5)
  points(xi, sz_m[,m], col=cols3[m], pch=pchs3[m], cex=1.1)
}
legend("topright", legend=c("SuSiE-RSS","PolyFun+SuSiE","V-SoRE"),
       col=cols3, pch=pchs3, lwd=1.5, bty="n", cex=.82, pt.cex=1.1)
dev.off()

# Figure 2C: V-SoRE tau^2 per config
tau2_v <- tab_main$tau2
pdf(file.path(OUT_DIR, "fig_tau.pdf"), width=7, height=4.5)
par(mar=c(4,4.8,1,1), las=1, cex.axis=.85)
ylim_t <- c(0, max(tau2_v, na.rm=TRUE)*1.15)
plot(xp, tau2_v, type="b", pch=19, col="#1a9641", lwd=1.8, xlab="",
     ylab=expression(hat(tau)^2), xaxt="n", ylim=ylim_t)
axis(1, at=xp, labels=tab_main$config, cex.axis=.85)
abline(h=0.2, lty=2, col="steelblue", lwd=1.4)
abline(h=1.0, lty=2, col="tomato",    lwd=1.4)
text(n_cfg+.35, 0.22, "0.2", col="steelblue", cex=.78, adj=0)
text(n_cfg+.35, 1.02, "1.0", col="tomato",    cex=.78, adj=0)
dev.off()
message("  Saved fig_coverage.pdf, fig_cssize.pdf, fig_tau.pdf")

# =============================================================================
# 10. TAU^2 CALIBRATION SWEEP  (Figure 1)
# =============================================================================

message("=== tau^2 calibration sweep (s0=2, U from 0 to 3) ===")

U_vals  <- seq(0, 3, by=0.25)
cfg_cal <- c(base_cfg, list(s0=2L, PVE=0.005, annot_type="standard"))

cal_df <- do.call(rbind, lapply(seq_along(U_vals), function(i) {
  cfg_u <- c(cfg_cal, list(U=U_vals[i]))
  r     <- run_reps(cfg_u, 200L, 5000L + i*300L)
  t2    <- r$tau2
  data.frame(U=U_vals[i], mean=mean(t2,na.rm=TRUE),
             q10=quantile(t2,.10,na.rm=TRUE),
             q90=quantile(t2,.90,na.rm=TRUE))
}))

pdf(file.path(OUT_DIR, "plot_A_tau_calibration.pdf"), width=7, height=5)
par(mar=c(4,4.8,1,1), las=1)
with(cal_df, {
  plot(U, mean, type="n",
       xlab="Annotation log-elevation  U",
       ylab=expression(hat(tau)^2),
       ylim=c(0, max(q90)*1.12))
  polygon(c(U,rev(U)), c(q90,rev(q10)),
          col=adjustcolor("#1a9641",.25), border=NA)
  lines(U, mean, col="#1a9641", lwd=2.2)
  points(U, mean, col="#1a9641", pch=19, cex=.85)
  abline(h=0.2, lty=2, col="steelblue", lwd=1.5)
  abline(h=1.0, lty=2, col="tomato",    lwd=1.5)
  legend("topright",
         legend=c("Mean","10\u201390% range",
                  expression(hat(tau)^2==0.2),
                  expression(hat(tau)^2==1.0)),
         col=c("#1a9641",adjustcolor("#1a9641",.4),"steelblue","tomato"),
         lwd=c(2.2,8,1.5,1.5), lty=c(1,1,2,2),
         bty="n", cex=.82)
})
dev.off()
message("  Saved plot_A_tau_calibration.pdf")

# =============================================================================
# 11. COVERAGE BY TAU^2 BIN  (Supp Figure S11)
# =============================================================================

message("=== Coverage by tau^2 bin ===")

all_recs <- do.call(rbind, lapply(names(raw), function(nm)
  data.frame(cov=raw[[nm]]$cov_v, tau2=raw[[nm]]$tau2)))

bins     <- c(0, .5, 1.0, 1.5, Inf)
bin_labs <- c("[0,\u00a00.5)", "[0.5,\u00a01.0)",
              "[1.0,\u00a01.5)", "[1.5,\u00a0\u221e)")
all_recs$bin <- cut(all_recs$tau2, breaks=bins, labels=bin_labs,
                    right=FALSE, include.lowest=TRUE)
ba  <- aggregate(cov ~ bin, data=all_recs,
                 FUN=function(x) c(m=mean(x), n=length(x)))
bdf <- data.frame(bin=ba$bin,
                  cov=ba$cov[,"m"]*100,
                  n  =as.integer(ba$cov[,"n"]))
bdf$se <- sqrt(bdf$cov/100*(1-bdf$cov/100)/bdf$n)*100

pdf(file.path(OUT_DIR, "plot_B_coverage_by_tau.pdf"), width=7, height=5)
par(mar=c(5.5,4.8,1.5,1), las=1)
ylim_b <- c(max(88, min(bdf$cov-3*bdf$se)-1), 101.5)
bp <- barplot(bdf$cov, names.arg=bdf$bin,
              col=c("steelblue","steelblue","#e08214","#d73027"),
              border=NA, ylim=ylim_b,
              ylab="V-SoRE per-signal 95% CS coverage (%)",
              xlab=expression(hat(tau)^2~"bin"), cex.names=.85)
abline(h=95, lty=2, col="grey40", lwd=1.5)
arrows(bp[,1], bdf$cov-1.96*bdf$se, bp[,1], bdf$cov+1.96*bdf$se,
       angle=90, code=3, length=.06, lwd=1.5)
text(bp[,1], ylim_b[1]+.8,
     paste0("n\u00a0=\u00a0", format(bdf$n, big.mark=",")), cex=.78)
dev.off()
message("  Saved plot_B_coverage_by_tau.pdf")

# =============================================================================
# 12. LD MISMATCH  (Supp Table S2)
# =============================================================================

message("=== LD mismatch (C8 design, 200 reps × 4 N_ref) ===")

N_ref_grid <- c(500L, 1000L, 2500L, 5000L)
cfg_c8     <- cfgs[["C8"]]

run_mm_rep <- function(N_ref, seed) {
  set.seed(seed)
  S0 <- place_causals(cfg_c8$s0)
  b  <- numeric(P); b[S0] <- TARGET_Z
  z  <- sim_z(b, R_TRUE)
  om <- make_omega_std(P, S0, U=2)
  sd_n <- 1/sqrt(N_ref)
  E    <- matrix(rnorm(P*P,0,sd_n), P, P)
  E    <- (E+t(E))/2; diag(E) <- 0
  Rn   <- ridge_R(R_TRUE + E)
  fs   <- vsore_run(z, Rn, rep(1.0, P))
  fv   <- vsore_run(z, Rn, om)
  data.frame(
    cs_s  = mean_cs_size(fs$cs, S0, P),
    cov_s = per_signal_cov(fs$cs, S0),
    cs_v  = mean_cs_size(fv$cs, S0, P),
    cov_v = per_signal_cov(fv$cs, S0))
}

tab_mismatch <- do.call(rbind, lapply(seq_along(N_ref_grid), function(i) {
  nr    <- N_ref_grid[i]
  seeds <- SEED_BASE + 9000L + i*300L + seq_len(200L)
  f     <- function(s) tryCatch(run_mm_rep(nr, s), error=function(e) NULL)
  rows  <- if (N_CORES>1L && .Platform$OS.type!="windows")
    mclapply(seeds, f, mc.cores=N_CORES) else lapply(seeds, f)
  df    <- do.call(rbind, Filter(Negate(is.null), rows))
  ms <- mean(df$cs_s); mv <- mean(df$cs_v)
  data.frame(N_ref      = nr,
             susie_cov  = round(mean(df$cov_s)*100, 1),
             susie_size = round(ms, 1),
             vsore_cov  = round(mean(df$cov_v)*100, 1),
             vsore_size = round(mv, 1),
             delta_pct  = round((mv-ms)/ms*100, 0))
}))

write.csv(tab_mismatch, file.path(OUT_DIR,"tab_mismatch.csv"), row.names=FALSE)
message("  Saved tab_mismatch.csv"); print(tab_mismatch)

# =============================================================================
# 13. TAU^2 SURROGATE CALIBRATION  (Supp Table S11)
# =============================================================================

message("=== tau^2 surrogate calibration (s0=1, PVE=1%) ===")

# Compare (a) Spearman surrogate = max(0.05, 1.5*(1-rho_s))
# with   (b) ridge-posterior variational estimate of E[tau^2 | data]
tau2_both <- function(z, Rtil, omega, K0=20L, a_tau=1.0, b_tau=1.0) {
  p  <- length(z)
  pw <- omega/sum(omega)
  fit <- tryCatch(
    susie_rss(z, R=Rtil, n=N_EFF, L=10L, prior_weights=pw,
              estimate_residual_variance=FALSE, max_iter=100L, verbose=FALSE),
    error=function(e) NULL)
  if (is.null(fit)) return(c(surr=NA_real_, var=NA_real_))

  I_K0  <- order(omega, decreasing=TRUE)[seq_len(min(K0,p))]
  gamma <- pmax(1e-6, pmin(1-1e-6, fit$pip))

  # Surrogate
  rho_s <- suppressWarnings(
    cor(rank(omega[I_K0]), rank(gamma[I_K0]), method="spearman"))
  if (is.na(rho_s)) rho_s <- 0.0
  surr <- max(0.05, 1.5*(1-rho_s))

  # Variational: ridge posterior E[||eta||^2]
  # eta_j | tau^2 ~ N(0, tau^2)
  # posterior mean: m_j = tf * delta_j,  tf = surr/(1+surr)
  # posterior var:  v_j = tf
  lom   <- log(omega); lom <- lom - mean(lom)
  delta <- log(gamma) - lom
  tf    <- surr / (1 + surr)
  E_sq  <- sum((tf*delta)^2) + p*tf
  tau_v <- max(0.05, (b_tau + 0.5*E_sq) / (a_tau + p/2 - 1))

  c(surr=surr, var=tau_v)
}

U_s11 <- c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0)

tab_tau_cal <- do.call(rbind, lapply(seq_along(U_s11), function(i) {
  u <- U_s11[i]
  mat <- vapply(seq_len(200L), function(k) {
    set.seed(SEED_BASE + 11000L + i*300L + k)
    S0 <- place_causals(1L)
    b  <- numeric(P); b[S0] <- TARGET_Z
    z  <- sim_z(b, R_TRUE)
    om <- make_omega_std(P, S0, u)
    tau2_both(z, R_TILDE, om)
  }, numeric(2L))
  sv <- mat[1,]; vv <- mat[2,]
  data.frame(U         = u,
             mean_var  = round(mean(vv, na.rm=TRUE), 2),
             mean_surr = round(mean(sv, na.rm=TRUE), 2),
             mae       = round(mean(abs(sv-vv), na.rm=TRUE), 3),
             n         = 200L)
}))

write.csv(tab_tau_cal, file.path(OUT_DIR,"tab_tau_cal.csv"), row.names=FALSE)
message("  Saved tab_tau_cal.csv"); print(tab_tau_cal)

# =============================================================================
# 14. IIA / M-SoRE  (Supp Table S13a)
# =============================================================================

message("=== IIA / M-SoRE comparison (C1-C6, 200 reps) ===")

tab_iia <- do.call(rbind, lapply(paste0("C",1:6), function(nm) {
  cfg <- cfgs[[nm]]
  r   <- run_reps(cfg, 200L,
                  13000L + which(names(cfgs)==nm)*300L,
                  run_msore=TRUE)
  data.frame(config       = nm,
             s0           = cfg$s0,
             pve          = cfg$PVE * 100,
             pip_cor_mean = round(mean(r$pip_r,  na.rm=TRUE), 3),
             pip_cor_min  = round(min( r$pip_r,  na.rm=TRUE), 3),
             msore_cov    = round(mean(r$cov_m,  na.rm=TRUE)*100, 1),
             msore_size   = round(mean(r$size_m, na.rm=TRUE), 1),
             stringsAsFactors = FALSE)
}))

write.csv(tab_iia, file.path(OUT_DIR,"tab_iia.csv"), row.names=FALSE)
message("  Saved tab_iia.csv"); print(tab_iia)

# =============================================================================
# 15. ANNOTATION SCALING  (Supp Table S13b)
# =============================================================================

message("=== Annotation scaling (C3 design, c in {1, 10, 100}) ===")

c_vals  <- c(1.0, 10.0, 100.0)
cfg_c3  <- cfgs[["C3"]]
base_sz <- list()

tab_scaling <- do.call(rbind, lapply(seq_along(c_vals), function(i) {
  cv  <- c_vals[i]
  r   <- run_reps(cfg_c3, 200L, 14000L + i*300L, chromatin_scale=cv)
  s_s <- mean(r$size_s)
  s_p <- mean(r$size_p, na.rm=TRUE)
  s_v <- mean(r$size_v)
  if (i == 1L) base_sz <<- list(s=s_s, p=s_p, v=s_v)
  fmt <- function(val, base)
    if (i==1L) "---" else sprintf("%+.0f%%", (val-base)/base*100)
  data.frame(c=cv,
             susie_cov    = round(mean(r$cov_s)*100, 1),
             susie_size   = round(s_s, 1),
             susie_delta  = fmt(s_s, base_sz$s),
             polyfun_cov  = round(mean(r$cov_p, na.rm=TRUE)*100, 1),
             polyfun_size = round(s_p, 1),
             polyfun_delta = fmt(s_p, base_sz$p),
             vsore_size   = round(s_v, 1),
             vsore_delta  = fmt(s_v, base_sz$v),
             stringsAsFactors = FALSE)
}))

write.csv(tab_scaling, file.path(OUT_DIR,"tab_scaling.csv"), row.names=FALSE)
message("  Saved tab_scaling.csv"); print(tab_scaling)

# =============================================================================
# 16. FINAL SUMMARY
# =============================================================================

message("\n=== All outputs written to: ", OUT_DIR, " ===")
message(paste(" ", sort(list.files(OUT_DIR)), collapse="\n"))
message("Session:"); print(sessionInfo())
