## =============================================================================
## Main simulation: reproduces Table 1 of the manuscript.
##
## Eight configurations (C1-C8) crossing:
##   - sparsity s_0 in {1, 2, 3}
##   - annotation regime: informative, partial, uninformative, anti-informative
##
## Each configuration is replicated R = 500 times. Output is a data frame with
## per-configuration mean and Monte Carlo SE for coverage and mean CS size,
## for both SuSiE-RSS and V-SoRE.
##
## Usage:
##   Rscript simulate_main.R [output_path]
## Default output: sim_results.rds in the working directory.
## =============================================================================

source("vsore.R")
source("weights.R")
source("evaluate.R")
source("ld_utils.R")

#' Generate one simulation replicate.
#'
#' @param p        Number of variants.
#' @param s0       Number of true causal variants.
#' @param U        Annotation elevation. U = 0 -> uninformative; U > 0 ->
#'                 causal variants get exp(U) extra weight; U < 0 ->
#'                 anti-informative.
#' @param effect   True effect-size magnitude (per causal variant).
#' @param n_eff   Effective sample size for SuSiE-RSS.
#' @param p_blocks Number of LD blocks (passed to build_ld).
#' @param rho      AR(1) parameter (passed to build_ld).
#' @param noise_sd Standard deviation of Gaussian noise added to anchor weights
#'                 (creates partial-information regime when 0 < U < 1).
#' @return One-row data frame with method, coverage, cs_size, hat_tau2,
#'         and other diagnostics.
generate_replicate <- function(p, s0, U, effect, n_eff,
                               p_blocks = 5, rho = 0.85, noise_sd = 0.5) {
  R <- build_ld(p = p, n_blocks = p_blocks, rho = rho)
  causal_idx <- sample.int(p, s0)
  b <- numeric(p); b[causal_idx] <- effect
  z <- sample_z(R, b)

  ## Anchor weights: log-normal noise + annotation elevation U on causal indices.
  log_omega <- stats::rnorm(p, sd = noise_sd)
  log_omega[causal_idx] <- log_omega[causal_idx] + U
  omega0 <- make_anchor_weights(log_omega)

  ## SuSiE-RSS baseline (uniform prior).
  fit_susie <- susieR::susie_rss(z = z, R = R, n = n_eff, L = 10)
  cs_susie <- susieR::susie_get_cs(fit_susie)$cs
  cs_susie <- cs_susie[lengths(cs_susie) > 0]
  ev_susie <- evaluate_cs(causal_idx, cs_susie)

  ## V-SoRE with annotation-derived anchors.
  res_vsore <- vsore(z = z, R = R, omega0 = omega0, n_eff = n_eff)
  ev_vsore  <- evaluate_cs(causal_idx, res_vsore$cs)

  rbind(
    data.frame(method   = "SuSiE-RSS",
               coverage = ev_susie["coverage"],
               cs_size  = ev_susie["cs_size"],
               n_cs     = ev_susie["n_cs"],
               hat_tau2 = NA_real_,
               n_iter   = NA_integer_,
               row.names = NULL),
    data.frame(method   = "V-SoRE",
               coverage = ev_vsore["coverage"],
               cs_size  = ev_vsore["cs_size"],
               n_cs     = ev_vsore["n_cs"],
               hat_tau2 = res_vsore$hat_tau2,
               n_iter   = res_vsore$n_iter,
               row.names = NULL)
  )
}


#' Configuration table for the main simulation.
make_configs <- function() {
  data.frame(
    config   = paste0("C", 1:8),
    s0       = c(1, 1, 1, 1, 2, 2, 1, 3),
    U        = c(2.0, 1.5, 1.0, 2.5, 2.0, 1.5, 0.0, 2.0),
    effect   = c(0.20, 0.20, 0.18, 0.22, 0.18, 0.18, 0.18, 0.18),
    p        = rep(500, 8),
    p_blocks = rep(5, 8),
    rho      = rep(0.85, 8),
    n_eff    = rep(50000, 8),
    noise_sd = rep(0.5, 8),
    stringsAsFactors = FALSE
  )
}


#' Run all configurations and aggregate per-config results.
#'
#' @param R_replicates Number of replicates per configuration. Default 500.
#' @param verbose      If TRUE, print progress per configuration.
#' @return data.frame with per-configuration aggregated metrics.
run_main_simulation <- function(R_replicates = 500L, verbose = TRUE) {
  configs <- make_configs()
  out <- list()

  for (k in seq_len(nrow(configs))) {
    cfg <- configs[k, ]
    if (verbose)
      message(sprintf("[%s] s0=%d, U=%.1f, effect=%.2f, R=%d",
                      cfg$config, cfg$s0, cfg$U, cfg$effect, R_replicates))

    rep_results <- vector("list", R_replicates)
    for (r in seq_len(R_replicates)) {
      rep_results[[r]] <- generate_replicate(
        p        = cfg$p,
        s0       = cfg$s0,
        U        = cfg$U,
        effect   = cfg$effect,
        n_eff    = cfg$n_eff,
        p_blocks = cfg$p_blocks,
        rho      = cfg$rho,
        noise_sd = cfg$noise_sd
      )
    }

    df <- do.call(rbind, rep_results)
    by_method <- split(df, df$method)

    cfg_out <- lapply(names(by_method), function(m) {
      sub <- by_method[[m]]
      ## evaluate_cs() rows are columns coverage, cs_size, n_cs in df.
      ev_list <- lapply(seq_len(nrow(sub)), function(i)
        c(coverage = sub$coverage[i], cs_size = sub$cs_size[i]))
      agg <- aggregate_replicates(ev_list)
      data.frame(config        = cfg$config,
                 method        = m,
                 coverage_mean = agg["coverage_mean"],
                 coverage_mcse = agg["coverage_mcse"],
                 cs_size_mean  = agg["cs_size_mean"],
                 cs_size_mcse  = agg["cs_size_mcse"],
                 mean_hat_tau2 = if (m == "V-SoRE") mean(sub$hat_tau2, na.rm = TRUE)
                                 else NA_real_,
                 mean_n_iter   = if (m == "V-SoRE") mean(sub$n_iter,   na.rm = TRUE)
                                 else NA_real_,
                 R_replicates  = R_replicates,
                 stringsAsFactors = FALSE)
    })
    out[[k]] <- do.call(rbind, cfg_out)
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}


## --- Entry point when run as a script -----------------------------------------
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  output_path <- if (length(args) >= 1L) args[1] else "sim_results.rds"

  set.seed(20260506)
  results <- run_main_simulation(R_replicates = 500L, verbose = TRUE)

  saveRDS(results, output_path)
  message("\nSaved results to: ", output_path)
  print(results, row.names = FALSE)
}
