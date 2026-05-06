## =============================================================================
## Supplementary simulations: Sections S12 and S13 of the manuscript.
##
## - S12: Reference-panel LD mismatch sensitivity. Uses make_mismatch_ld() to
##        replace in-sample LD with a noisy reference-panel approximation.
## - S13: PL vs sqrt-flattened (IIA) sensitivity. Compares vsore() and
##        vsore_flattened() at matched configurations.
##
## Usage:
##   Rscript simulate_supplement.R [output_path]
## Default output: supp_results.rds
## =============================================================================

source("vsore.R")
source("weights.R")
source("evaluate.R")
source("ld_utils.R")

#' Reference-panel mismatch experiment (Supp Section S12).
#'
#' Generates z-scores under R_true and runs both SuSiE-RSS and V-SoRE on a
#' mismatched LD matrix derived from R_true via make_mismatch_ld().
run_mismatch_experiment <- function(p = 500, s0 = 2, U = 2.0, effect = 0.18,
                                    n_eff = 50000, n_ref_grid = c(500, 1000, 2500, 5000),
                                    R_replicates = 200L, verbose = TRUE) {
  configs <- expand.grid(n_ref = n_ref_grid, stringsAsFactors = FALSE)
  out <- list()
  for (k in seq_len(nrow(configs))) {
    n_ref <- configs$n_ref[k]
    if (verbose) message(sprintf("[mismatch] n_ref = %d", n_ref))

    rep_metrics_susie <- vector("list", R_replicates)
    rep_metrics_vsore <- vector("list", R_replicates)

    for (r in seq_len(R_replicates)) {
      R_true <- build_ld(p = p, n_blocks = 5, rho = 0.85)
      causal_idx <- sample.int(p, s0)
      b <- numeric(p); b[causal_idx] <- effect
      z <- sample_z(R_true, b)

      log_omega <- stats::rnorm(p, sd = 0.5)
      log_omega[causal_idx] <- log_omega[causal_idx] + U
      omega0 <- make_anchor_weights(log_omega)

      R_mismatch <- make_mismatch_ld(R_true, n_ref = n_ref)

      fit_s <- susieR::susie_rss(z = z, R = R_mismatch, n = n_eff, L = 10)
      cs_s <- susieR::susie_get_cs(fit_s)$cs
      rep_metrics_susie[[r]] <- evaluate_cs(causal_idx, cs_s)

      res_v <- vsore(z, R = R_mismatch, omega0 = omega0, n_eff = n_eff)
      rep_metrics_vsore[[r]] <- evaluate_cs(causal_idx, res_v$cs)
    }

    agg_s <- aggregate_replicates(rep_metrics_susie)
    agg_v <- aggregate_replicates(rep_metrics_vsore)
    out[[k]] <- rbind(
      data.frame(experiment = "mismatch", n_ref = n_ref, method = "SuSiE-RSS",
                 t(as.list(agg_s)), stringsAsFactors = FALSE),
      data.frame(experiment = "mismatch", n_ref = n_ref, method = "V-SoRE",
                 t(as.list(agg_v)), stringsAsFactors = FALSE)
    )
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}


#' PL vs sqrt-flattened (IIA) sensitivity (Supp Section S13).
#'
#' Compares vsore() (Plackett-Luce on omega0) and vsore_flattened()
#' (Plackett-Luce on sqrt(omega0)) at the same configurations.
#' Reports PIP correlation between the two methods.
run_iia_experiment <- function(p = 500, s0 = 2, U = 2.0, effect = 0.18,
                               n_eff = 50000, R_replicates = 200L,
                               verbose = TRUE) {
  if (verbose) message(sprintf("[iia] p=%d, s0=%d, U=%.1f", p, s0, U))
  pip_cor <- numeric(R_replicates)
  cs_size_pl <- numeric(R_replicates)
  cs_size_fl <- numeric(R_replicates)

  for (r in seq_len(R_replicates)) {
    R <- build_ld(p = p, n_blocks = 5, rho = 0.85)
    causal_idx <- sample.int(p, s0)
    b <- numeric(p); b[causal_idx] <- effect
    z <- sample_z(R, b)

    log_omega <- stats::rnorm(p, sd = 0.5)
    log_omega[causal_idx] <- log_omega[causal_idx] + U
    omega0 <- make_anchor_weights(log_omega)

    res_pl <- vsore(z, R, omega0, n_eff = n_eff)
    res_fl <- vsore_flattened(z, R, omega0, n_eff = n_eff)

    pip_cor[r]    <- stats::cor(res_pl$pip, res_fl$pip)
    cs_size_pl[r] <- mean(lengths(res_pl$cs))
    cs_size_fl[r] <- mean(lengths(res_fl$cs))
  }
  data.frame(
    experiment        = "iia",
    R_replicates      = R_replicates,
    pip_cor_mean      = mean(pip_cor),
    pip_cor_min       = min(pip_cor),
    pip_cor_max       = max(pip_cor),
    cs_size_pl_mean   = mean(cs_size_pl, na.rm = TRUE),
    cs_size_fl_mean   = mean(cs_size_fl, na.rm = TRUE),
    stringsAsFactors  = FALSE
  )
}


## --- Entry point when run as a script -----------------------------------------
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  output_path <- if (length(args) >= 1L) args[1] else "supp_results.rds"

  set.seed(20260506)
  res_mismatch <- run_mismatch_experiment(R_replicates = 200L)
  res_iia      <- run_iia_experiment(R_replicates = 200L)

  results <- list(mismatch = res_mismatch, iia = res_iia)
  saveRDS(results, output_path)
  message("\nSaved supplementary results to: ", output_path)
  message("\nReference-panel mismatch:")
  print(res_mismatch, row.names = FALSE)
  message("\nIIA / sqrt-flattening sensitivity:")
  print(res_iia, row.names = FALSE)
}
