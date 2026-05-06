## =============================================================================
## Test suite for SoRE / V-SoRE.
##
## Run with:
##   library(testthat); test_file("test-vsore.R")
## =============================================================================

library(testthat)

source("vsore.R")
source("weights.R")
source("evaluate.R")
source("ld_utils.R")


## =============================================================================
## vsore() core
## =============================================================================

test_that("vsore returns the expected output structure", {
  set.seed(1)
  p <- 50
  R <- build_ld(p = p, n_blocks = 5, rho = 0.8)
  z <- stats::rnorm(p)
  w <- runif(p, 0.5, 2)
  res <- vsore(z, R, w, n_eff = 5000, L = 5, T_max = 3)
  expect_named(res, c("fit", "pip", "hat_tau2", "cs", "n_iter",
                      "converged", "bypassed"))
  expect_length(res$pip, p)
  expect_true(all(res$pip >= 0) && all(res$pip <= 1))
  expect_true(is.na(res$hat_tau2) || res$hat_tau2 > 0)
  expect_type(res$cs, "list")
  expect_false(res$bypassed)
})


test_that("vsore triggers the SuSiE bypass under uniform anchor weights", {
  skip_if_not_installed("susieR")
  set.seed(2)
  p <- 50
  R <- build_ld(p = p, n_blocks = 5, rho = 0.8)
  b <- numeric(p); b[25] <- 3
  z <- sample_z(R, b)

  ## Uniform weights: bypass should fire.
  res_uniform <- vsore(z, R, omega0 = rep(1, p), n_eff = 50000, L = 5)
  expect_true(res_uniform$bypassed)
  expect_true(is.na(res_uniform$hat_tau2))
  expect_equal(res_uniform$n_iter, 0L)

  ## Output should match SuSiE-RSS exactly under bypass.
  fit_s <- susieR::susie_rss(z = z, R = R, n = 50000, L = 5)
  pip_s <- susieR::susie_get_pip(fit_s)
  expect_equal(res_uniform$pip, pip_s, tolerance = 1e-12)
})


test_that("vsore enters the reweighting loop when weights are non-uniform", {
  skip_if_not_installed("susieR")
  set.seed(3)
  p <- 50
  R <- build_ld(p = p, n_blocks = 5, rho = 0.8)
  b <- numeric(p); b[25] <- 3
  z <- sample_z(R, b)

  w <- runif(p, 0.5, 2)  # non-uniform
  res <- vsore(z, R, w, n_eff = 50000, L = 5, T_max = 5)
  expect_false(res$bypassed)
  expect_gte(res$n_iter, 1L)
  expect_true(!is.na(res$hat_tau2))
  expect_gte(res$hat_tau2, 0.05)
  expect_lte(res$hat_tau2, 1.5)
})


test_that("vsore_flattened delegates to vsore with sqrt(omega0)", {
  skip_if_not_installed("susieR")
  set.seed(4)
  p <- 50
  R <- build_ld(p = p, n_blocks = 5, rho = 0.8)
  z <- stats::rnorm(p)
  w <- runif(p, 0.5, 2)

  res_a <- vsore_flattened(z, R, w, n_eff = 5000, L = 5, T_max = 3)
  set.seed(4); .junk_seed <- runif(1)  # advance RNG identically
  res_b <- vsore(z, R, sqrt(w), n_eff = 5000, L = 5, T_max = 3)
  expect_equal(res_a$pip, res_b$pip)
})


test_that("vsore validates its inputs", {
  set.seed(5)
  p <- 30
  R <- build_ld(p = p, n_blocks = 3, rho = 0.5)
  z <- stats::rnorm(p)
  w <- rep(1, p)
  expect_error(vsore(z, R[1:5, 1:5], w, n_eff = 1000), "p x p")
  expect_error(vsore(z, R, w[1:5], n_eff = 1000), "length p")
  expect_error(vsore(z, R, c(-1, w[-1]), n_eff = 1000), "strictly positive")
  expect_error(vsore(z, R, w, n_eff = -1), "positive")
})


## =============================================================================
## Spearman helper handles ties
## =============================================================================

test_that(".spearman_r returns 0 on degenerate inputs and matches cor() with ties", {
  expect_equal(.spearman_r(c(1, 1, 1, 1, 1), c(1, 2, 3, 4, 5)), 0)
  expect_equal(.spearman_r(c(1, 2, 3, 4, 5), c(1, 1, 1, 1, 1)), 0)
  expect_equal(.spearman_r(1, 1), 0)              # too short
  ## Matches base::cor() under ties.
  set.seed(7)
  x <- sample(c(0, 1, 2), 30, replace = TRUE)
  y <- sample(c(0, 1, 2), 30, replace = TRUE)
  expect_equal(.spearman_r(x, y),
               suppressWarnings(stats::cor(x, y, method = "spearman")))
})


## =============================================================================
## .is_uniform helper
## =============================================================================

test_that(".is_uniform handles ties, scaling, and edge cases", {
  expect_true(.is_uniform(rep(1, 10)))
  expect_true(.is_uniform(rep(1e-7, 10)))
  expect_true(.is_uniform(rep(1e7, 10)))
  expect_false(.is_uniform(c(1, 1.0001, 1)))
  expect_true(.is_uniform(numeric(0)))
  expect_true(.is_uniform(7))
})


## =============================================================================
## compute_neff and the variance-correction formula
## =============================================================================

test_that("compute_neff matches the formula for the AD application parameters", {
  ## Computed value with K = 0.05, n1 = 21,982, n0 = 41,944.
  ## This is the value implied by eq. (7) of the manuscript.
  neff <- compute_neff(n_cases = 21982, n_controls = 41944, prevalence = 0.05)
  expect_equal(neff, 54238.4, tolerance = 1)
})


test_that("compute_neff_unadjusted gives 4 n1 n0 / N", {
  neff <- compute_neff_unadjusted(n_cases = 21982, n_controls = 41944)
  expect_equal(neff, 4 * 21982 * 41944 / (21982 + 41944), tolerance = 1e-6)
})


test_that("compute_neff validates inputs", {
  expect_error(compute_neff(0, 100, 0.05), "positive")
  expect_error(compute_neff(100, 100, 0), "in \\(0, 1\\)")
  expect_error(compute_neff(100, 100, 1.0), "in \\(0, 1\\)")
})


## =============================================================================
## evaluate_cs() semantics (NA for empty CS list)
## =============================================================================

test_that("evaluate_cs returns correct coverage and size for non-empty CS", {
  cs <- list(c(10, 11, 12), c(30, 31))
  ev <- evaluate_cs(causal_idx = c(11, 30), cs_list = cs)
  expect_equal(ev[["coverage"]], 1.0)
  expect_equal(ev[["cs_size"]], 2.5)
  expect_equal(ev[["n_cs"]], 2)
})


test_that("evaluate_cs returns coverage = 0, cs_size = NA for empty CS list", {
  ev <- evaluate_cs(causal_idx = c(1, 2), cs_list = list())
  expect_equal(ev[["coverage"]], 0.0)
  expect_true(is.na(ev[["cs_size"]]))
  expect_equal(ev[["n_cs"]], 0)
})


test_that("evaluate_cs drops empty CS entries before computing coverage", {
  ev <- evaluate_cs(causal_idx = c(5),
                    cs_list = list(c(5, 6), integer(0), c(20)))
  expect_equal(ev[["coverage"]], 1.0)
  expect_equal(ev[["cs_size"]], 1.5)
  expect_equal(ev[["n_cs"]], 2)
})


## =============================================================================
## aggregate_replicates() Monte Carlo SEs
## =============================================================================

test_that("aggregate_replicates computes binomial MCSE for coverage", {
  ev_list <- replicate(500, c(coverage = 0.95, cs_size = 100), simplify = FALSE)
  agg <- aggregate_replicates(ev_list)
  ## Coverage SE for p = 0.95, R = 500: sqrt(0.95*0.05/500) = 0.00975
  expect_equal(agg[["coverage_mean"]], 0.95, tolerance = 1e-6)
  expect_equal(agg[["coverage_mcse"]], sqrt(0.95 * 0.05 / 500), tolerance = 1e-6)
})


test_that("aggregate_replicates handles NA cs_size gracefully", {
  ev_list <- list(c(coverage = 1, cs_size = 10),
                  c(coverage = 0, cs_size = NA),
                  c(coverage = 1, cs_size = 8))
  agg <- aggregate_replicates(ev_list)
  expect_equal(agg[["coverage_mean"]], 2/3, tolerance = 1e-6)
  expect_equal(agg[["cs_size_mean"]], 9, tolerance = 1e-6)
})


## =============================================================================
## Annotation weights
## =============================================================================

test_that("build_annot_score combines log-additive sources correctly", {
  p <- 100
  mat <- cbind(c(rep(1, 10), rep(0, 90)),
               c(rep(0, 50), rep(1, 50)))
  a <- build_annot_score(mat, weights = c(1, 2))
  expect_length(a, p)
  expect_equal(a[1],   1)  # source 1 only
  expect_equal(a[51],  2)  # source 2 only
  expect_equal(a[5],   1)  # source 1 only (not in source 2)
  expect_equal(a[100], 2)  # source 2 only
})


test_that("build_annot_score uses default log weights for K = 4", {
  p <- 5
  mat <- matrix(0, p, 4); mat[1, 1] <- 1; mat[2, 2] <- 1
  a <- build_annot_score(mat)
  expect_equal(a[1], log(50), tolerance = 1e-12)
  expect_equal(a[2], log(5),  tolerance = 1e-12)
  expect_equal(a[3], 0)
})


test_that("build_prior_gwas_term gives log(1 + |z|) and is non-negative", {
  z <- c(-3, -1, 0, 1, 3)
  v <- build_prior_gwas_term(z)
  expect_equal(v, log(1 + abs(z)))
  expect_true(all(v >= 0))
})


test_that("make_anchor_weights exponentiates a log-weight vector", {
  log_w <- c(0, log(2), log(5))
  expect_equal(make_anchor_weights(log_w), c(1, 2, 5))
})


## =============================================================================
## ld_utils
## =============================================================================

test_that("build_ld produces a positive-definite p x p matrix", {
  R <- build_ld(p = 50, n_blocks = 5, rho = 0.85)
  expect_equal(dim(R), c(50, 50))
  ev_min <- min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
  expect_gt(ev_min, 0)
})


test_that("build_ld errors when p is not divisible by n_blocks", {
  expect_error(build_ld(p = 50, n_blocks = 7), "divisible")
})


test_that("make_mismatch_ld returns a positive-definite correlation matrix", {
  set.seed(11)
  R_true <- build_ld(p = 50, n_blocks = 5, rho = 0.85)
  R_m <- make_mismatch_ld(R_true, n_ref = 500)
  expect_equal(dim(R_m), c(50, 50))
  expect_equal(diag(R_m), rep(1, 50), tolerance = 1e-6)
  ev_min <- min(eigen(R_m, symmetric = TRUE, only.values = TRUE)$values)
  expect_gt(ev_min, 0)
})


test_that("sample_z generates a length-p vector", {
  set.seed(13)
  R <- build_ld(p = 30, n_blocks = 3, rho = 0.5)
  b <- numeric(30); b[5] <- 1
  z <- sample_z(R, b)
  expect_length(z, 30L)
  expect_true(all(is.finite(z)))
})


test_that("summarise_vsore returns a data.frame with the expected columns", {
  set.seed(15)
  p <- 30
  R <- build_ld(p = p, n_blocks = 3, rho = 0.7)
  b <- numeric(p); b[10] <- 2
  z <- sample_z(R, b)
  w <- runif(p, 0.5, 2)
  res <- vsore(z, R, w, n_eff = 5000, L = 3, T_max = 3)
  df <- summarise_vsore(res)
  expect_true(is.data.frame(df))
  expect_named(df, c("cs_index", "cs_size", "lead_pip",
                     "lead_variant", "hat_tau2"))
})
