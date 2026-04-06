## tests/testthat/test-vsore.R

test_that("vsore returns valid output structure", {
  set.seed(1)
  p  <- 50
  R  <- build_ld(p = p, n_blocks = 5, rho = 0.8)
  z  <- rnorm(p)
  w  <- rep(1, p)

  res <- vsore(z, R, w, n_eff = 5000, L = 5, T_max = 3)

  expect_named(res, c("fit", "pip", "hat_tau2", "cs", "n_iter"))
  expect_length(res$pip, p)
  expect_true(all(res$pip >= 0) && all(res$pip <= 1))
  expect_true(res$hat_tau2 > 0)
  expect_type(res$cs, "list")
})

test_that("vsore with uniform weights approximates SuSiE-RSS", {
  skip_if_not_installed("susieR")
  set.seed(2)
  p     <- 50
  R     <- build_ld(p = p, n_blocks = 5, rho = 0.8)
  b     <- numeric(p); b[25] <- 3
  z     <- as.numeric(R %*% b + t(chol(R)) %*% rnorm(p))
  w     <- rep(1, p)

  res_v <- vsore(z, R, w, n_eff = 50000, L = 5, T_max = 5)
  fit_s <- susieR::susie_rss(z, R, n = 50000, L = 5)
  pip_s <- susieR::susie_get_pip(fit_s)

  cor_pips <- cor(res_v$pip, pip_s)
  expect_gt(cor_pips, 0.99)
})

test_that("vsore_mallows returns valid structure", {
  set.seed(3)
  p   <- 50
  R   <- build_ld(p = p, n_blocks = 5, rho = 0.8)
  z   <- rnorm(p)
  w   <- runif(p, 0.5, 2)

  res <- vsore_mallows(z, R, w, n_eff = 5000, L = 5, T_max = 3)
  expect_length(res$pip, p)
  expect_true(all(res$pip >= 0) && all(res$pip <= 1))
})

test_that("compute_neff matches Kunkle 2019 expected value", {
  neff <- compute_neff(n_cases = 21982, n_controls = 41944, prevalence = 0.05)
  expect_gt(neff, 40000)
  expect_lt(neff, 45000)
})

test_that("evaluate_cs returns correct coverage", {
  cs  <- list(c(10, 11, 12), c(30, 31))
  ev  <- evaluate_cs(causal_idx = c(11, 30), cs_list = cs)
  expect_equal(ev[["coverage"]], 1.0)
  expect_equal(ev[["cs_size"]], 2.5)
})

test_that("evaluate_cs handles empty cs_list", {
  ev <- evaluate_cs(causal_idx = c(1, 2), cs_list = list())
  expect_equal(ev[["coverage"]], 0.0)
  expect_equal(ev[["cs_size"]], 0.0)
})

test_that("build_annot_score combines sources correctly", {
  p   <- 100
  mat <- matrix(c(rep(1, 10), rep(0, 90),
                  rep(0, 50), rep(1, 50)), ncol = 2)
  a   <- build_annot_score(mat, weights = c(1, 2))
  expect_length(a, p)
  expect_equal(a[1],  1)   # only source 1
  expect_equal(a[51], 2)   # only source 2
  expect_equal(a[5],  1)   # source 1 only (not in source 2)
})
