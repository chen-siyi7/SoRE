test_that("pl_sample produces valid permutations", {
  set.seed(1)
  omega <- c(10, 5, 3, 1, 1)
  for (i in 1:10) {
    r <- pl_sample(omega)
    expect_equal(sort(r), 1:5)
  }
})

test_that("pl_sample concentrates mass on high-weight variants", {
  set.seed(1)
  omega <- c(100, 10, 1, 1, 1)
  draws <- pl_sample(omega, n = 2000)
  first_position <- draws[, 1]
  prop_top1 <- mean(first_position == 1)
  expect_gt(prop_top1, 0.7)
})

test_that("pl_logmass is invariant to scaling of omega", {
  r <- c(2, 4, 1, 3, 5)
  omega <- c(3, 1, 4, 1, 5)
  expect_equal(pl_logmass(r, omega), pl_logmass(r, 100 * omega))
})

test_that("pl_logmass is maximized at the rank ordering", {
  omega <- c(10, 5, 3, 1, 1)
  best_r <- order(omega, decreasing = TRUE)
  worst_r <- rev(best_r)
  expect_gt(pl_logmass(best_r, omega), pl_logmass(worst_r, omega))
})

test_that("build_anchor_weights returns positive weights", {
  set.seed(1)
  w <- build_anchor_weights(
    coding     = sample(1:3, 50, replace = TRUE),
    chrom_open = runif(50) < 0.1,
    conserved  = runif(50) < 0.1,
    z_prior    = rnorm(50)
  )
  expect_equal(length(w), 50)
  expect_true(all(w > 0))
})

test_that("build_anchor_weights with no args needs p", {
  expect_error(build_anchor_weights(), "Provide p")
  expect_equal(build_anchor_weights(p = 10), rep(1, 10))
})

test_that("credible_sets returns CS with cumulative PIP >= coverage", {
  set.seed(1)
  pips <- c(0.5, 0.3, 0.1, 0.05, 0.05)
  cs <- credible_sets(pips, coverage = 0.9)
  expect_gte(cs$coverage, 0.9)
  expect_equal(cs$lead_pip, 0.5)
})

test_that("credible_sets handles a single-variant signal", {
  pips <- c(0.99, 0.001, 0.001, 0.001)
  cs <- credible_sets(pips, coverage = 0.95)
  expect_equal(cs$size, 1)
  expect_equal(cs$cs, 1)
})

test_that("spearman_top_k returns Spearman correlation", {
  x <- c(1, 2, 3, 4, 5)
  y <- c(2, 1, 3, 4, 5)
  expect_lt(spearman_top_k(x, y), 1)
  expect_gt(spearman_top_k(x, y), 0)
})

test_that("tau2_diagnostic returns expected categories", {
  expect_equal(tau2_diagnostic(c(0.1, 0.5, 1.5, NA)),
               c("strong", "moderate", "uninformative", "no annotation"))
})
