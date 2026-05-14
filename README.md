# SoRE: Annotation-Ranked Bayesian Fine-Mapping of GWAS Signals

<!-- badges: start -->
[![R-CMD-check](https://github.com/chen-siyi7/SoRE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/chen-siyi7/SoRE/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

R package implementation of **SoRE** (Sum of Ranked Effects), a Bayesian
fine-mapping prior that encodes functional annotation through a
Plackett-Luce distribution on variant orderings, with a locus-specific
variance parameter that adaptively controls how closely the prior follows
the annotation ranking.

The package implements two interfaces:

- **`vsore()`** — empirical-Bayes variational implementation that wraps
  SuSiE-RSS in an outer reweighting loop. Recommended for genome-wide
  pipelines. Runs in under two seconds per locus at p = 500.
- **`sore_gibbs()`** — full blocked Gibbs sampler for the complete
  hierarchical posterior. Use when full uncertainty quantification is
  required.

## Installation

```r
# install.packages("remotes")
remotes::install_github("chen-siyi7/SoRE")
```

The package depends on `susieR (>= 0.12.0)`. Install it first if needed:

```r
install.packages("susieR")
```

## Quick start

```r
library(SoRE)

set.seed(1)
p <- 200; n <- 50000
R <- 0.9 ^ abs(outer(1:p, 1:p, "-"))
causal <- c(50, 150)
b <- numeric(p); b[causal] <- 5
z <- as.numeric(R %*% b + crossprod(chol(R), rnorm(p)))

# Build anchor weights (here: strong annotation at the causal positions)
omega <- rep(1, p); omega[causal] <- 100

# Fit V-SoRE
fit <- vsore(z, R, anchor_weights = omega, n = n)

# Inspect
fit$tau2                              # annotation-data agreement
sum(fit$pip > 0.5)                    # variants with PIP > 0.5
credible_sets(fit$pip, R = R, min_purity = 0.5)
```

See `vignette("intro-to-sore")` for a longer walkthrough.

## Method summary

SoRE places a Plackett-Luce prior on the latent variant ordering r and
maps each rank position to a prior inclusion probability through a logistic
function. The model is

```
omega_j   = tilde-omega_j * exp(eta_j),   eta_j ~ N(0, tau^2)
r | omega ~ PL(omega)
z_j | k, alpha, r ~ Bernoulli(expit(alpha (k - j)))
b_j | z_j ~ z_j N(0, W) + (1 - z_j) delta_0
zhat | b  ~ N(R b, R)
```

The key design choices:

1. **Scale invariance.** PL(r | c·omega) = PL(r | omega) for any c > 0, so
   the prior depends on the anchor weights only through the ordering they
   induce. Heterogeneous annotation sources do not need cross-source
   calibration.

2. **Adaptive trust.** The locus-specific tau^2 is estimated from the
   Spearman rank correlation between the anchor weights and the current
   PIPs. When annotation and data agree, tau^2 is small and the prior
   anchors tightly to the annotation. When they disagree, tau^2 grows and
   V-SoRE reverts to SuSiE-RSS.

3. **Theoretical guarantee.** Under an oracle annotation condition, the
   posterior contraction rate improves from O(s_0 log p / N) to
   O(s_0 log K_0 / N), a log(p / K_0) reduction.

## Reproducing the manuscript results

The `inst/scripts/` directory contains the analysis scripts:

- **`simulate_sore.R`** — reproduces the eight-configuration simulation
  study (Table 1, Figure 1 of the manuscript) at both N_ref = 10,000 and
  N_ref = 500. Calls FINEMAP, SuSiE-RSS, PolyFun+SuSiE, and V-SoRE.
- **`ad_pipeline.R`** — end-to-end Alzheimer's disease application
  driver. Reads Kunkle 2019 summary statistics, extracts the 17 loci,
  computes LD from 1000 Genomes Phase 3 EUR, and runs all three methods.
- **`make_ad_figures.R`** — generates Figures 2 and 3 from the saved
  results.

To run from a clone of this repository:

```bash
git clone https://github.com/chen-siyi7/SoRE.git
cd SoRE
R CMD INSTALL .
Rscript inst/scripts/simulate_sore.R  # ~3 hours on 24 cores
```

The AD application requires the Kunkle 2019 summary statistics (NIAGADS
access agreement) and the 1000 Genomes Phase 3 EUR genotypes. See
`inst/scripts/ad_pipeline.R` for the expected directory layout.

## Diagnostic interpretation

The estimated tau^2 is the per-locus diagnostic. The calibrated
thresholds are:

| tau^2 range | Interpretation |
|---|---|
| < 0.2     | Strong annotation-data agreement |
| 0.2 - 1.0 | Moderate agreement |
| > 1.0     | Data has overridden annotation; treat as effectively unannotated |

The helper `tau2_diagnostic()` returns these categories.

## Citation

If you use SoRE in your work, please cite:

> Chen, S. (2026). Bayesian Fine-Mapping with an Annotation Ranking Prior
> Using GWAS Summary Statistics. *Biometrics*.

```bibtex
@article{chen2026sore,
  title   = {Bayesian Fine-Mapping with an Annotation Ranking Prior Using GWAS Summary Statistics},
  author  = {Chen, Siyi},
  journal = {Biometrics},
  year    = {2026}
}
```

## License

MIT License. See `LICENSE` for details.

## Contact

Siyi Chen
School of Public Health
LSU Health Sciences Center New Orleans
sche11@lsuhsc.edu

Bug reports and feature requests:
[GitHub Issues](https://github.com/chen-siyi7/SoRE/issues)
