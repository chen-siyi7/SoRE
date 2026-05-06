# SoRE: Annotation-Ranked Bayesian Fine-Mapping

R implementation of the **SoRE** (Sum of Ranked Effects) prior and the **V-SoRE** empirical-Bayes reweighting algorithm for fine-mapping genome-wide association study (GWAS) summary statistics, accompanying the manuscript:

> Chen S. (2026) *Annotation-Ranked Bayesian Fine-Mapping of GWAS Signals with Data-Adaptive Reversion to Exchangeability.*

## Overview

V-SoRE wraps `susieR::susie_rss()` as an inner solver inside an outer empirical-Bayes loop that reweights variants by the agreement between annotation ranks and current posterior inclusion probabilities. The locus-specific dispersion parameter `hat_tau2` quantifies that agreement: small values (`< 0.2`) indicate strong agreement, intermediate values (`0.2`-`1`) moderate agreement, and values above `1` indicate that the data have substantially overridden the annotation prior.

Under uniform anchor weights, V-SoRE bypasses the reweighting loop and returns the SuSiE-RSS solution directly with `hat_tau2 = NA`, as required by the algorithmic-reduction guarantee in Proposition 3 of the manuscript.

## Installation

The code uses base R and depends on `susieR`. Optional helper dependencies are `Matrix` (for block-diagonal LD construction) and `testthat` (for the test suite).

```r
install.packages(c("susieR", "Matrix", "testthat"))
```

There is no need to install this code as a package; source the relevant files at the start of your script:

```r
source("vsore.R")
source("weights.R")
source("evaluate.R")
source("ld_utils.R")
```

## Quick example

```r
source("vsore.R")
source("weights.R")
source("ld_utils.R")

set.seed(1)
p <- 200
R <- build_ld(p, n_blocks = 5, rho = 0.85)

## Two true causal variants with effect size 0.18.
b <- numeric(p)
b[c(45, 130)] <- 0.18
z <- sample_z(R, b)

## Anchor weights with mild annotation elevation at the causal variants.
log_omega <- rnorm(p, sd = 0.5)
log_omega[c(45, 130)] <- log_omega[c(45, 130)] + 2.0
omega0 <- make_anchor_weights(log_omega)

res <- vsore(z, R, omega0, n_eff = 50000)
res$hat_tau2          # estimated annotation-data agreement
res$pip[c(45, 130)]   # PIPs at the true causal variants
res$cs                # 95% credible sets
```

## Reproducing the manuscript

The two simulation scripts reproduce the figures and tables of the paper:

```bash
Rscript simulate_main.R       # reproduces Table 1 (8 configs x 500 reps)
Rscript simulate_supplement.R # reproduces Sections S12 (mismatch) and S13 (IIA)
```

Both scripts save their output as RDS files. The main simulation runs in roughly 30 minutes on a single core; reduce `R_replicates` for a faster smoke test.

The Alzheimer's-disease application of Section 6 requires the following external data files, which are not redistributed with this code:

- **GWAS summary statistics:** Kunkle et al.\ 2019 Stage 1 (NIAGADS dataset NG00075).
- **Reference-panel LD:** 1000 Genomes Phase 3, European subpopulation, available through the International Genome Sample Resource.
- **Microglia ATAC-seq peaks:** Corces et al.\ 2020 (ENCODE accession ENCSR724KET).
- **Prior-GWAS z-scores:** Bellenguez et al.\ 2022 (GWAS Catalog GCST90027158).

A worked example using these inputs is described in Section 6.1 of the manuscript.

## Files

| File                   | Purpose                                                     |
| ---------------------- | ----------------------------------------------------------- |
| `vsore.R`              | Core V-SoRE algorithm with uniform-weight bypass.           |
| `weights.R`            | Annotation-weight construction and `compute_neff`.          |
| `evaluate.R`           | Credible-set coverage / size and replicate aggregation.     |
| `ld_utils.R`           | LD-matrix utilities and z-score sampler.                    |
| `simulate_main.R`      | Main simulation reproducing Table 1.                        |
| `simulate_supplement.R`| Supplementary simulations (mismatch, IIA).                  |
| `test-vsore.R`         | testthat suite covering all exported functions.             |

## Tests

Run the full test suite:

```r
library(testthat)
test_file("test-vsore.R")
```

All tests should pass on R `>= 4.1` with `susieR >= 0.12.0`.

## Notes on the effective sample size

The `compute_neff()` function implements eq. (7) of the manuscript, the liability-threshold-corrected effective sample size of Lee et al.\ 2012. With the AD application parameters (`n1 = 21,982`, `n0 = 41,944`, `K = 0.05`) it returns `54,238`. An unadjusted alternative `compute_neff_unadjusted()` returns `4 n1 n0 / N = 57,693` and is provided for comparison; it is not used in the manuscript.

If you obtain a different value when reproducing the AD analysis, please verify that the population-prevalence value `K` you have used matches the one assumed in the original analysis.

## License

This code is distributed under the MIT License. See `LICENSE` for the full text.

## Citation

If you use this code in published work, please cite:

> Chen S. (2026) *Annotation-Ranked Bayesian Fine-Mapping of GWAS Signals with Data-Adaptive Reversion to Exchangeability.* Submitted to *Statistics in Medicine*.

## Contact

Siyi Chen, School of Public Health, Louisiana State University Health Sciences Center New Orleans. Email: sche11@lsuhsc.edu.
