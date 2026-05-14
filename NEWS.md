# SoRE 0.1.0

Initial release accompanying the manuscript

> Chen, S. (2026). Bayesian Fine-Mapping with an Annotation Ranking Prior
> Using GWAS Summary Statistics. *Biometrics*.

Features:

* `vsore()` — empirical-Bayes variational fine-mapping with annotation ranking.
* `sore_gibbs()` — full blocked Gibbs sampler for the hierarchical posterior.
* `pl_sample()` / `pl_logmass()` — Plackett-Luce distribution utilities.
* `build_anchor_weights()` — anchor-weight construction from coding,
  chromatin, conservation, and prior-GWAS sources.
* `credible_sets()` — greedy credible-set construction with optional purity
  filter.
* `tau2_diagnostic()` — categorical interpretation of the V-SoRE
  diagnostic.
* `inst/scripts/` — reproducibility scripts for the manuscript simulations
  and Alzheimer's disease application.
