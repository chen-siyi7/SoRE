#!/usr/bin/env Rscript
# =============================================================================
# ad_pipeline.R
# End-to-end Alzheimer's disease application driver.
# Reproduces Table 2 and Figures 2-3 of:
#   Chen (2026). Annotation-Ranked Bayesian Fine-Mapping of GWAS Signals.
#
# Inputs (place in data/):
#   data/kunkle_2019_stage1.tsv.gz  - Kunkle Stage 1 AD GWAS summary statistics
#   data/bellenguez_2022.tsv.gz     - Bellenguez 2022 prior-GWAS summary stats
#   data/1kg_phase3_eur.bed/bim/fam - 1000 Genomes Phase 3 European LD panel
#   data/corces_2020_microglia.bed  - Corces 2020 microglia ATAC-seq peaks
#
# Outputs:
#   output/ad_results.rds   - SuSiE-RSS, PolyFun+SuSiE, V-SoRE results
#   output/table2.csv       - Per-locus CS-size, purity, lead PIP, tau^2
#   output/fig_ad_baseline.pdf
#   output/fig_ad_annotated.pdf
#
# Expected runtime: 30-60 minutes on a 4-core node.
# =============================================================================

suppressPackageStartupMessages({
  library(SoRE)
  library(susieR)
  library(BEDMatrix)
  library(data.table)
})

# ── Configuration ─────────────────────────────────────────────────────────────
data_dir   <- "data"
out_dir    <- "output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

LD_LAMBDA  <- 0.01
WINDOW_KB  <- 500
N_CASES    <- 21982
N_CTRLS    <- 41944
N_EFF      <- 4 * N_CASES * N_CTRLS / (N_CASES + N_CTRLS)

loci_file <- system.file("extdata", "kunkle_2019_loci.csv", package = "SoRE")
loci <- fread(loci_file)

# ── Per-locus loop ────────────────────────────────────────────────────────────
results <- vector("list", nrow(loci))

for (i in seq_len(nrow(loci))) {
  locus <- loci[i]
  cat(sprintf("[%2d/%2d] %s (chr%d:%d)\n", i, nrow(loci),
              locus$name, locus$chrom, locus$pos))

  # Load z-scores and LD for this locus (user supplies extract_locus_data())
  d <- extract_locus_data(
    chrom = locus$chrom, pos = locus$pos,
    window_kb = WINDOW_KB,
    gwas_file = file.path(data_dir, "kunkle_2019_stage1.tsv.gz"),
    bed_prefix = file.path(data_dir, "1kg_phase3_eur")
  )

  R <- (1 - LD_LAMBDA) * d$R + LD_LAMBDA * diag(nrow(d$R))

  # Anchor weights from microglia ATAC-seq + Bellenguez prior GWAS
  omega <- build_anchor_weights(
    chrom_open = d$in_chrom_peak,
    z_prior    = d$z_prior
  )

  # SuSiE-RSS baseline
  fit_susie <- susie_rss(z = d$z, R = R, n = N_EFF, L = 10)

  # PolyFun+SuSiE (fixed prior weights, no outer loop)
  fit_polyfun <- susie_rss(z = d$z, R = R, n = N_EFF, L = 10,
                           prior_weights = omega / sum(omega))

  # V-SoRE (outer reweighting loop)
  fit_vsore <- vsore(z = d$z, R = R, anchor_weights = omega, n = N_EFF)

  results[[i]] <- list(
    locus    = locus$name,
    chrom    = locus$chrom,
    n_var    = length(d$z),
    susie    = fit_susie,
    polyfun  = fit_polyfun,
    vsore    = fit_vsore
  )
}

# ── Save results ──────────────────────────────────────────────────────────────
saveRDS(results, file.path(out_dir, "ad_results.rds"))

# ── Summary table (Table 2 of the manuscript) ─────────────────────────────────
summary_table <- rbindlist(lapply(results, function(r) {
  cs_s  <- credible_sets(SoRE:::susie_pips(r$susie),   0.95, R = NULL)
  cs_p  <- credible_sets(SoRE:::susie_pips(r$polyfun), 0.95, R = NULL)
  cs_v  <- credible_sets(r$vsore$pip,                  0.95, R = NULL)
  data.table(locus = r$locus, chrom = r$chrom, n = r$n_var,
             susie_size = cs_s$size, polyfun_size = cs_p$size,
             vsore_size = cs_v$size,
             vsore_lead_pip = cs_v$lead_pip,
             tau2 = r$vsore$tau2)
}))
fwrite(summary_table, file.path(out_dir, "table2.csv"))

cat("\nDone. Results saved to", out_dir, "\n")
