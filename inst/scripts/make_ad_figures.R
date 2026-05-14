#!/usr/bin/env Rscript
# =============================================================================
# make_ad_figures.R
# Generates the three application figures for the main text:
#
#   fig_ad_baseline.pdf      -- Figure: baseline (V-SoRE = SuSiE-RSS)
#   plot_G_ad_tau_vs_delta.pdf -- Figure: tau^2 vs CS reduction at 13 annotated loci
#   fig_ad_annotated.pdf     -- Figure: annotated V-SoRE results (4 panels)
#
# All data are taken directly from Table 2 of the main paper.
# No external data files are required.
#
# Output directory: ./output/  (created if absent)
# =============================================================================

OUT_DIR <- "output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. TABLE 2 DATA  (transcribed verbatim from the paper)
# =============================================================================

ad <- data.frame(
  locus    = c("CR1","BIN1","chr6 locus","HLA","CD2AP",
               "EPHA1","CLU","CELF1","MS4A6A","PICALM",
               "SPI1","MS4A","SORL1",
               "SLC24A4","APOE","ABCA7","CD33"),
  chr      = c(1,2,6,6,6,7,8,10,11,11,11,11,11,14,19,19,19),
  n        = c(2159,2801,2697,13757,2610,1307,2860,2880,
               2140,2592,1396,1709,2214,3072,2776,3581,2452),
  susie_cs = c(18.0,75.7,1874.0,419.0,81.0,268.5,24.2,2225.4,
               76.0,34.0,104.0,88.0,1733.8,
               1274.8,2.8,2135.3,20.5),
  vsore_cs = c(16.0,19.9,199.0,142.0,55.0,18.0,9.5,1657.2,
               64.0,23.0,53.0,64.0,977.3,
               1274.8,2.8,2135.3,20.5),
  delta    = c(-11,-74,-89,-66,-32,-93,-61,-26,
               -16,-32,-49,-27,-44,
               0,0,0,0),
  purity   = c(0.91,0.00,0.00,0.01,0.26,0.00,0.00,0.00,
               0.57,0.78,0.18,0.32,0.00,
               0.00,1.00,0.00,0.77),
  lead_pip = c(0.81,0.99,0.83,0.30,0.50,0.64,0.86,0.73,
               0.70,0.78,0.52,0.47,0.94,
               0.68,1.00,0.78,1.00),
  tau2     = c(0.63,1.13,0.70,2.36,2.04,1.03,0.45,0.37,
               0.96,0.98,1.49,0.74,0.70,
               NA,NA,NA,NA),
  annotated = c(rep(TRUE,13), rep(FALSE,4)),
  stringsAsFactors = FALSE
)

# Colour scale by tau^2 category
tau2_col <- function(tau2) {
  ifelse(is.na(tau2),     "grey70",
  ifelse(tau2 <= 0.2,     "#2166ac",   # strong agreement  (blue)
  ifelse(tau2 <= 1.0,     "#4dac26",   # moderate          (green)
                          "#d73027"))) # caution           (red)
}
ad$col <- tau2_col(ad$tau2)

# Short locus labels for axis ticks
ad$label <- sub(" \\(chr.*\\)", "", ad$locus)

ann  <- ad[ad$annotated, ]       # 13 annotated loci
all17 <- ad                      # all 17

# =============================================================================
# 2. HELPER: multi-panel layout
# =============================================================================

panel_letter <- function(x, y, lab, cex=1.1)
  mtext(lab, side=3, adj=0, line=0.2, at=x, cex=cex, font=2)

# =============================================================================
# 3. fig_ad_baseline.pdf
#    Under uniform weights V-SoRE = SuSiE-RSS exactly.
#    Panel A: CS-size scatter (identity). B: bar chart of SuSiE CS size.
#    C: lead-variant PIP.
# =============================================================================

pdf(file.path(OUT_DIR,"fig_ad_baseline.pdf"), width=10, height=4)
layout(matrix(1:3, nrow=1), widths=c(1.1, 1.4, 1.4))

# --- Panel A: scatter (identity) -------------------------------------------
par(mar=c(4,4.5,2.5,1), las=1, cex.axis=0.82)
cs_range <- range(c(all17$susie_cs, all17$vsore_cs), na.rm=TRUE)
lims <- c(1, max(cs_range)*1.05)
plot(all17$susie_cs, all17$susie_cs, log="xy",
     xlim=lims, ylim=lims,
     xlab="SuSiE-RSS CS size", ylab="V-SoRE CS size",
     pch=21, bg="grey50", col="grey30", cex=1.0)
abline(0, 1, lty=2, col="grey50")
mtext("A", side=3, adj=0, line=0.4, font=2, cex=1.0)

# --- Panel B: bar chart of SuSiE CS size ------------------------------------
par(mar=c(6,4.5,2.5,0.5), las=2, cex.axis=0.68)
ord <- order(all17$susie_cs)
bp <- barplot(all17$susie_cs[ord],
              names.arg = all17$label[ord],
              col = "steelblue", border = NA,
              ylab = "SuSiE-RSS CS size", log = "y",
              ylim = c(1, max(all17$susie_cs)*1.5),
              cex.names = 0.65)
mtext("B", side=3, adj=0, line=0.4, font=2, cex=1.0)

# --- Panel C: lead-variant PIP ----------------------------------------------
par(mar=c(6,4,2.5,0.5), las=2, cex.axis=0.68)
barplot(all17$lead_pip[ord],
        names.arg = all17$label[ord],
        col = "steelblue", border = NA,
        ylab = "Lead-variant PIP",
        ylim = c(0, 1.05),
        cex.names = 0.65)
abline(h = 0.5, lty = 2, col = "grey50", lwd = 1.2)
mtext("C", side=3, adj=0, line=0.4, font=2, cex=1.0)

dev.off()
message("Saved fig_ad_baseline.pdf")

# =============================================================================
# 4. plot_G_ad_tau_vs_delta.pdf
#    Per-locus tau^2 vs CS reduction Delta for the 13 annotated loci.
#    Points grouped by diagnostic category.
# =============================================================================

pdf(file.path(OUT_DIR,"plot_G_ad_tau_vs_delta.pdf"), width=6.5, height=5)
par(mar=c(4.5,4.5,1,1), las=1)

cat_cols  <- c("#2166ac","#4dac26","#d73027")
cat_pchs  <- c(19L, 17L, 15L)
cat_labs  <- c(expression(hat(tau)^2 <= 0.2),
               expression(0.2 < hat(tau)^2 <= 1),
               expression(hat(tau)^2 > 1))

plot(ann$tau2, ann$delta,
     xlab = expression(hat(tau)^2),
     ylab = expression(paste("CS reduction ", Delta, " (%)")),
     xlim = c(0, max(ann$tau2)*1.08),
     ylim = c(min(ann$delta)*1.05, 0),
     type = "n")
abline(h=0, lty=2, col="grey60", lwd=1.1)
abline(v=c(0.2,1.0), lty=3, col="grey75", lwd=1.0)

# plot points by category
for (i in seq_along(ann$locus)) {
  t2 <- ann$tau2[i]; dl <- ann$delta[i]
  pch <- if (t2 <= 0.2) cat_pchs[1] else if (t2 <= 1.0) cat_pchs[2] else cat_pchs[3]
  col <- if (t2 <= 0.2) cat_cols[1] else if (t2 <= 1.0) cat_cols[2] else cat_cols[3]
  points(t2, dl, pch=pch, col=col, cex=1.3, lwd=1.5)
}

# label selected loci
label_idx <- c(which(ann$locus=="CLU"), which(ann$locus=="EPHA1"),
               which(ann$locus=="BIN1"), which(ann$locus=="CELF1"),
               which(ann$locus=="HLA"))
for (i in label_idx) {
  nudge_x <- 0.04; nudge_y <- 1.5
  text(ann$tau2[i] + nudge_x, ann$delta[i] + nudge_y,
       ann$label[i], cex=0.75, adj=0)
}

legend("bottomleft",
       legend = cat_labs,
       col    = cat_cols,
       pch    = cat_pchs,
       bty    = "n", cex = 0.85, pt.cex = 1.2)
dev.off()
message("Saved plot_G_ad_tau_vs_delta.pdf")

# =============================================================================
# 5. fig_ad_annotated.pdf
#    Panel A: CS-size scatter (SuSiE vs V-SoRE) coloured by tau^2.
#    Panel B: side-by-side CS size bars per locus.
#    Panel C: lead-variant PIP per locus.
#    Panel D: tau^2 per annotated locus with reference lines.
# =============================================================================

pdf(file.path(OUT_DIR,"fig_ad_annotated.pdf"), width=11, height=8.5)
layout(matrix(c(1,2,3,4), nrow=2, byrow=TRUE),
       widths=c(1.1, 1.6), heights=c(1,1))

# --- Panel A: scatter SuSiE vs V-SoRE --------------------------------------
par(mar=c(4.5,4.5,2.5,1), las=1, cex.axis=0.82)
lims2 <- c(1, max(all17$susie_cs)*1.05)
plot(all17$susie_cs, all17$vsore_cs, log="xy",
     xlim=lims2, ylim=lims2,
     xlab="SuSiE-RSS CS size", ylab="V-SoRE CS size",
     type="n")
abline(0,1, lty=2, col="grey55")
# unannotated (grey)
una <- all17[!all17$annotated, ]
points(una$susie_cs, una$vsore_cs, pch=21, bg="grey65", col="grey40", cex=1.0)
# annotated, coloured by tau^2
points(ann$susie_cs, ann$vsore_cs,
       pch=21, bg=ann$col, col="grey20", cex=1.25)
# label a few
for (nm in c("CLU","EPHA1","BIN1","CELF1","HLA")) {
  idx <- which(all17$locus == nm)
  text(all17$susie_cs[idx]*1.4, all17$vsore_cs[idx], nm, cex=0.72)
}
legend("topleft",
       legend=c(expression(hat(tau)^2<=0.2),
                expression(0.2<hat(tau)^2<=1),
                expression(hat(tau)^2>1),
                "No annotation"),
       pch=21,
       pt.bg=c("#2166ac","#4dac26","#d73027","grey65"),
       bty="n", cex=0.78, pt.cex=1.2)
mtext("A", side=3, adj=0, line=0.4, font=2, cex=1.0)

# --- Panel B: CS size bars per locus (SuSiE and V-SoRE side by side) --------
par(mar=c(6.5,4.5,2.5,0.5), las=2, cex.axis=0.65)
ord2 <- order(all17$susie_cs)
locs_ord <- all17$label[ord2]
n_loc <- nrow(all17)
x_pos <- seq_len(n_loc)
x_s   <- x_pos - 0.2
x_v   <- x_pos + 0.2

# log scale via manual placement
ys  <- log10(pmax(all17$susie_cs[ord2], 0.5))
yv  <- log10(pmax(all17$vsore_cs[ord2], 0.5))
ymax <- log10(max(all17$susie_cs)*1.5)

plot(NA, xlim=c(0.5, n_loc+0.5), ylim=c(0, ymax),
     xlab="", ylab="CS size (log\u2081\u2080 scale)",
     xaxt="n", yaxt="n")
axis(1, at=x_pos, labels=locs_ord, cex.axis=0.65, las=2)
y_breaks <- c(1,2,5,10,20,50,100,200,500,1000,2000)
axis(2, at=log10(y_breaks),
     labels=format(y_breaks, big.mark=","), cex.axis=0.7, las=1)
rect(x_s-0.18, 0, x_s+0.18, ys, col="steelblue", border=NA)
rect(x_v-0.18, 0, x_v+0.18, yv,
     col=ifelse(is.na(all17$tau2[ord2]),"grey60",
                tau2_col(all17$tau2[ord2])),
     border=NA)
legend("topleft",
       legend=c("SuSiE-RSS","V-SoRE"),
       fill=c("steelblue","#4dac26"), border=NA, bty="n", cex=0.78)
mtext("B", side=3, adj=0, line=0.4, font=2, cex=1.0)

# --- Panel C: lead-variant PIP per locus ------------------------------------
par(mar=c(6.5,4.5,2.5,0.5), las=2, cex.axis=0.65)
barplot(all17$lead_pip[ord2],
        names.arg=locs_ord,
        col=ifelse(is.na(all17$tau2[ord2]),"grey60",
                   tau2_col(all17$tau2[ord2])),
        border=NA,
        ylab="Lead-variant PIP",
        ylim=c(0,1.05),
        cex.names=0.65)
abline(h=0.5, lty=2, col="grey50", lwd=1.2)
mtext("C", side=3, adj=0, line=0.4, font=2, cex=1.0)

# --- Panel D: tau^2 per annotated locus -------------------------------------
par(mar=c(6.5,4.8,2.5,0.5), las=2, cex.axis=0.72)
ord3 <- order(ann$tau2)
bp_d <- barplot(ann$tau2[ord3],
        names.arg=ann$label[ord3],
        col=tau2_col(ann$tau2[ord3]),
        border=NA,
        ylab=expression(hat(tau)^2),
        ylim=c(0, max(ann$tau2)*1.15),
        cex.names=0.70)
abline(h=0.2, lty=2, col="steelblue", lwd=1.4)
abline(h=1.0, lty=2, col="tomato",    lwd=1.4)
text(max(bp_d)+0.4, 0.22, "0.2", col="steelblue", cex=0.75, adj=0)
text(max(bp_d)+0.4, 1.02, "1.0", col="tomato",    cex=0.75, adj=0)
legend("topleft",
       legend=c(expression(hat(tau)^2 <= 0.2),
                expression(0.2 < hat(tau)^2 <= 1),
                expression(hat(tau)^2 > 1)),
       fill=c("#2166ac","#4dac26","#d73027"),
       border=NA, bty="n", cex=0.78)
mtext("D", side=3, adj=0, line=0.4, font=2, cex=1.0)

dev.off()
message("Saved fig_ad_annotated.pdf")

message("All application figures written to: ", OUT_DIR)
