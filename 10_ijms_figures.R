# =============================================================================
# 10 - Final figures for the IJMS submission
#   Runs scripts 07 and 08 (which also run 02, 03 and 05) to rebuild all objects,
#   then draws Figure 1 (GENIE flow), redraws Figures 3-5 and S1 with labels matching
#   the manuscript text and copies Figure 2. The final article uses GENIE data only.
#   Figure numbers follow the IJMS manuscript:
#     1 GENIE flow, 2 frequency, 3 sensitivity, 4 lollipop, 5 co-alterations, S1 comparison genes
# =============================================================================
while (sink.number() > 0) sink()
PROJ <- path.expand("~/GENIE_BRCA2"); setwd(PROJ)
source("07_revision_D_second_review.R", encoding = "UTF-8")
source("08_revision_D_figure_fixes.R", encoding = "UTF-8")
while (sink.number() > 0) sink()
FI <- file.path(PROJ, "revizyon", "F_ijms"); dir.create(FI, showWarnings = FALSE)
unlink(file.path(FI, c("FigureS2_TCGA.png", "FigureS2_TCGA.tiff")))   # TCGA is not part of the final article
save_f <- function(g, name, w, h, dpi = 600) {
  ggsave(file.path(FI, paste0(name, ".png")),  g, width = w, height = h, dpi = dpi, bg = "white")
  ggsave(file.path(FI, paste0(name, ".tiff")), g, width = w, height = h, dpi = dpi, compression = "lzw", bg = "white")
}

# Figure 3: sensitivity and robustness (was D Figure 5)
fp[, lab := as.character(lab)]
fp[lab == "Metastasis samples", lab := "Metastatic samples"]
fp[, lab := factor(lab, levels = rev(unique(lab)))]
g5f <- g5 %+% fp
save_f(g5f, "Figure3_sensitivity", 10, 9)

# Figure 4: lollipop (was D Figure 3); legend text "at or near", labels kept apart
rlab <- sprintf("Frameshift at or near a mononucleotide repeat (≥%d nt)", REP)
lp[, rep := factor(fifelse(rep == "Other", "Other", rlab), c(rlab, "Other"))]
lab_lp[, hgvs := sub("_splice$", " splice", hgvs)]
setorder(lab_lp, ctx, pos)
# when two labels in a panel are close, put the first to the left and the second to the right
lab_lp[, `:=`(xl = as.numeric(pos), hj = 0.5)]
for (k in unique(lab_lp$ctx)) {
  ix <- which(lab_lp$ctx == k)
  for (j in seq_along(ix)[-1]) {
    a <- ix[j - 1]; b <- ix[j]
    if (lab_lp$pos[b] - lab_lp$pos[a] < 250 && lab_lp$hj[a] == 0.5) {
      lab_lp[a, `:=`(xl = pos + 15, hj = 1)]; lab_lp[b, `:=`(xl = pos - 15, hj = 0)]
    }
  }
}
g3f <- g3top %+% lp
g3f$layers[[3]] <- geom_text(data = lab_lp, aes(xl, N, label = hgvs, hjust = hj), vjust = -0.9, size = 2.5)
save_f(g3f / g3dom + plot_layout(heights = c(6, 1)), "Figure4_lollipop", 10, 7.5)

# Figure 5: co-alterations (was D Figure 4); clean axis labels, q with three decimals as in Table S10
g4f <- g4 + scale_x_log10(limits = c(0.05, 40), breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 10, 25),
                          labels = c("0.1", "0.25", "0.5", "1", "2", "4", "10", "25")) +
  labs(x = "Odds ratio (BRCA2-altered vs unaltered), log scale")
t4[, txt := g_co8[, sprintf("%.2f (%.2f-%.2f)   q = %.3f", OR, OR_lo, OR_hi, q)]]   # t4 was built from g_co8 in the same row order
g4tf <- ggplot(t4, aes(0, gen, label = txt)) + geom_text(hjust = 0, size = 3) + facet_wrap(~hdr) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0.02)) + labs(x = " ", y = NULL) + thm +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), legend.position = "none")
save_f(g4f + g4tf + plot_layout(widths = c(4, 1.25)), "Figure5_coalterations", 14, 5)

# Figure S1: comparison genes, ratio with one decimal
gs1f <- gs1
gs1f$layers[[3]] <- geom_text(data = g_cmp, aes(gene, pmax(hyper_pct, nonhyper_pct) + 7, label = sprintf("%.1fx", as.numeric(ratio))), inherit.aes = FALSE, size = 3.2)
save_f(gs1f, "FigureS1_comparison_genes", 9, 5, dpi = 300)

# Figure 1: GENIE flow (numbers read from the flow table written by script 02)
fw <- fread(file.path(PROJ, "makale", "Flow.csv")); fm <- function(v) format(v, big.mark = ",")
mb <- data.table(y = c(9, 7, 5, 3, 1), fill = c(rep("white", 4), "grey93"), face = c(rep("plain", 4), "bold"), lab = c(
  sprintf("Esophagogastric cancer in AACR Project GENIE v20.0\n%s samples; %s patients", fm(fw$samples[1]), fm(fw$patients[1])),
  sprintf("Adenocarcinoma histology\n%s samples; %s patients", fm(fw$samples[2]), fm(fw$patients[2])),
  sprintf("Tumor tissue samples\n%s samples; %s patients", fm(fw$samples[3]), fm(fw$patients[3])),
  sprintf("One sample per patient\n%s samples; %s patients", fm(fw$samples[4]), fm(fw$patients[4])),
  sprintf("Final study cohort\n%s patients", fm(fw$patients[5]))))
eb <- data.table(y = c(8, 6, 4, 2), lab = c(
  sprintf("Excluded: non-adenocarcinoma histology\n%s samples; %s patients", fm(fw$samples[1] - fw$samples[2]), fm(fw$patients[1] - fw$patients[2])),
  sprintf("Excluded: cell-free DNA samples\n%s samples; %s patients", fw$samples[2] - fw$samples[3], fw$patients[2] - fw$patients[3]),
  sprintf("Excluded: additional samples\nfrom the same patient\n%s samples", fw$samples[3] - fw$samples[4]),
  sprintf("Excluded: panel did not cover\nboth BRCA1 and BRCA2\n%s patients", fw$patients[4] - fw$patients[5])))
ar <- arrow(length = unit(0.18, "cm"), type = "closed")
g1 <- ggplot() +
  geom_segment(data = mb[-.N], aes(x = 2.7, xend = 2.7, y = y - 0.48, yend = y - 1.52), arrow = ar, linewidth = 0.5) +
  geom_segment(data = eb, aes(x = 2.7, xend = 5.75, y = y, yend = y), arrow = ar, linewidth = 0.5) +
  geom_rect(data = mb, aes(xmin = 0.1, xmax = 5.3, ymin = y - 0.48, ymax = y + 0.48, fill = fill), colour = "black", linewidth = 0.6) +
  geom_rect(data = eb, aes(xmin = 5.8, xmax = 10.1, ymin = y - 0.55, ymax = y + 0.55), fill = "white", colour = "black", linewidth = 0.6) +
  geom_text(data = mb, aes(2.7, y, label = lab, fontface = face), size = 4, lineheight = 1.05) +
  geom_text(data = eb, aes(7.95, y, label = lab), size = 3.7, lineheight = 1.05) +
  scale_fill_identity() + coord_cartesian(xlim = c(0, 10.2), ylim = c(0.4, 9.6), expand = FALSE) + theme_void()
save_f(g1, "Figure1_flow", 10, 9)

# Figure 2: same plot objects as script 02 (g2a, g2b), saved at 600 dpi like the other figures
save_f(g2a + g2b + plot_layout(widths = c(1.6, 1)), "Figure2_frequency", 11, 5)
writeLines(c(format(Sys.time()), capture.output(sessionInfo())), file.path(FI, "F_done.txt"))
cat("IJMS figures done\n")
