# =============================================================================
# REVISION C - Remaining items from the editorial comments
#   C1 Figure 1: samples and patients excluded at each step
#   C2 Figure 3: BRCA2 protein domains and HGVS (p.) labels
#   C3 Figure 4: third panel adjusted for continuous mutation burden (A8)
#   C4 OncoKB data version (at the time of download)
#   C5 TCGA: LOH at the BRCA2 locus (ABSOLUTE) and biallelic loss
#   C6 TCGA: sensitivity analysis with curated deleterious BRCA2 mutations (Knijnenburg)
# All numbers are computed in this R session.
# =============================================================================
while (sink.number() > 0) sink()
PROJ <- path.expand("~/GENIE_BRCA2"); setwd(PROJ)
RV <- file.path(PROJ, "revizyon"); TD <- file.path(RV, "tcga")

# GENIE objects and Figure 1-4 code (pt, bv, flow_ms, crude, nh, gene_lab, save_fig, thm, pal)
source("02_manuscript_tables_figures.R", encoding = "UTF-8")
while (sink.number() > 0) sink()
for (p in c("httr", "jsonlite", "patchwork")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(httr); library(jsonlite); library(patchwork)

sink(file.path(RV, "C_log.txt"), split = TRUE)
cat("Revision C date:", format(Sys.time()), "\n")
pfm <- function(p) ifelse(is.na(p), NA, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
save_rv <- function(g, name, w, h) {
  ggsave(file.path(RV, paste0(name, ".png")),  g, width = w, height = h, dpi = 300)
  ggsave(file.path(RV, paste0(name, ".tiff")), g, width = w, height = h, dpi = 300, compression = "lzw")
}

# ---- C1 Figure 1 ------------------------------------------------------------
cat("\n=== C1 FLOW AND EXCLUSIONS ===\n")
excl <- data.table(
  reason = c("Non-adenocarcinoma histology", "Cell-free DNA samples",
             "Additional samples from the same patient", "Panel did not cover both BRCA1 and BRCA2"),
  samples  = -diff(flow_ms$samples),
  patients = -diff(flow_ms$patients))
print(flow_ms); print(excl)
cat("\nExcluded histologies (OncoTree, number of samples):\n")
print(clin[!Oncotree.Code %in% names(ADENO_CODES), .N, by = Oncotree.Code][order(-N)])

fl <- copy(flow_ms)[, y := rev(seq_len(.N)) * 2]
fl[, lab := sprintf("%s\n%s samples, %s patients", step, format(samples, big.mark = ","), format(patients, big.mark = ","))]
fl[.N, lab := sprintf("%s\n%s patients (study cohort)", step, format(patients, big.mark = ","))]
ex <- copy(excl)[, y := fl$y[-nrow(fl)] - 1]
ex[, lab := c(sprintf("Excluded: non-adenocarcinoma histology\n%s samples, %s patients", format(samples[1], big.mark = ","), format(patients[1], big.mark = ",")),
             sprintf("Excluded: cell-free DNA samples\n%s samples, %s patients", samples[2], patients[2]),
             sprintf("Excluded: additional samples\nfrom the same patient\n%s samples", samples[3]),
             sprintf("Excluded: panel did not cover\nboth BRCA1 and BRCA2\n%s patients", patients[4]))]
g1 <- ggplot() +
  geom_segment(data = fl[-.N], aes(x = 0, xend = 0, y = y - 0.45, yend = y - 1.55), arrow = arrow(length = unit(0.18, "cm")), colour = "grey40") +
  geom_segment(data = ex, aes(x = 0, xend = 0.95, y = y, yend = y), arrow = arrow(length = unit(0.15, "cm")), colour = "grey40") +
  geom_label(data = fl, aes(0, y, label = lab), size = 3.1, label.padding = unit(0.45, "lines"), fill = "white") +
  geom_label(data = ex, aes(1.9, y, label = lab), size = 2.8, label.padding = unit(0.4, "lines"), fill = "grey95") +
  xlim(-1.3, 3) + theme_void()
save_rv(g1, "Figure1_flow", 8, 7)

# ---- C2 Figure 3 ------------------------------------------------------------
cat("\n=== C2 LOLLIPOP LABELS (HGVS) ===\n")
lp <- bv[!is.na(pos), .N, by = .(pos, cls, ctx = fifelse(hyper, "Hypermutated tumors", "Non-hypermutated tumors"))]
lp[, ctx := factor(ctx, c("Non-hypermutated tumors", "Hypermutated tumors"))]
lp[, cls := factor(cls, c("Frameshift", "Nonsense", "Splice", "Missense", "Other/ambiguous"))]
top_var <- bv[!is.na(pos), .N, by = .(pos, variant)][order(pos, -N)][, .(hgvs = paste0("p.", variant[1])), by = pos]
lab_lp <- lp[, .(N = sum(N)), by = .(pos, ctx)][N >= 3]
lab_lp <- merge(lab_lp, top_var, by = "pos")
print(lab_lp[order(ctx, -N)])
# Approximate domain boundaries, UniProt P51587
dom <- data.table(name = c("PALB2-binding", "BRC repeats 1-8", "DNA-binding domain", "C-terminal RAD51 binding / NLS"),
                  start = c(10, 1002, 2479, 3265), end = c(40, 2085, 3190, 3330))
cat("\nDomains:\n"); print(dom)
g3top <- ggplot(lp, aes(pos, N)) +
  geom_segment(aes(xend = pos, y = 0, yend = N), colour = "grey65", linewidth = 0.4) +
  geom_point(aes(fill = cls), shape = 21, size = 2.6, colour = "black", stroke = 0.3) +
  geom_text(data = lab_lp, aes(pos, N, label = hgvs), vjust = -0.9, size = 2.5) +
  facet_wrap(~ctx, ncol = 1) +
  scale_x_continuous(limits = c(0, 3418), breaks = seq(0, 3400, 500)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.4))) +
  scale_fill_manual(values = c(Frameshift = "#2C6E9B", Nonsense = "#1B3A57", Splice = "#E08E45",
                               Missense = "#4DAF4A", `Other/ambiguous` = "grey70"), name = "Mutation type") +
  labs(x = NULL, y = "Number of tumors") + thm + theme(axis.text.x = element_blank())
g3dom <- ggplot() +
  annotate("rect", xmin = 1, xmax = 3418, ymin = 0.4, ymax = 0.6, fill = "grey85") +
  geom_rect(data = dom, aes(xmin = start, xmax = end, ymin = 0.2, ymax = 0.8, fill = name), colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("#F2C14E", "#8E7CC3", "#E06666", "#76A5AF"), name = "Domain") +
  scale_x_continuous(limits = c(0, 3418), breaks = seq(0, 3400, 500)) +
  labs(x = "Position in the BRCA2 protein (amino acid 1 to 3418, NP_000050.3)", y = NULL) +
  thm + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(), legend.position = "bottom")
save_rv(g3top / g3dom + plot_layout(heights = c(6, 1)), "Figure3_lollipop", 10, 7)

# ---- C3 Figure 4 ------------------------------------------------------------
cat("\n=== C3 CO-ALTERATIONS, THREE PANELS ===\n")
a8 <- fread(file.path(RV, "A8_koalt_surekli.csv"))
f4 <- rbind(crude[, .(gen, set = "All tumors", OR = crude_or, lo = crude_lo, hi = crude_hi, q = crude_q)],
            nh[, .(gen, set = "Non-hypermutated tumors only", OR, lo = OR_alt, hi = OR_ust, q = q_BH)],
            a8[, .(gen = gene, set = "Adjusted for continuous mutation burden", OR, lo = OR_lo, hi = OR_hi, q)])
f4[, set := factor(set, c("All tumors", "Non-hypermutated tumors only", "Adjusted for continuous mutation burden"))]
f4 <- f4[is.finite(OR) & OR > 0 & is.finite(hi)]
f4[, gen := factor(gene_lab[gen], levels = rev(gene_lab))]
f4[, sig := fifelse(q < 0.05, "q < 0.05", "Not significant")]
print(f4[order(set, gen)])
g4 <- ggplot(f4, aes(OR, gen)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = pmax(lo, 0.05), xmax = pmin(hi, 40)), height = 0.25, colour = "grey30") +
  geom_point(aes(fill = sig), shape = 21, size = 2.8) +
  facet_wrap(~set, labeller = label_wrap_gen(28)) + scale_x_log10(limits = c(0.05, 40), breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 10, 25)) +
  scale_fill_manual(values = c(`q < 0.05` = pal[["red"]], `Not significant` = "white"), name = NULL) +
  labs(x = "Odds ratio (BRCA2-altered vs wild-type), log scale", y = NULL) + thm
save_rv(g4, "Figure4_coalterations", 12, 5)

# ---- C4 OncoKB version ------------------------------------------------------
cat("\n=== C4 OncoKB / cBioPortal VERSION ===\n")
for (u in c("https://www.oncokb.org/api/v1/info", "https://genie.cbioportal.org/api/info")) {
  r <- tryCatch(GET(u, timeout(30)), error = function(e) NULL)
  cat(u, "\n")
  if (!is.null(r) && status_code(r) == 200) print(str(fromJSON(content(r, "text", encoding = "UTF-8")))) else cat("  not reachable\n")
}

# ---- C5-C6 TCGA ------------------------------------------------------------
sink()
source(file.path(PROJ, "05_revision_B2_tcga_validation.R"), encoding = "UTF-8")   # cl
while (sink.number() > 0) sink()
sink(file.path(RV, "C_log.txt"), split = TRUE, append = TRUE)
pfm <- function(p) ifelse(is.na(p), NA, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

cat("\n=== C5 TCGA LOH AT THE BRCA2 LOCUS (ABSOLUTE) ===\n")
segf <- file.path(TD, "TCGA_mastercalls.abs_segtabs.fixed.txt")
if (!file.exists(segf)) GET("https://api.gdc.cancer.gov/data/0f4f5701-7b61-41ae-bda9-2805d1ca9781", write_disk(segf, overwrite = TRUE), timeout(1800))
seg <- fread(segf)
cat("ABSOLUTE segment columns:", paste(names(seg), collapse = ", "), "\n")
seg[, sampleId := substr(Sample, 1, 15)]
seg <- seg[sampleId %in% cl$sampleId & as.character(Chromosome) == "13" & Start <= 32973809 & End >= 32889611]
cat("Samples with an ABSOLUTE segment at the BRCA2 locus:", uniqueN(seg$sampleId), "/", nrow(cl), "\n")
loh_col <- intersect(c("LOH", "loh"), names(seg))
seg[, loh := if (length(loh_col)) get(loh_col) == 1 else FALSE]
if ("Modal_HSCN_1" %in% names(seg)) seg[, loh := loh | Modal_HSCN_1 == 0]
if ("Modal_Total_CN" %in% names(seg)) seg[, homdel_abs := Modal_Total_CN == 0] else seg[, homdel_abs := FALSE]
sl <- seg[, .(BRCA2_LOH = any(loh, na.rm = TRUE), BRCA2_homdel_abs = any(homdel_abs, na.rm = TRUE)), by = sampleId]
cl <- merge(cl, sl, by = "sampleId", all.x = TRUE)
cl[, abs_avail := !is.na(BRCA2_LOH)]

# C6 Curated deleterious mutations (Knijnenburg)
gn <- trimws(readLines(file.path(TD, "ddr", "Genes.tsv")))
smp <- fread(file.path(TD, "ddr", "Samples.tsv"), header = FALSE, col.names = c("sampleId", "ct"))
gm <- fread(file.path(TD, "ddr", "GeneMutations.tsv"), header = FALSE)
kn <- data.table(sampleId = smp$sampleId, KN_mut = as.numeric(unlist(gm[which(gn == "BRCA2")])))
cl <- merge(cl, kn, by = "sampleId", all.x = TRUE)
cl[, BRCA2_curated := KN_mut %in% 1]

cl[, biallelic := abs_avail & ((BRCA2_lof & (BRCA2_LOH | n_trunc >= 2)) | cna %in% -2 | BRCA2_homdel_abs)]
cl[, hgrp := fifelse(HYPER_true, "MSI or POLE", "Microsatellite stable")]

c5 <- cl[BRCA2_any == TRUE, .(
  BRCA2_altered = .N, ABSOLUTE_available = sum(abs_avail),
  LOH_at_locus = sprintf("%d/%d", sum(BRCA2_LOH, na.rm = TRUE), sum(abs_avail)),
  putative_LoF = sum(BRCA2_lof),
  LoF_with_LOH = sum(BRCA2_lof & BRCA2_LOH, na.rm = TRUE),
  biallelic = sprintf("%d/%d", sum(biallelic), sum(abs_avail)),
  curated_deleterious = sum(BRCA2_curated)), by = hgrp]
print(c5)
pl <- cl[BRCA2_any == TRUE & abs_avail]
p_loh <- fisher.test(table(pl$hgrp, pl$BRCA2_LOH))$p.value
p_bi <- fisher.test(table(pl$hgrp, pl$biallelic))$p.value
cat("LOH, MSI/POLE vs MSS (BRCA2-altered) Fisher P =", pfm(p_loh), "\n")
cat("Biallelic, MSI/POLE vs MSS (BRCA2-altered) Fisher P =", pfm(p_bi), "\n")
cat("\nLOH at the locus in BRCA2 wild-type tumors:\n")
print(cl[BRCA2_any == FALSE & abs_avail, .(N = .N, LOH = sum(BRCA2_LOH), pct = round(100 * mean(BRCA2_LOH), 1)), by = hgrp])

cat("\nTumors with likely biallelic loss:\n")
print(cl[biallelic == TRUE, .(sampleId, SUBTYPE, n_mut, n_trunc, cna, BRCA2_LOH, HRD_Score, HRD_high = HRD_Score >= 42, mutSig3)])

cat("\nHRD score: biallelic vs other BRCA2-altered vs WT (MSS only):\n")
cl[, bgrp := factor(fifelse(!BRCA2_any, "Wild-type", fifelse(biallelic, "Biallelic", "Other BRCA2")), c("Wild-type", "Other BRCA2", "Biallelic"))]
hb <- cl[HYPER_true == FALSE, .(N = .N, HRD_n = sum(!is.na(HRD_Score)), HRD_score = mq(HRD_Score),
                                HRD_high = sprintf("%d/%d", sum(HRD_Score >= 42, na.rm = TRUE), sum(!is.na(HRD_Score)))), by = bgrp][order(bgrp)]
print(hb)
x1 <- cl[HYPER_true == FALSE & bgrp == "Biallelic", HRD_Score]; x0 <- cl[HYPER_true == FALSE & bgrp == "Wild-type", HRD_Score]
p_hb <- if (sum(!is.na(x1)) >= 2) wilcox.test(x1, x0)$p.value else NA
cat("MSS biallelic vs WT HRD P =", pfm(p_hb), "\n")
cat("\nHRD score: biallelic vs other (all tumors):\n")
print(cl[BRCA2_any == TRUE, .(N = .N, HRD_score = mq(HRD_Score)), by = .(hgrp, biallelic)][order(hgrp, biallelic)])

cat("\n=== C6 CURATED DELETERIOUS BRCA2 (Knijnenburg 2018) ===\n")
c6 <- cbind(analysis = "Curated deleterious BRCA2 vs MSI or POLE", tab_or(cl$BRCA2_curated, cl$HYPER_true))
print(c6)
cg <- cl[, .(N = .N, curated = sum(BRCA2_curated), HRD_score = mq(HRD_Score[BRCA2_curated]),
             HRD_high = sprintf("%d/%d", sum(HRD_Score[BRCA2_curated] >= 42, na.rm = TRUE), sum(!is.na(HRD_Score[BRCA2_curated])))), by = hgrp]
print(cg)
cat("Tumors with a curated deleterious BRCA2 mutation:\n")
print(cl[BRCA2_curated == TRUE, .(sampleId, SUBTYPE, n_mut, n_trunc, cna, BRCA2_LOH, biallelic, HRD_Score, mutSig3)])

fwrite(excl, file.path(RV, "C1_cikarilanlar.csv"))
fwrite(c5, file.path(RV, "C5_tcga_biallelik.csv")); fwrite(hb, file.path(RV, "C5_tcga_biallelik_hrd.csv"))
fwrite(c6, file.path(RV, "C6_kurasyonlu_brca2.csv")); fwrite(cg, file.path(RV, "C6_kurasyonlu_hrd.csv"))
cat("\nDone.\n")
sink()
