# =============================================================================
# 09 - TCGA validation with the pathogenic BRCA2 definition only, and Figure 1
#   E1 Pathogenic BRCA2 (curated deleterious mutation or deep deletion) by MSI/HM-SNV
#   E2 HRD score and locus LOH: pathogenic vs no BRCA2 variant, within each stratum
#   E3 Patient-level table of pathogenic carriers
#   E4 Figure 1 (GENIE flow, black and white) and Figure S2 (TCGA flow and HRD)
# VUS-only tumors are kept out of the comparisons and only counted.
# =============================================================================
while (sink.number() > 0) sink()
PROJ <- path.expand("~/GENIE_BRCA2"); setwd(PROJ)
source("05_revision_B2_tcga_validation.R", encoding = "UTF-8")   # gives cl, mu, cn, st, rd_w
while (sink.number() > 0) sink()
library(data.table); library(ggplot2); library(patchwork)
RV <- file.path(PROJ, "revizyon"); TD <- file.path(RV, "tcga"); ED <- file.path(RV, "E"); dir.create(ED, showWarnings = FALSE)
sink(file.path(ED, "E_log.txt"), split = TRUE)
cat("Script 09 date:", format(Sys.time()), "\n")
pfm <- function(p) ifelse(is.na(p), NA, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
mq  <- function(x) sprintf("%.1f (%.1f-%.1f)", median(x, na.rm = TRUE), quantile(x, .25, na.rm = TRUE), quantile(x, .75, na.rm = TRUE))

# ---- Pathogenic definition --------------------------------------------------
gn  <- trimws(readLines(file.path(TD, "ddr", "Genes.tsv")))
smp <- fread(file.path(TD, "ddr", "Samples.tsv"), header = FALSE, col.names = c("sampleId", "ct"))
gm  <- fread(file.path(TD, "ddr", "GeneMutations.tsv"), header = FALSE)
cl  <- merge(cl, data.table(sampleId = smp$sampleId, KN = as.numeric(unlist(gm[which(gn == "BRCA2")]))), by = "sampleId", all.x = TRUE)
cl[, pathogenic := KN %in% 1 | cna %in% -2]
cl[, b2grp := factor(fifelse(pathogenic, "Pathogenic", fifelse(BRCA2_any, "VUS only", "No BRCA2 variant")),
                     levels = c("No BRCA2 variant", "VUS only", "Pathogenic"))]
cl[, strat := factor(fifelse(HYPER_true, "MSI/HM-SNV", "Non-MSI/HM-SNV"), levels = c("MSI/HM-SNV", "Non-MSI/HM-SNV"))]

# ---- Locus LOH (ABSOLUTE segment at the variant, or any BRCA2 segment) ------
seg <- fread(file.path(TD, "TCGA_mastercalls.abs_segtabs.fixed.txt"))
seg[, sampleId := substr(Sample, 1, 15)]
seg <- seg[sampleId %in% cl$sampleId & as.character(Chromosome) == "13" & Start <= 32973809 & End >= 32889611]
seg[, loh := (LOH == 1) | (Modal_HSCN_1 == 0)]
anyl <- seg[, .(any_loh = any(loh, na.rm = TRUE)), by = sampleId]
mv <- merge(mu[gene == "BRCA2" & mutationType != "Silent", .(sampleId, startPosition)], seg[, .(sampleId, Start, End, loh)], by = "sampleId", allow.cartesian = TRUE)
mv <- mv[startPosition >= Start & startPosition <= End, .(var_loh = any(loh, na.rm = TRUE)), by = sampleId]
cl <- merge(merge(cl, anyl, by = "sampleId", all.x = TRUE), mv, by = "sampleId", all.x = TRUE)
cl[, abs_avail := !is.na(any_loh)]
cl[, locus_loh := fifelse(!abs_avail, NA, fifelse(n_mut > 0 & !is.na(var_loh), var_loh, any_loh))]

# ---- E1 Frequency -----------------------------------------------------------
cat("\n=== E1 BRCA2 GROUPS BY STRATUM ===\n"); print(dcast(cl[, .N, by = .(b2grp, strat)], b2grp ~ strat, value.var = "N"))
t1 <- table(cl$HYPER_true, cl$pathogenic); f1 <- fisher.test(t1)
e1 <- data.table(cohort_N = nrow(cl), msi_n = sum(cl$HYPER_true),
  pathogenic = sprintf("%d (%.2f%%)", sum(cl$pathogenic), 100 * mean(cl$pathogenic)),
  in_msi = sprintf("%d/%d (%.1f%%)", t1[2, 2], sum(t1[2, ]), 100 * t1[2, 2] / sum(t1[2, ])),
  in_non_msi = sprintf("%d/%d (%.1f%%)", t1[1, 2], sum(t1[1, ]), 100 * t1[1, 2] / sum(t1[1, ])),
  carriers_in_msi = sprintf("%.1f%%", 100 * t1[2, 2] / sum(t1[, 2])),
  OR = sprintf("%.2f (%.2f-%.2f)", f1$estimate, f1$conf.int[1], f1$conf.int[2]), P = pfm(f1$p.value),
  pathogenic_truncating_or_deepdel = sum(cl$pathogenic & (cl$n_trunc > 0 | cl$cna %in% -2)),
  pathogenic_deep_deletion = sum(cl$pathogenic & cl$cna %in% -2),
  pathogenic_in_proxy_hyper = sprintf("%d/%d", sum(cl$pathogenic & cl$proxy), sum(cl$pathogenic)))
print(t(e1))

# ---- E2 HRD and locus LOH ---------------------------------------------------
cat("\n=== E2 HRD SCORE AND LOCUS LOH ===\n")
e2 <- cl[, .(N = .N, HRD_n = sum(!is.na(HRD_Score)), HRD = mq(HRD_Score),
             HRD42 = sprintf("%d/%d", sum(HRD_Score >= 42, na.rm = TRUE), sum(!is.na(HRD_Score))),
             locus_LOH = sprintf("%d/%d", sum(locus_loh, na.rm = TRUE), sum(!is.na(locus_loh)))), by = .(strat, b2grp)][order(strat, b2grp)]
print(e2)
wp <- function(s) { x <- cl[strat == s & b2grp == "Pathogenic", HRD_Score]; y <- cl[strat == s & b2grp == "No BRCA2 variant", HRD_Score]
  if (sum(!is.na(x)) < 2) NA else wilcox.test(x, y, exact = FALSE)$p.value }
lp <- function(s) { d <- cl[strat == s & b2grp != "VUS only" & !is.na(locus_loh)]; tryCatch(fisher.test(table(factor(d$b2grp == "Pathogenic", c(FALSE, TRUE)), factor(d$locus_loh, c(FALSE, TRUE))))$p.value, error = function(e) NA) }
e2p <- data.table(stratum = levels(cl$strat), P_HRD_pathogenic_vs_none = pfm(sapply(levels(cl$strat), wp)),
                  P_LOH_pathogenic_vs_none = pfm(sapply(levels(cl$strat), lp)))
print(e2p)
cat("Pathogenic MSI vs pathogenic non-MSI HRD, Wilcoxon P =",
    pfm(wilcox.test(cl[b2grp == "Pathogenic" & HYPER_true, HRD_Score], cl[b2grp == "Pathogenic" & !HYPER_true, HRD_Score], exact = FALSE)$p.value), "\n")

# ---- E3 Patient-level table -------------------------------------------------
lv <- mu[gene == "BRCA2" & mutationType != "Silent", .(variants = paste(proteinChange, collapse = ", ")), by = sampleId]
e3 <- merge(cl[pathogenic == TRUE, .(sampleId, site = fifelse(study == st[1], "Stomach", "Esophagus"),
  subtype = fifelse(is.na(SUBTYPE), "Missing", sub("^POLE$", "HM-SNV", sub("^(STAD|ESCA)_", "", SUBTYPE))),
  MANTIS = round(MSI_SCORE_MANTIS, 2), TMB = round(TMB_NONSYNONYMOUS, 1),
  event = fifelse(cna %in% -2, "Deep deletion", fifelse(n_trunc > 0, "Truncating", "Curated missense or other")),
  locus_LOH = locus_loh, HRD = HRD_Score)], lv, by = "sampleId", all.x = TRUE)
setorder(e3, subtype, sampleId); cat("\n=== E3 PATHOGENIC CARRIERS (n =", nrow(e3), ") ===\n"); print(e3)

# ---- E4 Figures -------------------------------------------------------------
box <- function(d) list(geom_rect(data = d, aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - 0.38, ymax = y + 0.38), fill = "white", colour = "black", linewidth = 0.4),
                        geom_text(data = d, aes(x, y, label = lab), size = 3, lineheight = 0.95))
arr <- arrow(length = unit(0.15, "cm"), type = "closed")
# Figure 1, GENIE
fw <- fread(file.path(PROJ, "makale", "Flow.csv")); fm <- function(v) format(v, big.mark = ",")
mainb <- data.table(x = 0, w = 2.6, y = c(9, 7, 5, 3, 1), lab = c(
  sprintf("Esophagogastric cancer samples\nin AACR Project GENIE v20.0\n%s samples, %s patients", fm(fw$samples[1]), fm(fw$patients[1])),
  sprintf("Adenocarcinoma histology\n%s samples, %s patients", fm(fw$samples[2]), fm(fw$patients[2])),
  sprintf("Tumor tissue samples\n%s samples, %s patients", fm(fw$samples[3]), fm(fw$patients[3])),
  sprintf("One sample per patient\n%s patients", fm(fw$patients[4])),
  sprintf("Final study cohort\n%s patients", fm(fw$patients[5]))))
exb <- data.table(x = 2.9, w = 2.4, y = c(8, 6, 4, 2), lab = c(
  sprintf("Excluded: non-adenocarcinoma\n%s samples, %s patients", fm(fw$samples[1] - fw$samples[2]), fm(fw$patients[1] - fw$patients[2])),
  sprintf("Excluded: cell-free DNA\n%s samples, %s patients", fw$samples[2] - fw$samples[3], fw$patients[2] - fw$patients[3]),
  sprintf("Excluded: additional samples\nfrom the same patient\n%s samples", fw$samples[3] - fw$samples[4]),
  sprintf("Excluded: panel did not cover\nboth BRCA1 and BRCA2\n%s patients", fw$patients[4] - fw$patients[5])))
g1 <- ggplot() + geom_segment(data = mainb[-.N], aes(x = 0, xend = 0, y = y - 0.38, yend = y - 1.62), arrow = arr) +
  geom_segment(data = exb, aes(x = 0, xend = x - w / 2, y = y, yend = y), arrow = arr) + box(mainb) + box(exb) +
  coord_cartesian(xlim = c(-1.4, 4.2), ylim = c(0.5, 9.5)) + theme_void()
ggsave(file.path(ED, "Figure1_flow.png"), g1, width = 7.5, height = 8, dpi = 600, bg = "white")
ggsave(file.path(ED, "Figure1_flow.tiff"), g1, width = 7.5, height = 8, dpi = 600, compression = "lzw", bg = "white")
# Figure S2, TCGA flow and HRD
cl0 <- rbindlist(lapply(st, function(s) { a <- rd_w(file.path(TD, paste0(s, "_clin_sample_long.csv")), "sampleId + patientId")
  b <- rd_w(file.path(TD, paste0(s, "_clin_patient_long.csv")), "patientId"); merge(a, b, by = "patientId", all.x = TRUE) }), fill = TRUE)
is_ad <- grepl("Adenocarcinoma|Signet Ring", cl0$CANCER_TYPE_DETAILED) & !grepl("ESCC", fifelse(is.na(cl0$SUBTYPE), "", cl0$SUBTYPE))
n_nomut <- sum(is_ad & is.na(suppressWarnings(as.numeric(cl0$MUTATION_COUNT))))
cat("\nTCGA flow:", nrow(cl0), "samples; non-adenocarcinoma", sum(!is_ad), "; no mutation data", n_nomut, "; cohort", nrow(cl),
    "; ABSOLUTE at BRCA2", sum(cl$abs_avail), "; HRD score", sum(!is.na(cl$HRD_Score)), "\n")
tb <- data.table(x = 0, w = 2.6, y = c(7, 5, 3), lab = c(sprintf("TCGA PanCancer Atlas stomach\nand esophageal samples\n%d samples", nrow(cl0)),
  sprintf("Adenocarcinomas with\nmutation data\n%d tumors", nrow(cl)),
  sprintf("With ABSOLUTE copy number: %d\nWith HRD score: %d", sum(cl$abs_avail), sum(!is.na(cl$HRD_Score)))))
te <- data.table(x = 2.9, w = 2.6, y = 6, lab = sprintf("Excluded:\nnot adenocarcinoma, %d samples\nno mutation data, %d samples", sum(!is_ad), n_nomut))
ga <- ggplot() + geom_segment(data = tb[-.N], aes(x = 0, xend = 0, y = y - 0.38, yend = y - 1.62), arrow = arr) +
  geom_segment(data = te, aes(x = 0, xend = x - w / 2, y = y, yend = y), arrow = arr) + box(tb) + box(te) +
  coord_cartesian(xlim = c(-1.4, 4.2), ylim = c(2.5, 7.5)) + labs(title = "A") + theme_void() + theme(plot.title = element_text(face = "bold"))
hd <- cl[!is.na(HRD_Score) & b2grp != "VUS only"]
nl <- hd[, .N, by = .(strat, b2grp)]
hd <- merge(hd, nl, by = c("strat", "b2grp")); hd[, xl := factor(sprintf("%s\n(n = %d)", fifelse(b2grp == "Pathogenic", "Pathogenic BRCA2", "No BRCA2 variant"), N))]
gb <- ggplot(hd, aes(xl, HRD_Score)) + geom_boxplot(outlier.shape = NA, width = 0.55, fill = "grey90") +
  geom_jitter(width = 0.12, size = 0.9, alpha = 0.6) + geom_hline(yintercept = 42, linetype = "dashed") +
  facet_wrap(~strat, scales = "free_x") + labs(x = NULL, y = "HRD score", title = "B") +
  theme_bw(base_size = 10) + theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
ggsave(file.path(ED, "FigureS2_TCGA.png"), ga | gb, width = 11, height = 4.8, dpi = 300, bg = "white")
ggsave(file.path(ED, "FigureS2_TCGA.tiff"), ga | gb, width = 11, height = 4.8, dpi = 300, compression = "lzw", bg = "white")

fwrite(e1, file.path(ED, "E1_frequency.csv")); fwrite(e2, file.path(ED, "E2_hrd_loh.csv")); fwrite(e2p, file.path(ED, "E2_tests.csv"))
fwrite(e3, file.path(ED, "E3_patient_level.csv"))
cat("\nDone.\n"); sink()
