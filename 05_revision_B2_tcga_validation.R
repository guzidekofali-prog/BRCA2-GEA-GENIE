# =============================================================================
# REVISION B2 - TCGA STAD + ESCA adenocarcinoma validation cohort
#   B2.1 Cohort: TCGA PanCancer Atlas adenocarcinomas with mutation data
#   B2.2 BRCA2 frequency by true MSI status (TCGA molecular subtype)
#   B2.3 HRD score, LOH fraction, signature 3 and second hits in four BRCA2 x MSI groups
#   B2.4 Accuracy of the GENIE hypermutation proxy against true MSI
# =============================================================================
while (sink.number() > 0) sink()
library(data.table); library(ggplot2)
if (!requireNamespace("logistf", quietly = TRUE)) install.packages("logistf"); library(logistf)
TD <- path.expand("~/GENIE_BRCA2/revizyon/tcga"); RV <- path.expand("~/GENIE_BRCA2/revizyon")
sink(file.path(RV, "B2_log.txt"), split = TRUE)
cat("TCGA validation date:", format(Sys.time()), "\n")
pfm <- function(p) ifelse(is.na(p), NA, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
mq <- function(x) sprintf("%.1f (%.1f-%.1f)", median(x, na.rm = TRUE), quantile(x, .25, na.rm = TRUE), quantile(x, .75, na.rm = TRUE))

# ---- B2.1 Cohort -----------------------------------------------------------
st <- c("stad_tcga_pan_can_atlas_2018", "esca_tcga_pan_can_atlas_2018")
rd_w <- function(f, id) { d <- fread(f, colClasses = "character"); dcast(d, as.formula(paste(id, "~ clinicalAttributeId")), value.var = "value") }
cl <- rbindlist(lapply(st, function(s) {
  cs <- rd_w(file.path(TD, paste0(s, "_clin_sample_long.csv")), "sampleId + patientId")
  cp <- rd_w(file.path(TD, paste0(s, "_clin_patient_long.csv")), "patientId")
  x <- merge(cs, cp, by = "patientId", all.x = TRUE); x[, study := s]; x }), fill = TRUE)
mu <- rbindlist(lapply(st, function(s) fread(file.path(TD, paste0(s, "_mutations.csv")))))
cn <- rbindlist(lapply(st, function(s) fread(file.path(TD, paste0(s, "_cna.csv")))))
num <- c("MSI_SCORE_MANTIS", "MSI_SENSOR_SCORE", "TMB_NONSYNONYMOUS", "MUTATION_COUNT", "FRACTION_GENOME_ALTERED", "ANEUPLOIDY_SCORE", "AGE")
for (v in num) cl[, (v) := suppressWarnings(as.numeric(get(v)))]
n0 <- nrow(cl)
cl <- cl[grepl("Adenocarcinoma|Signet Ring", CANCER_TYPE_DETAILED) & !grepl("ESCC", fifelse(is.na(SUBTYPE), "", SUBTYPE))]
cl <- cl[!is.na(MUTATION_COUNT)]
cat("TCGA STAD+ESCA samples:", n0, " -> adenocarcinomas with mutation data:", nrow(cl), "\n")
print(cl[, .N, by = .(study, SUBTYPE)][order(study, SUBTYPE)])

# Variables
cl[, MSI := grepl("MSI", SUBTYPE) | (is.na(SUBTYPE) & MSI_SCORE_MANTIS > 0.4 & !is.na(MSI_SCORE_MANTIS))]
cat("Subtype missing, called MSI by MANTIS >0.4:", cl[is.na(SUBTYPE) & MSI_SCORE_MANTIS > 0.4, .N], "\n")
cl[, POLE_sub := grepl("POLE", SUBTYPE)]
cl[, HYPER_true := MSI | POLE_sub]
cl[, MSI_mantis := MSI_SCORE_MANTIS > 0.4]
trunc_types <- c("Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Nonstop_Mutation")
gstat <- function(g) {
  m <- mu[gene == g & mutationType != "Silent", .(n_mut = .N, n_trunc = sum(mutationType %in% trunc_types)), by = sampleId]
  c2 <- cn[gene == g, .(cna = alteration[1]), by = sampleId]
  x <- merge(cl[, .(sampleId)], m, by = "sampleId", all.x = TRUE); x <- merge(x, c2, by = "sampleId", all.x = TRUE)
  x[is.na(n_mut), n_mut := 0L]; x[is.na(n_trunc), n_trunc := 0L]; x
}
b2 <- gstat("BRCA2")
cl <- merge(cl, b2, by = "sampleId", all.x = TRUE)
cl[, BRCA2_any := n_mut > 0 | cna %in% -2]
cl[, BRCA2_lof := n_trunc > 0 | cna %in% -2]
cl[, BRCA2_2hit := BRCA2_any & (n_mut >= 2 | cna %in% c(-1, -2))]
for (g in c("MLH1", "MSH2", "MSH6", "PMS2")) { x <- gstat(g); cl[x, on = "sampleId", (paste0(g, "_lof")) := (i.n_trunc > 0 | i.cna %in% -2)] }
POLE_HOT <- c("P286R", "V411L", "S459F", "A456P", "S297F", "F367S", "L424I", "M444K", "P436R", "M295R")
pole <- mu[gene == "POLE" & proteinChange %in% POLE_HOT, unique(sampleId)]
cl[, MMR_lof := MLH1_lof | MSH2_lof | MSH6_lof | PMS2_lof]; cl[is.na(MMR_lof), MMR_lof := FALSE]
cl[, POLE_hot := sampleId %in% pole]

# DDR scores (Knijnenburg 2018)
sc <- fread(file.path(TD, "ddr", "DDRscores.tsv"), header = FALSE)
setnames(sc, trimws(readLines(file.path(TD, "ddr", "Scores.tsv"))))
smp <- fread(file.path(TD, "ddr", "Samples.tsv"), header = FALSE, col.names = c("sampleId", "ct"))
sc <- cbind(smp, sc)
cl <- merge(cl, sc[, .(sampleId, HRD_Score, HRD_LOH, HRD_TAI, HRD_LST, LOH_frac_altered, mutSig3, purity)], by = "sampleId", all.x = TRUE)
cat("Samples with an HRD score:", sum(!is.na(cl$HRD_Score)), "/", nrow(cl), "\n")

# ---- B2.2 BRCA2 by MSI ------------------------------------------------------
tab_or <- function(y, x) { t <- table(factor(x, c(FALSE, TRUE)), factor(y, c(FALSE, TRUE))); f <- fisher.test(t)
  data.table(MSI_or_hyper = sprintf("%d/%d (%.1f%%)", t[2, 2], sum(t[2, ]), 100 * t[2, 2] / sum(t[2, ])),
             MSS = sprintf("%d/%d (%.1f%%)", t[1, 2], sum(t[1, ]), 100 * t[1, 2] / sum(t[1, ])),
             altered_in_MSI_group = sprintf("%.1f%%", 100 * t[2, 2] / sum(t[, 2])),
             OR = f$estimate, lo = f$conf.int[1], hi = f$conf.int[2], P = pfm(f$p.value)) }
b22 <- rbind(cbind(analysis = "Any BRCA2, true MSI (TCGA subtype)", tab_or(cl$BRCA2_any, cl$MSI)),
             cbind(analysis = "Any BRCA2, MSI or POLE subtype", tab_or(cl$BRCA2_any, cl$HYPER_true)),
             cbind(analysis = "Putative LoF BRCA2, true MSI", tab_or(cl$BRCA2_lof, cl$MSI)),
             cbind(analysis = "Any BRCA2, MANTIS >0.4", tab_or(cl$BRCA2_any[!is.na(cl$MSI_mantis)], cl$MSI_mantis[!is.na(cl$MSI_mantis)])))
cat("\n=== B2.2 TCGA: BRCA2 and true MSI ===\n"); print(b22)

# ---- B2.3 Four groups: HRD score, LOH, signature 3, second hits ---------------------
cl[, grp := factor(fifelse(BRCA2_any, fifelse(HYPER_true, "BRCA2-altered / MSI or POLE", "BRCA2-altered / MSS"),
                           fifelse(HYPER_true, "BRCA2 wild-type / MSI or POLE", "BRCA2 wild-type / MSS")),
                   levels = c("BRCA2 wild-type / MSS", "BRCA2-altered / MSS", "BRCA2 wild-type / MSI or POLE", "BRCA2-altered / MSI or POLE"))]
cl[, HRD_high := HRD_Score >= 42]
g4 <- cl[, .(N = .N, HRD_score = mq(HRD_Score), HRD_high_42 = sprintf("%d/%d (%.0f%%)", sum(HRD_high, na.rm = TRUE), sum(!is.na(HRD_high)), 100 * mean(HRD_high, na.rm = TRUE)),
             LOH_fraction = mq(LOH_frac_altered), signature3 = mq(mutSig3),
             truncating_or_deepdel = sum(BRCA2_lof), second_hit = sum(BRCA2_2hit)), by = grp][order(grp)]
cat("\n=== B2.3 FOUR GROUPS ===\n"); print(g4)
pw <- function(v, a, b) { x <- cl[grp == a][[v]]; y <- cl[grp == b][[v]]; if (sum(!is.na(x)) < 2) return(NA); wilcox.test(x, y)$p.value }
cmpt <- rbindlist(lapply(c("HRD_Score", "LOH_frac_altered", "mutSig3"), function(v) data.table(measure = v,
  P_altered_vs_WT_in_MSS = pfm(pw(v, "BRCA2-altered / MSS", "BRCA2 wild-type / MSS")),
  P_altered_vs_WT_in_MSI = pfm(pw(v, "BRCA2-altered / MSI or POLE", "BRCA2 wild-type / MSI or POLE")),
  P_altered_MSS_vs_altered_MSI = pfm(pw(v, "BRCA2-altered / MSS", "BRCA2-altered / MSI or POLE")),
  P_kruskal_4groups = pfm(kruskal.test(cl[[v]] ~ cl$grp)$p.value))))
print(cmpt)
# HRD score and putative LoF BRCA2 in MSS tumors
cl[, grp_lof := factor(fifelse(!BRCA2_any, "WT", fifelse(BRCA2_lof, "LoF", "Other BRCA2")), levels = c("WT", "Other BRCA2", "LoF"))]
mss_lof <- cl[HYPER_true == FALSE, .(N = .N, HRD_score = mq(HRD_Score), HRD_high = sprintf("%d/%d", sum(HRD_high, na.rm = TRUE), sum(!is.na(HRD_high))),
                                     second_hit = sum(BRCA2_2hit)), by = grp_lof][order(grp_lof)]
cat("\nBRCA2 type and HRD in MSS tumors:\n"); print(mss_lof)
p_lof <- if (cl[HYPER_true == FALSE & grp_lof == "LoF", sum(!is.na(HRD_Score))] >= 2)
  wilcox.test(cl[HYPER_true == FALSE & grp_lof == "LoF", HRD_Score], cl[HYPER_true == FALSE & grp_lof == "WT", HRD_Score])$p.value else NA
cat("MSS: LoF BRCA2 vs WT HRD score P =", pfm(p_lof), "\n")
# List of BRCA2-altered tumors
lst <- cl[BRCA2_any == TRUE, .(sampleId, SUBTYPE, MSI_SCORE_MANTIS, TMB_NONSYNONYMOUS, n_mut, n_trunc, cna, BRCA2_2hit, HRD_Score, mutSig3)]
lst <- merge(lst, mu[gene == "BRCA2", .(proteinChange = paste(proteinChange, collapse = ", "), types = paste(mutationType, collapse = ", ")), by = sampleId], by = "sampleId", all.x = TRUE)
print(lst[order(SUBTYPE)])

# ---- B2.4 Accuracy of the proxy --------------------------------------------
cl[, top5 := MUTATION_COUNT > quantile(MUTATION_COUNT, 0.95, na.rm = TRUE)]
cl[, proxy := top5 | MMR_lof | POLE_hot]
dx <- function(pred, truth, lab) { tp <- sum(pred & truth); fp <- sum(pred & !truth); fn <- sum(!pred & truth); tn <- sum(!pred & !truth)
  ci <- function(x, n) sprintf("%.1f%% (%.1f-%.1f)", 100 * x / n, 100 * binom.test(x, n)$conf.int[1], 100 * binom.test(x, n)$conf.int[2])
  po <- (tp + tn) / (tp + fp + fn + tn); pe <- ((tp + fp) * (tp + fn) + (fn + tn) * (fp + tn)) / (tp + fp + fn + tn)^2
  data.table(reference = lab, TP = tp, FP = fp, FN = fn, TN = tn, sensitivity = ci(tp, tp + fn), specificity = ci(tn, tn + fp),
             PPV = ci(tp, tp + fp), NPV = ci(tn, tn + fn), kappa = round((po - pe) / (1 - pe), 2)) }
for (q in c(0.85, 0.90)) cl[, (paste0("top", round(100 * (1 - q)))) := MUTATION_COUNT > quantile(MUTATION_COUNT, q, na.rm = TRUE)]
cl[, proxy15 := top15 | MMR_lof | POLE_hot]
cat("Share of MSI or POLE tumors in TCGA:", round(100 * mean(cl$HYPER_true), 1), "%\n")
acc <- rbind(dx(cl$proxy15, cl$HYPER_true, "Proxy with top 15% vs MSI or POLE subtype"), dx(cl$top10 | cl$MMR_lof | cl$POLE_hot, cl$HYPER_true, "Proxy with top 10% vs MSI or POLE subtype"),
dx(cl$proxy, cl$HYPER_true, "Proxy vs MSI or POLE subtype"), dx(cl$proxy, cl$MSI, "Proxy vs MSI subtype"),
             dx(cl$top5, cl$HYPER_true, "Mutation count top 5% alone vs MSI or POLE subtype"))
cat("\n=== B2.4 ACCURACY OF THE PROXY (TCGA) ===\n"); print(acc)
cat("\nMUTATION_COUNT and MANTIS of false negatives (MSI tumors missed by the proxy):\n")
print(cl[HYPER_true & !proxy, .(sampleId, SUBTYPE, MUTATION_COUNT, MSI_SCORE_MANTIS, BRCA2_any)])

fwrite(b22, file.path(RV, "B2_2_brca2_msi.csv")); fwrite(g4, file.path(RV, "B2_3_dort_grup.csv"))
fwrite(cmpt, file.path(RV, "B2_3_karsilastirma.csv")); fwrite(mss_lof, file.path(RV, "B2_3_mss_lof.csv"))
fwrite(lst, file.path(RV, "B2_3_brca2_tcga_liste.csv")); fwrite(acc, file.path(RV, "B2_4_vekil_dogruluk.csv"))

# ---- Figure 7 -----------------------------------------------------------------
pd <- melt(cl[!is.na(HRD_Score), .(grp, `HRD score` = HRD_Score, `Signature 3 weight` = mutSig3)], id.vars = "grp")
gl <- cl[!is.na(HRD_Score), .N, by = grp]; lv <- setNames(paste0(gsub(" / ", "\n", levels(cl$grp)), "\n(n = ", gl$N[match(levels(cl$grp), gl$grp)], ")"), levels(cl$grp))
pd[, grp := factor(lv[as.character(grp)], levels = lv)]
g7 <- ggplot(pd, aes(grp, value, fill = grepl("altered", grp))) + geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) + facet_wrap(~variable, scales = "free_y") +
  scale_fill_manual(values = c("grey85", "#C0392B"), guide = "none") + labs(x = NULL, y = NULL) +
  theme_bw(base_size = 10)
ggsave(file.path(RV, "Figure7_TCGA_HRD.png"), g7, width = 10, height = 5, dpi = 300)
ggsave(file.path(RV, "Figure7_TCGA_HRD.tiff"), g7, width = 10, height = 5, dpi = 300, compression = "lzw")
cat("\nDone.\n"); sink()
