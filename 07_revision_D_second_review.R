# =============================================================================
# REVISION D - Second editorial review
#   D1  GENIE: P80 and P85 cutoffs added to the sensitivity analyses
#   D2  GENIE: putative LoF defined without hypermutation; split by hypermutation
#   D3  GENIE: Table 3 tested n/N; Table 4 Benjamini-Hochberg q values
#   D4  GENIE: sensitivity analysis including poorly differentiated carcinomas (SPDAC, EPDCA)
#   D5  BRCA2 frameshifts at mononucleotide repeats (coding sequence NM_000059.4), GENIE and TCGA
#   D6  TCGA: cohort flow, missing HRD scores, POLE subtype vs POLE hotspot
#   D7  TCGA: HRD score ~ BRCA2 x MSI/POLE interaction (bootstrap CI), subtype analysis
#   D8  TCGA: proxy accuracy by BRCA2 status (differential misclassification)
#   D9  TCGA: allelic status tiers with variant-matched ABSOLUTE segments; patient-level table
#   D10 GENIE vs TCGA harmonized BRCA2 definitions (broad, pathogenic, truncating)
#   D11 TCGA: absolute signature 3 contribution
#   D12 Figures: 1 (GENIE + TCGA flow), 3 (repeats), 4 (adjusted OR panel), 5 (P80/P85,
#       numeric OR), 6 (TCGA HRD with 42 line, N and P), 7 (TCGA patient-level oncoprint),
#       S1 (comparison genes with 95% CI)
# All numbers are computed in this R session.
# =============================================================================
while (sink.number() > 0) sink()
PROJ <- path.expand("~/GENIE_BRCA2"); setwd(PROJ)
RV <- file.path(PROJ, "revizyon"); TD <- file.path(RV, "tcga"); FD <- file.path(RV, "D_figures")
dir.create(FD, showWarnings = FALSE)

# ---- GENIE objects -----------------------------------------------------------
source("02_manuscript_tables_figures.R", encoding = "UTF-8")
while (sink.number() > 0) sink()
m_crude <- copy(crude); m_nh <- copy(co_strat[tabaka == "Hipermutasyonsuz"]); m_flow <- copy(flow_ms); m_gene_lab <- gene_lab
source("03_revision_A_genie_sensitivity.R", encoding = "UTF-8")
while (sink.number() > 0) sink()
g_b2 <- copy(b2); g_rob <- copy(rob); g_s2602 <- copy(s2602); g_cmp <- copy(cmp); g_co8 <- copy(co8); g_lab_en <- lab_en
for (p in c("httr", "jsonlite", "patchwork")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(httr); library(jsonlite); library(patchwork); library(splines)

sink(file.path(RV, "D_log.txt"), split = TRUE)
cat("Revision D date:", format(Sys.time()), "\n")
pfm <- function(p) ifelse(is.na(p), NA, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
ci3 <- function(o, l, u) sprintf("%.2f (%.2f-%.2f)", o, l, u)
save_d <- function(g, name, w, h) {
  ggsave(file.path(FD, paste0(name, ".png")),  g, width = w, height = h, dpi = 300)
  ggsave(file.path(FD, paste0(name, ".tiff")), g, width = w, height = h, dpi = 300, compression = "lzw")
}

# ---- D1 Sensitivity with P80 and P85 ----------------------------------------
cat("\n=== D1 SENSITIVITY ANALYSES INCLUDING P80 AND P85 ===\n")
pt[, h_p80 := hyper_def(mc_x, 0.80)]
pt[, h_p85 := hyper_def(mc_x, 0.85)]
defs <- c(hyper = "Main: P95 (total count) or MMR/POLE",
          h_p80 = "P80 (BRCA2-free count) or MMR/POLE", h_p85 = "P85 (BRCA2-free count) or MMR/POLE",
          h_p90 = "P90 (BRCA2-free count) or MMR/POLE", h_p95 = "P95 (BRCA2-free count) or MMR/POLE",
          h_p975 = "P97.5 (BRCA2-free count) or MMR/POLE", h_p99 = "P99 (BRCA2-free count) or MMR/POLE",
          h_mmr = "MMR/POLE alteration only", h_mc95 = "Mutation count P95 only (BRCA2-free count)")
pt[, lof := trunc | BRCA2_homdel]
d1 <- rbindlist(lapply(names(defs), function(h) {
  s <- summ(pt, h, defs[[h]])
  x <- sum(pt$lof & !pt[[h]]); r <- prop_ci(x, nrow(pt))
  y <- sum(pt$lof & pt[[h]])
  cbind(s, lof_nonhyper_n = x, lof_nonhyper_pct = sprintf("%.2f (%.2f-%.2f)", r[1], r[2], r[3]), lof_hyper_n = y)
}))
print(d1)

# ---- D2 Putative LoF independent of hypermutation ---------------------------
cat("\n=== D2 PUTATIVE LoF (TRUNCATING OR DEEP DELETION) ===\n")
d2 <- data.table(
  category = c("Any BRCA2 alteration", "Putative LoF BRCA2 alteration", "  in hypermutated tumors", "  in non-hypermutated tumors",
               "Recurrent frameshift (any)", "Recurrent frameshift as the only BRCA2 change"),
  n = c(sum(pt$BRCA2_alt), sum(pt$lof), sum(pt$lof & pt$hyper), sum(pt$lof & !pt$hyper),
        sum(g_b2$rec_fs), sum(pt$only_recurrent_fs & pt$BRCA2_alt)))
d2[, pct := sapply(n, function(x) { r <- prop_ci(x, nrow(pt)); sprintf("%.2f (%.2f-%.2f)", r[1], r[2], r[3]) })]
print(d2)
cat("Recurrent frameshift (any) by hypermutation:", sum(g_b2$rec_fs & g_b2$hyper), "hypermutated,", sum(g_b2$rec_fs & !g_b2$hyper), "non-hypermutated\n")
cat("Mutation count source: GENIE clinical field Mutation.Count (all mutations reported for the sample). Zero counts: ",
    sum(pt$mut_count == 0, na.rm = TRUE), "\n")

# ---- D3 Table 3 n/N and Table 4 q values ------------------------------------
cat("\n=== D3a TABLE 3 TESTED n/N (ALL TUMORS) ===\n")
d3a <- rbindlist(lapply(co_genes, function(g) {
  dd <- if (g == "ERBB2") pt[ERBB2_cna_prof == TRUE] else pt[get(paste0(g, "_mut_prof")) == TRUE]
  tb <- table(factor(dd[[paste0(g, "_alt")]], c(FALSE, TRUE)), factor(dd$BRCA2_alt, c(FALSE, TRUE)))
  data.table(gene = g, tested = nrow(dd),
             in_BRCA2_altered = sprintf("%d/%d (%.1f%%)", tb[2, 2], sum(tb[, 2]), 100 * tb[2, 2] / sum(tb[, 2])),
             in_BRCA2_wt = sprintf("%d/%d (%.1f%%)", tb[2, 1], sum(tb[, 1]), 100 * tb[2, 1] / sum(tb[, 1])))
}))
print(d3a)
cat("\n=== D3b TABLE 4 WITH BENJAMINI-HOCHBERG q ===\n")
nh4 <- sum(g_b2$hyper); nn4 <- sum(!g_b2$hyper)
vars4 <- c("trunc_hd", "fs", "non", "spl", "mis_only", "rec_fs", "multi", "BRCA2_homdel", "br1_pal")
d3b <- rbindlist(lapply(vars4, function(v) {
  a <- sum(g_b2[[v]] & g_b2$hyper); b <- sum(g_b2[[v]] & !g_b2$hyper)
  data.table(feature = v, hyper = sprintf("%d (%.1f)", a, 100 * a / nh4), nonhyper = sprintf("%d (%.1f)", b, 100 * b / nn4),
             p = fisher.test(matrix(c(a, nh4 - a, b, nn4 - b), 2))$p.value)
}))
d3b[, q := p.adjust(p, "BH")]; d3b[, `:=`(P = pfm(p), Q = pfm(q))]
print(d3b)

# ---- D4 Including poorly differentiated carcinomas ---------------------------
cat("\n=== D4 COHORT INCLUDING SPDAC AND EPDCA ===\n")
build <- function(codes) {
  d <- merge(clin, alt[, !c("Study.ID", "Patient.ID"), with = FALSE], by = "Sample.ID")
  keep <- c("Sample.ID", grep(paste0("^(", paste(c(MMR_GENES, "POLE"), collapse = "|"), ")"), names(mmr), value = TRUE))
  d <- merge(d, mmr[, ..keep], by = "Sample.ID", all.x = TRUE)
  d <- d[Oncotree.Code %in% names(codes) & Sample.class == "Tumor"]
  d[, site := factor(codes[Oncotree.Code], levels = c("Stomach", "GEJ", "Esophagus", "Esophagogastric NOS"))]
  for (g in c("BRCA1", "BRCA2")) {
    d[, paste0(g, "_mut_prof") := !is_np(get(paste0(g, "..MUT")))]
    d[, paste0(g, "_cna_prof") := !is_np(get(paste0(g, "..HOMDEL")))]
  }
  d[, BRCA2_alt := is_drv(BRCA2..MUT) | is_drv(BRCA2..HOMDEL)]
  d[, prio_primary := Sample.Type == "Primary"]; d[, seq_year := suppressWarnings(as.numeric(Sequence.Year))]
  setorder(d, Patient.ID, -BRCA2_mut_prof, -BRCA2_cna_prof, -prio_primary, seq_year)
  p <- d[, .SD[1], by = Patient.ID][BRCA2_mut_prof == TRUE & BRCA1_mut_prof == TRUE]
  p[, mut_count := suppressWarnings(as.numeric(Mutation.Count))]
  p[, top := mut_count > quantile(mut_count, 0.95, na.rm = TRUE), by = Sequence.Assay.ID]; p[is.na(top), top := FALSE]
  for (g in MMR_GENES) { mc <- paste0(g, "..MUT"); hc <- paste0(g, "..HOMDEL")
    p[, paste0(g, "_a") := is_drv(get(mc)) | (hc %in% names(p) & is_drv(get(hc)))] }
  p[, mmr_pole := MLH1_a | MSH2_a | MSH6_a | PMS2_a | grepl(paste0("\\b(", paste(POLE_HOT, collapse = "|"), ")\\b"), POLE..MUT)]
  p[, hyper := top | mmr_pole]
  ar <- as.character(p$Age.at.Which.Sequencing.was.Reported..Years.); ar[ar == ">89"] <- "90"
  p[, age_grp := cut(suppressWarnings(as.numeric(ar)), c(-Inf, 49, 64, 74, Inf))]
  p[, sex := fifelse(Sex %in% c("Male", "Female"), Sex, NA_character_)]
  tc <- p[, .N, by = Center][N >= MIN_CENTER_N, Center]
  p[, center_g := factor(fifelse(Center %in% tc, Center, "Other"))]
  p[, center_g := relevel(center_g, ref = names(sort(table(center_g), decreasing = TRUE))[1])]
  p
}
chk <- build(ADENO_CODES)
cat("Check with original codes: patients", nrow(chk), " BRCA2", sum(chk$BRCA2_alt), " hypermutated", sum(chk$hyper), "\n")
ext <- build(c(ADENO_CODES, SPDAC = "Stomach", EPDCA = "Esophagus"))
cat("Extended cohort: patients", nrow(ext), " (added", nrow(ext) - nrow(chk), ")  BRCA2", sum(ext$BRCA2_alt),
    " hypermutated", sum(ext$hyper), "\n")
de <- ext[!is.na(age_grp) & !is.na(sex)]
me <- logistf(BRCA2_alt ~ hyper + age_grp + sex + site + center_g, data = de); ie <- which(names(coef(me)) == "hyperTRUE")
d4 <- data.table(cohort = "Including SPDAC and EPDCA", N = nrow(ext), BRCA2 = sum(ext$BRCA2_alt),
  BRCA2_in_hyper = sprintf("%d/%d (%.1f%%)", sum(ext$BRCA2_alt & ext$hyper), sum(ext$hyper), 100 * mean(ext$BRCA2_alt[ext$hyper])),
  BRCA2_in_nonhyper = sprintf("%d/%d (%.2f%%)", sum(ext$BRCA2_alt & !ext$hyper), sum(!ext$hyper), 100 * mean(ext$BRCA2_alt[!ext$hyper])),
  carriers_hyper = sprintf("%.1f%%", 100 * mean(ext$hyper[ext$BRCA2_alt])),
  OR = ci3(exp(coef(me)[ie]), exp(me$ci.lower[ie]), exp(me$ci.upper[ie])), P = pfm(me$prob[ie]))
print(d4)
cat("Added patients by code:\n"); print(ext[!Patient.ID %in% chk$Patient.ID, .N, by = Oncotree.Code])

# ---- D5 Frameshifts at mononucleotide repeats ---------------------------------
cat("\n=== D5 BRCA2 FRAMESHIFTS AT MONONUCLEOTIDE REPEATS ===\n")
cds <- tryCatch({
  r <- GET("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NM_000059.4&rettype=fasta_cds_na&retmode=text", timeout(60))
  x <- strsplit(content(r, "text", encoding = "UTF-8"), "\n")[[1]]; paste(x[-1], collapse = "") }, error = function(e) NA)
cat("BRCA2 CDS length (nt):", nchar(cds), " codons:", nchar(cds) / 3, "\n")
cs <- strsplit(cds, "")[[1]]
max_run <- function(pos, flank = 6) {
  s <- max(1, 3 * pos - 2 - flank); e <- min(length(cs), 3 * pos + flank)
  r <- rle(cs[s:e]); i <- which.max(r$lengths); c(len = r$lengths[i], base = r$values[i])
}
REP <- 5
bg <- sapply(1:3418, function(p) as.numeric(max_run(p)["len"]))
cat("Background: share of BRCA2 codons with a run of >=", REP, "nt within +/-6 nt:", round(100 * mean(bg >= REP), 1), "%\n")
gfs <- bv[cls == "Frameshift" & !is.na(pos) & pos <= 3418]
gfs[, `:=`(run = as.numeric(sapply(pos, function(p) max_run(p)["len"])), base = sapply(pos, function(p) max_run(p)["base"]))]
gfs[, at_rep := run >= REP]
rep_tab <- function(d, hv, lab) {
  a <- sum(d$at_rep & d[[hv]]); n1 <- sum(d[[hv]]); b <- sum(d$at_rep & !d[[hv]]); n0 <- sum(!d[[hv]])
  ft <- fisher.test(matrix(c(a, n1 - a, b, n0 - b), 2))
  bt1 <- binom.test(a, n1, mean(bg >= REP))$p.value; bt0 <- binom.test(b, n0, mean(bg >= REP))$p.value
  data.table(cohort = lab, hyper_at_repeat = sprintf("%d/%d (%.1f%%)", a, n1, 100 * a / n1),
             nonhyper_at_repeat = sprintf("%d/%d (%.1f%%)", b, n0, 100 * b / max(1, n0)),
             background = sprintf("%.1f%%", 100 * mean(bg >= REP)),
             P_hyper_vs_nonhyper = pfm(ft$p.value), P_hyper_vs_background = pfm(bt1), P_nonhyper_vs_background = pfm(bt0))
}
print(rep_tab(gfs, "hyper", "GENIE frameshift variants"))
cat("\nRecurrent GENIE frameshift codons and repeat context:\n")
d5b <- gfs[, .(n = .N, hyper = sum(hyper), run = run[1], base = base[1], context = paste(cs[max(1, 3 * pos[1] - 8):min(length(cs), 3 * pos[1] + 6)], collapse = "")),
           by = codon_key][n >= 3][order(-n)]
print(d5b)

# ---- TCGA objects -------------------------------------------------------------
sink()
source("05_revision_B2_tcga_validation.R", encoding = "UTF-8")
while (sink.number() > 0) sink()
sink(file.path(RV, "D_log.txt"), split = TRUE, append = TRUE)
pfm <- function(p) ifelse(is.na(p), NA, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

# D5 for TCGA has to be recomputed now that cl/mu exist
mu_b <- mu[gene == "BRCA2" & mutationType %in% c("Frame_Shift_Del", "Frame_Shift_Ins")]
mu_b[, pos := suppressWarnings(as.numeric(sub("^[A-Z*X](\\d+).*", "\\1", proteinChange)))]
mu_b <- merge(mu_b[!is.na(pos) & pos <= 3418], cl[, .(sampleId, HYPER_true)], by = "sampleId")
mu_b[, run := as.numeric(sapply(pos, function(p) max_run(p)["len"]))]; mu_b[, at_rep := run >= REP]
d5 <- rbind(rep_tab(gfs, "hyper", "GENIE frameshift variants"), rep_tab(mu_b, "HYPER_true", "TCGA frameshift variants"))
cat("\n=== D5 (final) FRAMESHIFTS AT MONONUCLEOTIDE REPEATS, GENIE AND TCGA ===\n"); print(d5)

# ---- D6 TCGA flow, missing HRD, POLE -----------------------------------------
cat("\n=== D6a TCGA FLOW ===\n")
cl0 <- rbindlist(lapply(st, function(s) {
  cs0 <- rd_w(file.path(TD, paste0(s, "_clin_sample_long.csv")), "sampleId + patientId")
  cp0 <- rd_w(file.path(TD, paste0(s, "_clin_patient_long.csv")), "patientId")
  x <- merge(cs0, cp0, by = "patientId", all.x = TRUE); x[, study := s]; x }), fill = TRUE)
cl0[, MUTATION_COUNT := suppressWarnings(as.numeric(MUTATION_COUNT))]
is_ad <- grepl("Adenocarcinoma|Signet Ring", cl0$CANCER_TYPE_DETAILED) & !grepl("ESCC", fifelse(is.na(cl0$SUBTYPE), "", cl0$SUBTYPE))
tflow <- data.table(step = c("TCGA PanCancer Atlas stomach and esophageal samples", "Excluded: not adenocarcinoma",
                             "Excluded: no mutation data", "Adenocarcinomas with mutation data", "With HRD score", "With ABSOLUTE segment at BRCA2"),
                    n = c(nrow(cl0), sum(!is_ad), sum(is_ad & is.na(cl0$MUTATION_COUNT)), nrow(cl),
                          sum(!is.na(cl$HRD_Score)), NA))
print(tflow)
cat("\n=== D6b MISSING HRD SCORE ===\n")
cl[, hrd_avail := !is.na(HRD_Score)]
print(cl[, .(N = .N, with_HRD = sum(hrd_avail), pct = round(100 * mean(hrd_avail), 1)), by = grp][order(grp)])
print(cl[, .(N = .N, with_HRD = sum(hrd_avail)), by = SUBTYPE][order(SUBTYPE)])
cat("Missing HRD vs BRCA2 status, Fisher P =", pfm(fisher.test(table(cl$hrd_avail, cl$BRCA2_any))$p.value),
    "; vs MSI/POLE, Fisher P =", pfm(fisher.test(table(cl$hrd_avail, cl$HYPER_true))$p.value), "\n")
cat("\n=== D6c POLE SUBTYPE vs POLE MUTATIONS ===\n")
pm <- mu[gene == "POLE", .(POLE_variants = paste(proteinChange, collapse = ", ")), by = sampleId]
d6c <- merge(cl[POLE_sub == TRUE, .(sampleId, SUBTYPE, MUTATION_COUNT, POLE_hot, BRCA2_any)], pm, by = "sampleId", all.x = TRUE)
print(d6c)
cat("POLE subtype tumors with a POLE exonuclease hotspot:", sum(d6c$POLE_hot), "/", nrow(d6c), "\n")

# ---- D7 Interaction model ------------------------------------------------------
cat("\n=== D7 HRD SCORE ~ BRCA2 x MSI/POLE ===\n")
dm <- cl[!is.na(HRD_Score)]
dm[, `:=`(BRCA2 = as.integer(BRCA2_any), MSIPOLE = as.integer(HYPER_true), esca = as.integer(study == "esca_tcga_pan_can_atlas_2018"))]
f1 <- HRD_Score ~ BRCA2 * MSIPOLE + esca
f2 <- HRD_Score ~ BRCA2 * MSIPOLE + esca + purity + AGE
boot_lm <- function(f, d, B = 2000) {
  set.seed(2026); m <- lm(f, data = d)
  bs <- t(replicate(B, { i <- sample(nrow(d), replace = TRUE); cf <- coef(lm(f, data = d[i])); cf[names(coef(m))] }))
  data.table(term = names(coef(m)), estimate = coef(m), model_lo = confint(m)[, 1], model_hi = confint(m)[, 2],
             boot_lo = apply(bs, 2, quantile, 0.025, na.rm = TRUE), boot_hi = apply(bs, 2, quantile, 0.975, na.rm = TRUE),
             P = summary(m)$coefficients[, 4], n = nobs(m))
}
d7a <- boot_lm(f1, dm); d7b <- boot_lm(f2, dm[!is.na(purity) & !is.na(AGE)])
cat("Model 1 (BRCA2 x MSI/POLE + site):\n"); print(d7a)
cat("Model 2 (+ purity + age):\n"); print(d7b)
cat("\nHRD-high (>=42) logistic, Firth, BRCA2 x MSI/POLE + site:\n")
dm[, hrd_high := as.integer(HRD_Score >= 42)]
lf <- logistf(hrd_high ~ BRCA2 * MSIPOLE + esca, data = dm)
d7c <- data.table(term = names(coef(lf)), OR = exp(coef(lf)), lo = exp(lf$ci.lower), hi = exp(lf$ci.upper), P = lf$prob)
print(d7c)
cat("\nHRD by molecular subtype (non-MSI/POLE tumors):\n")
dm[, sub := fifelse(is.na(SUBTYPE), "Missing", sub("^(STAD|ESCA)_", "", SUBTYPE))]
d7d <- dm[MSIPOLE == 0, .(N = .N, BRCA2_altered = sum(BRCA2), HRD = mq(HRD_Score),
                          HRD_altered = mq(HRD_Score[BRCA2 == 1]), HRD_wt = mq(HRD_Score[BRCA2 == 0])), by = sub][order(sub)]
print(d7d)
d7e <- boot_lm(HRD_Score ~ BRCA2 + sub + esca, dm[MSIPOLE == 0])
cat("Non-MSI/POLE tumors, HRD ~ BRCA2 + subtype + site:\n"); print(d7e)

# ---- D8 Proxy accuracy by BRCA2 status -----------------------------------------
cat("\n=== D8 PROXY SENSITIVITY BY BRCA2 STATUS ===\n")
cl[, proxy10 := top10 | MMR_lof | POLE_hot]
d8 <- rbindlist(lapply(list(c("proxy", "P95"), c("proxy10", "P90"), c("proxy15", "P85")), function(z) rbindlist(lapply(c(TRUE, FALSE), function(b) {
  d <- cl[BRCA2_any == b]; tp <- sum(d[[z[1]]] & d$HYPER_true); fn <- sum(!d[[z[1]]] & d$HYPER_true)
  tn <- sum(!d[[z[1]]] & !d$HYPER_true); fp <- sum(d[[z[1]]] & !d$HYPER_true)
  data.table(proxy = z[2], BRCA2 = ifelse(b, "BRCA2-altered", "BRCA2 wild-type"), TP = tp, FN = fn, FP = fp, TN = tn,
             sensitivity = sprintf("%.1f%%", 100 * tp / max(1, tp + fn)), specificity = sprintf("%.1f%%", 100 * tn / max(1, tn + fp)))
}))))
print(d8)
cat("Sensitivity difference (P95), Fisher P =", pfm(fisher.test(matrix(c(d8[proxy == "P95"]$TP, d8[proxy == "P95"]$FN), 2))$p.value), "\n")

# ---- D9 Allelic status tiers ---------------------------------------------------
cat("\n=== D9 ALLELIC STATUS (variant-matched ABSOLUTE segment) ===\n")
seg <- fread(file.path(TD, "TCGA_mastercalls.abs_segtabs.fixed.txt"))
seg[, sampleId := substr(Sample, 1, 15)]
seg <- seg[sampleId %in% cl$sampleId & as.character(Chromosome) == "13" & Start <= 32973809 & End >= 32889611]
seg[, loh := (LOH == 1) | (Modal_HSCN_1 == 0)]
nseg <- seg[, .(n_seg = .N, any_loh = any(loh, na.rm = TRUE), all_loh = all(loh, na.rm = TRUE), homdel_abs = any(Modal_Total_CN == 0, na.rm = TRUE)), by = sampleId]
cat("BRCA2-altered tumors with >1 overlapping segment:", nseg[sampleId %in% cl[BRCA2_any == TRUE, sampleId] & n_seg > 1, .N],
    "; with discordant LOH calls:", nseg[sampleId %in% cl[BRCA2_any == TRUE, sampleId] & n_seg > 1 & any_loh != all_loh, .N], "\n")
mv <- mu[gene == "BRCA2" & mutationType != "Silent", .(sampleId, proteinChange, mutationType, startPosition)]
mv <- merge(mv, seg[, .(sampleId, Start, End, loh)], by = "sampleId", allow.cartesian = TRUE, all.x = TRUE)
mv <- mv[is.na(Start) | (startPosition >= Start & startPosition <= End)]
mvl <- mv[, .(loh_at_variant = any(loh, na.rm = TRUE)), by = sampleId]
gn <- trimws(readLines(file.path(TD, "ddr", "Genes.tsv")))
smp <- fread(file.path(TD, "ddr", "Samples.tsv"), header = FALSE, col.names = c("sampleId", "ct"))
gm <- fread(file.path(TD, "ddr", "GeneMutations.tsv"), header = FALSE)
kn <- data.table(sampleId = smp$sampleId, KN_mut = as.numeric(unlist(gm[which(gn == "BRCA2")])))
cl <- merge(cl, kn, by = "sampleId", all.x = TRUE); cl[, curated := KN_mut %in% 1]
cl <- merge(cl, nseg, by = "sampleId", all.x = TRUE); cl <- merge(cl, mvl, by = "sampleId", all.x = TRUE)
cl[is.na(loh_at_variant), loh_at_variant := FALSE]
cl[, locus_loh := fifelse(n_mut > 0, loh_at_variant, any_loh %in% TRUE)]
cl[, tier := fifelse(!BRCA2_any, NA_character_,
             fifelse(cna %in% -2 | homdel_abs %in% TRUE, "1 Homozygous deletion",
             fifelse((n_trunc > 0 | curated) & locus_loh, "2 Pathogenic mutation with locus LOH",
             fifelse(n_trunc >= 2, "3 Multiple truncating hits, phase unresolved",
             fifelse(n_trunc > 0 | curated, "4 Single truncating or curated event", "5 Non-truncating variant(s) only")))))]
d9 <- dcast(cl[BRCA2_any == TRUE, .N, by = .(tier, grp = fifelse(HYPER_true, "MSI or POLE", "Non-MSI/POLE"))], tier ~ grp, value.var = "N", fill = 0)
print(d9)
cat("Locus LOH among BRCA2-altered: MSI/POLE", cl[BRCA2_any & HYPER_true, sprintf("%d/%d", sum(locus_loh), .N)],
    "; non-MSI/POLE", cl[BRCA2_any & !HYPER_true, sprintf("%d/%d", sum(locus_loh), .N)], "\n")
cat("Locus LOH among BRCA2 wild-type (any segment): MSI/POLE", cl[!BRCA2_any & HYPER_true & !is.na(n_seg), sprintf("%d/%d", sum(any_loh), .N)],
    "; non-MSI/POLE", cl[!BRCA2_any & !HYPER_true & !is.na(n_seg), sprintf("%d/%d", sum(any_loh), .N)], "\n")
cat("\nHRD by tier:\n"); print(cl[BRCA2_any == TRUE, .(N = .N, HRD = mq(HRD_Score), HRD_high = sprintf("%d/%d", sum(HRD_Score >= 42, na.rm = TRUE), sum(!is.na(HRD_Score)))), by = .(tier, MSIPOLE = HYPER_true)][order(tier, MSIPOLE)])

# ---- D11 Absolute signature 3 ---------------------------------------------------
cl[, sbs3_abs := mutSig3 * MUTATION_COUNT]
cat("\n=== D11 SIGNATURE 3: RELATIVE WEIGHT AND ESTIMATED NUMBER OF MUTATIONS ===\n")
print(cl[!is.na(mutSig3), .(N = .N, sbs3_rel = mq(mutSig3), sbs3_abs = mq(sbs3_abs)), by = grp][order(grp)])
pw2 <- function(v, a, b) { x <- cl[grp == a][[v]]; y <- cl[grp == b][[v]]; wilcox.test(x, y)$p.value }
cat("Estimated SBS3 mutations, altered vs WT: non-MSI/POLE P =", pfm(pw2("sbs3_abs", "BRCA2-altered / MSS", "BRCA2 wild-type / MSS")),
    "; MSI/POLE P =", pfm(pw2("sbs3_abs", "BRCA2-altered / MSI or POLE", "BRCA2 wild-type / MSI or POLE")), "\n")

# patient-level table
lv2 <- mu[gene == "BRCA2" & mutationType != "Silent", .(variants = paste(proteinChange, collapse = ", "),
          types = paste(unique(sub("_Mutation|_Site", "", mutationType)), collapse = ", ")), by = sampleId]
d9p <- merge(cl[BRCA2_any == TRUE, .(sampleId, site = fifelse(study == "stad_tcga_pan_can_atlas_2018", "Stomach", "Esophagus"),
  subtype = fifelse(is.na(SUBTYPE), "Missing", sub("^(STAD|ESCA)_", "", SUBTYPE)), MSI_POLE = HYPER_true,
  MANTIS = round(MSI_SCORE_MANTIS, 2), TMB = round(TMB_NONSYNONYMOUS, 1), n_trunc, curated, deep_del = cna %in% -2,
  locus_LOH = locus_loh, tier, HRD = HRD_Score, SBS3 = round(mutSig3, 3), SBS3_abs = round(sbs3_abs, 1))],
  lv2, by = "sampleId", all.x = TRUE)
setorder(d9p, -MSI_POLE, tier, sampleId)
cat("\nPatient-level table (n =", nrow(d9p), "):\n"); print(d9p)

# ---- D10 Harmonized definitions ---------------------------------------------------
cat("\n=== D10 HARMONIZED BRCA2 DEFINITIONS ===\n")
hdef <- function(y, h, lab, coh, hl) { t <- table(factor(h, c(FALSE, TRUE)), factor(y, c(FALSE, TRUE))); f <- fisher.test(t)
  data.table(cohort = coh, definition = lab, hyper_definition = hl, n = sum(y), N = length(y), pct = round(100 * mean(y), 2),
             in_hyper = sprintf("%d/%d (%.1f%%)", t[2, 2], sum(t[2, ]), 100 * t[2, 2] / sum(t[2, ])),
             in_nonhyper = sprintf("%d/%d (%.2f%%)", t[1, 2], sum(t[1, ]), 100 * t[1, 2] / sum(t[1, ])),
             carriers_hyper = sprintf("%.1f%%", 100 * t[2, 2] / sum(t[, 2])), OR = ci3(f$estimate, f$conf.int[1], f$conf.int[2])) }
cl[, pathogenic := curated | cna %in% -2]
d10 <- rbind(
  hdef(pt$BRCA2_alt, pt$hyper, "Pathogenic (OncoKB oncogenic or likely oncogenic)", "GENIE", "Panel proxy"),
  hdef(pt$lof, pt$hyper, "Truncating or deep deletion", "GENIE", "Panel proxy"),
  hdef(cl$BRCA2_any, cl$HYPER_true, "Broad (any non-silent variant or deep deletion)", "TCGA", "MSI or POLE subtype"),
  hdef(cl$pathogenic, cl$HYPER_true, "Pathogenic (curated deleterious mutation or deep deletion)", "TCGA", "MSI or POLE subtype"),
  hdef(cl$BRCA2_lof, cl$HYPER_true, "Truncating or deep deletion", "TCGA", "MSI or POLE subtype"),
  hdef(cl$pathogenic, cl$proxy, "Pathogenic (curated deleterious mutation or deep deletion)", "TCGA", "Panel-style proxy (P95)"),
  hdef(cl$BRCA2_lof, cl$proxy, "Truncating or deep deletion", "TCGA", "Panel-style proxy (P95)"))
print(d10)

# ---- D12 Figures --------------------------------------------------------------------
cat("\n=== D12 FIGURES ===\n")
pal <- c(blue = "#2C6E9B", red = "#C0392B", grey = "grey80")

# Figure 1: GENIE (A) and TCGA (B)
fl <- copy(m_flow)[, y := rev(seq_len(.N)) * 2]
fl[, lab := sprintf("%s\n%s samples, %s patients", step, format(samples, big.mark = ","), format(patients, big.mark = ","))]
fl[.N, lab := sprintf("%s\n%s patients (study cohort)", step, format(patients, big.mark = ","))]
ex <- data.table(samples = -diff(m_flow$samples), patients = -diff(m_flow$patients))[, y := fl$y[-nrow(fl)] - 1]
ex[, lab := c(sprintf("Excluded: non-adenocarcinoma histology\n%s samples, %s patients", format(samples[1], big.mark = ","), format(patients[1], big.mark = ",")),
             sprintf("Excluded: cell-free DNA samples\n%s samples, %s patients", samples[2], patients[2]),
             sprintf("Excluded: additional samples\nfrom the same patient\n%s samples", samples[3]),
             sprintf("Excluded: panel did not cover\nboth BRCA1 and BRCA2\n%s patients", patients[4]))]
flowplot <- function(fl, ex, ttl) ggplot() +
  geom_segment(data = fl[-.N], aes(x = 0, xend = 0, y = y - 0.45, yend = y - 1.55), arrow = arrow(length = unit(0.18, "cm")), colour = "grey40") +
  geom_segment(data = ex, aes(x = 0, xend = 0.95, y = y, yend = y), arrow = arrow(length = unit(0.15, "cm")), colour = "grey40") +
  geom_label(data = fl, aes(0, y, label = lab), size = 2.9, label.padding = unit(0.4, "lines"), fill = "white") +
  geom_label(data = ex, aes(1.9, y, label = lab), size = 2.6, label.padding = unit(0.35, "lines"), fill = "grey95") +
  xlim(-1.4, 3.1) + labs(title = ttl) + theme_void() + theme(plot.title = element_text(face = "bold", size = 12))
nabs <- cl[!is.na(n_seg), .N]
tfl <- data.table(lab = c(sprintf("TCGA PanCancer Atlas stomach and\nesophageal samples\n%d samples", tflow$n[1]),
                          sprintf("Adenocarcinomas with mutation data\n%d tumors (validation cohort)", nrow(cl)),
                          sprintf("Allele-specific copy number at BRCA2\n%d tumors", nabs),
                          sprintf("Genomic HRD score available\n%d tumors", sum(!is.na(cl$HRD_Score)))), y = c(8, 6, 4, 2))
tex <- data.table(lab = c(sprintf("Excluded: not adenocarcinoma\n%d samples", tflow$n[2]), sprintf("Excluded: no mutation data\n%d samples", tflow$n[3])), y = c(7.2, 6.8))
tex <- tex[1:2][, y := c(7.3, 6.7)]
tflow_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = 0, y = 7.55, yend = 6.45), arrow = arrow(length = unit(0.18, "cm")), colour = "grey40") +
  geom_segment(aes(x = 0, xend = 0, y = 5.55, yend = 4.45), arrow = arrow(length = unit(0.18, "cm")), colour = "grey40") +
  geom_segment(aes(x = 0, xend = 0, y = 3.55, yend = 2.45), arrow = arrow(length = unit(0.18, "cm")), colour = "grey40") +
  geom_segment(data = tex, aes(x = 0, xend = 0.95, y = y, yend = y), arrow = arrow(length = unit(0.15, "cm")), colour = "grey40") +
  geom_label(data = tfl, aes(0, y, label = lab), size = 2.9, label.padding = unit(0.4, "lines"), fill = "white") +
  geom_label(data = tex, aes(1.9, y, label = lab), size = 2.6, label.padding = unit(0.35, "lines"), fill = "grey95") +
  xlim(-1.4, 3.1) + ylim(1, 11) + labs(title = "B. Validation cohort (TCGA)") + theme_void() + theme(plot.title = element_text(face = "bold", size = 12))
save_d(flowplot(fl, ex, "A. Discovery cohort (AACR Project GENIE)") | tflow_plot, "Figure1_flow", 14, 7.5)

# Figure 3: lollipop with repeat context
rep_pos <- unique(gfs[at_rep == TRUE, pos])
lp <- bv[!is.na(pos), .N, by = .(pos, cls, ctx = fifelse(hyper, "Hypermutated tumors", "Non-hypermutated tumors"))]
lp[, ctx := factor(ctx, c("Non-hypermutated tumors", "Hypermutated tumors"))]
lp[, cls := factor(cls, c("Frameshift", "Nonsense", "Splice", "Missense", "Other/ambiguous"))]
lp[, rep := factor(fifelse(cls == "Frameshift" & pos %in% rep_pos, sprintf("Frameshift at a mononucleotide repeat (>=%d nt)", REP), "Other"),
                   c(sprintf("Frameshift at a mononucleotide repeat (>=%d nt)", REP), "Other"))]
top_var <- bv[!is.na(pos), .N, by = .(pos, variant)][order(pos, -N)][, .(hgvs = paste0("p.", variant[1])), by = pos]
lab_lp <- merge(lp[, .(N = sum(N)), by = .(pos, ctx)][N >= 3], top_var, by = "pos")
dom <- data.table(name = c("PALB2-binding", "BRC repeats 1-8", "DNA-binding domain", "C-terminal RAD51 binding / NLS"),
                  start = c(10, 1002, 2479, 3265), end = c(40, 2085, 3190, 3330))
g3top <- ggplot(lp, aes(pos, N)) +
  geom_segment(aes(xend = pos, y = 0, yend = N), colour = "grey65", linewidth = 0.4) +
  geom_point(aes(fill = cls, shape = rep), size = 2.6, colour = "black", stroke = 0.3) +
  geom_text(data = lab_lp, aes(pos, N, label = hgvs), vjust = -0.9, size = 2.5) +
  facet_wrap(~ctx, ncol = 1) +
  scale_x_continuous(limits = c(0, 3418), breaks = seq(0, 3400, 500)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.4))) +
  scale_shape_manual(values = c(24, 21), name = NULL) +
  scale_fill_manual(values = c(Frameshift = "#2C6E9B", Nonsense = "#1B3A57", Splice = "#E08E45",
                               Missense = "#4DAF4A", `Other/ambiguous` = "grey70"), name = "Mutation type") +
  guides(fill = guide_legend(override.aes = list(shape = 21), order = 1), shape = guide_legend(order = 2)) +
  labs(x = NULL, y = "Number of tumors") + thm + theme(axis.text.x = element_blank(), legend.box = "vertical")
g3dom <- ggplot() + annotate("rect", xmin = 1, xmax = 3418, ymin = 0.4, ymax = 0.6, fill = "grey85") +
  geom_rect(data = dom, aes(xmin = start, xmax = end, ymin = 0.2, ymax = 0.8, fill = name), colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("#F2C14E", "#8E7CC3", "#E06666", "#76A5AF"), name = "Domain") +
  scale_x_continuous(limits = c(0, 3418), breaks = seq(0, 3400, 500)) +
  labs(x = "Position in the BRCA2 protein (amino acid 1 to 3418, NP_000050.3)", y = NULL) +
  thm + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(), legend.position = "bottom")
save_d(g3top / g3dom + plot_layout(heights = c(6, 1)), "Figure3_lollipop", 10, 7.5)

# Figure 4: three panels + adjusted OR text
f4 <- rbind(m_crude[, .(gen, set = "All tumors", OR = crude_or, lo = crude_lo, hi = crude_hi, q = crude_q)],
            m_nh[, .(gen, set = "Non-hypermutated tumors only", OR, lo = OR_alt, hi = OR_ust, q = q_BH)],
            g_co8[, .(gen = gene, set = "Adjusted for continuous mutation burden", OR, lo = OR_lo, hi = OR_hi, q)])
f4[, set := factor(set, c("All tumors", "Non-hypermutated tumors only", "Adjusted for continuous mutation burden"))]
f4 <- f4[is.finite(OR) & OR > 0 & is.finite(hi)]
f4[, gen := factor(m_gene_lab[gen], levels = rev(m_gene_lab))]
f4[, sig := fifelse(q < 0.05, "q < 0.05", "Not significant")]
g4 <- ggplot(f4, aes(OR, gen)) + geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = pmax(lo, 0.05), xmax = pmin(hi, 40)), height = 0.25, colour = "grey30") +
  geom_point(aes(fill = sig), shape = 21, size = 2.8) +
  facet_wrap(~set, labeller = label_wrap_gen(26)) + scale_x_log10(limits = c(0.05, 40), breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 10, 25)) +
  scale_fill_manual(values = c(`q < 0.05` = pal[["red"]], `Not significant` = "white"), name = NULL) +
  labs(x = "Odds ratio (BRCA2-altered vs wild-type), log scale", y = NULL) + thm
t4 <- g_co8[, .(gen = factor(m_gene_lab[gene], levels = rev(m_gene_lab)), txt = sprintf("%.2f (%.2f-%.2f)   q = %.2f", OR, OR_lo, OR_hi, q))]
g4t <- ggplot(t4, aes(0, gen, label = txt)) + geom_text(hjust = 0, size = 3) + xlim(0, 1) +
  labs(title = "Adjusted OR (95% CI)") + theme_void() + theme(plot.title = element_text(size = 9, hjust = 0.1), plot.margin = margin(32, 5, 26, 0))
save_d(g4 + g4t + plot_layout(widths = c(4, 1.1)), "Figure4_coalterations", 14, 5)

# Figure 5: sensitivity with P80/P85 and numeric OR
fp <- rbind(d1[, .(grp = "Hypermutation definition", lab = analiz, OR, OR_alt, OR_ust)],
            g_rob[, .(grp = "Subgroup and robustness", lab = analiz, OR, OR_alt, OR_ust)],
            g_s2602[, .(grp = "Subgroup and robustness", lab = analiz, OR, OR_alt, OR_ust)],
            data.table(grp = "Subgroup and robustness", lab = "Including poorly differentiated carcinomas",
                       OR = exp(coef(me)[ie]), OR_alt = exp(me$ci.lower[ie]), OR_ust = exp(me$ci.upper[ie])))
fp[, lab := ifelse(lab %in% names(g_lab_en), g_lab_en[lab], sub("Çıkarılan merkez: ", "Without ", lab))]
fp[lab == "Main: P95 (total count) or MMR/POLE" & grp == "Subgroup and robustness", lab := "All patients (main definition)"]
fp[lab == "All patients", lab := "All patients (main definition)"]
fp <- fp[!(grp == "Subgroup and robustness" & duplicated(paste(grp, lab)))]
fp[, lab := factor(lab, levels = rev(unique(lab)))]
fp[, txt := sprintf("%.2f (%.2f-%.2f)", OR, OR_alt, OR_ust)]
g5 <- ggplot(fp, aes(OR, lab)) + geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = OR_alt, xmax = OR_ust), height = 0.2) + geom_point(size = 2.2, colour = pal[["red"]]) +
  geom_text(aes(x = 110, label = txt), hjust = 0, size = 2.8) +
  scale_x_log10(limits = c(1, 900), breaks = c(1, 3, 10, 30, 100)) + facet_grid(grp ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Adjusted odds ratio for BRCA2 alteration in hypermutated tumors (log scale); right: OR (95% CI)", y = NULL) +
  theme_bw(base_size = 10) + theme(strip.text.y = element_text(angle = 0))
save_d(g5, "Figure5_sensitivity", 10, 9)

# Figure 6: TCGA HRD with 42 line, analysed N and P values
hd <- cl[!is.na(HRD_Score)]
gl <- hd[, .N, by = grp]
lvn <- setNames(paste0(gsub(" / ", "\n", gsub("MSS", "non-MSI/POLE", levels(cl$grp))), "\n(n = ", gl$N[match(levels(cl$grp), gl$grp)], ")"), levels(cl$grp))
hd[, grpl := factor(lvn[as.character(grp)], levels = lvn)]
pv <- function(a, b) pfm(wilcox.test(hd[grp == a, HRD_Score], hd[grp == b, HRD_Score])$p.value)
p_mss <- pv("BRCA2-altered / MSS", "BRCA2 wild-type / MSS"); p_msi <- pv("BRCA2-altered / MSI or POLE", "BRCA2 wild-type / MSI or POLE")
p_alt <- pv("BRCA2-altered / MSS", "BRCA2-altered / MSI or POLE")
pint <- d7a[term == "BRCA2:MSIPOLE", P]
br <- data.table(x = c(1, 3, 2), xend = c(2, 4, 4), y = c(102, 102, 112), lab = paste("P =", c(p_mss, p_msi, p_alt)))
g6 <- ggplot(hd, aes(grpl, HRD_Score, fill = BRCA2_any)) + geom_hline(yintercept = 42, linetype = 2, colour = "grey40") +
  geom_boxplot(outlier.shape = NA, width = 0.6) + geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
  geom_segment(data = br, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE) +
  geom_text(data = br, aes(x = (x + xend) / 2, y = y + 3.5, label = lab), inherit.aes = FALSE, size = 3) +
  annotate("text", x = 4.45, y = 44.5, label = "HRD score 42", size = 2.8, hjust = 1, colour = "grey30") +
  scale_fill_manual(values = c("grey85", pal[["red"]]), guide = "none") +
  scale_y_continuous(limits = c(0, 118), breaks = seq(0, 100, 25)) +
  labs(x = NULL, y = "Genomic HRD score", caption = sprintf("BRCA2 x MSI/POLE interaction in a linear model adjusted for site: P = %s", pfm(pint))) +
  theme_bw(base_size = 10)
save_d(g6, "Figure6_TCGA_HRD", 8, 5.5)

# Figure 7: TCGA patient-level oncoprint
op <- copy(d9p); op[, id := factor(sampleId, levels = sampleId)]
op[, grp := factor(fifelse(MSI_POLE, "MSI or POLE", "Non-MSI/POLE"), c("MSI or POLE", "Non-MSI/POLE"))]
op[, ev := fifelse(deep_del, "Deep deletion", fifelse(n_trunc > 0, "Truncating", "Non-truncating only"))]
tiles <- rbind(op[, .(id, grp, row = "Molecular subtype", val = subtype)],
               op[, .(id, grp, row = "BRCA2 event", val = ev)],
               op[, .(id, grp, row = "Curated deleterious (Knijnenburg)", val = fifelse(curated, "Yes", "No"))],
               op[, .(id, grp, row = "LOH at BRCA2 locus", val = fifelse(locus_LOH, "Yes", "No"))],
               op[, .(id, grp, row = "Allelic status", val = sub("^\\d ", "", tier))],
               op[, .(id, grp, row = "HRD score 42 or higher", val = fifelse(is.na(HRD), "Not available", fifelse(HRD >= 42, "Yes", "No")))])
tiles[, row := factor(row, rev(c("Molecular subtype", "BRCA2 event", "Curated deleterious (Knijnenburg)", "LOH at BRCA2 locus", "Allelic status", "HRD score 42 or higher")))]
cols <- c(MSI = "#C0392B", POLE = "#E08E45", CIN = "#2C6E9B", GS = "#76A5AF", EBV = "#8E7CC3", Missing = "grey70",
          `Deep deletion` = "#1B3A57", Truncating = "#2C6E9B", `Non-truncating only` = "#A9CCE3", Yes = "#C0392B", No = "grey88",
          `Not available` = "white", `Homozygous deletion` = "#1B3A57", `Pathogenic mutation with locus LOH` = "#C0392B",
          `Multiple truncating hits, phase unresolved` = "#E08E45", `Single truncating or curated event` = "#F2C14E",
          `Non-truncating variant(s) only` = "grey80")
g7a <- ggplot(op, aes(id, HRD, fill = grp)) + geom_col(width = 0.8) + geom_hline(yintercept = 42, linetype = 2, colour = "grey40") +
  facet_grid(~grp, scales = "free_x", space = "free_x") + scale_fill_manual(values = c("#C0392B", "#2C6E9B"), guide = "none") +
  labs(x = NULL, y = "HRD score") + theme_bw(base_size = 9) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
g7b <- ggplot(tiles, aes(id, row, fill = val)) + geom_tile(colour = "white", linewidth = 0.4) +
  facet_grid(~grp, scales = "free_x", space = "free_x") + scale_fill_manual(values = cols, name = NULL) +
  labs(x = "BRCA2-altered TCGA tumors (n = 43)", y = NULL) + theme_minimal(base_size = 9) +
  theme(axis.text.x = element_blank(), panel.grid = element_blank(), strip.text = element_blank(), legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 4))
save_d(g7a / g7b + plot_layout(heights = c(1, 1.6)), "Figure7_TCGA_oncoprint", 11, 7)

# Figure S1: comparison genes with 95% CI
cmp_ci <- rbindlist(lapply(g_cmp$gene, function(g) {
  d <- pt[get(paste0(g, "_mut_prof")) == TRUE]; y <- d[[paste0(g, "_mut")]]
  b1 <- binom.test(sum(y & d$hyper), sum(d$hyper)); b0 <- binom.test(sum(y & !d$hyper), sum(!d$hyper))
  data.table(gene = g, grp = c("Hypermutated", "Non-hypermutated"), pct = 100 * c(b1$estimate, b0$estimate),
             lo = 100 * c(b1$conf.int[1], b0$conf.int[1]), hi = 100 * c(b1$conf.int[2], b0$conf.int[2]))
}))
cmp_ci[, gene := factor(gene, levels = g_cmp[order(-ratio), gene])]
gs1 <- ggplot(cmp_ci, aes(gene, pct, fill = grp)) + geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(0.8), width = 0.25) +
  geom_text(data = g_cmp, aes(gene, pmax(hyper_pct, nonhyper_pct) + 7, label = paste0(ratio, "x")), inherit.aes = FALSE, size = 3.2) +
  scale_fill_manual(values = c(pal[["red"]], pal[["blue"]])) + labs(x = NULL, y = "Tumors with a driver mutation (%), with 95% CI", fill = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "top")
save_d(gs1, "FigureS1_comparative_genes", 9, 5)
print(cmp_ci)

# ---- Save ------------------------------------------------------------------------
fwrite(d1, file.path(RV, "D1_sensitivity_P80_P85.csv")); fwrite(d2, file.path(RV, "D2_putative_lof.csv"))
fwrite(d3a, file.path(RV, "D3a_table3_nN.csv")); fwrite(d3b, file.path(RV, "D3b_table4_q.csv"))
fwrite(d4, file.path(RV, "D4_spdac_epdca.csv")); fwrite(d5, file.path(RV, "D5_repeats.csv")); fwrite(d5b, file.path(RV, "D5b_repeat_codons.csv"))
fwrite(tflow, file.path(RV, "D6a_tcga_flow.csv")); fwrite(d6c, file.path(RV, "D6c_pole.csv"))
fwrite(rbind(cbind(model = "M1", d7a), cbind(model = "M2", d7b), cbind(model = "nonMSI_subtype", d7e)), file.path(RV, "D7_interaction.csv"))
fwrite(d7c, file.path(RV, "D7c_hrd_high_firth.csv")); fwrite(d7d, file.path(RV, "D7d_subtype.csv"))
fwrite(d8, file.path(RV, "D8_proxy_by_brca2.csv")); fwrite(d9, file.path(RV, "D9_tiers.csv")); fwrite(d9p, file.path(RV, "D9_patient_level.csv"))
fwrite(d10, file.path(RV, "D10_harmonized.csv")); fwrite(cmp_ci, file.path(RV, "D12_figS1_ci.csv"))
cat("\nDone.\n")
sink()
