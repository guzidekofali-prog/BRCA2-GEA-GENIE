# =============================================================================
# BRCA2 in gastric, gastroesophageal junction and esophageal adenocarcinoma - AACR GENIE v20.0-public
# Cohort building and main analysis
#
# This script builds the study cohort (one sample per patient, panels covering
# BRCA1 and BRCA2) and runs the main analyses:
#   A. BRCA2 variant class (truncating / missense driver / ambiguous),
#      recurrent frameshifts, protein position
#   B. Hypermutation proxy: top 5% mutation count within panel, or MMR/POLE driver
#      alteration (from veri/alterations_mmr.tsv), combined definition
#   C. Putative loss-of-function BRCA2 (called "HR-relevant" in object names) =
#      truncating mutation or deep deletion in a non-hypermutated tumor
#   D. Prevalence by age group, sex, tumor site and sample type
#   E. Co-alterations: Fisher tests by hypermutation stratum and Mantel-Haenszel OR
#   F. Firth logistic regression (rare events)
#   G. Between-center heterogeneity (random-effects meta-analysis, I2)
#   H. Overlap of core HR genes (BRCA1/2, PALB2, ATM)
#   I. Exploratory figures: lollipop, oncoprint, age group prevalence, center forest plot
# =============================================================================

# ---- 0. SETTINGS -------------------------------------------------------------
PROJ      <- "~/GENIE_BRCA2"
DATA_DIR  <- file.path(PROJ, "veri")
CLIN_FILE <- "genie_public_clinical_data (2).tsv"
ALT_FILE  <- "alterations_across_samples.tsv"
MMR_FILE  <- "alterations_mmr.tsv"          # MLH1, MSH2, MSH6, PMS2, POLE (optional)
OUT_DIR   <- file.path(PROJ, "cikti_v2")
MIN_CENTER_N <- 100                         # smallest center size kept in the heterogeneity analysis

ADENO_CODES <- c(STAD = "Stomach", DSTAD = "Stomach", ISTAD = "Stomach", TSTAD = "Stomach",
                 MSTAD = "Stomach", PSTAD = "Stomach", USTAD = "Stomach", SSRCC = "Stomach",
                 GEJ = "GEJ", ESCA = "Esophagus", EGC = "Esophagogastric NOS")
GENES <- c("BRCA1", "BRCA2", "PALB2", "ATM", "ERBB2", "TP53", "KRAS", "CDH1", "ARID1A", "PIK3CA")
MMR_GENES <- c("MLH1", "MSH2", "MSH6", "PMS2")
# Pathogenic POLE exonuclease domain hotspots
POLE_HOT <- c("P286R", "V411L", "S459F", "A456P", "S297F", "F367S", "L424I", "M444K", "P436R", "M295R")

# ---- 1. PACKAGES -----------------------------------------------------------
pkgs <- c("data.table", "ggplot2", "logistf", "metafor")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
invisible(lapply(pkgs, library, character.only = TRUE))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
sink(file.path(OUT_DIR, "00_konsol_log_v2.txt"), split = TRUE)
cat("Analysis date:", format(Sys.time()), "\n")

rd <- function(f) { d <- fread(file.path(DATA_DIR, f), sep = "\t", quote = "", na.strings = c("", "NA"))
                    setnames(d, make.names(names(d))); d }
is_np  <- function(x) x == "not profiled" | is.na(x)
is_drv <- function(x) grepl("(driver)", x, fixed = TRUE)
fmt_or <- function(o, l, u) sprintf("%.2f (%.2f-%.2f)", o, l, u)

# ---- 2. DATA AND COHORT --------------------------------------
clin <- rd(CLIN_FILE); alt <- rd(ALT_FILE)
dt <- merge(clin, alt[, !c("Study.ID", "Patient.ID"), with = FALSE], by = "Sample.ID")

has_mmr <- file.exists(file.path(DATA_DIR, MMR_FILE))
if (has_mmr) {
  mmr <- rd(MMR_FILE)
  keep <- c("Sample.ID", grep(paste0("^(", paste(c(MMR_GENES, "POLE"), collapse = "|"), ")"), names(mmr), value = TRUE))
  dt <- merge(dt, mmr[, ..keep], by = "Sample.ID", all.x = TRUE)
  cat("MMR/POLE file found and merged.\n")
} else cat("WARNING: no MMR/POLE file; only the top 5% mutation count rule will be used.\n")

dt <- dt[Oncotree.Code %in% names(ADENO_CODES) & Sample.class == "Tumor"]
dt[, site := factor(ADENO_CODES[Oncotree.Code], levels = c("Stomach", "GEJ", "Esophagus", "Esophagogastric NOS"))]

for (g in GENES) {
  dt[, paste0(g, "_mut_prof") := !is_np(get(paste0(g, "..MUT")))]
  dt[, paste0(g, "_cna_prof") := !is_np(get(paste0(g, "..HOMDEL")))]
  dt[, paste0(g, "_mut")      := is_drv(get(paste0(g, "..MUT")))]
  dt[, paste0(g, "_homdel")   := is_drv(get(paste0(g, "..HOMDEL")))]
  dt[, paste0(g, "_amp")      := is_drv(get(paste0(g, "..AMP")))]
}
for (g in c("BRCA1", "BRCA2", "PALB2", "ATM", "TP53", "CDH1", "ARID1A"))
  dt[, paste0(g, "_alt") := get(paste0(g, "_mut")) | get(paste0(g, "_homdel"))]
dt[, ERBB2_alt := ERBB2_amp]; dt[, KRAS_alt := KRAS_mut | KRAS_amp]; dt[, PIK3CA_alt := PIK3CA_mut]

dt[, prio_primary := Sample.Type == "Primary"]
dt[, seq_year := suppressWarnings(as.numeric(Sequence.Year))]
setorder(dt, Patient.ID, -BRCA2_mut_prof, -BRCA2_cna_prof, -prio_primary, seq_year)
pt <- dt[, .SD[1], by = Patient.ID][BRCA2_mut_prof == TRUE & BRCA1_mut_prof == TRUE]
cat("Study cohort:", nrow(pt), "patients\n")

age_raw <- as.character(pt$Age.at.Which.Sequencing.was.Reported..Years.); age_raw[age_raw == ">89"] <- "90"
pt[, age := suppressWarnings(as.numeric(age_raw))]
pt[, age_grp := cut(age, c(-Inf, 49, 64, 74, Inf), labels = c("<50", "50-64", "65-74", ">=75"))]
pt[, sex := fifelse(Sex %in% c("Male", "Female"), Sex, NA_character_)]
pt[, sample_type := fifelse(Sample.Type %in% c("Primary", "Metastasis"), Sample.Type, "Other/Unknown")]
pt[, mut_count := suppressWarnings(as.numeric(Mutation.Count))]

# ---- 3. HYPERMUTATION DEFINITIONS -------------------------------------------
pt[, H1_top5 := mut_count > quantile(mut_count, 0.95, na.rm = TRUE), by = Sequence.Assay.ID]
pt[is.na(H1_top5), H1_top5 := FALSE]
if (has_mmr) {
  for (g in MMR_GENES) {
    mc <- paste0(g, "..MUT"); hc <- paste0(g, "..HOMDEL")
    pt[, paste0(g, "_alt") := is_drv(get(mc)) | (hc %in% names(pt) & is_drv(get(hc)))]
  }
  pt[, MMR_alt := MLH1_alt | MSH2_alt | MSH6_alt | PMS2_alt]
  pt[, POLE_hot := grepl(paste0("\\b(", paste(POLE_HOT, collapse = "|"), ")\\b"), POLE..MUT)]
  pt[, MMR_prof := !is_np(MLH1..MUT)]
  pt[, H2_mmr_pole := MMR_alt | POLE_hot]
  pt[, hyper := H1_top5 | H2_mmr_pole]
} else {
  pt[, H2_mmr_pole := NA]; pt[, hyper := H1_top5]
}
pt[, hyper_f := factor(fifelse(hyper, "Hypermutated", "Non-hypermutated"),
                       levels = c("Non-hypermutated", "Hypermutated"))]

# ---- 4. BRCA2 VARIANT CLASSIFICATION -----------------------------------------
split_var <- function(x) { v <- trimws(unlist(strsplit(x, ","))); v <- v[!v %in% c("no alteration", "not profiled", "")]
                           sub("\\s*\\(driver\\)", "", v[grepl("(driver)", v, fixed = TRUE)]) }
var_class <- function(v) fifelse(grepl("fs", v), "Frameshift",
                         fifelse(grepl("splice", v, ignore.case = TRUE), "Splice",
                         fifelse(grepl("^[A-Z]\\d+\\*$", v), "Nonsense",
                         fifelse(grepl("^[A-Z]\\d+[A-Z]$", v), "Missense", "Other/ambiguous"))))
bv <- pt[BRCA2_mut == TRUE, .(variant = split_var(BRCA2..MUT)), by = .(Patient.ID, hyper)]
bv[, cls := var_class(variant)]
bv[, pos := suppressWarnings(as.numeric(sub("^[A-Z*X](\\d+).*", "\\1", variant)))]
bv[, codon_key := sub("^([A-Z*X]\\d+).*", "\\1", variant)]
rec <- bv[, .N, by = codon_key][N >= 3, codon_key]
bv[, recurrent := codon_key %in% rec]
bv[, homopolymer_fs := recurrent & cls == "Frameshift"]

pt_cls <- bv[, .(trunc = any(cls %in% c("Frameshift", "Nonsense", "Splice")),
                 missense = any(cls == "Missense"),
                 n_var = .N,
                 only_recurrent_fs = all(homopolymer_fs)), by = Patient.ID]
pt <- merge(pt, pt_cls, by = "Patient.ID", all.x = TRUE)
for (v in c("trunc", "missense", "only_recurrent_fs")) pt[is.na(get(v)), (v) := FALSE]
pt[is.na(n_var), n_var := 0L]
pt[, biallelic_proxy := BRCA2_alt & (n_var >= 2 | (BRCA2_mut & BRCA2_homdel))]
# Putative LoF BRCA2: truncating mutation or deep deletion in a non-hypermutated tumor
pt[, BRCA2_hr := (trunc | BRCA2_homdel) & !hyper]
pt[, BRCA2_grp := factor(fifelse(!BRCA2_alt, "WT",
                         fifelse(hyper, "BRCA2 in hypermutated", "BRCA2 in non-hypermutated")),
                         levels = c("WT", "BRCA2 in non-hypermutated", "BRCA2 in hypermutated"))]

var_tab <- bv[, .N, by = .(hiper = hyper, sinif = cls)][order(hiper, -N)]
rec_tab <- bv[recurrent == TRUE, .(n = .N, hipermutasyonlu = sum(hyper), sinif = cls[1]), by = codon_key][order(-n)]

# ---- 5. PREVALENCE TABLES -------------------------------------------------
prop_ci <- function(x, n) { if (n == 0) return(c(NA, NA, NA)); ci <- binom.test(x, n)$conf.int; c(100 * x / n, 100 * ci) }
prev_by <- function(var, outcome) {
  pt[!is.na(get(var)), {
    r <- prop_ci(sum(get(outcome)), .N)
    .(n_alt = sum(get(outcome)), n = .N, yuzde = round(r[1], 2), ci_alt = round(r[2], 2), ci_ust = round(r[3], 2))
  }, by = var][, .(degisken = var, kategori = as.character(get(var)), sonuc = outcome, n_alt, n, yuzde, ci_alt, ci_ust)]
}
pt[, tum := "Tüm hastalar"]
prev2 <- rbindlist(lapply(c("BRCA2_alt", "BRCA2_hr"), function(o)
  rbindlist(lapply(c("tum", "site", "age_grp", "sex", "sample_type", "hyper_f"), prev_by, outcome = o))))

ozet <- data.table(
  olcu = c("BRCA2 altere (herhangi)", "  hipermutasyonlu tümörde", "  hipermutasyonsuz tümörde",
           "  truncating veya derin delesyon", "  sadece driver missense",
           "  sadece tekrarlayan homopolimer frameshift", "  bialelik proxy",
           "HR açısından anlamlı BRCA2 (truncating/HOMDEL, hipermutasyonsuz)"),
  n = c(sum(pt$BRCA2_alt), sum(pt$BRCA2_alt & pt$hyper), sum(pt$BRCA2_alt & !pt$hyper),
        sum(pt$BRCA2_alt & (pt$trunc | pt$BRCA2_homdel)), sum(pt$BRCA2_alt & pt$missense & !pt$trunc & !pt$BRCA2_homdel),
        sum(pt$BRCA2_alt & pt$only_recurrent_fs), sum(pt$biallelic_proxy), sum(pt$BRCA2_hr)))
ozet[, yuzde_kohort := round(100 * n / nrow(pt), 2)]
ozet[, yuzde_brca2 := round(100 * n / sum(pt$BRCA2_alt), 1)]
if (has_mmr) {
  ozet <- rbind(ozet, data.table(
    olcu = c("MMR geni driver alterasyonu (kohort)", "POLE ekzonükleaz hotspot (kohort)",
             "BRCA2 altere hastalarda MMR/POLE alterasyonu"),
    n = c(sum(pt$MMR_alt, na.rm = TRUE), sum(pt$POLE_hot, na.rm = TRUE), sum(pt$BRCA2_alt & pt$H2_mmr_pole, na.rm = TRUE)),
    yuzde_kohort = NA, yuzde_brca2 = NA))
  ozet[olcu == "BRCA2 altere hastalarda MMR/POLE alterasyonu", yuzde_brca2 := round(100 * n / sum(pt$BRCA2_alt), 1)]
}

# ---- 6. CO-ALTERATIONS: STRATIFIED FISHER + MANTEL-HAENSZEL -------------------
co_genes <- c("ERBB2", "TP53", "KRAS", "CDH1", "ARID1A", "PIK3CA", "ATM", "PALB2", "BRCA1")
co_strat <- rbindlist(lapply(co_genes, function(g) {
  dd <- if (g == "ERBB2") pt[ERBB2_cna_prof == TRUE] else pt[get(paste0(g, "_mut_prof")) == TRUE]
  one <- function(d, lab) {
    tb <- table(factor(d[[paste0(g, "_alt")]], c(FALSE, TRUE)), factor(d$BRCA2_alt, c(FALSE, TRUE)))
    ft <- fisher.test(tb)
    data.table(gen = g, tabaka = lab, n = nrow(d),
               brca2_alt_ile = sprintf("%d/%d (%.1f%%)", tb[2, 2], sum(tb[, 2]), 100 * tb[2, 2] / max(1, sum(tb[, 2]))),
               brca2_wt_ile  = sprintf("%d/%d (%.1f%%)", tb[2, 1], sum(tb[, 1]), 100 * tb[2, 1] / sum(tb[, 1])),
               OR = unname(ft$estimate), OR_alt = ft$conf.int[1], OR_ust = ft$conf.int[2], p = ft$p.value)
  }
  arr <- table(factor(dd[[paste0(g, "_alt")]], c(FALSE, TRUE)), factor(dd$BRCA2_alt, c(FALSE, TRUE)), dd$hyper_f)
  mh <- tryCatch(mantelhaen.test(arr, exact = TRUE), error = function(e) NULL)
  rbind(one(dd[hyper == FALSE], "Hipermutasyonsuz"), one(dd[hyper == TRUE], "Hipermutasyonlu"),
        data.table(gen = g, tabaka = "MH (hipermutasyona göre düzeltilmiş)", n = nrow(dd),
                   brca2_alt_ile = NA, brca2_wt_ile = NA,
                   OR = if (is.null(mh)) NA else unname(mh$estimate),
                   OR_alt = if (is.null(mh)) NA else mh$conf.int[1],
                   OR_ust = if (is.null(mh)) NA else mh$conf.int[2],
                   p = if (is.null(mh)) NA else mh$p.value))
}))
co_strat[, q_BH := p.adjust(p, "BH"), by = tabaka]

# ---- 7. FIRTH LOGISTIC REGRESSION -------------------------------------------
top_c <- pt[, .N, by = Center][N >= MIN_CENTER_N, Center]
pt[, center_g := factor(fifelse(Center %in% top_c, Center, "Other"))]
pt[, center_g := relevel(center_g, ref = names(sort(table(center_g), decreasing = TRUE))[1])]
mdat <- pt[!is.na(age_grp) & !is.na(sex)]
firth_tab <- function(fml, lab) {
  m <- logistf(fml, data = mdat)
  data.table(model = lab, terim = names(coef(m)), OR = exp(coef(m)),
             OR_alt = exp(m$ci.lower), OR_ust = exp(m$ci.upper), p = m$prob,
             n = m$n, olay = sum(model.response(model.frame(fml, mdat))))[!grepl("Intercept|center_g", terim)]
}
firth <- rbind(
  firth_tab(BRCA2_alt ~ age_grp + sex + site + hyper + center_g, "M1: BRCA2 herhangi"),
  firth_tab(BRCA2_hr  ~ age_grp + sex + site + sample_type + center_g, "M2: HR-anlamlı BRCA2 (hipermutasyonsuz)")
)

# ---- 8. CENTER HETEROGENEITY (random effects, logit-transformed proportions) -------
cen <- pt[, .(x_any = sum(BRCA2_alt), x_hr = sum(BRCA2_hr), n = .N,
              hyper_pct = round(100 * mean(hyper), 1)), by = Center][n >= MIN_CENTER_N][order(-n)]
meta_one <- function(xcol) {
  es <- escalc(measure = "PLO", xi = cen[[xcol]], ni = cen$n)
  m <- rma(yi, vi, data = es, method = "REML")
  data.table(sonuc = xcol, havuz_yuzde = round(100 * plogis(m$b[1]), 2),
             ci_alt = round(100 * plogis(m$ci.lb), 2), ci_ust = round(100 * plogis(m$ci.ub), 2),
             I2 = round(m$I2, 1), tau2 = signif(m$tau2, 3), Q_p = signif(m$QEp, 3), k = m$k)
}
meta <- rbind(meta_one("x_any"), meta_one("x_hr"))

# ---- 9. OVERLAP OF CORE HR GENES ----------------------------------------
hrr <- pt[, .N, by = .(BRCA1 = BRCA1_alt, BRCA2 = BRCA2_alt, PALB2 = PALB2_alt, ATM = ATM_alt)][order(-N)]
hrr_any <- pt[, .(herhangi_HRR = sum(BRCA1_alt | BRCA2_alt | PALB2_alt | ATM_alt),
                  hipermutasyonsuz_HRR = sum((BRCA1_alt | BRCA2_alt | PALB2_alt | ATM_alt) & !hyper), n = .N)]

# ---- 10. OUTPUTS ----------------------------------------------------------
fwrite(ozet,      file.path(OUT_DIR, "v2_01_brca2_ozet.csv"))
fwrite(var_tab,   file.path(OUT_DIR, "v2_02_varyant_sinifi.csv"))
fwrite(rec_tab,   file.path(OUT_DIR, "v2_03_tekrarlayan_varyantlar.csv"))
fwrite(prev2,     file.path(OUT_DIR, "v2_04_prevalans_alt_gruplar.csv"))
fwrite(co_strat,  file.path(OUT_DIR, "v2_05_ko_alterasyon_tabakali.csv"))
fwrite(firth,     file.path(OUT_DIR, "v2_06_firth_lojistik.csv"))
fwrite(cen,       file.path(OUT_DIR, "v2_07_merkez_prevalans.csv"))
fwrite(meta,      file.path(OUT_DIR, "v2_08_merkez_meta_analiz.csv"))
fwrite(hrr,       file.path(OUT_DIR, "v2_09_hrr_ortusme.csv"))
fwrite(hrr_any,   file.path(OUT_DIR, "v2_10_hrr_toplam.csv"))
fwrite(bv,        file.path(OUT_DIR, "v2_11_brca2_varyant_listesi.csv"))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "v2_12_sessionInfo.txt"))

th <- theme_minimal(base_size = 12)
# F3: BRCA2 lollipop (protein position)
lp <- bv[!is.na(pos), .N, by = .(pos, cls, hiper = fifelse(hyper, "Hypermutated", "Non-hypermutated"))]
g3 <- ggplot(lp, aes(pos, N, colour = cls)) +
  geom_segment(aes(xend = pos, y = 0, yend = N), colour = "grey70") + geom_point(size = 2.4) +
  facet_wrap(~hiper, ncol = 1) + scale_x_continuous(limits = c(0, 3418), breaks = seq(0, 3400, 500)) +
  labs(x = "BRCA2 amino acid position (1-3418)", y = "Number of tumors", colour = "Variant class") + th
ggsave(file.path(OUT_DIR, "F3_brca2_lollipop.png"), g3, width = 9, height = 5, dpi = 300)

# F4: BRCA2 by age group (any and putative LoF)
f4 <- prev2[degisken == "age_grp"]
f4[, sonuc := factor(sonuc, c("BRCA2_alt", "BRCA2_hr"), c("Any BRCA2 alteration", "HR-relevant BRCA2 (non-hypermutated)"))]
f4[, kategori := factor(kategori, c("<50", "50-64", "65-74", ">=75"))]
g4 <- ggplot(f4, aes(kategori, yuzde, fill = sonuc)) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_errorbar(aes(ymin = ci_alt, ymax = ci_ust), position = position_dodge(0.75), width = 0.2) +
  scale_fill_manual(values = c("#2C6E9B", "#E08E45")) +
  labs(x = "Age group (years)", y = "Prevalence (%)", fill = NULL) + th + theme(legend.position = "top")
ggsave(file.path(OUT_DIR, "F4_yas_grubu_prevalans.png"), g4, width = 7, height = 4.5, dpi = 300)

# F5: Oncoprint-like plot (BRCA2-altered patients)
op_genes <- c("BRCA2", "BRCA1", "PALB2", "ATM", "ERBB2", "TP53", "KRAS", "CDH1", "ARID1A", "PIK3CA")
op <- pt[BRCA2_alt == TRUE, c("Patient.ID", "hyper", paste0(op_genes, "_alt")), with = FALSE]
setorderv(op, c("hyper", paste0(op_genes, "_alt")), order = -1)
op[, idx := .I]
opl <- melt(op, id.vars = c("Patient.ID", "idx", "hyper"), variable.name = "gen", value.name = "alt")
opl[, gen := factor(sub("_alt", "", gen), levels = rev(op_genes))]
opl[, durum := fifelse(alt, fifelse(hyper, "Altered (hypermutated)", "Altered (non-hypermutated)"), "None")]
g5 <- ggplot(opl, aes(idx, gen, fill = durum)) + geom_tile(colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = c("None" = "grey92", "Altered (non-hypermutated)" = "#2C6E9B",
                               "Altered (hypermutated)" = "#C0392B")) +
  labs(x = sprintf("BRCA2-altered tumors (n = %d)", nrow(op)), y = NULL, fill = NULL) +
  th + theme(axis.text.x = element_blank(), panel.grid = element_blank(), legend.position = "top")
ggsave(file.path(OUT_DIR, "F5_oncoprint.png"), g5, width = 10, height = 4, dpi = 300)

# F6: Prevalence by center, forest plot
f6 <- cen[, { r <- t(mapply(prop_ci, x_any, n)); .(Center, n, yuzde = r[, 1], lo = r[, 2], hi = r[, 3]) }]
f6[, Center := factor(Center, levels = Center[order(yuzde)])]
g6 <- ggplot(f6, aes(yuzde, Center)) + geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2) +
  geom_point(aes(size = n), colour = "#2C6E9B") +
  geom_vline(xintercept = meta[sonuc == "x_any", havuz_yuzde], linetype = 2) +
  labs(x = "BRCA2 alteration prevalence (%)", y = NULL, size = "N") + th
ggsave(file.path(OUT_DIR, "F6_merkez_forest.png"), g6, width = 7, height = 5, dpi = 300)

# ---- 11. CONSOLE SUMMARY ------------------------------------------------------
cat("\n=== BRCA2 SUMMARY ===\n"); print(ozet)
cat("\n=== VARIANT CLASS ===\n"); print(var_tab)
cat("\n=== RECURRENT VARIANTS (>=3 tumors) ===\n"); print(rec_tab)
cat("\n=== PREVALENCE BY SUBGROUP ===\n"); print(prev2)
cat("\n=== CO-ALTERATIONS (stratified + MH) ===\n"); print(co_strat)
cat("\n=== FIRTH LOGISTIC ===\n"); print(firth)
cat("\n=== CENTER META-ANALYSIS ===\n"); print(cen); print(meta)
cat("\n=== HR GENE OVERLAP ===\n"); print(hrr); print(hrr_any)
cat("\nOutputs:", normalizePath(OUT_DIR), "\n")
sink()
