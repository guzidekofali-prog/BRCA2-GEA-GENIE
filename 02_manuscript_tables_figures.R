# =============================================================================
# MANUSCRIPT TABLES AND FIGURES
# BRCA2 alterations in gastroesophageal adenocarcinoma (AACR GENIE v20.0)
#
# Runs script 01 first (same cohort and definitions), then writes the manuscript
# tables (Word + CSV) and Figures 1-4 (PNG + TIFF).
# All numbers are computed in this R session.
# =============================================================================

while (sink.number() > 0) sink()
PROJ <- "~/GENIE_BRCA2"
setwd(PROJ)
source("01_genie_cohort_main_analysis.R", encoding = "UTF-8")   # pt, dt, clin, co_strat, firth, meta, cen ...

MS_DIR <- file.path(PROJ, "makale")
dir.create(MS_DIR, showWarnings = FALSE)

for (p in c("flextable", "officer", "patchwork")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(flextable); library(officer); library(patchwork)

sink(file.path(MS_DIR, "00_makale_log.txt"), split = TRUE)
cat("Tables and figures date:", format(Sys.time()), "\n")

# Helpers ---------------------------------------------------------------
npct  <- function(x, n) sprintf("%d (%.1f)", x, 100 * x / n)
medq  <- function(x) sprintf("%.0f (%.0f-%.0f)", median(x, na.rm = TRUE),
                             quantile(x, .25, na.rm = TRUE), quantile(x, .75, na.rm = TRUE))
pfmt  <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
ci_txt <- function(o, l, u) sprintf("%.2f (%.2f-%.2f)", o, l, u)
site_en <- c("Stomach" = "Stomach", "GEJ" = "Gastroesophageal junction",
             "Esophagus" = "Esophagus", "Esophagogastric NOS" = "Esophagogastric, not specified")

pt[, grp := factor(fifelse(BRCA2_alt, "BRCA2-altered", "BRCA2 wild-type"),
                   levels = c("BRCA2 wild-type", "BRCA2-altered"))]
pt[, site_en := factor(site_en[as.character(site)], levels = site_en)]
pt[, sex_en := factor(sex, levels = c("Male", "Female"))]
pt[, stype_en := factor(fifelse(sample_type == "Other/Unknown", "Other or unknown", sample_type),
                        levels = c("Primary", "Metastasis", "Other or unknown"))]
pt[, age_en := factor(age_grp, levels = c("<50", "50-64", "65-74", ">=75"),
                      labels = c("<50", "50-64", "65-74", "75 or older"))]

# ---- FLOW (for Figure 1) ---------------------------------------------------
adeno_s <- clin[Oncotree.Code %in% names(ADENO_CODES)]
flow_ms <- data.table(
  step = c("Esophagogastric cancer samples in GENIE v20.0",
           "Adenocarcinoma histology",
           "Tumor tissue samples (cfDNA removed)",
           "One sample per patient",
           "BRCA1 and BRCA2 covered by the sequencing panel"),
  samples  = c(nrow(clin), nrow(adeno_s), nrow(dt), uniqueN(dt$Patient.ID), nrow(pt)),
  patients = c(uniqueN(clin$Patient.ID), uniqueN(adeno_s$Patient.ID), uniqueN(dt$Patient.ID),
               uniqueN(dt$Patient.ID), nrow(pt)))
cat("\n=== FLOW ===\n"); print(flow_ms)

# ---- TABLE 1: Patient characteristics -------------------------------------------
cat_rows <- function(var, label) {
  d <- pt[!is.na(get(var))]
  tb <- table(d[[var]], d$grp)
  p <- tryCatch(chisq.test(tb)$p.value, warning = function(w) fisher.test(tb, simulate.p.value = TRUE, B = 1e5)$p.value)
  lev <- rownames(tb)
  rbind(data.table(Characteristic = label, Overall = "", WT = "", ALT = "", P = pfmt(p)),
        data.table(Characteristic = paste0("   ", lev),
                   Overall = mapply(npct, rowSums(tb), sum(tb)),
                   WT  = mapply(npct, tb[, 1], sum(tb[, 1])),
                   ALT = mapply(npct, tb[, 2], sum(tb[, 2])), P = ""))
}
pt[, hyper_en := factor(fifelse(hyper, "Yes", "No"), levels = c("No", "Yes"))]
pt[, mmr_en := factor(fifelse(H2_mmr_pole, "Yes", "No"), levels = c("No", "Yes"))]
t1 <- rbind(
  data.table(Characteristic = "Number of patients", Overall = as.character(nrow(pt)),
             WT = as.character(sum(pt$grp == "BRCA2 wild-type")),
             ALT = as.character(sum(pt$grp == "BRCA2-altered")), P = ""),
  data.table(Characteristic = "Age at sequencing, median (IQR), years",
             Overall = medq(pt$age), WT = medq(pt[grp == "BRCA2 wild-type", age]),
             ALT = medq(pt[grp == "BRCA2-altered", age]),
             P = pfmt(wilcox.test(age ~ grp, data = pt)$p.value)),
  cat_rows("age_en", "Age group, n (%)"),
  cat_rows("sex_en", "Sex, n (%)"),
  cat_rows("site_en", "Tumor site, n (%)"),
  cat_rows("stype_en", "Sample type, n (%)"),
  data.table(Characteristic = "Mutation count, median (IQR)", Overall = medq(pt$mut_count),
             WT = medq(pt[grp == "BRCA2 wild-type", mut_count]), ALT = medq(pt[grp == "BRCA2-altered", mut_count]),
             P = pfmt(wilcox.test(mut_count ~ grp, data = pt)$p.value)),
  cat_rows("mmr_en", "MMR gene or POLE hotspot alteration, n (%)"),
  cat_rows("hyper_en", "Hypermutated tumor, n (%)")
)
setnames(t1, c("WT", "ALT", "P"), c("BRCA2 wild-type", "BRCA2-altered", "P value"))
cat("\n=== TABLE 1 ===\n"); print(t1)

# ---- TABLE 2: Features of BRCA2 alterations ---------------------------
nb <- sum(pt$BRCA2_alt)
prev_line <- function(lab, x, n) { ci <- binom.test(x, n)$conf.int
  data.table(Item = lab, n = sprintf("%d/%d", x, n),
             Value = sprintf("%.2f%% (%.2f-%.2f)", 100 * x / n, 100 * ci[1], 100 * ci[2])) }
cls_n <- function(cc) sum(pt$BRCA2_alt & pt$Patient.ID %in% bv[cls %in% cc, Patient.ID])
t2 <- rbind(
  data.table(Item = "A. Frequency, % (95% CI)", n = "", Value = ""),
  prev_line("   Any BRCA2 alteration, all patients", nb, nrow(pt)),
  rbindlist(lapply(levels(pt$site_en), function(s) prev_line(paste0("   ", s), sum(pt[site_en == s, BRCA2_alt]), pt[site_en == s, .N]))),
  prev_line("   In hypermutated tumors", sum(pt[hyper == TRUE, BRCA2_alt]), pt[hyper == TRUE, .N]),
  prev_line("   In non-hypermutated tumors", sum(pt[hyper == FALSE, BRCA2_alt]), pt[hyper == FALSE, .N]),
  prev_line("   HR-relevant BRCA2 alteration*, all patients", sum(pt$BRCA2_hr), nrow(pt)),
  data.table(Item = "B. Type of alteration among BRCA2-altered patients, n (%)", n = "", Value = ""),
  data.table(Item = c("   Truncating mutation or deep deletion", "   Driver missense mutation only",
                      "   Recurrent frameshift at a repeat site only", "   Two or more hits (possible biallelic loss)"),
             n = "",
             Value = c(npct(sum(pt$BRCA2_alt & (pt$trunc | pt$BRCA2_homdel)), nb),
                       npct(sum(pt$BRCA2_alt & pt$missense & !pt$trunc & !pt$BRCA2_homdel), nb),
                       npct(sum(pt$BRCA2_alt & pt$only_recurrent_fs), nb),
                       npct(sum(pt$biallelic_proxy), nb))),
  data.table(Item = "C. Tumor context among BRCA2-altered patients, n (%)", n = "", Value = ""),
  data.table(Item = c("   Hypermutated tumor", "   MMR gene or POLE hotspot alteration"), n = "",
             Value = c(npct(sum(pt$BRCA2_alt & pt$hyper), nb), npct(sum(pt$BRCA2_alt & pt$H2_mmr_pole), nb)))
)
setnames(t2, c("Item", "n", "Value"), c("Measure", "n/N", "Value"))
cat("\n=== TABLE 2 ===\n"); print(t2)

# ---- TABLE 3: Co-alterations ---------------------------------------------
crude <- rbindlist(lapply(co_genes, function(g) {
  dd <- if (g == "ERBB2") pt[ERBB2_cna_prof == TRUE] else pt[get(paste0(g, "_mut_prof")) == TRUE]
  tb <- table(factor(dd[[paste0(g, "_alt")]], c(FALSE, TRUE)), factor(dd$BRCA2_alt, c(FALSE, TRUE)))
  ft <- fisher.test(tb)
  data.table(gen = g, alt_pct = sprintf("%.1f", 100 * tb[2, 2] / sum(tb[, 2])),
             wt_pct = sprintf("%.1f", 100 * tb[2, 1] / sum(tb[, 1])),
             crude_or = unname(ft$estimate), crude_lo = ft$conf.int[1], crude_hi = ft$conf.int[2], crude_p = ft$p.value)
}))
crude[, crude_q := p.adjust(crude_p, "BH")]
nh <- co_strat[tabaka == "Hipermutasyonsuz"]
mh <- co_strat[grepl("^MH", tabaka)]
gene_lab <- c(ERBB2 = "ERBB2 amplification", TP53 = "TP53", KRAS = "KRAS", CDH1 = "CDH1",
              ARID1A = "ARID1A", PIK3CA = "PIK3CA", ATM = "ATM", PALB2 = "PALB2", BRCA1 = "BRCA1")
t3 <- data.table(
  Gene = gene_lab[crude$gen],
  `All tumors: BRCA2-altered vs wild-type, %` = paste(crude$alt_pct, "vs", crude$wt_pct),
  `All tumors: OR (95% CI)` = ci_txt(crude$crude_or, crude$crude_lo, crude$crude_hi),
  `All tumors: q` = pfmt(crude$crude_q),
  `Non-hypermutated: BRCA2-altered vs wild-type` = paste(nh$brca2_alt_ile, "vs", nh$brca2_wt_ile),
  `Non-hypermutated: OR (95% CI)` = ci_txt(nh$OR, nh$OR_alt, nh$OR_ust),
  `Non-hypermutated: q` = pfmt(nh$q_BH),
  `Adjusted for hypermutation: OR (95% CI)` = ci_txt(mh$OR, mh$OR_alt, mh$OR_ust),
  `Adjusted: q` = pfmt(mh$q_BH))
cat("\n=== TABLE 3 ===\n"); print(t3)

# ---- SUPPLEMENTARY TABLES -----------------------------------------------------------
term_lab <- c("age_grp50-64" = "Age 50-64 (vs <50)", "age_grp65-74" = "Age 65-74 (vs <50)",
              "age_grp>=75" = "Age 75 or older (vs <50)", "sexMale" = "Male sex",
              "siteGEJ" = "Gastroesophageal junction (vs stomach)", "siteEsophagus" = "Esophagus (vs stomach)",
              "siteEsophagogastric NOS" = "Esophagogastric, not specified (vs stomach)",
              "hyperTRUE" = "Hypermutated tumor", "sample_typePrimary" = "Primary tumor sample (vs metastasis)",
              "sample_typeOther/Unknown" = "Other or unknown sample (vs metastasis)")
s1 <- firth[, .(Model = fifelse(grepl("^M1", model), "Any BRCA2 alteration", "HR-relevant BRCA2 alteration"),
                Variable = term_lab[terim], `OR (95% CI)` = ci_txt(OR, OR_alt, OR_ust), `P value` = pfmt(p))]
s2 <- cen[, .(Center, Patients = n, `Any BRCA2, n` = x_any, `HR-relevant BRCA2, n` = x_hr,
              `Hypermutated, %` = hyper_pct)]
s2b <- meta[, .(Outcome = fifelse(sonuc == "x_any", "Any BRCA2", "HR-relevant BRCA2"),
                `Pooled %, (95% CI)` = sprintf("%.2f (%.2f-%.2f)", havuz_yuzde, ci_alt, ci_ust),
                `I2, %` = I2, `Q test P` = pfmt(Q_p), Centers = k)]
cat("\n=== TABLE S1 ===\n"); print(s1); cat("\n=== TABLE S2 ===\n"); print(s2); print(s2b)

# ---- WORD OUTPUT ----------------------------------------------------------
ft_make <- function(d, title, foot = NULL) {
  f <- flextable(d) |> theme_booktabs() |> fontsize(size = 9, part = "all") |> font(fontname = "Arial", part = "all") |>
    bold(part = "header") |> align(j = -1, align = "center", part = "all") |> autofit() |>
    set_caption(caption = title)
  if (!is.null(foot)) f <- add_footer_lines(f, foot) |> fontsize(size = 8, part = "footer")
  f
}
doc <- read_docx() |>
  body_add_flextable(ft_make(t1, "Table 1. Characteristics of the patients by BRCA2 status",
    "IQR, interquartile range; MMR, mismatch repair. Hypermutated: mutation count in the top 5% within the same sequencing panel, or a driver alteration in MLH1, MSH2, MSH6 or PMS2, or a POLE exonuclease hotspot mutation.")) |>
  body_add_break() |>
  body_add_flextable(ft_make(t2, "Table 2. Frequency and type of BRCA2 alterations",
    "*HR-relevant: truncating mutation (frameshift, nonsense or splice site) or deep deletion in a tumor that is not hypermutated. CI, confidence interval.")) |>
  body_add_break() |>
  body_add_flextable(ft_make(t3, "Table 3. Genes altered together with BRCA2",
    "OR, odds ratio from Fisher's exact test; adjusted OR from the Mantel-Haenszel test stratified by hypermutation status. q, P value corrected for multiple testing (Benjamini-Hochberg). ERBB2 was analyzed only in samples tested for copy number changes.")) |>
  body_add_break() |>
  body_add_flextable(ft_make(s1, "Supplementary Table S1. Factors linked to BRCA2 alteration (Firth logistic regression)",
    "Both models were also adjusted for sequencing center.")) |>
  body_add_break() |>
  body_add_flextable(ft_make(s2, "Supplementary Table S2. BRCA2 alterations by sequencing center (centers with 100 or more patients)")) |>
  body_add_flextable(ft_make(s2b, "Pooled estimate (random-effects model)"))
print(doc, target = file.path(MS_DIR, "Tables.docx"))
fwrite(t1, file.path(MS_DIR, "Table1.csv")); fwrite(t2, file.path(MS_DIR, "Table2.csv"))
fwrite(t3, file.path(MS_DIR, "Table3.csv")); fwrite(s1, file.path(MS_DIR, "TableS1.csv"))
fwrite(s2, file.path(MS_DIR, "TableS2.csv")); fwrite(flow_ms, file.path(MS_DIR, "Flow.csv"))

# ---- FIGURES --------------------------------------------------------------
pal <- c(blue = "#2C6E9B", red = "#C0392B", grey = "grey80", orange = "#E08E45")
thm <- theme_classic(base_size = 11, base_family = "sans") +
  theme(plot.title = element_text(face = "bold", size = 12), legend.position = "top")
save_fig <- function(g, name, w, h) {
  ggsave(file.path(MS_DIR, paste0(name, ".png")),  g, width = w, height = h, dpi = 300)
  ggsave(file.path(MS_DIR, paste0(name, ".tiff")), g, width = w, height = h, dpi = 300, compression = "lzw")
}

# Figure 1: Flow diagram (the final version with exclusion boxes is made in script 06)
fl <- copy(flow_ms)[, y := rev(seq_len(.N))]
fl[, lab := sprintf("%s\n%s samples, %s patients", step, format(samples, big.mark = ","), format(patients, big.mark = ","))]
fl[.N, lab := sprintf("%s\n%s patients (study cohort)", step, format(patients, big.mark = ","))]
g1 <- ggplot(fl, aes(0, y)) +
  geom_segment(data = fl[-.N], aes(x = 0, xend = 0, y = y - 0.25, yend = y - 0.75),
               arrow = arrow(length = unit(0.18, "cm")), colour = "grey40") +
  geom_label(aes(label = lab), size = 3.3, label.padding = unit(0.45, "lines"), fill = "white") +
  xlim(-1, 1) + theme_void()
save_fig(g1, "Figure1_flow", 6, 6)

# Figure 2: Frequency (A: by tumor site, colored by hypermutation; B: by hypermutation status)
f2a <- pt[, .(n = .N, hyper_brca2 = sum(BRCA2_alt & hyper), nonhyper_brca2 = sum(BRCA2_alt & !hyper)), by = site_en]
f2a <- melt(f2a, id.vars = c("site_en", "n"), variable.name = "ctx", value.name = "x")
f2a[, pct := 100 * x / n]
f2a[, ctx := factor(ctx, c("hyper_brca2", "nonhyper_brca2"),
                    c("BRCA2 alteration in a hypermutated tumor", "BRCA2 alteration in a non-hypermutated tumor"))]
tot2a <- f2a[, .(pct = sum(pct), x = sum(x), n = n[1]), by = site_en]
g2a <- ggplot(f2a, aes(site_en, pct, fill = ctx)) + geom_col(width = 0.65) +
  geom_text(data = tot2a, aes(site_en, pct, label = sprintf("%.1f%%\n(%d/%d)", pct, x, n)), inherit.aes = FALSE,
            vjust = -0.3, size = 3) +
  scale_fill_manual(values = c(pal[["red"]], pal[["blue"]]), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  labs(title = "A. BRCA2 alterations by tumor site", x = NULL, y = "Patients with a BRCA2 alteration (%)") + thm +
  guides(fill = guide_legend(ncol = 1))
f2b <- pt[, { ci <- binom.test(sum(BRCA2_alt), .N)$conf.int
              .(pct = 100 * mean(BRCA2_alt), lo = 100 * ci[1], hi = 100 * ci[2], x = sum(BRCA2_alt), n = .N) }, by = hyper_f]
g2b <- ggplot(f2b, aes(hyper_f, pct, fill = hyper_f)) + geom_col(width = 0.55) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15) +
  geom_text(aes(y = hi, label = sprintf("%.1f%%\n(%d/%d)", pct, x, n)), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c(pal[["blue"]], pal[["red"]]), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(title = "B. BRCA2 alterations by hypermutation status", x = NULL, y = "Patients with a BRCA2 alteration (%)") + thm
save_fig(g2a + g2b + plot_layout(widths = c(1.6, 1)), "Figure2_frequency", 11, 5)

# Figure 3: Lollipop (the final version with domains is made in script 06)
lp <- bv[!is.na(pos), .N, by = .(pos, cls, ctx = fifelse(hyper, "Hypermutated tumors", "Non-hypermutated tumors"))]
lp[, ctx := factor(ctx, c("Non-hypermutated tumors", "Hypermutated tumors"))]
lp[, cls := factor(cls, c("Frameshift", "Nonsense", "Splice", "Missense", "Other/ambiguous"))]
lab_lp <- lp[N >= 3, .(N = sum(N)), by = .(pos, ctx)]
lab_lp <- merge(lab_lp, bv[, .(codon_key = codon_key[1]), by = pos], by = "pos")
g3 <- ggplot(lp, aes(pos, N)) +
  geom_segment(aes(xend = pos, y = 0, yend = N), colour = "grey65", linewidth = 0.4) +
  geom_point(aes(fill = cls), shape = 21, size = 2.6, colour = "black", stroke = 0.3) +
  geom_text(data = lab_lp, aes(pos, N, label = codon_key), vjust = -0.9, size = 2.8) +
  facet_wrap(~ctx, ncol = 1) +
  scale_x_continuous(limits = c(0, 3418), breaks = seq(0, 3400, 500)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
  scale_fill_manual(values = c(Frameshift = "#2C6E9B", Nonsense = "#1B3A57", Splice = "#E08E45",
                               Missense = "#4DAF4A", `Other/ambiguous` = "grey70"), name = "Mutation type") +
  labs(x = "Position in the BRCA2 protein (amino acid 1 to 3418)", y = "Number of tumors") + thm
save_fig(g3, "Figure3_lollipop", 10, 6)

# Figure 4: Co-alteration forest plot (the final three-panel version is made in script 06)
f4 <- rbind(crude[, .(gen, set = "All tumors", OR = crude_or, lo = crude_lo, hi = crude_hi, q = crude_q)],
            nh[, .(gen, set = "Non-hypermutated tumors only", OR, lo = OR_alt, hi = OR_ust, q = q_BH)])
f4 <- f4[is.finite(OR) & OR > 0 & is.finite(hi)]
f4[, gen := factor(gene_lab[gen], levels = rev(gene_lab))]
f4[, sig := fifelse(q < 0.05, "q < 0.05", "Not significant")]
g4 <- ggplot(f4, aes(OR, gen)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = pmax(lo, 0.05), xmax = pmin(hi, 40)), height = 0.25, colour = "grey30") +
  geom_point(aes(fill = sig), shape = 21, size = 2.8) +
  facet_wrap(~set) + scale_x_log10(limits = c(0.05, 40), breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 10, 25)) +
  scale_fill_manual(values = c(`q < 0.05` = pal[["red"]], `Not significant` = "white"), name = NULL) +
  labs(x = "Odds ratio (BRCA2-altered vs wild-type), log scale", y = NULL) + thm
save_fig(g4, "Figure4_coalterations", 10, 5)

cat("\nFiles:", normalizePath(MS_DIR), "\n"); print(list.files(MS_DIR))
sink()
