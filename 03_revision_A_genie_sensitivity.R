# =============================================================================
# REVISION ANALYSES - PART A (GENIE v20.0)
# Additional analyses done in response to editorial comments.
# Script 01 is run first so that the cohort and definitions are the same.
#   A0 Missing data and model N / number of events
#   A1 Circularity: mutation count after removing BRCA2 variants
#   A2 Sensitivity analyses of the hypermutation definition (P90-P99, MMR/POLE only,
#      mutation count only, continuous mutation burden)
#   A3 Analysis without the ambiguous *2602 / X2602_splice variant
#   A4 Robustness: leave-one-center-out, MSK only, primary, metastasis, mutations only,
#      copy number profiled samples, alternative sample selection
#   A5 Prevalence by sequencing panel
#   A6 BRCA2 variant spectrum in hypermutated vs non-hypermutated tumors (per patient)
#   A7 Comparison genes (BRCA1, PALB2, ATM, ARID1A, TP53 ...)
#   A8 Co-alterations: Firth models adjusted for continuous mutation burden
#   A9 Hypermutation by tumor site
# =============================================================================
while (sink.number() > 0) sink()
PROJ <- "~/GENIE_BRCA2"; setwd(PROJ)
source("01_genie_cohort_main_analysis.R", encoding = "UTF-8")
library(splines)
RV <- file.path(PROJ, "revizyon"); dir.create(RV, showWarnings = FALSE)
sink(file.path(RV, "A_log.txt"), split = TRUE)
cat("Revision A date:", format(Sys.time()), "\n")
f2 <- function(x) sprintf("%.2f", x)
pfm <- function(p) ifelse(is.na(p), NA, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

# ---- helper: number and class of variants for one gene -----------------------
gene_vars <- function(d, g) {
  col <- paste0(g, "..MUT")
  v <- d[, .(variant = split_var(get(col))), by = Patient.ID]
  if (nrow(v)) { v[, cls := var_class(variant)]; v[, codon_key := sub("^([A-Z*X]\\d+).*", "\\1", variant)] }
  v
}
nvar_gene <- function(d, g) { v <- gene_vars(d, g); n <- v[, .N, by = Patient.ID]
  out <- n$N[match(d$Patient.ID, n$Patient.ID)]; out[is.na(out)] <- 0L; out }

# ---- A0. Missing data ---------------------------------------------------------
miss <- data.table(degisken = c("Yaş", "Cinsiyet", "Mutasyon sayısı"),
  eksik_n = c(sum(is.na(pt$age)), sum(is.na(pt$sex)), sum(is.na(pt$mut_count))),
  eksik_BRCA2_altere = c(sum(is.na(pt$age) & pt$BRCA2_alt), sum(is.na(pt$sex) & pt$BRCA2_alt), sum(is.na(pt$mut_count) & pt$BRCA2_alt)))
cat("\n=== A0 MISSING DATA ===\n"); print(miss)
cat("\nFirth model N and number of events:\n"); print(unique(firth[, .(model, n, olay)]))

# ---- A1. Circularity: mutation count without BRCA2 variants -----------
pt[, nv_BRCA2 := nvar_gene(pt, "BRCA2")]
pt[, mc_x := pmax(mut_count - nv_BRCA2, 0)]
hyper_def <- function(mc, q, use_mmr = TRUE, use_mc = TRUE) {
  top <- if (use_mc) ave(mc, pt$Sequence.Assay.ID, FUN = function(z) z > quantile(z, q, na.rm = TRUE)) else rep(FALSE, length(mc))
  top <- as.logical(top); top[is.na(top)] <- FALSE
  if (use_mmr) top | pt$H2_mmr_pole else top
}
pt[, hyper_x := hyper_def(mc_x, 0.95)]
cat("\n=== A1 CIRCULARITY ===\n")
cat("Hypermutated, original definition:", sum(pt$hyper), " with BRCA2-free count:", sum(pt$hyper_x), "\n")
cat("Patients who changed class:", sum(pt$hyper != pt$hyper_x), "\n")
print(table(orijinal = pt$hyper, brca2_cikarilmis = pt$hyper_x))

# ---- Shared model function --------------------------------------------------
mdl_base <- c("age_grp", "sex", "site")
fit_or <- function(d, hv, outcome = "BRCA2_alt", center = TRUE) {
  d <- d[!is.na(age_grp) & !is.na(sex)]
  d[, H := get(hv)]; d[, center_g := droplevels(center_g)]
  cv <- c(mdl_base, if (center && uniqueN(d$center_g) > 1) "center_g")
  fml <- as.formula(paste(outcome, "~ H +", paste(cv, collapse = " + ")))
  m <- tryCatch(logistf(fml, data = d), error = function(e) NULL)
  if (is.null(m)) return(c(NA, NA, NA, NA))
  i <- which(names(coef(m)) == "HTRUE")
  c(exp(coef(m)[i]), exp(m$ci.lower[i]), exp(m$ci.upper[i]), m$prob[i])
}
summ <- function(d, hv, lab, outcome = "BRCA2_alt", center = TRUE) {
  H <- d[[hv]]; Y <- d[[outcome]]
  r <- fit_or(d, hv, outcome, center)
  data.table(analiz = lab, N = nrow(d), n_hiper = sum(H), hiper_yuzde = round(100 * mean(H), 1),
             brca2_n = sum(Y), brca2_yuzde = round(100 * mean(Y), 2),
             hiperde_brca2 = sprintf("%d/%d (%.1f%%)", sum(Y & H), sum(H), 100 * sum(Y & H) / max(1, sum(H))),
             hipersizde_brca2 = sprintf("%d/%d (%.2f%%)", sum(Y & !H), sum(!H), 100 * sum(Y & !H) / sum(!H)),
             brca2_icinde_hiper = sprintf("%.1f%%", 100 * sum(Y & H) / max(1, sum(Y))),
             OR = r[1], OR_alt = r[2], OR_ust = r[3], p = r[4])
}

# ---- A2. Sensitivity analyses of the hypermutation definition --------------------------
pt[, h_p90  := hyper_def(mc_x, 0.90)]
pt[, h_p95  := hyper_def(mc_x, 0.95)]
pt[, h_p975 := hyper_def(mc_x, 0.975)]
pt[, h_p99  := hyper_def(mc_x, 0.99)]
pt[, h_mmr  := H2_mmr_pole]
pt[, h_mc95 := hyper_def(mc_x, 0.95, use_mmr = FALSE)]
sens <- rbind(
  summ(pt, "hyper",  "Ana analiz: P95 (toplam mutasyon) veya MMR/POLE"),
  summ(pt, "h_p90",  "P90 (BRCA2 hariç) veya MMR/POLE"),
  summ(pt, "h_p95",  "P95 (BRCA2 hariç) veya MMR/POLE"),
  summ(pt, "h_p975", "P97.5 (BRCA2 hariç) veya MMR/POLE"),
  summ(pt, "h_p99",  "P99 (BRCA2 hariç) veya MMR/POLE"),
  summ(pt, "h_mmr",  "Yalnız MMR/POLE alterasyonu"),
  summ(pt, "h_mc95", "Yalnız mutasyon sayısı P95 (BRCA2 hariç)"))
# Continuous mutation burden (log mutation count without BRCA2, standardized within panel)
pt[, lmc_z := ave(log1p(mc_x), Sequence.Assay.ID, FUN = function(z) (z - mean(z, na.rm = TRUE)) / sd(z, na.rm = TRUE))]
pt[!is.finite(lmc_z), lmc_z := 0]
dc <- pt[!is.na(age_grp) & !is.na(sex)]
mc_cont <- logistf(BRCA2_alt ~ lmc_z + age_grp + sex + site + center_g, data = dc)
i <- which(names(coef(mc_cont)) == "lmc_z")
cont_row <- data.table(analiz = "Sürekli mutasyon yükü (1 SD artış, panel içi, BRCA2 hariç)", N = nrow(dc),
  OR = exp(coef(mc_cont)[i]), OR_alt = exp(mc_cont$ci.lower[i]), OR_ust = exp(mc_cont$ci.upper[i]), p = mc_cont$prob[i])
cat("\n=== A2 SENSITIVITY ANALYSES ===\n"); print(sens); print(cont_row)
# Putative LoF BRCA2 in non-hypermutated tumors under each definition
hr_sens <- rbindlist(lapply(c("hyper", "h_p90", "h_p95", "h_p975", "h_p99", "h_mmr", "h_mc95"), function(h) {
  x <- sum((pt$trunc | pt$BRCA2_homdel) & !pt[[h]]); r <- prop_ci(x, nrow(pt))
  data.table(tanim = h, putatif_LoF_hipersiz_n = x, yuzde = round(r[1], 2), ci_alt = round(r[2], 2), ci_ust = round(r[3], 2))
}))
cat("\nPutative LoF BRCA2 in non-hypermutated tumors, by definition:\n"); print(hr_sens)

# ---- A3. Without the *2602 variant ------------------------------------------------
amb <- bv[grepl("2602", variant), unique(Patient.ID)]
only_amb <- bv[, .(only = all(grepl("2602", variant))), by = Patient.ID][only == TRUE, Patient.ID]
amb_info <- pt[Patient.ID %in% amb, .(Patient.ID, Center, Sequence.Assay.ID, hyper, MMR_alt, BRCA2..MUT)]
cat("\n=== A3 *2602 VARIANT ===\n"); print(amb_info)
cat("Patients with this variant only:", length(only_amb), "\n")
pt[, BRCA2_alt_no2602 := BRCA2_alt & !(Patient.ID %in% only_amb & !BRCA2_homdel)]
s2602 <- summ(pt, "hyper", "*2602 varyantı hariç", outcome = "BRCA2_alt_no2602")
print(s2602)

# ---- A4. Robustness subgroups ---------------------------------------------
pt[, BRCA2_mutonly := BRCA2_mut]
alt_sel <- {   # Samples covering BRCA1 and BRCA2 first, then one sample per patient
  d2 <- dt[BRCA2_mut_prof == TRUE & BRCA1_mut_prof == TRUE]
  setorder(d2, Patient.ID, -BRCA2_cna_prof, -prio_primary, seq_year)
  d2[, .SD[1], by = Patient.ID]$Sample.ID }
cat("\nPatients with alternative sample selection:", length(alt_sel), " (main cohort:", nrow(pt), ")\n")
cat("Patients with a different sample than in the main cohort:", sum(!pt$Sample.ID %in% alt_sel), "\n")
rob <- list(summ(pt, "hyper", "Tüm kohort"))
for (cc in cen$Center) rob[[length(rob) + 1]] <- summ(pt[Center != cc], "hyper", paste("Çıkarılan merkez:", cc))
rob[[length(rob) + 1]] <- summ(pt[Center == "MSK"], "hyper", "Yalnız MSK (tümör-normal eşleşmeli)", center = FALSE)
rob[[length(rob) + 1]] <- summ(pt[sample_type == "Primary"], "hyper", "Yalnız primer tümör örnekleri")
rob[[length(rob) + 1]] <- summ(pt[sample_type == "Metastasis"], "hyper", "Yalnız metastaz örnekleri")
rob[[length(rob) + 1]] <- summ(pt, "hyper", "Yalnız BRCA2 mutasyonu (derin delesyon hariç)", outcome = "BRCA2_mutonly")
rob[[length(rob) + 1]] <- summ(pt[BRCA2_cna_prof == TRUE], "hyper", "Yalnız kopya sayısı profilli örnekler")
rob[[length(rob) + 1]] <- summ(pt[Sample.ID %in% alt_sel], "hyper", "Alternatif örnek seçimi")
rob <- rbindlist(rob)
cat("\n=== A4 ROBUSTNESS ===\n"); print(rob)

# ---- A5. By panel -------------------------------------------------------
pan <- pt[, .(N = .N, hiper_yuzde = round(100 * mean(hyper), 1), brca2_n = sum(BRCA2_alt),
              brca2_yuzde = round(100 * mean(BRCA2_alt), 2), putatif_LoF_hipersiz = sum(BRCA2_hr),
              CNA_profilli = round(100 * mean(BRCA2_cna_prof), 0), merkez = paste(unique(Center), collapse = "/")),
          by = Sequence.Assay.ID][order(-N)]
cat("\n=== A5 BY PANEL (>=50 patients) ===\n"); print(pan[N >= 50])

# ---- A6. Variant spectrum (per patient) -------------------------------------
pv <- bv[, .(fs = any(cls == "Frameshift"), non = any(cls == "Nonsense"), spl = any(cls == "Splice"),
             mis_only = all(cls == "Missense"), rec_fs = any(homopolymer_fs), multi = .N >= 2), by = Patient.ID]
b2 <- merge(pt[BRCA2_alt == TRUE, .(Patient.ID, hyper, trunc, BRCA2_homdel, BRCA1_alt, PALB2_alt)], pv, by = "Patient.ID", all.x = TRUE)
for (v in c("fs", "non", "spl", "mis_only", "rec_fs", "multi")) b2[is.na(get(v)), (v) := FALSE]
b2[, trunc_hd := trunc | BRCA2_homdel]; b2[, br1_pal := BRCA1_alt | PALB2_alt]
spec_rows <- list(trunc_hd = "Truncating mutation or deep deletion", fs = "Frameshift", non = "Nonsense",
  spl = "Splice site", mis_only = "Driver missense only", rec_fs = "Recurrent frameshift (codon in >=3 tumors)",
  multi = "Two or more BRCA2 mutations (phase unknown)", BRCA2_homdel = "Deep deletion", br1_pal = "BRCA1 or PALB2 co-alteration")
nh <- sum(b2$hyper); nn <- sum(!b2$hyper)
spec <- rbindlist(lapply(names(spec_rows), function(v) {
  a <- sum(b2[[v]] & b2$hyper); b <- sum(b2[[v]] & !b2$hyper)
  data.table(feature = spec_rows[[v]],
             hypermutated = sprintf("%d (%.1f)", a, 100 * a / nh), non_hypermutated = sprintf("%d (%.1f)", b, 100 * b / nn),
             P = pfm(fisher.test(matrix(c(a, nh - a, b, nn - b), 2))$p.value))
}))
cat(sprintf("\n=== A6 VARIANT SPECTRUM (hypermutated n=%d, non-hypermutated n=%d) ===\n", nh, nn)); print(spec)

# ---- A7. Comparison genes ---------------------------------------------
cmp_genes <- c("BRCA2", "BRCA1", "PALB2", "ATM", "ARID1A", "TP53", "CDH1", "PIK3CA", "KRAS")
cmp <- rbindlist(lapply(cmp_genes, function(g) {
  d <- pt[get(paste0(g, "_mut_prof")) == TRUE]
  y <- d[[paste0(g, "_mut")]]                       # mutations only (for comparability)
  a <- sum(y & d$hyper); b <- sum(y & !d$hyper); nh <- sum(d$hyper); nn <- sum(!d$hyper)
  ft <- fisher.test(matrix(c(a, nh - a, b, nn - b), 2))
  v <- gene_vars(d[y == TRUE], g)
  vh <- merge(v, d[, .(Patient.ID, hyper)], by = "Patient.ID")
  rec <- vh[, .N, by = codon_key][N >= 3, codon_key]
  data.table(gene = g, N = nrow(d),
    hyper_pct = round(100 * a / nh, 2), nonhyper_pct = round(100 * b / nn, 2), ratio = round((a / nh) / (b / nn), 1),
    OR = ft$estimate, OR_lo = ft$conf.int[1], OR_hi = ft$conf.int[2],
    carriers_hyper_pct = round(100 * a / max(1, a + b), 1),
    frameshift_pct_of_variants = round(100 * mean(vh$cls == "Frameshift"), 1),
    recurrent_fs_pct_of_variants = round(100 * mean(vh$cls == "Frameshift" & vh$codon_key %in% rec), 1),
    hyper_fs_pct = round(100 * mean(vh[hyper == TRUE]$cls == "Frameshift"), 1),
    nonhyper_fs_pct = round(100 * mean(vh[hyper == FALSE]$cls == "Frameshift"), 1))
}))
cat("\n=== A7 COMPARISON GENES (mutations only) ===\n"); print(cmp)

# ---- A8. Co-alterations: Firth models with continuous mutation burden --------------------
co8 <- rbindlist(lapply(co_genes, function(g) {
  d <- if (g == "ERBB2") pt[ERBB2_cna_prof == TRUE] else pt[get(paste0(g, "_mut_prof")) == TRUE]
  d <- d[!is.na(mut_count)]
  d[, mc_g := pmax(mut_count - nv_BRCA2 - (if (g == "ERBB2") 0 else nvar_gene(d, g)), 0)]
  d[, lz := ave(log1p(mc_g), Sequence.Assay.ID, FUN = function(z) (z - mean(z)) / ifelse(sd(z) > 0, sd(z), 1))]
  d[, G := get(paste0(g, "_alt"))]
  m <- tryCatch(logistf(BRCA2_alt ~ G + ns(lz, df = 3) + site + center_g, data = d), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  i <- which(names(coef(m)) == "GTRUE")
  data.table(gene = g, N = nrow(d), BRCA2_events = sum(d$BRCA2_alt), co_alt_in_BRCA2 = sum(d$G & d$BRCA2_alt),
             OR = exp(coef(m)[i]), OR_lo = exp(m$ci.lower[i]), OR_hi = exp(m$ci.upper[i]), p = m$prob[i])
}))
co8[, q := p.adjust(p, "BH")]
cat("\n=== A8 CO-ALTERATIONS: continuous mutation burden (spline, without BRCA2 and the gene studied) + site + center ===\n"); print(co8)

# ---- A9. Hypermutation by tumor site ---------------------------------------
site_h <- pt[, .(N = .N, hiper_n = sum(hyper), hiper_yuzde = round(100 * mean(hyper), 1),
                 brca2_hiperde = sum(BRCA2_alt & hyper), brca2_hipersiz = sum(BRCA2_alt & !hyper)), by = site]
cat("\n=== A9 HYPERMUTATION BY TUMOR SITE ===\n"); print(site_h)

# ---- Save ------------------------------------------------------------------
fwrite(miss, file.path(RV, "A0_eksik_veri.csv")); fwrite(sens, file.path(RV, "A2_duyarlilik.csv"))
fwrite(cont_row, file.path(RV, "A2_surekli_yuk.csv")); fwrite(hr_sens, file.path(RV, "A2_putatifLoF_tanim.csv"))
fwrite(amb_info, file.path(RV, "A3_2602_hastalar.csv")); fwrite(s2602, file.path(RV, "A3_2602_haric.csv"))
fwrite(rob, file.path(RV, "A4_saglamlik.csv")); fwrite(pan, file.path(RV, "A5_panel.csv"))
fwrite(spec, file.path(RV, "A6_varyant_spektrumu.csv")); fwrite(cmp, file.path(RV, "A7_karsilastirmali_genler.csv"))
fwrite(co8, file.path(RV, "A8_koalt_surekli.csv")); fwrite(site_h, file.path(RV, "A9_yerlesim_hiper.csv"))

# ---- Figure 5: sensitivity and robustness forest plot -----------------------------------
fp <- rbind(sens[, .(grp = "Hypermutation definition", analiz, OR, OR_alt, OR_ust)],
            rob[, .(grp = "Subgroup / robustness", analiz, OR, OR_alt, OR_ust)], s2602[, .(grp = "Subgroup / robustness", analiz, OR, OR_alt, OR_ust)])
lab_en <- c("Ana analiz: P95 (toplam mutasyon) veya MMR/POLE" = "Main: P95 (total count) or MMR/POLE",
  "P90 (BRCA2 hariç) veya MMR/POLE" = "P90 (BRCA2 excluded) or MMR/POLE", "P95 (BRCA2 hariç) veya MMR/POLE" = "P95 (BRCA2 excluded) or MMR/POLE",
  "P97.5 (BRCA2 hariç) veya MMR/POLE" = "P97.5 (BRCA2 excluded) or MMR/POLE", "P99 (BRCA2 hariç) veya MMR/POLE" = "P99 (BRCA2 excluded) or MMR/POLE",
  "Yalnız MMR/POLE alterasyonu" = "MMR/POLE alteration only", "Yalnız mutasyon sayısı P95 (BRCA2 hariç)" = "Mutation count P95 only (BRCA2 excluded)",
  "Tüm kohort" = "All patients", "Yalnız MSK (tümör-normal eşleşmeli)" = "MSK only (tumor-normal)",
  "Yalnız primer tümör örnekleri" = "Primary tumor samples", "Yalnız metastaz örnekleri" = "Metastasis samples",
  "Yalnız BRCA2 mutasyonu (derin delesyon hariç)" = "BRCA2 mutations only", "Yalnız kopya sayısı profilli örnekler" = "Copy number profiled samples",
  "Alternatif örnek seçimi" = "Alternative sample selection", "*2602 varyantı hariç" = "Ambiguous *2602 variant excluded")
fp[, lab := ifelse(analiz %in% names(lab_en), lab_en[analiz], sub("Çıkarılan merkez: ", "Without ", analiz))]
fp[, lab := factor(lab, levels = rev(unique(lab)))]
g <- ggplot(fp, aes(OR, lab)) + geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = OR_alt, xmax = OR_ust), height = 0.2) + geom_point(size = 2.2, colour = "#C0392B") +
  scale_x_log10() + facet_grid(grp ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Adjusted odds ratio for BRCA2 alteration in hypermutated tumors (log scale)", y = NULL) +
  theme_bw(base_size = 11) + theme(strip.text.y = element_text(angle = 0))
ggsave(file.path(RV, "Figure5_sensitivity.png"), g, width = 9, height = 8, dpi = 300)
ggsave(file.path(RV, "Figure5_sensitivity.tiff"), g, width = 9, height = 8, dpi = 300, compression = "lzw")

# Figure 6: comparison genes
cl <- melt(cmp[, .(gene, Hypermutated = hyper_pct, `Non-hypermutated` = nonhyper_pct)], id.vars = "gene")
cl[, gene := factor(gene, levels = cmp[order(-ratio), gene])]
g7 <- ggplot(cl, aes(gene, value, fill = variable)) + geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_text(data = cmp, aes(gene, pmax(hyper_pct, nonhyper_pct) + 3, label = paste0(ratio, "x")), inherit.aes = FALSE, size = 3.2) +
  scale_fill_manual(values = c("#C0392B", "#2C6E9B")) + labs(x = NULL, y = "Tumors with a driver mutation (%)", fill = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(RV, "Figure6_comparative_genes.png"), g7, width = 9, height = 5, dpi = 300)
ggsave(file.path(RV, "Figure6_comparative_genes.tiff"), g7, width = 9, height = 5, dpi = 300, compression = "lzw")
cat("\nDone. Outputs:", normalizePath(RV), "\n")
sink()
