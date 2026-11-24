# =============================================================================
# REVISION B1 - Download TCGA STAD + ESCA (PanCancer Atlas) data
# Sources: cBioPortal public API (clinical, mutation and copy number data) and
# the GDC PanCan-DDR-2018 data resources (Knijnenburg et al. 2018).
# Raw files are saved in ~/GENIE_BRCA2/revizyon/tcga/.
# =============================================================================
while (sink.number() > 0) sink()
for (p in c("httr", "jsonlite", "data.table")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(httr); library(jsonlite); library(data.table)
TD <- path.expand("~/GENIE_BRCA2/revizyon/tcga"); dir.create(TD, showWarnings = FALSE, recursive = TRUE)
sink(file.path(TD, "B1_log.txt"), split = TRUE)
cat("Download date:", format(Sys.time()), "\n")
API <- "https://www.cbioportal.org/api"
studies <- c("stad_tcga_pan_can_atlas_2018", "esca_tcga_pan_can_atlas_2018")
genes <- c(BRCA2 = 675, BRCA1 = 672, PALB2 = 79728, ATM = 472, MLH1 = 4292, MSH2 = 4436,
           MSH6 = 2956, PMS2 = 5395, POLE = 5426, ARID1A = 8289, TP53 = 7157)
getj <- function(url) { r <- GET(url, timeout(300)); stop_for_status(r); as.data.table(fromJSON(content(r, "text", encoding = "UTF-8"))) }
postj <- function(url, body) { r <- POST(url, body = body, encode = "json", timeout(300)); stop_for_status(r)
  as.data.table(fromJSON(content(r, "text", encoding = "UTF-8"))) }
for (s in studies) {
  cat("\n=====", s, "=====\n")
  cs <- getj(sprintf("%s/studies/%s/clinical-data?clinicalDataType=SAMPLE&projection=SUMMARY&pageSize=10000000", API, s))
  cp <- getj(sprintf("%s/studies/%s/clinical-data?clinicalDataType=PATIENT&projection=SUMMARY&pageSize=10000000", API, s))
  fwrite(cs, file.path(TD, paste0(s, "_clin_sample_long.csv"))); fwrite(cp, file.path(TD, paste0(s, "_clin_patient_long.csv")))
  cat("Sample-level clinical attributes:\n"); print(sort(unique(cs$clinicalAttributeId)))
  cat("Patient-level clinical attributes:\n"); print(sort(unique(cp$clinicalAttributeId)))
  mp <- getj(sprintf("%s/studies/%s/molecular-profiles", API, s)); cat("Molecular profiles:\n"); print(mp$molecularProfileId)
  sl <- getj(sprintf("%s/studies/%s/sample-lists", API, s)); cat("Sample lists:\n"); print(sl$sampleListId)
  mu <- postj(sprintf("%s/molecular-profiles/%s_mutations/mutations/fetch?projection=SUMMARY", API, s),
              list(sampleListId = paste0(s, "_sequenced"), entrezGeneIds = unname(genes)))
  keep <- intersect(c("sampleId", "patientId", "entrezGeneId", "proteinChange", "mutationType", "variantType",
                      "tumorAltCount", "tumorRefCount", "chr", "startPosition", "referenceAllele", "variantAllele"), names(mu))
  mu <- mu[, ..keep]; mu[, gene := names(genes)[match(entrezGeneId, genes)]]
  fwrite(mu, file.path(TD, paste0(s, "_mutations.csv"))); cat("Mutation rows:", nrow(mu), "\n"); print(mu[, .N, by = gene])
  cn <- tryCatch(postj(sprintf("%s/molecular-profiles/%s_gistic/discrete-copy-number/fetch?discreteCopyNumberEventType=ALL&projection=SUMMARY", API, s),
              list(sampleListId = paste0(s, "_cna"), entrezGeneIds = unname(genes))), error = function(e) { cat("CNA error:", conditionMessage(e), "\n"); NULL })
  if (!is.null(cn)) { cn <- cn[, intersect(c("sampleId", "patientId", "entrezGeneId", "alteration"), names(cn)), with = FALSE]
    cn[, gene := names(genes)[match(entrezGeneId, genes)]]; fwrite(cn, file.path(TD, paste0(s, "_cna.csv")))
    cat("CNA rows:", nrow(cn), "\n"); print(cn[gene == "BRCA2", .N, by = alteration]) }
}
# GDC PanCan-DDR-2018 data resources
z <- file.path(TD, "DDR_data_resources.zip")
ok <- tryCatch({ download.file("https://api.gdc.cancer.gov/data/d9fdd899-bffe-4325-b0ce-58cbf9924c0a", z, mode = "wb", quiet = TRUE); TRUE },
               error = function(e) { cat("GDC download error:", conditionMessage(e), "\n"); FALSE })
if (ok) {
  cat("\nZIP contents:\n"); print(unzip(z, list = TRUE))
  unzip(z, exdir = file.path(TD, "ddr"))
  for (f in list.files(file.path(TD, "ddr"), recursive = TRUE, full.names = TRUE)) {
    cat("\n---", basename(f), "---\n")
    x <- tryCatch(fread(f, nrows = 3), error = function(e) NULL)
    if (!is.null(x)) { cat("Columns:", paste(head(names(x), 60), collapse = " | "), "\n") }
  }
}
cat("\nDone.\n"); sink()
