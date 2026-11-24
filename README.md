# BRCA2 and hypermutation in gastroesophageal adenocarcinoma

R code for the study:

> Kofalı Ayakdaş G. Genomic context of BRCA2 alterations in gastroesophageal adenocarcinoma: an integrated analysis of AACR Project GENIE and TCGA.

The code reproduces every number, table and figure in the manuscript and in Additional file 1. It uses only public data from AACR Project GENIE (release v20.0-public) and the TCGA PanCancer Atlas. No patient-level data are stored in this repository.

Author: Güzide Kofalı Ayakdaş, University of Health Sciences Türkiye, İzmir City Hospital, Department of Medical Oncology, İzmir, Türkiye. ORCID: 0000-0001-6112-0127

## Scripts

| Script | What it does | Main outputs |
|---|---|---|
| `01_genie_cohort_main_analysis.R` | Builds the GENIE cohort (one tumor sample per patient, panels covering BRCA1 and BRCA2), defines hypermutation and putative loss-of-function (LoF) BRCA2, and runs the main analyses: frequencies, co-alterations, Firth models and the center meta-analysis. | `cikti_v2/` |
| `02_manuscript_tables_figures.R` | Runs script 01, then writes Tables 1-3, Tables S1-S2 and Figures 1-4 (first versions). | `makale/` |
| `03_revision_A_genie_sensitivity.R` | Runs script 01, then the sensitivity and robustness analyses: missing data, circularity check, other hypermutation cutoffs, continuous mutation burden, leave-one-center-out, panels, variant spectrum, comparison genes, co-alteration models adjusted for continuous burden and hypermutation by site. Writes Figures 5 and 6. | `revizyon/` |
| `04_revision_B1_tcga_download.R` | Downloads TCGA stomach and esophageal PanCancer Atlas data (clinical, mutations, copy number) from the cBioPortal API and the DNA damage response data resources (Knijnenburg et al. 2018) from the GDC. | `revizyon/tcga/` |
| `05_revision_B2_tcga_validation.R` | TCGA validation: BRCA2 by true MSI status, HRD score, LOH fraction and signature 3 in four groups, and the accuracy of the GENIE hypermutation proxy. Writes Figure 7. | `revizyon/` |
| `06_revision_C_figures_biallelic.R` | Runs scripts 02 and 05, reads the OncoKB data version and downloads the TCGA ABSOLUTE allele-specific copy number file (needed by script 07). | `revizyon/` |
| `07_revision_D_second_review.R` | Runs scripts 02, 03 and 05, then: P80 and P85 cutoffs, putative LoF by hypermutation, Table 3 n/N, Table 4 q values, a cohort including poorly differentiated carcinomas, the sequence context of BRCA2 frameshifts (NM_000059.4 from NCBI), TCGA cohort flow, missing HRD scores, POLE subtype check, the HRD interaction models with bootstrap, proxy accuracy by BRCA2 status, allelic status tiers with variant-matched ABSOLUTE segments, the patient-level TCGA table, harmonized BRCA2 definitions and absolute signature 3. Writes the final Figures 1, 3, 4, 5, 6, 7 and S1. | `revizyon/`, `revizyon/D_figures/` |
| `08_revision_D_figure_fixes.R` | Small layout fixes for Figures 4, 6 and 7. Run in the same R session right after script 07. | `revizyon/D_figures/` |

Where each manuscript item comes from:

| Manuscript item | Script |
|---|---|
| Table 1, Table 2 (part A), Table S1, Table S2 | 02 |
| Table 2 (parts B and C), Table 3 (n/N), Table 4 (q values) | 03 and 07 |
| Table 5 | 05 and 07 |
| Table S3 (sensitivity), Table S4 (robustness), Table S5 (panels) | 03 and 07 |
| Table S6 (comparison genes), Table S7 (continuous-burden models), Table S8 (site) | 03 |
| Table S9 (proxy accuracy) | 05 and 07 |
| Table S10 (frameshift sequence context), Table S11 (HRD models), Table S12 (harmonized definitions), Table S13 (TCGA patient-level data) | 07 |
| Figure 2 | 02 |
| Figures 1, 3, 4, 5, 6, 7 and S1 | 07 and 08 |

The tables in the manuscript were formatted by hand from the CSV files and log files that these scripts write. All numbers come from those outputs.

## Folder layout

All scripts expect a project folder at `~/GENIE_BRCA2`. To use another location, change `PROJ` (and `TD`/`RV` in scripts 04 and 05) at the top of each script.

```
~/GENIE_BRCA2/
  01_genie_cohort_main_analysis.R ... 06_revision_C_figures_biallelic.R
  veri/                      GENIE input files (you download these, see below)
  cikti_v2/                  outputs of script 01
  makale/                    outputs of script 02
  revizyon/                  outputs of scripts 03, 05 and 06
  revizyon/tcga/             TCGA files downloaded by scripts 04 and 06
```

## Input data

### AACR Project GENIE v20.0-public

Access to GENIE requires agreeing to the AACR Project GENIE Terms of Access. Because the terms do not allow the data to be shared again, the input files are not included here. They were downloaded in November 2026 from the GENIE cBioPortal (https://genie.cbioportal.org), with all esophagogastric cancer samples of the v20.0-public release selected:

1. `genie_public_clinical_data (2).tsv`: the clinical data table for the selected samples (7,837 samples). Rename the file or change `CLIN_FILE` in script 01 if your copy has another name.
2. `alterations_across_samples.tsv`: the alteration table for BRCA1, BRCA2, PALB2, ATM, ERBB2, TP53, KRAS, CDH1, ARID1A and PIK3CA, with OncoKB and hotspot driver annotation switched on, so that driver alterations are marked "(driver)".
3. `alterations_mmr.tsv`: the same table for MLH1, MSH2, MSH6, PMS2 and POLE.

Place the three files in `~/GENIE_BRCA2/veri/`. Driver labels in cBioPortal come from OncoKB; the OncoKB data version at the time of download was v7.7 (29 October 2026), as recorded by script 06.

### TCGA PanCancer Atlas

No login is needed. Script 04 downloads the TCGA data from the public cBioPortal API (https://www.cbioportal.org/api) and the DNA damage response data resources from the GDC. Script 06 downloads `TCGA_mastercalls.abs_segtabs.fixed.txt` (ABSOLUTE allele-specific copy number) from the GDC PanCancer Atlas page. These files are large and the download can take several minutes. Script 07 also downloads the BRCA2 coding sequence (NM_000059.4) from NCBI.

## How to run

Run the scripts in order in R (version 4.6.2 was used):

```r
setwd("~/GENIE_BRCA2")
source("02_manuscript_tables_figures.R", encoding = "UTF-8")    # also runs 01
source("03_revision_A_genie_sensitivity.R", encoding = "UTF-8") # also runs 01
source("04_revision_B1_tcga_download.R", encoding = "UTF-8")
source("05_revision_B2_tcga_validation.R", encoding = "UTF-8")
source("06_revision_C_figures_biallelic.R", encoding = "UTF-8") # also runs 02 and 05
source("07_revision_D_second_review.R", encoding = "UTF-8")     # also runs 02, 03 and 05
source("08_revision_D_figure_fixes.R", encoding = "UTF-8")      # same session as 07
```

Missing packages are installed automatically from CRAN. The packages used are data.table, ggplot2, logistf, metafor, splines, flextable, officer, patchwork, httr and jsonlite. The exact versions are listed in `sessionInfo.txt`.

Each script writes a log file with all printed results (`cikti_v2/00_konsol_log_v2.txt`, `makale/00_makale_log.txt`, `revizyon/A_log.txt`, `revizyon/tcga/B1_log.txt`, `revizyon/B2_log.txt`, `revizyon/C_log.txt`, `revizyon/D_log.txt`).

## Language note

Comments and console headings are in English. Some internal labels, column names and output file names are in Turkish, because the analysis was first written in Turkish. They are kept unchanged so that the code gives exactly the same output as the analysis reported in the paper. Common terms:

| Turkish | English |
|---|---|
| veri / cikti / makale / revizyon | data / output / manuscript / revision |
| hipermutasyonlu / hipermutasyonsuz | hypermutated / non-hypermutated |
| tabaka, olay, degisken, eksik | stratum, events, variable, missing |
| yalnız ..., ... hariç | ... only, ... excluded |
| Çıkarılan merkez: X | Center X left out |
| duyarlilik / saglamlik | sensitivity / robustness |
| varyant_spektrumu / karsilastirmali_genler | variant spectrum / comparison genes |
| koalt_surekli | co-alterations, continuous burden model |
| yerlesim_hiper | hypermutation by site |
| dort_grup / vekil_dogruluk | four groups / proxy accuracy |
| cikarilanlar / biallelik / kurasyonlu | exclusions / biallelic / curated |

In object names, `HR` or "HR-anlamlı" refers to what the paper calls a putative LoF BRCA2 alteration (a truncating mutation or deep deletion in a non-hypermutated tumor).

## License

MIT License (see `LICENSE`).
