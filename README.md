# BRCA2 alterations in gastroesophageal adenocarcinoma (AACR Project GENIE)

R code for the study:

> Kofalı Ayakdaş G, Kutlu Ö. Genomic Context of BRCA2 Alterations in Gastroesophageal Adenocarcinoma: A Large-Scale Analysis of AACR Project GENIE. Submitted to the International Journal of Molecular Sciences.

The code reproduces every number, table and figure in the article and its Supplementary Materials. It uses only public, de-identified data from AACR Project GENIE (release v20.0-public). No patient-level data are stored in this repository.

Authors:

- Güzide Kofalı Ayakdaş (corresponding author, guzidekofali@gmail.com), ORCID 0000-0001-6112-0127
- Özlem Kutlu, ORCID 0000-0001-5050-3333

University of Health Sciences Türkiye, İzmir City Hospital, Department of Medical Oncology, İzmir, Türkiye.

## Scripts used in the article

| Script | What it does | Main outputs |
|---|---|---|
| `01_genie_cohort_main_analysis.R` | Builds the GENIE cohort (one tumor sample per patient, panels covering BRCA1 and BRCA2), defines hypermutation and putative loss-of-function (LoF) BRCA2, and runs the main analyses: frequencies, co-alterations, Firth models and the center meta-analysis. | `cikti_v2/` |
| `02_manuscript_tables_figures.R` | Runs script 01, then writes the first versions of the main tables and figures. | `makale/` |
| `03_revision_A_genie_sensitivity.R` | Runs script 01, then the sensitivity and robustness analyses: missing data, circularity check, other hypermutation cutoffs, continuous mutation burden, leave-one-center-out, panels, variant spectrum, comparison genes, co-alteration models adjusted for continuous burden and hypermutation by site. | `revizyon/` |
| `07_revision_D_second_review.R` | Runs scripts 02 and 03, then: P80 and P85 cutoffs, putative LoF by hypermutation, co-alteration n/N and q values, a cohort including poorly differentiated carcinomas and the sequence context of BRCA2 frameshifts (NM_000059.4 from NCBI). | `revizyon/`, `revizyon/D_figures/` |
| `08_revision_D_figure_fixes.R` | Small layout fixes. Run in the same R session right after script 07. | `revizyon/D_figures/` |
| `10_ijms_figures.R` | Runs scripts 07 and 08, then draws all final figures with the numbering of the article. | `revizyon/F_ijms/` |

## Exploratory scripts (not used in the article)

Scripts 04, 05, 06 and 09, and the TCGA parts of script 07, test the findings in TCGA PanCancer Atlas data. This work was done during an earlier draft. The submitted article is based on GENIE data only, so none of these TCGA results appear in it. The scripts are kept for transparency.

| Script | What it does |
|---|---|
| `04_revision_B1_tcga_download.R` | Downloads TCGA stomach and esophageal data from the cBioPortal API and DNA damage response resources from the GDC. |
| `05_revision_B2_tcga_validation.R` | TCGA checks by MSI status, HRD score, LOH and signature 3. |
| `06_revision_C_figures_biallelic.R` | Reads the OncoKB data version and downloads the TCGA ABSOLUTE copy number file. |
| `09_tcga_pathogenic_check_and_figure1.R` | TCGA check with a pathogenic BRCA2 definition. |

Script 07 sources script 05, so running 07 also runs the TCGA steps. Their outputs are simply not used.

## Where each item in the article comes from

| Article item | Script |
|---|---|
| Table 1 (patient characteristics) | 02 |
| Table 2 (frequency and type of BRCA2 alterations) | 02, 03 and 07 |
| Table 3 (variant spectrum) | 03 and 07 |
| Table S1 (hypermutation by site) | 03 |
| Table S2 (Firth regression) | 02 |
| Tables S3 (sensitivity), S4 (robustness), S5 (panels) | 03 and 07 |
| Table S6 (centers with at least 100 patients) | 01 and 02 |
| Table S7 (frameshift sequence context) | 07 |
| Table S8 (comparison genes) | 03 |
| Table S9 (co-alterations) | 01, 03 and 07 |
| Table S10 (continuous-burden models) | 03 |
| Figures 1 to 5 and Figure S1 | 10 |

The tables were formatted by hand from the CSV and log files that these scripts write. All numbers come from those outputs.

## Folder layout

All scripts expect a project folder at `~/GENIE_BRCA2`. To use another location, change `PROJ` at the top of each script.

```
~/GENIE_BRCA2/
  01_genie_cohort_main_analysis.R ... 10_ijms_figures.R
  veri/              GENIE input files (you download these, see below)
  cikti_v2/          outputs of script 01
  makale/            outputs of script 02
  revizyon/          outputs of scripts 03 and 07
  revizyon/F_ijms/   final figures (script 10)
```

The `figures/` folder in this repository holds copies of the final figures.

## Input data

Access to GENIE requires agreeing to the AACR Project GENIE Terms of Access. The terms do not allow the data to be shared again, so the input files are not included here. They were downloaded in November 2026 from the GENIE cBioPortal (https://genie.cbioportal.org), with all esophagogastric cancer samples of the v20.0-public release selected:

1. `genie_public_clinical_data (2).tsv`: clinical data for the selected samples (7,837 samples). Rename the file or change `CLIN_FILE` in script 01 if your copy has another name.
2. `alterations_across_samples.tsv`: alteration table for BRCA1, BRCA2, PALB2, ATM, ERBB2, TP53, KRAS, CDH1, ARID1A and PIK3CA, with OncoKB and hotspot driver annotation switched on, so that driver alterations are marked "(driver)".
3. `alterations_mmr.tsv`: the same table for MLH1, MSH2, MSH6, PMS2 and POLE.

Place the three files in `~/GENIE_BRCA2/veri/`. Driver labels in cBioPortal come from OncoKB (data version v7.7, 29 October 2026, at the time of download). Script 07 also downloads the BRCA2 coding sequence (NM_000059.4) from NCBI.

## How to run

R version 4.6.2 was used. The quickest way to rebuild everything in the article:

```r
setwd("~/GENIE_BRCA2")
source("10_ijms_figures.R", encoding = "UTF-8")   # runs 07 and 08, which run 01, 02, 03 and 05
```

Because script 07 also runs the exploratory TCGA steps, it needs the TCGA files downloaded by scripts 04 and 06. If you have not run them, run them once first:

```r
source("04_revision_B1_tcga_download.R", encoding = "UTF-8")
source("06_revision_C_figures_biallelic.R", encoding = "UTF-8")
```

Missing packages are installed automatically from CRAN. The packages used are data.table, ggplot2, logistf, metafor, splines, flextable, officer, patchwork, httr and jsonlite. The exact versions are listed in `sessionInfo.txt`.

Each script writes a log file with all printed results (for example `cikti_v2/00_konsol_log_v2.txt`, `makale/00_makale_log.txt`, `revizyon/A_log.txt`, `revizyon/D_log.txt`).

## Language note

Comments and console headings are in English. Some internal labels, column names and output file names are in Turkish, because the analysis was first written in Turkish. They are kept as they are so that the code gives exactly the same output as reported in the article. Common terms:

| Turkish | English |
|---|---|
| veri / cikti / makale / revizyon | data / output / manuscript / revision |
| hipermutasyonlu / hipermutasyonsuz | hypermutated / non-hypermutated |
| tabaka, olay, degisken, eksik | stratum, events, variable, missing |
| duyarlilik / saglamlik | sensitivity / robustness |
| varyant_spektrumu / karsilastirmali_genler | variant spectrum / comparison genes |
| koalt_surekli | co-alterations, continuous burden model |
| yerlesim_hiper | hypermutation by site |

In object names, `HR` or "HR-anlamlı" refers to what the article calls a putative LoF BRCA2 alteration (a truncating mutation or deep deletion in a non-hypermutated tumor).

## License

MIT License (see `LICENSE`).
