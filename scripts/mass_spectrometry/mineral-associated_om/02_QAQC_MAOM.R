# -----------------------------------------------------------------------------#
# 02_QAQC_MAOM.R
#
# Purpose
#   QA/QC for MAOM (FTICR-MS): screens replicate outliers by per-sample
#   peak count using a median +/- k*MAD rule, generates pre- and
#   post-filter diagnostic plots, and exports the QAQC'd datasets
#   (presence/absence, intensities, sample-average and formula-level
#   molecular properties).
#
# Inputs
#   2023_MONet_FTICR-MS_mineral-associated_om_Intensities.csv
#   2023_MONet_FTICR-MS_mineral-associated_om_Presence-Absence.csv
#   2023_MONet_FTICR-MS_mineral-associated_om_Avg_Molecular-Properties.csv
#   2023_MONet_FTICR-MS_mineral-associated_om_Molecular-Properties.csv
#   (all written by 01_processing_MAOM.R)
#
# Outputs
#   Data (written to OD):
#     - 2023_FTICR-MS_mineral-associated_om_PA_MADfiltered.csv
#     - 2023_FTICR-MS_mineral-associated_om_Intensity_MADfiltered.csv
#     - 2023_FTICR-MS_mineral-associated_om_Sample_Properties_MADfiltered.csv
#     - 2023_FTICR-MS_mineral-associated_om_ChemProperties_MADfiltered.csv
#   Figures (written to PD):
#     - preQAQC_mineral-associated_om_boxplot.pdf / preQAQC_mineral-associated_om_histogram.pdf
#     - mineral-associated_om_MAD_boxplot.pdf
#
# Notes
#   - Outlier flagging uses k = 3 (median +/- 3*MAD) within
#     Site x Sampling.Event x Core.Section strata; see
#     mad_filter_samples() in fticr_qaqc_outlier_functions.R.
#   - After sample filtering, formulas with zero detections across the
#     retained samples are dropped from the presence/absence and
#     intensity tables (and the formula-property table is subset to
#     match).
#
# Created:  2026-01-25
# By: X. Takver
# -----------------------------------------------------------------------------#

## ----------------------------------------- ##
#             Import Packages ----
## ----------------------------------------- ##

# Core wrangling + plotting
librarian::shelf(tidyverse, here)

# Theme + PDF device export
librarian::shelf(tinythemes, Cairo)

## ----------------------------------------- ##
#             Source Functions ----
## ----------------------------------------- ##

# Source helper functions
source(here('R', 'functions', 'paths_and_directories.R'))
source(here('R', 'functions', 'fticr_qaqc_outlier_functions.R'))
source(here('R', 'functions', 'fonts.R'))          # e.g., register_fonts()
source(here('R', 'functions', 'plot_theme.R'))      # e.g., set_plot_theme()

## ----------------------------------------- ##
#           Directory Definitions ----
## ----------------------------------------- ##

# Base data directory and output locations
DD <- here('data', 'mass_spectrometry', 'processed', 'mineral-associated_om')

OD <- file.path(DD, 'QAQC')
PD <- here('output', 'figures', 'mass_spectrometry', 'mineral-associated_om', 'qaqc')

# Ensure output directories exist (idempotent)
ensure_dir(OD)
ensure_dir(PD)

## ----------------------------------------- ##
#             Set Plot Theme ----
## ----------------------------------------- ##

# Set global ggplot2 theme for consistent manuscript figure styling
set_plot_theme()

## ----------------------------------------- ##
#                    MAOM ----
## ----------------------------------------- ##

## ----------------------------------------- ##
#                Import Data ----
## ----------------------------------------- ##

# Read MAOM inputs (presence/absence drives peak-count QAQC; other tables are filtered downstream)
df_intensities <- read_csv(file.path(DD, '2023_MONet_FTICR-MS_mineral-associated_om_Intensities.csv'))
df_PA <- read_csv(file.path(DD, '2023_MONet_FTICR-MS_mineral-associated_om_Presence-Absence.csv'))
df_avg_mol_prop <- read_csv(file.path(DD, '2023_MONet_FTICR-MS_mineral-associated_om_Avg_Molecular-Properties.csv'))
df_mol_prop <- read_csv(file.path(DD, '2023_MONet_FTICR-MS_mineral-associated_om_Molecular-Properties.csv'))

## ----------------------------------------- ##
#             preQAQC Analysis ----
## ----------------------------------------- ##

# Summarise per-sample peak counts and parse sample metadata
df_peaks <- make_peak_counts(df_PA)

## ----------------------------------------- ##
#           PreQAQC Diagnostics ----
## ----------------------------------------- ##

# Boxplots of peak counts (used to visualise replicate-level outliers)
boxplt_all <- ggplot(df_peaks) +
  geom_boxplot(aes(x = Site, y = peak_no, fill = Core.Section)) +
  scale_fill_manual(values = c('#638B66FF', '#D9AF6BFF')) +
  facet_wrap(~Sampling.Event) +
  labs(y = 'Peak #', fill = 'Core Section', title = 'MAOM Outlier Boxplot w/ all data')

boxplt_all

CairoPDF(file.path(PD, 'preQAQC_mineral-associated_om_boxplot.pdf'), width = 6.5, height = 4)
print(boxplt_all)
dev.off()

# Histogram distributions of peak counts by site/event/core section
df_hist <- df_peaks |>
  mutate(Sample.ID = as_factor(paste(Site, Core.Section, sep = ' ')))

plt_hist <- ggplot(df_hist) +
  geom_histogram(aes(x = peak_no, fill = Site)) +
  scale_fill_manual(values = c('#79A7ACFF', '#AF6458FF')) +
  facet_wrap(Sampling.Event ~ Sample.ID) +
  labs(title = 'MAOM Outlier Distribution w/ all data')

plt_hist

CairoPDF(file.path(PD, 'preQAQC_mineral-associated_om_histogram.pdf'), width = 15, height = 8.5)
print(plt_hist)
dev.off()

rm(boxplt_all, plt_hist, df_hist)

## ----------------------------------------- ##
#               MAD Filtering ----
## ----------------------------------------- ##

# Identify outlier samples using median ± k*MAD within Site x Sampling.Event x Core.Section strata
mad_res <- mad_filter_samples(df_peaks, k = 3)

# Retained samples after MAD filtering
df_MAD <- mad_res$df_kept
kept_ids <- df_MAD$SampleID

## ----------------------------------------- ##
#           Post-filter Diagnostic ----
## ----------------------------------------- ##

# Visualise peak counts after filtering
boxplt_MAD <- ggplot(df_MAD) +
  geom_boxplot(aes(x = Site, y = peak_no, fill = Core.Section)) +
  scale_fill_manual(values = c('#638B66FF', '#D9AF6BFF')) +
  facet_wrap(~Sampling.Event) +
  labs(y = 'Peak #', fill = 'Core Section', title = 'Outlier Boxplot w/ MAD Filtering')

boxplt_MAD

CairoPDF(file.path(PD, 'mineral-associated_om_MAD_boxplot.pdf'), width = 15, height = 8.5)
print(boxplt_MAD)
dev.off()

rm(boxplt_MAD)

## ----------------------------------------- ##
#             Filter Data ----
## ----------------------------------------- ##

# Presence/absence: keep retained samples and remove formulas with zero detections across retained set
df_PA_MAD <- df_PA |>
  select(molecular_formula, any_of(kept_ids)) |>
  mutate(all.detections = rowSums(across(-molecular_formula))) |>
  filter(all.detections > 0) |>
  select(-all.detections)

# Intensities: keep retained samples and remove formulas with zero detections across retained set
df_intensity_MAD <- df_intensities |>
  select(molecular_formula, any_of(kept_ids)) |>
  mutate(all.detections = rowSums(across(-molecular_formula))) |>
  filter(all.detections > 0) |>
  select(-all.detections)

# Sample-average molecular properties: subset to retained samples
df_AvgMol_MAD <- df_avg_mol_prop |>
  filter(Sample %in% kept_ids)

# Formula-level properties: keep only formulas present after PA filtering
df_MolProp_MAD <- df_mol_prop |>
  filter(molecular_formula %in% df_PA_MAD$molecular_formula)

## ----------------------------------------- ##
#               Save Output ----
## ----------------------------------------- ##

write_csv(df_PA_MAD, file.path(OD, '2023_FTICR-MS_mineral-associated_om_PA_MADfiltered.csv'))
write_csv(df_intensity_MAD, file.path(OD, '2023_FTICR-MS_mineral-associated_om_Intensity_MADfiltered.csv'))
write_csv(df_AvgMol_MAD, file.path(OD, '2023_FTICR-MS_mineral-associated_om_Sample_Properties_MADfiltered.csv'))
write_csv(df_MolProp_MAD, file.path(OD, '2023_FTICR-MS_mineral-associated_om_ChemProperties_MADfiltered.csv'))
