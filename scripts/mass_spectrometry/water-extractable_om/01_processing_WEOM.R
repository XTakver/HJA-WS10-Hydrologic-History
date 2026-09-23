# -----------------------------------------------------------------------------#
# 01_processing_WEOM.R
#
# Purpose
#   Import the published FTICR-MS peak table (water-extractable organic
#   matter; WEOM) and generate analysis-ready outputs.
#
# Inputs
#   HJA_WS10_2023_WEOM_FTICR-MS.csv
#   (published, already-standardized table: molecular formulas, elemental
#   counts, and derived properties -- HC, OC, DBE, ai_mod, NOSC, GFE,
#   el_form, class -- arrive pre-computed, and sample columns are already
#   in the standardized WEOM_<Site><Pit>_<Event>_<Top|Btm>_<Rep>_<mRep>
#   format)
#
# Outputs (written to OD)
#   1) Intensities table:
#        - Peak heights for each molecular formula x sample
#   2) Presence-absence table:
#        - Binary detection (1 = detected, 0 = not detected) for each formula x sample
#   3) Molecular properties table:
#        - Per-formula properties (elemental counts, H/C, O/C, AImod, NOSC, etc.)
#   4) Sample-average molecular properties table:
#        - Per-sample intensity-weighted (or PA-weighted) means of key properties
#        - Per-sample composition summaries (percent elemental groups and compound classes)
#   5) Raw table:
#        - Full published table, property column names unchanged from
#          the published file
#
# Notes
#   - AImod, NOSC, and GFE are published pre-computed (Koch & Dittmar, 2016;
#     LaRowe & Van Cappellen, 2011); the NA/Inf cleanup on AImod is kept as
#     a safety check on the published values, not a recalculation.
#   - Percent metrics for elemental groups and compound classes are based on
#     detection (presence/absence), not intensity weighting.
#
# Created:  2026-01-25
# Author:   X. Takver
# -----------------------------------------------------------------------------#

## ----------------------------------------- ##
#             Package Import ----
## ----------------------------------------- ##

# Data Processing Packages
librarian::shelf(tidyverse, stringr, here)

## ----------------------------------------- ##
#             Source Functions ----
## ----------------------------------------- ##

# Source helper functions
source(here('R', 'functions', 'paths_and_directories.R'))
source(here('R', 'functions', 'fticr_sample_means_function.R'))

## ----------------------------------------- ##
#           Directory Navigation ----
## ----------------------------------------- ##

# Makes quick navigation directory commands
# DD: directory containing the published FTICR-MS table
# OD: directory for processed outputs written by this script

DD <- here('data', 'mass_spectrometry', 'raw')
OD <- here('data', 'mass_spectrometry', 'processed', 'water-extractable_om')

# Ensure output directories exist (idempotent)
ensure_dir(OD)

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

# Column names are kept exactly as published (molecular_formula, mass,
# HC, OC, ai_mod, el_form, etc.); all property values are pre-computed.
df_raw <- read_csv(file.path(DD, 'HJA_WS10_2023_WEOM_FTICR-MS.csv')) |>
  mutate(ai_mod = ifelse(is.na(ai_mod), 0, ai_mod),
         ai_mod = ifelse(ai_mod == 'Inf', 0, ai_mod),
         ai_mod = ifelse(ai_mod == '-Inf', 0, ai_mod))

## ----------------------------------------- ##
#               Data Separation ----
## ----------------------------------------- ##

## ----------------------------------------- ##
##           Intensities Data ----
## ----------------------------------------- ##

# Create intensities-only matrix: rows = molecular_formula, columns = sample spectra.
# Sample columns are already standardized in the published file, so no
# splitting/renaming across sampling periods is needed here.
df_intensities <- df_raw |>
  column_to_rownames('molecular_formula') |>
  select(-c(mass:el_form)) |>
  rownames_to_column('molecular_formula')

## ----------------------------------------- ##
##          Presence-Absence Data ----
## ----------------------------------------- ##

# Convert intensities to detection (0/1).
# Missing values are treated as non-detections.
df_PA <- df_raw |>
  column_to_rownames('molecular_formula') |>
  select(-c(mass:el_form)) |>
  mutate(across(everything(), ~replace_na(., 0)),
         across(everything(), ~as.numeric(.>0))) |>
  rownames_to_column('molecular_formula')

## ----------------------------------------- ##
##          Mol Property Data ----
## ----------------------------------------- ##

# Retain per-formula molecular properties used for downstream summarisation and plotting.
df_mol_prop <- df_raw |>
  select(molecular_formula, mass:el_form)

## ----------------------------------------- ##
##         Avg Mol Property Data ----
## ----------------------------------------- ##

# Calculate Sample Average Molecular Properties
# (wmean() and calc_fticr_means() are sourced from
# fticr_sample_means_function.R above)

# Long-form table linking each detected peak (or intensity) to its molecular properties.
df_long <- df_PA |>
  pivot_longer(cols = -molecular_formula,
               names_to = 'Sample',
               values_to = 'w') |>
  left_join(df_mol_prop, by = 'molecular_formula') |>
  mutate(present = w != 0)

# Compute per-sample summary statistics.
# Percent metrics are computed over the set of detected formulas (presence/absence).

fticr_means <- calc_fticr_means(
  df_long,
  el_forms = c('CHO', 'CHON', 'CHOS', 'CHOP', 'CHONS', 'CHONP', 'CHOSP', 'CHONSP')
)

# Parse experimental metadata from the standardized sample identifiers.
# Expected format: WEOM_<HS|LS><Pit>_<June|December>_<Top|Btm>_R<Rep>_mR<mRep>
df_avg_mol_prop <- fticr_means |>
  mutate(Hillslope = case_when(str_detect(Sample, 'HS') ~ 'High Storage',
                               str_detect(Sample, 'LS') ~ 'Low Storage',
                               TRUE ~ NA_character_),
         Sampling.Event = str_extract(Sample, 'June|December'),
         Rep  = as.integer(str_extract(Sample, '(?<=_R)\\d+')),
         mRep = as.integer(str_extract(Sample, '(?<=_mR)\\d+')),
         Core.Section = case_when(str_detect(Sample, 'Top') ~ 'Top',
                                  str_detect(Sample, 'Btm') ~ 'Bottom',
                                  TRUE ~ NA_character_),
         .after = Sample)

## ----------------------------------------- ##
##              Export Data ----
## ----------------------------------------- ##

# Export df_intensities
write_csv(df_intensities, file.path(OD, '2023_MONet_FTICR-MS_WEOM_Intensities.csv'))

# Save Presence Absence Data
write_csv(df_PA, file.path(OD, '2023_MONet_FTICR-MS_WEOM_Presence-Absence.csv'))

# Save molecular property data
write_csv(df_mol_prop, file.path(OD, '2023_MONet_FTICR-MS_WEOM_Molecular-Properties.csv'))

# Save Average Molecular Properties
write_csv(df_avg_mol_prop, file.path(OD, '2023_MONet_FTICR-MS_WEOM_Avg_Molecular-Properties.csv'))

# Save Raw Data with updated column names
write_csv(df_raw, file.path(OD, '2023_MONet_FTICR-MS_WEOM_Raw_Updated_ColNames.csv'))
