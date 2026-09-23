# -----------------------------------------------------------------------------#
# 04_unique_formula_subsetting_MAOM.R
#
# Purpose
#   Identifies MAOM molecular formulas that are unique to a single
#   sample (Site x Sampling.Event combination) within each core
#   section, and recomputes sample-average molecular properties from
#   that unique-formula subset.
#
# Inputs
#   2023_FTICR-MS_mineral-associated_om_aggregated.csv
#   2023_FTICR-MS_mineral-associated_om_ChemProperties_aggregated.csv
#   (both written by 03_replicate_aggregation_MAOM.R)
#
# Outputs (written to OD)
#   - 2023_FTICR-MS_mineral-associated_om_Presence-Absence_Unique.csv
#   - 2023_FTICR-MS_mineral-associated_om_SampleAvg_Unique.csv
#
# Notes
#   - "Unique" = detected in exactly one sample-level column within a
#     given core section (Top/Bottom); uniqueness is pooled across all
#     site x event combinations present for that section, not assessed
#     separately per event.
#
# Created:  2026-01-26
# Author:   X. Takver
# -----------------------------------------------------------------------------#

## ----------------------------------------- ##
#             Import Packages ----
## ----------------------------------------- ##

# Data Processing
librarian::shelf(tidyverse, here)

## ----------------------------------------- ##
#             Source Functions ----
## ----------------------------------------- ##

# Source helper functions
source(here('R', 'functions', 'paths_and_directories.R'))
source(here('R', 'functions', 'fticr_sample_means_function.R'))

## ----------------------------------------- ##
#             Directory Setup ----
## ----------------------------------------- ##

# Inputs are replicate-aggregated (sample-level) MAOM outputs.
DD <- here('data', 'mass_spectrometry', 'processed', 'mineral-associated_om', 'replicate_aggregation')
OD <- file.path(dirname(DD), 'unique_formulas')

# Ensure output directories exist (idempotent)
ensure_dir(OD)

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

hcl_pa <- read_csv(file.path(DD, '2023_FTICR-MS_mineral-associated_om_aggregated.csv'))
hcl_formulas <- read_csv(file.path(DD, '2023_FTICR-MS_mineral-associated_om_ChemProperties_aggregated.csv'))

## ----------------------------------------- ##
#            Unique Formulas ----
## ----------------------------------------- ##

# Define 'unique' as formulas detected in exactly one sample within each core section (Top/Bottom).
hcl_unique <- hcl_pa |> 
  pivot_longer(cols = -molecular_formula,
               names_to = 'Sample.ID',
               values_to = 'PA') |> 
  separate(Sample.ID, into = c('SiteID', 'Sampling.Event', 'Core.Section')) |> 
  group_by(molecular_formula, Core.Section) |> 
  filter(sum(PA) == 1) |> 
  ungroup() |> 
  unite('SampleID', SiteID, Sampling.Event, Core.Section, sep = '_') 

# Reconstitute wide presence/absence table for unique formulas only.
hcl_unique_wide <- hcl_unique |> 
  pivot_wider(names_from = SampleID,
              values_from = PA, 
              values_fill = 0)

# Subset molecular properties to formulas retained in the unique set.
hcl_unique_formulas <- hcl_unique_wide |> 
  select(molecular_formula) |> 
  left_join(hcl_formulas)

## ----------------------------------------- ##
#       Recalculate Avg Properties ----
## ----------------------------------------- ##

# Join sample-level detections to molecular properties to compute per-sample means and composition metrics.
df_long <- hcl_unique_wide |>
  pivot_longer(cols = -molecular_formula,
               names_to = 'Sample',
               values_to = 'w') |>
  left_join(hcl_unique_formulas, by = 'molecular_formula') |>
  mutate(present = w != 0)

# Recalculate average properties based on sample-level molecular formulas.
# el_forms is limited to the P-free groups, since phosphorus is not
# present in the published MAOM data.
df_hcl_AvgProp <- 
  calc_fticr_means(df_long,
                   el_forms = c('CHO', 'CHON', 'CHOS', 'CHONS'))

## ----------------------------------------- ##
#           Export Dataframes ----
## ----------------------------------------- ##

write_csv(hcl_unique_wide, file.path(OD, '2023_FTICR-MS_mineral-associated_om_Presence-Absence_Unique.csv'))
write_csv(df_hcl_AvgProp, file.path(OD, '2023_FTICR-MS_mineral-associated_om_SampleAvg_Unique.csv'))
