# ------------------------------------------------------------ #
# 01_processing_Nitrogen_Phosphorus.R
#
# Purpose
#   Process extractable inorganic nitrogen (NH4-N, NO3-N) and
#   bicarbonate-extractable phosphorus for WS10 MONet biogeochemistry.
#
# Inputs
#   HJA_WS10_2023_Inorganic_Nitrogen_and_Phosphorus.csv
#
# Outputs (written to OD)
#   - 2023_Inorganic_Nitrogen_and_Phosphorus_processed.csv
#
# Created:  2026-02-20
# Author: X. Takver
# ------------------------------------------------------------ #

## ----------------------------------------- ##
#             Package Import ----
## ----------------------------------------- ##

librarian::shelf(tidyverse, here)

## ----------------------------------------- ##
#             Source Functions ----
## ----------------------------------------- ##

# Directory helpers (ensure_dir)
source(here('R', 'functions', 'paths_and_directories.R'))

## ----------------------------------------- ##
#        Directory helpers + paths ----
## ----------------------------------------- ##

# DD: raw/published biogeochem inputs
# OD: processed products written by this script
DD <- here('data', 'biogeochemistry', 'raw')
OD <- here('data', 'biogeochemistry', 'processed')

ensure_dir(OD)

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

df_NP <- read_csv(file.path(DD, 'HJA_WS10_2023_Inorganic_Nitrogen_and_Phosphorus.csv'))

## ----------------------------------------- ##
#            Factor Conversion ----
## ----------------------------------------- ##

# Sets grouping/plotting order for downstream models and figures
df_NP <- df_NP |>
  mutate(
    hillslope = as_factor(hillslope),
    sampling_month = fct_relevel(as_factor(sampling_month), 'June', 'December'),
    core_section = fct_relevel(as_factor(core_section), 'Top', 'Bottom')
  )

## ----------------------------------------- ##
#                  Export ----
## ----------------------------------------- ##

write_csv(df_NP, file.path(OD, '2023_Inorganic_Nitrogen_and_Phosphorus_processed.csv'))
