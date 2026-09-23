# ------------------------------------------------------------ #
# 01_processing_Ions_Nutrients.R
#
# Purpose
#   Process WS10 MONet soil ion and macro/micronutrient concentrations
#   (ppm and meq/100g) into long-format tables for plotting, and
#   compute derived organic sulfur and base saturation.
#
# Inputs
#   HJA_WS10_2023_Soil_Ion_and_Nutrients.csv
#   (published, already-combined file; sample identifiers arrive
#   pre-populated, so no raw import/binding is done in this script)
#
# Outputs (written to OD)
#   - 2023_Soil_Ion_and_Nutrients_processed.csv (wide, as published)
#   - 2023_Soil_Ion_meq_long.csv
#   - 2023_Soil_Ion_ppm_long.csv        (includes derived organic S)
#   - 2023_Soil_Ion_base_saturation.csv (derived: total bases / CEC x 100)
#
# Notes
#   - Base saturation = (total_bases_meq / cec_meq) * 100.
#   - Organic S (ppm) = total_S_ppm - SO4S_ppm.
#   - K is published in both ppm and meq; K_ppm is excluded from the
#     ppm long table so K appears only once, via K_meq in the major
#     cations panel.
#   - These ion/nutrient measurements were collected once per
#     hillslope x core section, in June only (see manuscript
#     Methods 2.5.2); there is no December counterpart.
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

source(here('R', 'functions', 'paths_and_directories.R'))  # ensure_dir()

## ----------------------------------------- ##
#             Directories ----
## ----------------------------------------- ##

DD <- here('data', 'biogeochemistry', 'raw')
OD <- here('data', 'biogeochemistry', 'processed')

ensure_dir(OD)

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

df_ions <- read_csv(file.path(DD, 'HJA_WS10_2023_Soil_Ion_and_Nutrients.csv'))

## ----------------------------------------- ##
#            Factor Conversion ----
## ----------------------------------------- ##

# Sets grouping/plotting order for downstream figures
df_ions <- df_ions |>
  mutate(
    hillslope = as_factor(hillslope),
    core_section = fct_relevel(as_factor(core_section), 'Top', 'Bottom')
  )

## ----------------------------------------- ##
#      Identifier Columns ----
## ----------------------------------------- ##

# Carried through both long-format pivots below
id_cols <- c('doe_proposal_id', 'doe_sample_name', 'site_id', 'pit_id',
             'hillslope', 'core_section', 'sampling_month')

## ----------------------------------------- ##
#        Pivot meq fields to long ----
## ----------------------------------------- ##

meq_cols <- names(df_ions)[str_detect(names(df_ions), regex('_meq$', ignore_case = TRUE))]

df_meq_long <- df_ions |>
  select(all_of(id_cols), all_of(meq_cols)) |>
  pivot_longer(
    cols = all_of(meq_cols),
    names_to = 'ion_raw',
    values_to = 'meq'
  ) |>
  mutate(ion = str_remove(ion_raw, regex('_meq$', ignore_case = TRUE))) |>
  select(-ion_raw)

## ----------------------------------------- ##
#      Base Saturation (total bases & CEC) ----
## ----------------------------------------- ##

df_cec <- df_meq_long |>
  filter(ion %in% c('total_bases', 'cec')) |>
  pivot_wider(names_from = ion, values_from = meq)

df_base_saturation <- df_cec |>
  mutate(
    base_saturation = if_else(!is.na(total_bases) & !is.na(cec) & cec != 0,
                              (total_bases / cec) * 100,
                              NA_real_)
  ) |>
  select(all_of(id_cols), total_bases, cec, base_saturation)

## ----------------------------------------- ##
#        Pivot ppm fields to long ----
## ----------------------------------------- ##

ppm_cols <- names(df_ions)[str_detect(names(df_ions), regex('_ppm$', ignore_case = TRUE))]

# Organic S is not published directly: total_S includes both
# organic and sulfate-bound S, so subtract SO4-S to isolate it
df_ppm_long <- df_ions |>
  select(all_of(id_cols), all_of(ppm_cols)) |>
  mutate(organic_S_ppm = total_S_ppm - SO4S_ppm) |>
  pivot_longer(
    cols = matches('_ppm$', ignore.case = TRUE),
    names_to = 'ion_raw',
    values_to = 'ppm'
  ) |>
  mutate(
    ion = str_remove(ion_raw, regex('_ppm$', ignore_case = TRUE)),
    ion = case_when(ion == 'organic_S' ~ 'Org-S',
                    TRUE ~ ion)
  ) |>
  # K is reported in both ppm and meq; only the meq (exchangeable)
  # form belongs in the major-cations panel, so K_ppm is dropped
  # here to keep this table limited to micronutrients/trace elements
  filter(!ion %in% c('total_S', 'K')) |>
  select(-ion_raw)

## ----------------------------------------- ##
#                  Export ----
## ----------------------------------------- ##

write_csv(df_ions, file.path(OD, '2023_Soil_Ion_and_Nutrients_processed.csv'))
write_csv(df_meq_long, file.path(OD, '2023_Soil_Ion_meq_long.csv'))
write_csv(df_ppm_long, file.path(OD, '2023_Soil_Ion_ppm_long.csv'))
write_csv(df_base_saturation, file.path(OD, '2023_Soil_Ion_base_saturation.csv'))
