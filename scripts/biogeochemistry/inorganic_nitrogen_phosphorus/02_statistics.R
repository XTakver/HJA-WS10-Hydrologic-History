# ------------------------------------------------------------ #
# 02_statistics_Nitrogen_Phosphorus.R
#
# Purpose
#   Fit simple linear models for extractable inorganic nitrogen (NH4-N, NO3-N),
#   bicarbonate-extractable phosphorus (P_bicarb), and derived ratios, then
#   export model summaries and emmeans seasonal contrasts within hillslope x depth.
#
# Inputs
#   data/biogeochemistry/processed/2023_Inorganic_Nitrogen_and_Phosphorus_processed.csv
#
# Outputs (written to SD)
#   - nh4n_lm_summary.txt
#   - no3n_lm_summary.txt
#   - p_bicarb_lm_summary.txt
#   - nh4_no3_ratio_lm_summary.txt
#   - np_ratio_lm_summary.txt
#   - *_emmeans_sampling_month_by_hillslope_core_section.csv
#   - *_pairwise_sampling_month_by_hillslope_core_section.csv
#
# Notes
#   - Ratios:
#       * NH4:NO3 uses a small adjustment for NO3 = 0 to avoid division-by-zero.
#       * N:P uses Total_N (NH4 + NO3) divided by P_bicarb.
#
# Created:  2026-02-20
# Author:   X. Takver
# ------------------------------------------------------------ #

## ----------------------------------------- ##
#             Package Import ----
## ----------------------------------------- ##

# Wrangling + exports
librarian::shelf(tidyverse, here)

# Post-hoc contrasts
librarian::shelf(emmeans)

## ----------------------------------------- ##
#             Source Functions ----
## ----------------------------------------- ##

# ensure_dir()
source(here('R', 'functions', 'paths_and_directories.R'))

## ----------------------------------------- ##
#           Directories + paths ----
## ----------------------------------------- ##

# DD: processed products from 01_processing_Nitrogen_Phosphorus.R
DD <- here('data', 'biogeochemistry', 'processed')

# SD: stats outputs written by this script
SD <- here('output', 'statistics', 'biogeochemistry')
ensure_dir(SD)

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

df_NP <- read_csv(file.path(DD, '2023_Inorganic_Nitrogen_and_Phosphorus_processed.csv'))

## ----------------------------------------- ##
#             Data Preparation ----
## ----------------------------------------- ##

# Minimal guardrails: re-establish factor levels used by the models,
# since the CSV round-trip stores hillslope/sampling_month/core_section
# as plain character columns
df_NP <- df_NP |>
  mutate(
    hillslope = as_factor(hillslope),
    sampling_month = fct_relevel(as_factor(sampling_month), 'June', 'December'),
    core_section = fct_relevel(as_factor(core_section), 'Top', 'Bottom')
  )

# Derived variables for ratios and totals
df_NP_final <- df_NP |>
  mutate(
    no3_adj = if_else(is.na(NO3N), NA_real_, if_else(NO3N == 0, 0.0005, NO3N)),
    nh4_no3_ratio = NH4N / no3_adj,
    total_n = NH4N + NO3N,
    np_ratio = total_n / P_bicarb
  )

np_response_lookup <- tibble(
  response = c('NH4N', 'NO3N', 'P_bicarb', 'nh4_no3_ratio', 'np_ratio'),
  parameter = c('NH4-N', 'NO3-N', 'Bicarbonate-extractable P',
                'NH4:NO3', 'N:P'),
  units = c('mg N kg soil^-1', 'mg N kg soil^-1', 'mg P kg soil^-1',
            'unitless', 'unitless')
)

## ----------------------------------------- ##
#             Summary Statistics ----
## ----------------------------------------- ##

np_raw_summary <- df_NP_final |>
  pivot_longer(
    cols = all_of(np_response_lookup$response),
    names_to = 'response',
    values_to = 'value'
  ) |>
  left_join(np_response_lookup, by = 'response') |>
  group_by(parameter, response, units, hillslope, core_section, sampling_month) |>
  summarise(
    core_section_n = 1,
    measurement_n = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    .groups = 'drop'
  ) |>
  arrange(parameter, hillslope, core_section, sampling_month)

write_csv(
  np_raw_summary,
  file.path(SD, 'Nitrogen_Phosphorus_raw_summary_mean_sd_n.csv')
)

## ----------------------------------------- ##
#           Statistical Analysis ----
## ----------------------------------------- ##

## ------------------------------- ##
##             Linear Models ----

nh4_lm <- lm(NH4N ~ hillslope * sampling_month * core_section, data = df_NP_final)
no3_lm <- lm(NO3N ~ hillslope * sampling_month * core_section, data = df_NP_final)
p_lm   <- lm(P_bicarb ~ hillslope * sampling_month * core_section, data = df_NP_final)

# Ratio models (guard against missing values introduced by no3_adj or P_bicarb)
nh4no3_lm <- lm(nh4_no3_ratio ~ hillslope * sampling_month * core_section, data = df_NP_final)
np_lm     <- lm(np_ratio ~ hillslope * sampling_month * core_section, data = df_NP_final)

## ------------------------------- ##
##        Post-hoc Comparisons ----

# Seasonal differences (sampling_month) within each hillslope x core_section stratum.
nh4_emm <- emmeans(nh4_lm, ~ sampling_month | hillslope * core_section)
no3_emm <- emmeans(no3_lm, ~ sampling_month | hillslope * core_section)
p_emm   <- emmeans(p_lm,   ~ sampling_month | hillslope * core_section)

nh4no3_emm <- emmeans(nh4no3_lm, ~ sampling_month | hillslope * core_section)
np_emm     <- emmeans(np_lm,     ~ sampling_month | hillslope * core_section)

nh4_pairs    <- pairs(nh4_emm)
no3_pairs    <- pairs(no3_emm)
p_pairs      <- pairs(p_emm)
nh4no3_pairs <- pairs(nh4no3_emm)
np_pairs     <- pairs(np_emm)

## ----------------------------------------- ##
#                  Export ----
## ----------------------------------------- ##

# Model summaries (plain text for quick manuscript reporting)
writeLines(capture.output(summary(nh4_lm)), file.path(SD, 'nh4n_lm_summary.txt'))
writeLines(capture.output(summary(no3_lm)), file.path(SD, 'no3n_lm_summary.txt'))
writeLines(capture.output(summary(p_lm)),   file.path(SD, 'p_bicarb_lm_summary.txt'))
writeLines(capture.output(summary(nh4no3_lm)), file.path(SD, 'nh4_no3_ratio_lm_summary.txt'))
writeLines(capture.output(summary(np_lm)),     file.path(SD, 'np_ratio_lm_summary.txt'))

# Estimated marginal means (CSV)
write_csv(as.data.frame(nh4_emm), file.path(SD, 'nh4n_emmeans_sampling_month_by_hillslope_core_section.csv'))
write_csv(as.data.frame(no3_emm), file.path(SD, 'no3n_emmeans_sampling_month_by_hillslope_core_section.csv'))
write_csv(as.data.frame(p_emm),   file.path(SD, 'p_bicarb_emmeans_sampling_month_by_hillslope_core_section.csv'))
write_csv(as.data.frame(nh4no3_emm), file.path(SD, 'nh4_no3_ratio_emmeans_sampling_month_by_hillslope_core_section.csv'))
write_csv(as.data.frame(np_emm),     file.path(SD, 'np_ratio_emmeans_sampling_month_by_hillslope_core_section.csv'))

# Pairwise seasonal contrasts (CSV)
write_csv(as.data.frame(nh4_pairs), file.path(SD, 'nh4n_pairwise_sampling_month_by_hillslope_core_section.csv'))
write_csv(as.data.frame(no3_pairs), file.path(SD, 'no3n_pairwise_sampling_month_by_hillslope_core_section.csv'))
write_csv(as.data.frame(p_pairs),   file.path(SD, 'p_bicarb_pairwise_sampling_month_by_hillslope_core_section.csv'))
write_csv(as.data.frame(nh4no3_pairs), file.path(SD, 'nh4_no3_ratio_pairwise_sampling_month_by_hillslope_core_section.csv'))
write_csv(as.data.frame(np_pairs),     file.path(SD, 'np_ratio_pairwise_sampling_month_by_hillslope_core_section.csv'))

## ----------------------------------------- ##
#                 Cleanup ----
## ----------------------------------------- ##

rm(df_NP, df_NP_final)
