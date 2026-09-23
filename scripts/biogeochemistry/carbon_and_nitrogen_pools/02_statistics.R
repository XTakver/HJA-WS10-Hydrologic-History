# ------------------------------------------------------------ #
#
# Purpose
#   Fit simple linear models for WEOM C, MAOM C, and microbial
#   biomass C (MBC), and export emmeans estimates + pairwise
#   seasonal contrasts within hillslope x core_section strata.
#
# Inputs (from processing script)
#   2023_Carbon_and_Nitrogen_pools_processed.csv
#
# Outputs (written to SD)
#   - *_lm_summary.txt
#   - *_emmeans_sampling_month_by_hillslope_core_section.csv
#   - *_pairwise_sampling_month_by_hillslope_core_section.csv
#   - WEOM_MAOM_mBioC_model_objects.rds
#
# Notes
#   - Model formula: response ~ hillslope * sampling_month * core_section
#
# Created: 2026-02-19
# Author: X. Takver
# ------------------------------------------------------------ #

## ----------------------------------------- ##
#             Package Import ----
## ----------------------------------------- ##

# Core wrangling + model outputs
librarian::shelf(tidyverse, here, emmeans, glue, broom)

## ----------------------------------------- ##
#             Source Functions ----
## ----------------------------------------- ##

# Project-wide helpers (ensure_dir)
source(here('R', 'functions', 'paths_and_directories.R'))

## ----------------------------------------- ##
#           Directory Navigation ----
## ----------------------------------------- ##

# OD: processed biogeochemistry products
# SD: statistics outputs for manuscript tables/figures
OD <- here('data', 'biogeochemistry', 'processed')
SD <- here('output', 'statistics', 'biogeochemistry')

ensure_dir(SD)

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

df_OCarbon <- read_csv(file.path(OD, '2023_Carbon_and_Nitrogen_pools_processed.csv'))

# Minimal guardrails: re-establish factor levels used by the models,
# since the CSV round-trip stores hillslope/sampling_month/core_section
# as plain character columns
df_OCarbon <- df_OCarbon |>
  mutate(
    hillslope = as_factor(hillslope),
    sampling_month = fct_relevel(as_factor(sampling_month), 'June', 'December'),
    core_section = fct_relevel(as_factor(core_section), 'top', 'bottom')
  )

# Responses modeled in this script
response_lookup <- tibble(
  response = c('weom_c', 'maom_c', 'mbc'),
  parameter = c('WEOM-C', 'MAOM-C', 'Microbial biomass C'),
  units = c('mg C kg soil^-1', 'mg C kg soil^-1', 'mg C kg soil^-1')
)

## ----------------------------------------- ##
#           Statistical Analysis ----
## ----------------------------------------- ##

## ----------------------------------------- ##
#          Raw Data Summary Table ----
## ----------------------------------------- ##

# This directly addresses the reviewer request:
# mean, sd, and n for each sampling_month × hillslope × core_section group.

raw_summary <- df_OCarbon |>
  pivot_longer(
    cols = all_of(response_lookup$response),
    names_to = 'response',
    values_to = 'value'
  ) |>
  left_join(response_lookup, by = 'response') |>
  group_by(parameter, response, units, hillslope, core_section, sampling_month) |>
  summarise(
    n = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    se = sd / sqrt(n),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = 'drop'
  ) |>
  arrange(parameter, hillslope, core_section, sampling_month)

write_csv(
  raw_summary,
  file.path(SD, 'WEOM_MAOM_mBioC_raw_summary_mean_sd_n.csv')
)

# Helper: fit model + emmeans + seasonal contrasts, using a consistent structure
fit_lm_emm <- function(df, response) {
  fml <- as.formula(glue('{response} ~ hillslope * sampling_month * core_section'))
  mod <- lm(fml, data = df)
  
  # Seasonal means within each hillslope × depth stratum
  emm <- emmeans(mod, ~ sampling_month | hillslope * core_section)
  
  # Pairwise seasonal contrasts (December − June by default ordering of sampling_month)
  prs <- pairs(emm)
  
  list(model = mod, emmeans = emm, pairs = prs)
}

## ----------------------------------------- ##
##              Linear Models ----

# Fit models for each response (mg/kg)
weom_res <- fit_lm_emm(df_OCarbon, 'weom_c')
maom_res <- fit_lm_emm(df_OCarbon, 'maom_c')
mbio_res <- fit_lm_emm(df_OCarbon, 'mbc')

## ----------------------------------------- ##
#                  Export ----
## ----------------------------------------- ##

# --- Model summaries (plain text) ---
writeLines(capture.output(summary(weom_res$model)),
           file.path(SD, 'weom_c_lm_summary.txt'))
writeLines(capture.output(summary(maom_res$model)),
           file.path(SD, 'maom_c_lm_summary.txt'))
writeLines(capture.output(summary(mbio_res$model)),
           file.path(SD, 'mbio_c_lm_summary.txt'))

# --- EMMs (CSV) ---
write_csv(as_tibble(as.data.frame(weom_res$emmeans)),
          file.path(SD, 'weom_c_emmeans_sampling_month_by_hillslope_core_section.csv'))
write_csv(as_tibble(as.data.frame(maom_res$emmeans)),
          file.path(SD, 'maom_c_emmeans_sampling_month_by_hillslope_core_section.csv'))
write_csv(as_tibble(as.data.frame(mbio_res$emmeans)),
          file.path(SD, 'mbio_c_emmeans_sampling_month_by_hillslope_core_section.csv'))

# --- Pairwise seasonal contrasts (CSV) ---
write_csv(as_tibble(as.data.frame(weom_res$pairs)),
          file.path(SD, 'weom_c_pairwise_sampling_month_by_hillslope_core_section.csv'))
write_csv(as_tibble(as.data.frame(maom_res$pairs)),
          file.path(SD, 'maom_c_pairwise_sampling_month_by_hillslope_core_section.csv'))
write_csv(as_tibble(as.data.frame(mbio_res$pairs)),
          file.path(SD, 'mbio_c_pairwise_sampling_month_by_hillslope_core_section.csv'))

# --- Save model objects for downstream plotting/reporting ---
saveRDS(
  list(
    weom_lm = weom_res$model,
    maom_lm = maom_res$model,
    mbio_lm = mbio_res$model,
    weom_emm = weom_res$emmeans,
    maom_emm = maom_res$emmeans,
    mbio_emm = mbio_res$emmeans,
    weom_pairs = weom_res$pairs,
    maom_pairs = maom_res$pairs,
    mbio_pairs = mbio_res$pairs
  ),
  file.path(SD, 'WEOM_MAOM_mBioC_model_objects.rds')
)

## ----------------------------------------- ##
#                 Cleanup ----
## ----------------------------------------- ##

rm(df_OCarbon, weom_res, maom_res, mbio_res)
