# -----------------------------------------------------------------------------#
# Purpose
#   Import WS10 soil texture (particle size fraction; PSF) data, aggregate
#   replicate runs by particle-size class, apply basic QC screening, and export
#   processed tables for downstream analysis and plotting.
#
# Inputs (read from DD)
#   - SP03903_v1.csv
#
# Outputs (written to OD)
#   - 2023_WS10-HJA_PSF_Processed.csv
#       Wide-format PSF table with replicate runs and summary statistics.
#   - 2023_WS10-HJA_PSF_Processed-long.csv
#       Long-format PSF summary (mean, SD, SE) by pit_id, depth, and size class.
#
# Created: 01/26/2026
# Author: X. Takver
# -----------------------------------------------------------------------------#

## ----------------------------------------- ##
#             Import Packages ----
## ----------------------------------------- ##

# Data processing packages
librarian::shelf(tidyverse, glue, here)

## ----------------------------------------- ##
#             Source Functions ----
## ----------------------------------------- ##

# Directory helpers (ensure_dir)
source(here('R', 'functions', 'paths_and_directories.R'))

## ----------------------------------------- ##
#           Directory Creation ----
## ----------------------------------------- ##

# DD: directory containing raw soil texture inputs
# OD: directory for processed outputs written by this script
DD <- here('data', 'soil_texture', 'raw')
OD <- here('data', 'soil_texture', 'processed')

# Ensure output directory exists (idempotent)
ensure_dir(OD)

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

# Import raw particle size fraction dataframe
df_texture <- read_csv(file.path(DD, 'SP03903_v1.csv'))

## ----------------------------------------- ##
#              Data Processing ----
## ----------------------------------------- ##

## ----------------------------------------- ##
##        Calculate Particle Size ----

# Process and aggregate raw particle-size fraction data:
#   - Create depth range label
#   - Sum % within each particle-size class for each instrument replicate
#   - Pivot replicates to wide format (one column per replicate run)

df_texture_processed <- df_texture |> 
  mutate(depth_range_cm = glue('{top_depth}-{bottom_depth}'),
         .after = bottom_depth) |> 
  group_by(site_id, depth_range_cm, analytical_rep, technical_rep, psize_class, horizon) |> 
  summarise(psf = sum(percent_in_range),
            .groups = 'drop') |> 
  pivot_wider(names_from = c(analytical_rep, technical_rep),
              values_from = psf,
              names_sep = '.') |> 
  ungroup()

## ----------------------------------------- ##
##      Compute Summary Statistics ----

# Compute replicate mean, SD, and CV for QC screening.
# QC thresholds:
#   - CV <= 0.20: QC_Pass
#   - 0.20 < CV <= 0.30: QC_Review
#   - CV > 0.30: QC_Fail
df_psf_QC <- 
  df_texture_processed |> 
  rowwise() |> 
  mutate(psf_mean = mean(c_across(-c(site_id:horizon)), na.rm = TRUE),
         psf_sd = sd(c_across(-c(site_id:horizon)), na.rm = TRUE),
         psf_cv = psf_sd/psf_mean,
         psf_n_runs = sum(!is.na(c_across(-c(site_id:horizon, psf_mean, psf_sd, psf_cv))))) |> 
  mutate(cv_flag = case_when(psf_cv <= .2 ~ 'QC_Pass',
                             psf_cv > .2 & psf_cv <= .3 ~ 'QC_Review',
                             psf_cv > .3 ~ 'QC_Fail'))

# Exclude QC failures and gravel class, then recompute summary statistics
df_psf_clean <- 
  df_psf_QC |> 
  filter(cv_flag != 'QC_Fail',
         psize_class != 'gravel') |> 
  select(-c(psf_mean:cv_flag)) |> 
  rowwise() |> 
  mutate(psf_mean = mean(c_across(-c(site_id:horizon)), na.rm = TRUE),
         psf_sd = sd(c_across(-c(site_id:horizon)), na.rm = TRUE),
         psf_cv = psf_sd/psf_mean,
         psf_n_runs = sum(!is.na(c_across(-c(site_id:horizon, psf_mean, psf_sd, psf_cv))))) |> 
  ungroup()

## ----------------------------------------- ##
##            Long Format Data----

# Convert replicate columns to long format for plotting and compute mid-depth
texture_long <- 
  df_psf_clean |> 
  select(!(psf_mean:psf_n_runs)) |> 
  pivot_longer(cols = c(`1.1`:`3.3`),
               names_to = 'instrument_replicate',
               values_to = 'psize_percent') |> 
  separate(depth_range_cm, into = c('bottom_depth_cm', 'top_depth_cm'), sep = '-') |> 
  mutate(bottom_depth_cm = parse_number(bottom_depth_cm),
         top_depth_cm = parse_number(top_depth_cm)) |> 
  mutate(mid_depth_cm = ((top_depth_cm - bottom_depth_cm)/2) + bottom_depth_cm, .after = 'top_depth_cm')

# Summary statistics by pit, depth, and particle-size class
df_long_means <- 
  texture_long |> 
  group_by(site_id, mid_depth_cm, psize_class) |> 
  summarise(mean_psf = mean(psize_percent, na.rm = TRUE),
            sd_psf = sd(psize_percent, na.rm = TRUE),
            total_reps = sum(!is.na(psize_percent)),
            se_psf   = sd_psf / sqrt(total_reps),
            .groups = 'drop') 

## ----------------------------------------- ##
##            Export PSF Data ----

# Export wide data (replicates + summary statistics)
# Export long data (summaries by depth and particle-size class)
write_csv(df_psf_clean, file.path(OD, '2023_WS10-HJA_PSF_Processed.csv'))
write_csv(df_long_means, file.path(OD, '2023_WS10-HJA_PSF_Processed-long.csv'))


