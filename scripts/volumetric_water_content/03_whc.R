# -----------------------------------------------------------------------------#
# 03_whc.R
#
# Purpose
#   Import smoothed VWC time series, convert to water holding content (WHC; % of
#   series maximum by site_id x depth), linearly interpolate gaps, and quantify
#   time spent within an "optimal" WHC window (50-70%) by hillslope and depth.
#
# Inputs (read from DD)
#   - smoothed_vwc.csv
#   - HR_Calculations_June-Aug.csv  (imported for reference; not used below)
#
# Outputs (written to SD)
#   - WHC_Optimal_Days_Comparison.csv
#       Wide-format table of mean optimal days by hillslope x depth + HS/LS ratio.
#   Also written to OD:
#   - whc.csv (interpolated per-timestamp WHC series)
#
# Created:  2025-11-25
# Author:   X. Takver
# -----------------------------------------------------------------------------#

## ----------------------------------------- ##
#             Package Import ----
## ----------------------------------------- ##

# Data Analysis
librarian::shelf(tidyverse, pracma, quantmod, zoo, mgcv, here)

## ----------------------------------------- ##
##            Source functions ----
## ----------------------------------------- ##

# Project-wide helpers (consistent with manuscript scripts)
#   - ensure_dir(), font registration, set_plot_theme
source(here('R', 'functions', 'paths_and_directories.R')) 
source(here('R', 'functions', 'fonts.R'))                  
source(here('R', 'functions', 'plot_theme.R'))             

# VWC utilities (e.g., shared VWC helpers across scripts)
source(here('R', 'functions', 'vwc_processing_functions.R'))

## ----------------------------------------- ##
#           Directory Creation ----
## ----------------------------------------- ##

# DD: processed VWC inputs
# OD: processed outputs (not written by this script; retained for convention)
# SD: statistics outputs written by this script
DD <- here('data', 'volumetric_water_content', 'processed')
OD <- here('data', 'volumetric_water_content', 'processed')
SD <- here('output', 'statistics', 'volumetric_water_content', 'HR')


ensure_dir(OD)
ensure_dir(SD)

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

# Smoothed VWC (Sentek); site_id already matches the published/manuscript convention
df_vwc <- read_csv(file.path(DD, 'smoothed_vwc.csv'))

# HR calculations (imported for reference / potential downstream comparisons)
df_HR <- read_csv(file.path(DD, 'HR_Calculations_June-Aug.csv'))

## ----------------------------------------- ##
#             Data Processing ----
## ----------------------------------------- ##

## ----------------------------------------- ##
##               Convert to WHC ----

# Compute series maximum VWC by site_id x depth, then express VWC as % of max.
df_max <- df_vwc |> 
  group_by(site_id, depth) |> 
  summarize(max_vwc = max(Smooth_VWC_Avg, na.rm = TRUE), .groups = 'drop')

df_whc <- df_vwc |> 
  left_join(df_max) |> 
  mutate(whc = (Smooth_VWC_Avg/max_vwc * 100))

## ----------------------------------------- ##
##           Interpolate WHC gaps ----

# Linear interpolation across gaps (used for duration-based summaries).
df_whc_interp <- df_whc |>
  arrange(site_id, depth, datetime) |>
  group_by(site_id, depth) |>
  mutate(whc_interp = na.approx(whc, x = as.numeric(datetime),
                                na.rm = FALSE, maxgap = Inf))

## ----------------------------------------- ##
##           Export WHC Dataframe ----

write_csv(df_whc_interp, file.path(OD, 'whc.csv'))

## ----------------------------------------- ##
#               Calculations ----
## ----------------------------------------- ##

## ----------------------------------------- ##
#               Optimal WHC ----

# "Optimal" defined as WHC between 50-70% (inclusive). Assumes 10-minute timestep.
whc_summary <- df_whc_interp |>
  mutate(optimal = whc_interp >= 50 & whc_interp <= 70 ) |> 
  group_by(site_id, depth) |> 
  summarize(optimal_timestamps = sum(optimal, na.rm = TRUE)) |> 
  mutate(optimal_minutes = optimal_timestamps * 10,
         optimal_hours = optimal_timestamps/6,
         optimal_days = optimal_hours/24) |> 
  ungroup()


# Summarize optimal time by hillslope and depth (mean +/- SD/SE across pits).
whc_analysis <- whc_summary |> 
  mutate(hillslope = case_when(str_detect(site_id, 'LS') ~ 'Low Storage',
                               str_detect(site_id, 'HS') ~ 'High Storage'),
         .after = site_id) |> 
  group_by(hillslope, depth) |> 
  summarize(mean_opt_days = mean(optimal_days, na.rm = TRUE),
            sd_days_opt = sd(optimal_days, na.rm = TRUE),
            n_sites = n(),
            se_days_opt = sd_days_opt/sqrt(n_sites),
            .groups = 'drop')

# Range summaries (depth <= 60 cm); retained as exploratory outputs in object space.
opt_days_range <- whc_analysis |> 
  select(hillslope, depth, mean_opt_days) |> 
  filter(depth <= 60) |> 
  group_by(hillslope) |> 
  summarize(max_days = max(mean_opt_days),
            min_days = min(mean_opt_days),
            avg_days = (max_days + min_days)/2)

opt_days_range_depth <- whc_analysis |> 
  select(hillslope, depth, mean_opt_days) |> 
  filter(depth <= 60) |> 
  group_by(hillslope) |> 
  summarize(max_days = max(mean_opt_days),
            min_days = min(mean_opt_days),
            avg_days = (max_days + min_days)/2)

# Wide comparison table and HS/LS ratio by depth.
opt_days_comparison <- whc_analysis |> 
  select(hillslope, depth, mean_opt_days) |> 
  pivot_wider(names_from = hillslope, 
              values_from = mean_opt_days) |> 
  mutate(HS_vs_LS = `High Storage`/`Low Storage`)

# Output Optimal days
write_csv(opt_days_comparison, file.path(SD, 'WHC_Optimal_Days_Comparison.csv'))
