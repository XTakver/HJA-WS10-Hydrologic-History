# -----------------------------------------------------------------------------#
# 02_HR_and_predicted_vwc.R
#
# Purpose
#   Import smoothed VWC time series, calculate hydraulic redistribution (HR) and
#   no-HR predicted VWC, and summarize end-of-season minima and time spent in an
#   optimal WHC range (50-70%) by hillslope and depth.
#
# Inputs (read from DD)
#   - smoothed_vwc.csv
#
# Outputs (written to OD / SD)
#   - HR_Calculations_June-Aug.csv
#       Per-pit, per-depth HR calculations and predicted no-HR VWC time series.
#   - HR_hillslope_summary.csv
#       Hillslope-aggregated HR summary metrics.
#   - predicted_v_actual_vwc_end-of-season.csv
#       End-of-season actual vs predicted minima and percent additional water.
#   - Optimal_Time_Deficit_Pred.v.Actual.csv
#       Mean optimal-time deficit (predicted minus observed) by hillslope x depth.
#
# Created:  2026-02-10
# Author:   X. Takver
# -----------------------------------------------------------------------------#

## ----------------------------------------- ##
#             Package Import ----
## ----------------------------------------- ##

# Data Analysis
librarian::shelf(tidyverse, pracma, quantmod, zoo, mgcv, here)

# Data visualization
librarian::shelf(Cairo, tinythemes)

## ----------------------------------------- ##
##            Source functions ----
## ----------------------------------------- ##

# Project-wide helpers:
#   - ensure_dir(), font registration, set_plot_theme()
source(here('R', 'functions', 'paths_and_directories.R')) 
source(here('R', 'functions', 'fonts.R'))                  
source(here('R', 'functions', 'plot_theme.R'))             

# VWC/HR functions:
#   - calc_vwc_hr(), plus additional VWC utilities used below
source(here('R', 'functions', 'vwc_processing_functions.R'))

## ----------------------------------------- ##
#           Directory Creation ----
## ----------------------------------------- ##


# DD: processed VWC inputs
# OD: processed outputs written by this script
# SD: statistics outputs (summaries/tables)
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

## ----------------------------------------- ##
#                 HR Calcs ----
## ----------------------------------------- ##
###           Set Variables ----

# Analysis window for HR calculations
start_time <- ymd_hms("2024-06-17 00:00:00")
end_time   <- ymd_hms("2024-08-15 23:59:59")

# Target pits and depths (cm)
pits   <- c("LS1", "LS2", "HS1", "HS2")
depths <- c(10, 20, 30)

# Parameter grid (site_id x depth) for consistent iteration
param_grid <- tidyr::crossing(site_id = pits, depth = depths)

# Run HR calculation per series using calc_vwc_hr()
# Note: threshold/minpeakdistance/window control event detection and smoothing.
df_VWC_calcs <- param_grid %>%
  mutate(
    res = purrr::map2(
      site_id, depth,
      ~ calc_vwc_hr(
        pit_id     = .x,
        depth_cm   = .y,
        start_time = start_time,
        end_time   = end_time,
        df_vwc     = df_vwc,
        threshold       = 0.0003,
        minpeakdistance = 100,
        window          = 6,
        smooth_col      = "Smooth_VWC_Avg"
      )
    )
  ) %>%
  select(-site_id, -depth, res) %>%
  tidyr::unnest(res)

# Hillslope-level HR summaries (mean recharge/drawdown; observed vs predicted minima)
HR_summary <- df_VWC_calcs |> 
  mutate(hillslope = case_when(str_detect(site_id, 'HS') ~ 'High Storage',
                               str_detect(site_id, 'LS') ~ 'Low Storage'),
         .after = site_id) |> 
  group_by(hillslope) |> 
  summarise(average_HR = mean(recharge, na.rm = TRUE),
            average_drawdown = mean(drawdown, na.rm = TRUE),
            min_vwc_actual = min(Smooth_VWC_Avg, na.rm = TRUE),
            min_vwc_predicted = min(predicted_vwc_noHR, na.rm = TRUE))

# Depth-resolved "additional water" estimate based on end-of-season minima
Additional_water <- df_VWC_calcs |> 
  mutate(hillslope = case_when(str_detect(site_id, 'HS') ~ 'High Storage',
                               str_detect(site_id, 'LS') ~ 'Low Storage'),
         .after = site_id) |> 
  group_by(hillslope, depth) |> 
  summarise(min_vwc_actual = min(Smooth_VWC_Avg, na.rm = TRUE),
            min_vwc_predicted = min(predicted_vwc_noHR, na.rm = TRUE),
            percent_additional = min_vwc_actual/min_vwc_predicted * 100)


## ----------------------------------------- ##
#           Export HR Calculations ----
## ----------------------------------------- ##

write_csv(df_VWC_calcs, file.path(OD, 'HR_Calculations_June-Aug.csv'))

# Export Summaries
write_csv(HR_summary, file.path(SD, 'HR_hillslope_summary.csv'))
write_csv(Additional_water, file.path(SD, 'predicted_v_actual_vwc_end-of-season.csv'))

## ----------------------------------------- ##
##      Calculate Time in Optimal  ----
## ----------------------------------------- ##

# Series maximum VWC by site_id x depth (used to convert VWC to %WHC).
# See header Notes: this reconstructs a step (df_max_summary) the
# original script referenced but never defined, matching 03_whc.R's
# equivalent (unfiltered, full-year) computation. site_id/depth are
# cast to factor here to match calc_vwc_hr()'s output in df_VWC_calcs
# (which factors both internally), avoiding a factor-vs-numeric join
# mismatch below.
df_max_summary <- df_vwc |>
  group_by(site_id, depth) |>
  summarise(max_vwc = max(Smooth_VWC_Avg, na.rm = TRUE), .groups = 'drop')

df_max <- df_max_summary |> 
  filter(depth <= 30) |> 
  mutate(site_id = as_factor(site_id),
         depth = as_factor(depth))

# Convert VWC to WHC (% of max), for observed and predicted-noHR series
df_whc <- df_VWC_calcs |> 
  left_join(df_max) |> 
  mutate(whc_OG = (Smooth_VWC_Avg/max_vwc * 100),
         whc_predicted = (predicted_vwc_noHR/max_vwc * 100),
         hillslope = case_when(str_detect(site_id, 'HS') ~ 'High Storage',
                               str_detect(site_id, 'LS') ~ 'Low Storage'))

# Interpolate across gaps in WHC series (linear interpolation in time)
df_whc_interp <- df_whc |>
  arrange(site_id, depth, datetime) |>
  group_by(site_id, depth) |>
  mutate(whc_interp = na.approx(whc_OG, x = as.numeric(datetime),
                                na.rm = FALSE, maxgap = Inf),
         whc_predicted_interp = na.approx(whc_predicted, x = as.numeric(datetime),
                                          na.rm = FALSE, maxgap = Inf),
         hillslope = hillslope)

## ----------------------------------------- ##
#               Optimal WHC ----

# Time in "optimal" WHC range (50-70%); assumes 10-min timestep (x10 minutes)
whc_summary <- df_whc_interp |>
  mutate(optimal_OG = whc_interp >= 50 & whc_interp <= 70,
         optimal_predicted = whc_predicted_interp >= 50 & whc_predicted_interp <= 70) |> 
  group_by(site_id, depth, hillslope) |> 
  summarize(optimal_time_OG = sum(optimal_OG, na.rm = TRUE),
            optimal_time_pred = sum(optimal_predicted, na.rm = TRUE)) |> 
  mutate(optimal_OG_minutes = optimal_time_OG * 10,
         optimal_OG_hours = optimal_time_OG/6,
         optimal_OG_days = optimal_OG_hours/24,
         optimal_pred_minutes = optimal_time_pred * 10,
         optimal_pred_hours = optimal_time_pred/6,
         optimal_pred_days = optimal_pred_hours/24) |> 
  ungroup()


# Hillslope x depth summaries of optimal time (observed vs predicted-noHR)
whc_analysis <- whc_summary |> 
  mutate(hillslope = case_when(str_detect(site_id, 'LS') ~ 'Low Storage',
                               str_detect(site_id, 'HS') ~ 'High Storage'),
         .after = site_id) |> 
  group_by(hillslope, depth) |> 
  summarize(mean_opt_days_OG = mean(optimal_OG_days, na.rm = TRUE),
            mean_opt_days_pred = mean(optimal_pred_days, na.rm = TRUE),
            sd_days_opt_OG = sd(optimal_OG_days, na.rm = TRUE),
            sd_days_opt_pred = sd(optimal_pred_days, na.rm = TRUE),
            n_sites = n(),
            se_days_opt_OG = sd_days_opt_OG/sqrt(n_sites),
            se_days_opt_pred = sd_days_opt_pred/sqrt(n_sites),
            .groups = 'drop')

# Cross-hillslope ratios (HS / LS) for observed and predicted series
OG_opt_days_comparison <- whc_analysis |> 
  select(hillslope, depth, mean_opt_days_OG) |> 
  pivot_wider(names_from = hillslope, 
              values_from = mean_opt_days_OG) |> 
  mutate(HS_vs_LS = `High Storage`/`Low Storage`)

Pred_opt_days_comparison <- whc_analysis |> 
  select(hillslope, depth, mean_opt_days_pred) |> 
  pivot_wider(names_from = hillslope, 
              values_from = mean_opt_days_pred) |> 
  mutate(HS_vs_LS = `High Storage`/`Low Storage`)

# Observed vs predicted optimal-time deficit (predicted minus observed)
OGvPred_optimal <- whc_analysis |> 
  select(hillslope, depth, mean_opt_days_OG, mean_opt_days_pred) |> 
  mutate(Optimal_Deficit_noHR = mean_opt_days_pred - mean_opt_days_OG)

# Export optimal-time deficit table
write_csv(OGvPred_optimal, file.path(SD, 'Optimal_Time_Deficit_Pred.v.Actual.csv'))
