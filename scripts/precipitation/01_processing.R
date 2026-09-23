# ------------------------------------------------------------ #
# 01_processing.R
#
# Purpose
#   Import and tidy PRIMET (primary meteorological station)
#   precipitation records from the H.J. Andrews LTER, and compute
#   monthly precipitation totals by water year (2022-2025) for
#   Figure S1, plus a historical water-year total comparison
#   supporting the precipitation summary in Supplement S1.
#
# Inputs
#   primet_precip.csv (public repository export; MATLAB-derived)
#
# Outputs (written to OD)
#   - PRIMET_monthly_precip_WY2022_2025.csv          (feeds 03_figure-S1.R)
#   - PRIMET_WaterYear_Totals_HistoricalComparison.csv
#
# Notes
#   - Time is converted from MATLAB serial datenum to America/Los_Angeles.
#   - Water year = October through September (e.g., WY2023 = Oct 2022 - Sep 2023).
#   - The historical comparison excludes water years 2013 and 2026
#     (partial years at the start/end of the record) and compares
#     2023/2024 against the mean/SD of all remaining complete water years.
#
# Created:  2026-06-05
# Author:   X. Takver
# ------------------------------------------------------------ #

## ----------------------------------------- ##
#             Import Packages ----
## ----------------------------------------- ##

librarian::shelf(tidyverse, here, janitor)

## ----------------------------------------- ##
#             Source Functions ----
## ----------------------------------------- ##

source(here('R', 'functions', 'paths_and_directories.R'))  # ensure_dir()

## ----------------------------------------- ##
#          Directory Creation ----
## ----------------------------------------- ##

DD <- here('data', 'precipitation', 'raw')
OD <- here('data', 'precipitation', 'processed')

ensure_dir(OD)

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

# Read precip. Skip lines 1, 3, and 5; assign row 2 as the header
primet1 <- read_csv(file.path(DD, 'primet_precip.csv'), skip = 1) |>
  filter(!row_number() %in% c(2, 3)) |>
  row_to_names(row_number = 1) |>
  clean_names()

# Convert from MATLAB time to local time
tidy_precip1 <- primet1 |>
  mutate(date = as.numeric(date),
         datetime = as_datetime(round((date - 719529) * 86400), tz = 'UTC'),
         datetime = force_tz(datetime, 'America/Los_Angeles'),
         .after = date) |>
  select(site, datetime, precip_tot_100_0_01:swe_check)

# Prep for plotting: keep only unflagged precipitation totals
precip_plotting <- tidy_precip1 |>
  mutate(month = month(datetime),
         year = year(datetime),
         .after = datetime) |>
  mutate(precip_tot_100_0_01 = as.numeric(precip_tot_100_0_01)) |>
  filter(is.na(flag_precip_tot_100_0_01))

## ----------------------------------------- ##
#          Monthly Totals by Water Year ----
## ----------------------------------------- ##

# Water year = Oct-Sep; water_month reindexes Oct=1 ... Sep=12
precip_monthly_wy <- precip_plotting |>
  mutate(
    water_year = year(datetime) + if_else(month(datetime) >= 10, 1, 0),
    water_month = if_else(month(datetime) >= 10, month(datetime) - 9, month(datetime) + 3),
    month_label = month(datetime, label = TRUE, abbr = TRUE),
    month_label = fct_relevel(
      month_label,
      'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar',
      'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep'
    )
  ) |>
  group_by(water_year, water_month, month_label) |>
  summarise(
    total_monthly = sum(precip_tot_100_0_01, na.rm = TRUE),
    .groups = 'drop'
  )

# Figure S1: monthly totals for water years 2022-2025
precip_22_to_25 <- precip_monthly_wy |>
  filter(water_year >= 2022, water_year <= 2025) |>
  mutate(water_year = as_factor(water_year))

write_csv(precip_22_to_25, file.path(OD, 'PRIMET_monthly_precip_WY2022_2025.csv'))

## ----------------------------------------- ##
#      Historical Water-Year Comparison ----
## ----------------------------------------- ##

# Full-record water-year totals, excluding partial years at the
# start/end of the record (2013, 2026)
wy_totals <- precip_monthly_wy |>
  group_by(water_year) |>
  summarise(
    wy_total = sum(total_monthly, na.rm = TRUE),
    .groups = 'drop'
  ) |>
  filter(
    !water_year %in% c(2013, 2026),
    !is.na(water_year)
  )

target_years <- c(2023, 2024)

historical_wy <- wy_totals |>
  filter(!water_year %in% target_years)

# Compare 2023/2024 totals against the historical mean/SD, expressed
# as an anomaly, z-score, and percentile
target_wy_comparison <- wy_totals |>
  filter(water_year %in% target_years) |>
  mutate(
    historical_mean = mean(historical_wy$wy_total, na.rm = TRUE),
    historical_sd = sd(historical_wy$wy_total, na.rm = TRUE),
    anomaly = wy_total - historical_mean,
    percent_anomaly = 100 * anomaly / historical_mean,
    z_score = anomaly / historical_sd,
    percentile = map_dbl(
      wy_total,
      ~ mean(historical_wy$wy_total <= .x, na.rm = TRUE) * 100
    ),
    interpretation = case_when(
      percentile <= 10 ~ 'very dry',
      percentile <= 25 ~ 'dry',
      percentile < 75  ~ 'near normal',
      percentile < 90  ~ 'wet',
      TRUE             ~ 'very wet'
    )
  )

write_csv(target_wy_comparison, file.path(OD, 'PRIMET_WaterYear_Totals_HistoricalComparison.csv'))
