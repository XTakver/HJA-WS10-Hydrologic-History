# -----------------------------------------------------------------------------#
# Aggregate Size Distribution (ASD): Geometric Mean Diameter
#
# Description: Imports raw aggregate size distribution data and calculates the
# geometric mean aggregate diameter per replicate and summary statistics across
# replicates.
#
# Outputs:
#   - processed_asd_geometric_mean.csv: replicate-level geometric mean diameter
#   - processed_asd_GeoMean_SummaryStats.csv: grouped mean, SD, SE, and n
#
# Created: 2026-01-23
# Author: X. Takver
# -----------------------------------------------------------------------------#

## ----------------------------------------- ##
#             Import Packages ----
## ----------------------------------------- ##

# Data processing packages
librarian::shelf(tidyverse, here, glue)

## ----------------------------------------- ##
#             Source Functions ----
## ----------------------------------------- ##

# Directory helpers (ensure_dir)
source(here('R', 'functions', 'paths_and_directories.R'))

## ----------------------------------------- ##
#           Directory Navigation ----
## ----------------------------------------- ##

# DD: directory containing raw ASD data
# OD: directory for processed outputs written by this script
DD <- here('data', 'aggregate_size_distribution', 'raw')
OD <- here('data', 'aggregate_size_distribution', 'processed')

# Ensure output directories exist (idempotent)
ensure_dir(OD)

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

# Raw ASD data table:
#   - read_csv() imports replicate-level mass fractions by aggregate size class
#     (plus identifiers and depth metadata) used to calculate geometric means
df_rawASD <- read_csv(file.path(DD, 'SP03901_v1.csv'))

## ----------------------------------------- ##
#               Geometric Mean ----
## ----------------------------------------- ##

# Prepare analysis table:
#   - retain identifiers + mass fractions
#   - compute mid-depth for grouping and reporting
df_mass_frac <- df_rawASD |> 
  select(sample_id:instrument_rep, fraction_gt_4.76_mm:fraction_lt_0.21_mm, large_agg_mean) |> 
  mutate(mid_depth = (top_depth + bottom_depth)/2,
         .after = bottom_depth)

# Calculate geometric mean diameter:
#   - mean bin diameters are geometric means of adjacent sieve bounds
#   - mean_diam is a mass-fraction-weighted geometric mean across included bins
df_GeoMean <- df_mass_frac |> 
  mutate(mean_f4.76 = exp((log(large_agg_mean) + log(4.76))/2),
         mean_f2 = exp((log(4.76) + log(2))/2),
         mean_f1 = exp((log(2) + log(1))/2),
         mean_f0.21 = exp((log(1) + log(0.21))/2),
         sum_mass_frac_forMean = rowSums(across(fraction_gt_4.76_mm:fraction_0.21_to_1_mm)),
         sum_mass_frac_all = rowSums(across(fraction_gt_4.76_mm:fraction_lt_0.21_mm))) |> 
  mutate(mean_diam = exp(((fraction_gt_4.76_mm*log(mean_f4.76)) + (fraction_2_to_4.76_mm*log(mean_f2)) +
                           (fraction_1_to_2_mm*log(mean_f1)) + (fraction_0.21_to_1_mm*log(mean_f0.21)))/sum_mass_frac_forMean))

##     Calculate summary statistics ----

# Summarize replicate mean_diam within each hillslope/site/depth group.
# Replicates with missing mean_diam are excluded from summary calculations.
df_GeoMean_stats <- df_GeoMean |> 
  filter(!is.na(mean_diam)) |> 
  group_by(hillslope, pit_id, site_id, top_depth, bottom_depth, mid_depth) |> 
  summarise(final_mean_diam = mean(mean_diam, na.rm = TRUE),
            sd_diam = sd(mean_diam, na.rm = TRUE),
            count = n(),
            se_diam = sd_diam/sqrt(n()), .groups = 'drop')

## ----------------------------------------- ##
#                 Export DFs ----

# Export replicate-level and summarized geometric mean diameter tables
write_csv(df_GeoMean, file.path(OD, 'processed_asd_geometric_mean.csv'))
write_csv(df_GeoMean_stats, file.path(OD, 'processed_asd_GeoMean_SummaryStats.csv'))

