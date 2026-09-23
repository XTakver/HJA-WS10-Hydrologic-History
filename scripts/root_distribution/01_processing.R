# -----------------------------------------------------------------------------#
# Purpose
#   Import ImageJ root abundance scoring data, calculate percent root abundance
#   by depth, and export a tidy summary for downstream plotting/statistics.
#
# Input (from DD)
#   - 2023_WS10-HJA_ImageJ_Root_Abundance.csv
#
# Output (to OD)
#   - 2023_Root_Abundance_Percent_Processed.csv
#       Long-format table of mean percent root abundance by depth and hillslope,
#       split by root type (Fine, Coarse, Total).
#
# Notes
#   - Percent abundance computed from root counts standardized by width (see script).
#
# Created: 01/26/2026
# Author: X. Takver
# -----------------------------------------------------------------------------#

## ----------------------------------------- ##
#             Package Import ----
## ----------------------------------------- ##

librarian::shelf(tidyverse, here)

## ----------------------------------------- ##
#             Source Functions ----
## ----------------------------------------- ##

# Helper for creating output directories (idempotent)
source(here('R', 'functions', 'paths_and_directories.R'))

## ----------------------------------------- ##
#             Directory Setup ----
## ----------------------------------------- ##

# DD: directory containing raw root abundance data
# OD: directory for processed outputs written by this script
DD <- here('data', 'root_distribution', 'raw') 
OD <- here('data', 'root_distribution', 'processed') 

# Ensure output directories exist (idempotent)
ensure_dir(OD)

## ----------------------------------------- ##
#                Import Data ----
## ----------------------------------------- ##

df_root <- read_csv(file.path(DD, 'SP03902_v1.csv'))

## ----------------------------------------- ##
#             Data Processing ----
## ----------------------------------------- ##

# Convert fractions to percent and create a depth label for summaries/plots
df_root_processed <- 
  df_root |> 
  mutate(percent_fine = fine_fraction * 100,
         percent_coarse = coarse_fraction * 100,
         root_percent = root_fraction * 100) 

# Summarize mean percent abundance by depth and hillslope
df_root_summary <- 
  df_root_processed |> 
  group_by(bottom_depth_cm, hillslope) |> 
  summarise(avg_fine = mean(percent_fine),
            avg_coarse = mean(percent_coarse),
            avg_root = mean(root_percent),
            .groups = 'drop')

# Convert to long format for plotting and consistent downstream joins
df_root_long <- 
  df_root_summary |> 
  pivot_longer(cols = c(avg_fine, avg_coarse, avg_root),
               names_to = 'root_type',
               values_to = 'root_abundance_percent') |> 
  mutate(root_type = as_factor(case_when(root_type == 'avg_fine' ~ 'Fine',
                                         root_type == 'avg_coarse' ~ 'Coarse',
                                         root_type == 'avg_root' ~ 'Both'))) |> 
  rename('depth_cm' = bottom_depth_cm)

## ----------------------------------------- ##
#               Export Data ----
## ----------------------------------------- ##

# Export long dataframe for plotting and analysis
write_csv(df_root_long, file.path(OD, '2023_Root_Abundance_Percent_Processed.csv'))
