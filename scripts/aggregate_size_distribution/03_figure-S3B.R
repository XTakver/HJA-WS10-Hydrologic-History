# -----------------------------------------------------------------------------#
# Aggregate Size Distribution (ASD): Depth Profile Plot
#
# Purpose: Import processed ASD geometric mean diameter summary statistics and
# generate a depth-profile plot (mean ± SE) by site. The figure is exported as a
# PDF for manuscript or supplemental use.
#
# Input:
#   - Data/aggregate_size_distribution/Processed/processed_asd_GeoMean_SummaryStats.csv
#
# Output:
#   - Plots/Soil_Properties/Fig.S3B_ASD_Depth_Profile.pdf
#
# Created: 2026-01-23
# Author: X. Takver
# -----------------------------------------------------------------------------#

## ----------------------------------------- ##
#             Import Packages ----
## ----------------------------------------- ##

# Data processing packages
librarian::shelf(tidyverse, here)

# Data Visualization packages
librarian::shelf(tinythemes, Cairo)

## ----------------------------------------- ##
#             Source Functions ----
## ----------------------------------------- ##

# Directory helpers (ensure_dir)
source(here('R', 'functions', 'paths_and_directories.R'))
source(here('R', 'functions', 'fonts.R'))
source(here('R', 'functions', 'plot_theme.R'))

## ----------------------------------------- ##
#          Directory Creation ----
## ----------------------------------------- ##

# DD: directory containing processed ASD datasets
# PD: directory for plots written by this script
DD <- here('data', 'aggregate_size_distribution', 'processed')
PD <- here('output', 'figures', 'soil_properties')

# Ensure output directories exist (idempotent)
ensure_dir(PD)

## ----------------------------------------- ##
#             Set Plot Theme ----
## ----------------------------------------- ##

# Set global ggplot2 theme for consistent manuscript figure styling
set_plot_theme()

## ----------------------------------------- ##
#                 Data Import ----
## ----------------------------------------- ##

# Geometric mean diameter summary:
#   - rename() standardizes the mean diameter column name for plotting scripts
#   - arrange() orders rows for consistent line drawing by site and depth
df_GeoMean_Sum <- read_csv(file.path(DD, 'processed_asd_GeoMean_SummaryStats.csv')) |> 
  rename(mean_diam = 'final_mean_diam') |> 
  arrange(site_id, mid_depth)

## ----------------------------------------- ##
#                 Plot Data ----
## ----------------------------------------- ##

## ----------------------------------------- ##
##                All Depths ----

# Depth profile of mean aggregate diameter (± SE) by site.
# y-axis is reversed so increasing depth plots downward.
plt.all <- 
ggplot(df_GeoMean_Sum, aes(x = mean_diam, y = mid_depth, color = site_id)) +
  geom_errorbar(aes(xmin = mean_diam - se_diam, xmax = mean_diam + se_diam),
    orientation = 'y',  width = 2.5, linewidth = 0.5, alpha = 0.5) +
  geom_point(size = 1.5) +
  geom_line(orientation = 'y', linewidth = 0.75, alpha = 0.5) +
  scale_color_manual(values = c('HS1' = '#759295',
                                'HS2' = '#A9D9DE',
                                'LS1' = '#7E2E27',
                                'LS2' = '#C55648')) +
  scale_y_reverse() +
  labs(subtitle = 'Aggregate Size Distribution',
       x = 'Mean Aggregate Diameter (mm)',
       y = 'Depth (cm)',
       color = '')

plt.all

###     Save plots     ----

# Export manuscript figure panel as PDF
CairoPDF(file.path(PD, 'Fig.S3B_ASD_Depth_Profile.pdf'), width = 4.5, height = 3.5)
print(plt.all)
dev.off()




