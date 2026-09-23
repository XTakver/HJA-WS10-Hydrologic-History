# ------------------------------------------------------------------------------
# 01_processing_VWC.R
#
# VWC smoothing (Sentek) | reads published VWC data, interpolates/smooths by
# site_id x depth, and exports a single processed table + a PDF of per-series
# diagnostic plots.
#
# Inputs
#   data/volumetric_water_content/raw/SP03904_v2.csv
#
# Outputs
#   data/volumetric_water_content/processed/smoothed_vwc.csv
#   output/figures/volumetric_water_content/smoothing_analysis/VWC_Smoothed_Plots.pdf
#
# Notes
#   - Missing readings are coded in the published file as the literal
#     sentinel -9999 (not NA); these are converted to NA immediately on
#     import, before any interpolation/smoothing, since na.approx()
#     would otherwise treat -9999 as a real (and wildly wrong) reading.
#   - The published file stacks three readings per timestamp
#     (value_type: avg/min/max); only 'avg' is used here and in all
#     downstream VWC analyses, matching the original single-reading-
#     per-timestamp design this pipeline was built around.
#   - Smoothing window (2024-06-15 to 2024-10-31) matches the drydown-
#     to-rewetting period reported in the manuscript.
#
# Created:  2026-02-05
# Author:   X. Takver
# ------------------------------------------------------------------------------

## ----------------------------------------- ##
##                 Packages ----
## ----------------------------------------- ##

# Data processing 
librarian::shelf(tidyverse, here, zoo, pracma, quantmod)

# Plot export
librarian::shelf(Cairo, tinythemes)

## ----------------------------------------- ##
##            Source functions ----
## ----------------------------------------- ##

# Project-wide helpers (consistent with manuscript scripts)
#   - ensure_dir(), font registration, set_plot_theme
source(here('R', 'functions', 'paths_and_directories.R')) 
source(here('R', 'functions', 'fonts.R'))                  
source(here('R', 'functions', 'plot_theme.R'))             

# VWC processing functions
source(here('R', 'functions', 'vwc_processing_functions.R'))

## ----------------------------------------- ##
##           Paths / directories ----
## ----------------------------------------- ##

# DD: directory containing raw vwc inputs
# OD: directory for processed outputs written by this script
DD <- here('data', 'volumetric_water_content', 'raw')
OD <- here('data', 'volumetric_water_content', 'processed')

# Plot output directory
PD <- here('output', 'figures', 'volumetric_water_content', 'smoothing_analysis')

# Ensure directories exist
ensure_dir(OD)
ensure_dir(PD)

## ----------------------------------------- ##
##              Plot theme ----
## ----------------------------------------- ##

# Apply a consistent theme only once (plots produced by interp_and_smooth inherit it)
set_plot_theme()

## ----------------------------------------- ##
##               Data import ----
## ----------------------------------------- ##

# Import Sentek data. Keep only average readings (value_type == 'avg'):
# the published file also includes 'min'/'max' readings per timestamp,
# which are not used here or in any downstream VWC analysis. Missing
# readings are coded as the literal value -9999 and are converted to
# NA before any interpolation/smoothing.


df_sntk <- read_csv(file.path(DD, 'SP03904_v2.csv')) %>%
  filter(value_type == 'avg') %>%
  mutate(vwc = na_if(vwc, -9999)) %>%
  mutate(datetime = as_datetime(datetime)) %>%
  select(site_id, depth, datetime, vwc)

## ----------------------------------------- ##
##         Smoothing: down to 90 cm ----
## ----------------------------------------- ##

# Smoothing across the entire date range gets rid of minor fluctuations (HR) in the data.

# list depths
depths <- df_sntk %>% 
  pull(depth) %>% 
  unique()

# list sites (one per hillslope x pit; see published site_id/hillslope/pit_id columns)
site_ids <- df_sntk %>% 
  pull(site_id) %>% 
  unique()

## ----------------------------------------- ##
##              Smoothing ----

# Store smoothed data and smoothed plots as lists
smoothed_data <- list()
smoothed_plots <- list()

for (site in site_ids) {
  for (depth in depths) {
    key <- paste(site, depth, sep = '_')
    result <- interp_and_smooth(df_sntk, site, '2024-06-15', '2024-10-31', depth)
    
    smoothed_data[[key]] <- result$data
    smoothed_plots[[key]] <- result$plot
  }
}

## ----------------------------------------- ##
##               Save Output ----

# Combine list output to a single table. 
combined_df <- bind_rows(smoothed_data) %>% 
  mutate(datetime = format(datetime, '%Y-%m-%d %H:%M:%S'))

write_csv(combined_df, file.path(OD, 'smoothed_vwc.csv'))

# Diagnostics: all series plots in a single PDF
CairoPDF(file.path(PD, 'VWC_Smoothed_Plots.pdf'), width = 12, height = 10)

for (plot in smoothed_plots) {
  print(plot)  
}

dev.off()
