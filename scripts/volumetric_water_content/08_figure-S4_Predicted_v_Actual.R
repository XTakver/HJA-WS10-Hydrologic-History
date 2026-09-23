# ------------------------------------------------------------------------------
# Predicted vs. observed VWC (LS-1, 10 cm) -- Figure S4
# ------------------------------------------------------------------------------
# Purpose:
#   - Read HR calculation output and subset to LS-1 at 10 cm depth
#   - Restrict to a short July 2024 window used in the manuscript
#   - Plot predicted (no-HR) vs actual smoothed VWC time series
#   - Export a small, publication-style PDF panel for the paper
#
# Input:
#   - data/volumetric_water_content/processed/HR_Calculations_June-Aug.csv
#
# Output:
#   - output/figures/volumetric_water_content/Fig.S4-Actual-v-Predicted_VWC (CairoPDF)
#
# Created:  2026-01-27
# Author:   X. Takver
# ------------------------------------------------------------------------------

## ----------------------------------------- ##
#             Package Import ----
## ----------------------------------------- ##

# Data Analysis
librarian::shelf(tidyverse, zoo, here, glue)

# Data visualization
librarian::shelf(tinythemes, Cairo)

## ----------------------------------------- ##
##            Source functions ----
## ----------------------------------------- ##

# Project helpers:
#   - ensure_dir(): idempotent directory creation
#   - fonts.R / plot_theme.R: consistent manuscript styling
source(here('R', 'functions', 'paths_and_directories.R')) 
source(here('R', 'functions', 'fonts.R'))                  
source(here('R', 'functions', 'plot_theme.R'))       

## ----------------------------------------- ##
#             Set Plot Theme ----
## ----------------------------------------- ##

# Set global ggplot2 theme for consistent manuscript figure styling
set_plot_theme()

## ----------------------------------------- ##
#           Directory Creation ----
## ----------------------------------------- ##

# DD: processed VWC/HR outputs
# PD: figure output directory for this script
DD <- here('data', 'volumetric_water_content', 'processed')
PD <- here('output', 'figures', 'volumetric_water_content')

# Ensure output directory exists
ensure_dir(PD)

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

# HR calculations include:
#   - predicted_vwc_noHR (modeled drawdown without HR)
#   - Smooth_VWC_Avg     (observed/smoothed VWC)
df_hr <- read_csv(file.path(DD, 'HR_Calculations_June-Aug.csv'))

## ----------------------------------------- ##
#              Data Processing ----
## ----------------------------------------- ##

# Subset to the panel location (LS1, 10 cm) and parse timestamps for plotting
df_ls1 <- df_hr |> 
  filter(site_id == 'LS1',
         depth == 10) |> 
  mutate(datetime = as_datetime(datetime)) 

## ----------------------------------------- ##
#             Subset for paper ----
## ----------------------------------------- ##

# Restrict to the manuscript window (mid-July 2024)
df_subset <- df_ls1 |> 
  filter(datetime < '2024-07-19 23:59:00' & datetime > '2024-07-15 00:00:00')

# Plot predicted (no-HR) vs observed VWC; y-limits zoom to emphasize divergence
plt.subset <- 
  ggplot(df_subset) +
  geom_line(aes(x = datetime, y = predicted_vwc_noHR), 
            color = '#4B2E5A', linetype = 'longdash',
            linewidth = 0.35) +
  geom_line(aes(x = datetime, y = Smooth_VWC_Avg), 
            color = '#7E2E27',
            linewidth = 0.35) +
  coord_cartesian(ylim = c(15, 18)) +
  labs(y = 'Volumetric Water Content (%)',
       x = 'Date',
       subtitle = 'Predicted vs. Actual HR LS-1') +
  theme(plot.subtitle = element_text(margin = margin(b = 5)))

plt.subset

## ----------------------------------------- ##
#                Export Plot ----
## ----------------------------------------- ##

# Save as PDF with Cairo for consistent font/line rendering
CairoPDF(file.path(PD, 'Fig.S4-Actual-v-Predicted_VWC'), width = 3, height = 2.25)
print(plt.subset)
dev.off()


