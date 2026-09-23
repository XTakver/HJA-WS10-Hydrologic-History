# ------------------------------------------------------------------------------
# VWC and WHC time series by depth -- Figure S2A (WHC) and Figure S2B (VWC)
# ------------------------------------------------------------------------------
# Purpose:
#   - Read smoothed VWC, compute WHC (% of each sensor's max VWC)
#   - Plot depth-faceted time series for WHC and VWC, and export to PDF
#
# Input:
#   - data/volumetric_water_content/processed/smoothed_vwc.csv
#
# Outputs:
#   - Fig.S2A (WHC by depth)
#   - Fig.S2B (VWC by depth)
#
# Created:  2026-02-10
# Author:   X. Takver
# ------------------------------------------------------------------------------

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
#   - ensure_dir(): idempotent directory creation
#   - fonts.R / plot_theme.R: consistent manuscript styling
source(here('R', 'functions', 'paths_and_directories.R')) 
source(here('R', 'functions', 'fonts.R'))                  
source(here('R', 'functions', 'plot_theme.R'))             

## ----------------------------------------- ##
#           Directory Creation ----
## ----------------------------------------- ##

# DD: processed VWC data
# PD: figure output directory for this script
DD <- here('data', 'volumetric_water_content', 'processed')
PD <- here('output', 'figures', 'volumetric_water_content')

# Ensure output directory exists
ensure_dir(PD)

## ----------------------------------------- ##
#             Set Plot Theme ----
## ----------------------------------------- ##

# Set global ggplot2 theme for consistent manuscript figure styling
set_plot_theme()

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

# Read smoothed VWC time series (site_id already matches manuscript labels).
df_vwc <- read_csv(file.path(DD, 'smoothed_vwc.csv'))

## ----------------------------------------- ##
#              Compute WHC (%) ----
## ----------------------------------------- ##

# For each site_id x depth, compute the maximum observed smoothed VWC.
# This maximum is used to normalize VWC into WHC (% of max).
df_max <- df_vwc |> 
  group_by(site_id, depth) |> 
  summarize(max_vwc = max(Smooth_VWC_Avg, na.rm = TRUE))

# Join max values back to the full time series and compute WHC
# (WHC = Smooth_VWC_Avg / max_vwc * 100)
df_whc <- df_vwc |> 
  left_join(df_max) |> 
  mutate(whc = (Smooth_VWC_Avg/max_vwc * 100))


# Linear interpolation across time within each site_id x depth (draft choice).
# na.rm = FALSE preserves leading/trailing NAs; maxgap = Inf allows large gaps.
df_whc_interp <- df_whc |>
  arrange(site_id, depth, datetime) |>
  group_by(site_id, depth) |>
  mutate(whc_interp = na.approx(whc, x = as.numeric(datetime),
                                na.rm = FALSE, maxgap = Inf))

###               WHC ----

# Define full time range for the shaded reference band
x_rng <- range(df_whc_interp$datetime, na.rm = TRUE)

# WHC time series by depth (facets); shaded band marks 50-70% WHC reference range
plt.all.whc <- 
  ggplot(df_whc_interp) +
  geom_line(aes(x = datetime, y = whc, color = site_id),
            linewidth = 0.35) +
  facet_wrap(~depth, ncol = 1) +
  scale_color_manual(values = c('HS1' = '#759295',
                                'HS2' = '#A9D9DE',
                                'LS1' = '#7E2E27',
                                'LS2' = '#C55648')) +
  annotate('rect', 
           xmin = x_rng[1], xmax = x_rng[2], 
           ymin = 50, ymax = 70,
           fill = '#a3b88c', alpha = 0.33) +
  labs(subtitle = 'Water Holding Content: June-October 2024',
       x = 'Date',
       y = 'Water Holding Content (%)',
       color = 'Site ID') +
  theme(panel.spacing = unit(0.1, 'lines'),
        legend.position = 'bottom',
        legend.margin = margin(t = -12.5))

plt.all.whc


# Export WHC panel (Fig. S2A)
CairoPDF(file.path(PD, 'Fig.S2A'), width = 3.25, height = 7.5)
print(plt.all.whc)
dev.off()

###               VWC ----

# VWC time series by depth (facets); uses smoothed VWC directly
plt.all.vwc <- 
  ggplot(df_vwc) +
  geom_line(aes(x = datetime, y = Smooth_VWC_Avg, color = site_id),
            linewidth = 0.35) +
  facet_wrap(~depth, ncol = 1) +
  scale_color_manual(values = c('HS1' = '#759295',
                                'HS2' = '#A9D9DE',
                                'LS1' = '#7E2E27',
                                'LS2' = '#C55648')) +
  labs(subtitle = 'Volumetric Content: June-October 2024',
       x = 'Date',
       y = 'Volumetric Water Content (%)',
       color = 'Site ID') +
  theme(panel.spacing = unit(0.1, 'lines'),
        legend.position = 'bottom',
        legend.margin = margin(t = -12.5))

plt.all.vwc


# Export VWC panel (Fig. S2B)
CairoPDF(file.path(PD, 'Fig.S2B'), width = 3.25, height = 7.5)
print(plt.all.vwc)
dev.off()
