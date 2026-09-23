# -----------------------------------------------------------------------------#
# Purpose
#   Plot mean particle-size fraction (sand, silt, clay) by depth for each
#   site, faceted by particle-size class, to visualize texture depth
#   profiles across hillslopes.
#
# Inputs (read from DD)
#   - 2023_WS10-HJA_PSF_Processed-long.csv
#
# Outputs (written to PD)
#   - Fig.S16A_Texture_Depth_Profile.pdf
#
# Created: 01/23/2026
# Author: X. Takver
# -----------------------------------------------------------------------------#

## ----------------------------------------- ##
#             Import Packages ----
## ----------------------------------------- ##

# Data Processing
librarian::shelf(tidyverse, here)

# Data Visualization
librarian::shelf(tinythemes, Cairo, grid)

## ----------------------------------------- ##
#          Import Helper Functions ----
## ----------------------------------------- ##

# Directory helpers (ensure_dir)
source(here('R', 'functions', 'paths_and_directories.R'))

# Plot theme helper functions
source(here('R', 'functions', 'fonts.R'))
source(here('R', 'functions', 'plot_theme.R'))

## ----------------------------------------- ##
#           Directory Creation ----
## ----------------------------------------- ##

# Data Directories
DD <- here('data', 'soil_texture', 'processed')

# Plotting Directories
PD <- here('output', 'figures', 'soil_properties')

# Ensure plot directory exists
ensure_dir(PD)

## ----------------------------------------- ##
#             Set Plot Theme ----
## ----------------------------------------- ##

# Set global ggplot2 theme for consistent manuscript figure styling
set_plot_theme()

## ----------------------------------------- ##
#                 Data Import ----
## ----------------------------------------- ##

# Import processed texture data (long format: one row per site_id x depth x PS class)
df_texture <- read_csv(here(DD, '2023_WS10-HJA_PSF_Processed-long.csv')) |> 
  mutate(across(c(site_id, psize_class), as_factor)) |> 
  # Standardize PS class labels and control facet order (Sand, Silt, Clay)
  mutate(psize_class = fct_recode(psize_class,
                               Sand = 'sand', Silt = 'silt', Clay = 'clay')) |> 
  mutate(psize_class = fct_relevel(psize_class, 'Sand', 'Silt', 'Clay'))

## ----------------------------------------- ##
#               Depth Plots ----
## ----------------------------------------- ##  

# Plot mean particle size fraction by depth; horizontal error bars show +/- 1 SE
plt.depthProfile <- 
  ggplot(df_texture, aes(x = mean_psf, y = mid_depth_cm, color = site_id)) +
  geom_errorbarh(aes(xmin = mean_psf - se_psf,
                     xmax = mean_psf + se_psf),
                 height = 2, linewidth = 0.5) +
  geom_point(size = 1.5) +
  # Reverse y-axis so depth increases downward (soil profile convention)
  scale_y_reverse() +
  facet_wrap(~psize_class) +
  # Map site IDs to consistent color palette used across figures
  scale_color_manual(values = c('HS1' = '#759295',
                                'HS2' = '#A9D9DE',
                                'LS1' = '#7E2E27',
                                'LS2' = '#C55648')) +
  labs(title = 'Particle Size Fraction with Depth',
       x = 'Mean PSF (%)',
       y = 'Depth (cm)',
       color = 'Site ID')

# NOTE: geom_errorbarh's `height` aesthetic is internally translated to
# `width`; this produces an innocuous console message and can be ignored.
plt.depthProfile

## ----------------------------------------- ##
##           Save Depth Profiles ----

# Export depth profile figure
CairoPDF(here(PD, 'Fig.S16A_Texture_Depth_Profile.pdf'), width = 6.5, height = 3.5)
print(plt.depthProfile)
dev.off()


