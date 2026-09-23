# -----------------------------------------------------------------------------#
# Purpose
#   Plot depth profiles of root abundance (fine vs coarse) by hillslope using the
#   processed root abundance dataset produced in the corresponding processing script.
#
# Inputs (read from DD)
#   - 2023_Root_Abundance_Percent_Processed.csv
#
# Outputs (written to PD)
#   - Fig.S2_Root_Abundance_with_Depth.pdf
#
# Created: 01/26/2026
# Author: X. Takver
# -----------------------------------------------------------------------------#


## ----------------------------------------- ##
#             Import Packages ----
## ----------------------------------------- ##

# Data processing and file paths
librarian::shelf(tidyverse, here)

# Plot styling and PDF export
librarian::shelf(tinythemes, Cairo)

## ----------------------------------------- ##
#             Source Functions ----
## ----------------------------------------- ##

# Helper for creating output directories (idempotent)
source(here('R', 'functions', 'paths_and_directories.R'))
source(here('R', 'functions', 'fonts.R'))
source(here('R', 'functions', 'plot_theme.R'))

## ----------------------------------------- ##
#             Directory Setup ----
## ----------------------------------------- ##

# DD: directory containing raw root abundance data
# PD: directory for rendered plots written by this script
DD <- here('data', 'root_distribution', 'processed') 
PD <- here('output', 'figures', 'vegetation') 

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

# Import processed root abundance (processed long format)
df_roots <- read_csv(file.path(DD, '2023_Root_Abundance_Percent_Processed.csv'))

## ----------------------------------------- ##
#               Plot Data ----
## ----------------------------------------- ##

# Set facet order to control panel placement in the final figure
df_roots <- df_roots |>
  filter(root_type != 'Both') |> 
  mutate(root_type = factor(root_type, levels = c('Fine', 'Coarse')))

# Render plot
plt.root.abundances <- 
  ggplot(df_roots) + 
  geom_point(aes(x = root_abundance_percent, y = depth_cm, color = hillslope), size = 0.9) +
  scale_y_reverse() +
  scale_color_manual(values = c('#7E2E27', '#759295')) +
  facet_wrap(~root_type) +
  coord_cartesian(xlim = c(0, 100)) +
  labs(x = 'Abundance (%)',
       y = 'Depth (cm)',
       color = 'Hillslope') +
  theme(panel.spacing = unit(1, 'lines'),
        plot.title = element_blank(),
        plot.subtitle = element_blank())

plt.root.abundances

## ----------------------------------------- ##
#               Export Plot(s) ----
## ----------------------------------------- ##

# Export as vector PDF for manuscripts/supplement
CairoPDF(file.path(PD, 'Fig.S2_Root_Abundance_with_Depth.pdf'), width = 5, height = 5)
print(plt.root.abundances)
dev.off()
