# ------------------------------------------------------------ #
# 04_figure-S14B.R
#
# Purpose
#   Plot per-year PERMANOVA R^2 (hillslope effect) trends for the
#   full WS10 watershed vegetation dataset and the subset co-located
#   with our study sites, with fitted linear trends and confidence
#   intervals. Produces manuscript Figure S14, panel B.
#
# Inputs
#   output/statistics/vegetation/Veg_Permanova_byYear.csv
#   output/statistics/vegetation/coLoc_Permanova_byYear.csv
#   (both written by 02_statistics_Vegetation_Community.R)
#
# Outputs
#   output/figures/vegetation/Fig.S14B.pdf
#
# Created:  2025-11-24
# Author:   X. Takver
# ------------------------------------------------------------ #

## ----------------------------------------- ##
#             Package Import ----
## ----------------------------------------- ##

librarian::shelf(tidyverse, here, Cairo)

## ----------------------------------------- ##
#             Source Functions ----
## ----------------------------------------- ##

source(here('R', 'functions', 'paths_and_directories.R'))  # ensure_dir()
source(here('R', 'functions', 'fonts.R'))                  # e.g., register_fonts()
source(here('R', 'functions', 'plot_theme.R'))              # e.g., set_plot_theme()

## ----------------------------------------- ##
#          Directory Navigation ----
## ----------------------------------------- ##

# Data Directory Variables
DD <- here('output', 'statistics', 'vegetation')

# Plot directory variables
PD <- here('output', 'figures', 'vegetation')

ensure_dir(PD)


## ----------------------------------------- ##
#             Set Plot Theme ----
## ----------------------------------------- ##

# Set global ggplot2 theme for consistent manuscript figure styling
set_plot_theme()

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

veg_permanova <- read_csv(file.path(DD, 'Veg_Permanova_byYear.csv'))
coLoc_permanova <- read_csv(file.path(DD, 'coLoc_Permanova_byYear.csv'))

## ----------------------------------------- ##
#                 Plot ----
## ----------------------------------------- ##

# 1. Add dataset labels and combine
veg_combined <- veg_permanova |>
  mutate(dataset = 'Watershed')

coloc_combined <- coLoc_permanova |>
  mutate(dataset = 'Co-located')

all_permanova <- bind_rows(veg_combined, coloc_combined)

# 2. Plot both with a shared legend
plt.both <-
  ggplot(all_permanova, aes(x = year, y = R2, color = dataset)) +
  geom_point(size = 1) +
  geom_smooth(method = 'lm', se = TRUE, alpha = 0.3, linewidth = 0.5) +
  scale_color_manual(name   = 'Site type',
    values = c('Watershed'  = '#9B8357FF', 'Co-located' = '#BD582CFF')) +
  labs(x = 'Years Since Harvest', y = expression(R^2))

plt.both

## ----------------------------------------- ##
#               Export Plots ----
## ----------------------------------------- ##

CairoPDF(file.path(PD, 'Fig.S14B.pdf'), width = 5, height = 2)
print(plt.both)
dev.off()


