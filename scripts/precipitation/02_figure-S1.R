# ------------------------------------------------------------ #
# 02_figure-S1.R
#
# Purpose
#   Plot monthly PRIMET precipitation totals by water year
#   (2022-2025). Produces manuscript Figure S1.
#
# Inputs
#   data/precipitation/processed/PRIMET_monthly_precip_WY2022_2025.csv
#
# Outputs
#   output/figures/precipitation/Fig.S1.pdf
#
# Created:  2026-06-05
# Author:   X. Takver
# ------------------------------------------------------------ #

## ----------------------------------------- ##
#             Import Packages ----
## ----------------------------------------- ##

librarian::shelf(tidyverse, here, Cairo)

## ----------------------------------------- ##
#             Source Functions ----
## ----------------------------------------- ##

source(here('R', 'functions', 'paths_and_directories.R'))  # ensure_dir()
source(here('R', 'functions', 'fonts.R'))                  # e.g., register_fonts()
source(here('R', 'functions', 'plot_theme.R'))              # e.g., set_plot_theme()

## ----------------------------------------- ##
#          Directory Creation ----
## ----------------------------------------- ##

DD <- here('data', 'precipitation', 'processed')
PD <- here('output', 'figures', 'precipitation')

ensure_dir(PD)

## ----------------------------------------- ##
#             Set Plot Theme ----
## ----------------------------------------- ##

# Set global ggplot2 theme for consistent manuscript figure styling
set_plot_theme()

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

precip_22_to_25 <- read_csv(file.path(DD, 'PRIMET_monthly_precip_WY2022_2025.csv')) |>
  mutate(
    water_year = as_factor(water_year),
    month_label = fct_relevel(
      month_label,
      'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar',
      'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep'
    )
  )

## ----------------------------------------- ##
#                   Plot ----
## ----------------------------------------- ##

# Fixed y-axis breaks
n_intervals <- 5
monthly_ymax <- 600
monthly_breaks <- seq(0, monthly_ymax, length.out = n_intervals + 1)

plt.monthly <-
  ggplot(precip_22_to_25, aes(x = month_label, y = total_monthly, fill = water_year)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c('#9BC7D9', '#5F9AA8', '#2F6F73', '#123C45')) +
  scale_y_continuous(limits = c(0, monthly_ymax),
                     breaks = monthly_breaks,
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL,
       y = 'Monthly precipitation (mm)',
       fill = 'Water year') +
  theme(plot.margin = margin(b = 0, t = 2.5, l = 7.5, r = 0))

plt.monthly

## ----------------------------------------- ##
#               Export Plot ----
## ----------------------------------------- ##

CairoPDF(file.path(PD, 'Fig.S1.pdf'), width = 6.5, height = 3)
print(plt.monthly)
dev.off()
