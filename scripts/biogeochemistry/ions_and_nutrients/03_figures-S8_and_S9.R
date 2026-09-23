# ------------------------------------------------------------ #
# 03_figures-S8_and_S9.R
#
# Purpose
#   Plot WS10 MONet ions (meq/100g) and micro/macronutrients (ppm),
#   plus base saturation, using processed outputs from processing script.
#   Produces manuscript Figure S8 (panels A-D) and Figure S9.
#
# Inputs (processed)
#   2023_Soil_Ion_meq_long.csv
#   2023_Soil_Ion_ppm_long.csv
#   2023_Soil_Ion_base_saturation.csv
#
# Outputs (figures)
#   output/figures/biogeochemistry/
#     - Fig.S8A.pdf (major cations/anions, fixed axes)
#     - Fig.S8B.pdf (micronutrients/trace elements, fixed axes)
#     - Fig.S8C.pdf (major cations/anions, faceted by ion)
#     - Fig.S8D.pdf (micronutrients/trace elements, faceted by ion)
#     - Fig.S9.pdf  (base saturation)
#
# Notes
#   - Uses shared font + theme helpers (fonts.R, plot_theme.R).
#   - This dataset is June-only, one sample per hillslope x core
#     section (see 01_processing_Ions_Nutrients.R), so all panels
#     show a single time point.
#
# Created:  2026-02-20
# Author:   X. Takver
# ------------------------------------------------------------ #

## ----------------------------------------- ##
#             Package Import ----
## ----------------------------------------- ##

librarian::shelf(tidyverse, here, Cairo, wesanderson, tinythemes, glue)

## ----------------------------------------- ##
#             Source Functions ----
## ----------------------------------------- ##

source(here('R', 'functions', 'paths_and_directories.R'))  # ensure_dir()
source(here('R', 'functions', 'fonts.R'))                  # e.g., register_fonts()
source(here('R', 'functions', 'plot_theme.R'))             # e.g., set_plot_theme()

## ----------------------------------------- ##
#             Directories ----
## ----------------------------------------- ##

DD <- here('data', 'biogeochemistry', 'processed')
PD <- here('output', 'figures', 'biogeochemistry')

ensure_dir(PD)

## ----------------------------------------- ##
#             Set Plot Theme ----
## ----------------------------------------- ##

# Set global ggplot2 theme for consistent manuscript figure styling
set_plot_theme()

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

df_meq <- read_csv(file.path(DD, '2023_Soil_Ion_meq_long.csv'), show_col_types = FALSE)
df_ppm <- read_csv(file.path(DD, '2023_Soil_Ion_ppm_long.csv'), show_col_types = FALSE)
df_bs  <- read_csv(file.path(DD, '2023_Soil_Ion_base_saturation.csv'),  show_col_types = FALSE)

# Keep site_depth ordering consistent across figures
site_levels <- c('HS 0-10cm', 'HS 20-30cm', 'LS 0-10cm', 'LS 20-30cm')

# site_depth combines hillslope with a depth label derived from
# core_section (Top/Bottom), matching the published figure legends
add_site_depth <- function(df) {
  df |>
    mutate(depth_label = case_when(core_section == 'Top' ~ '0-10cm',
                                   core_section == 'Bottom' ~ '20-30cm')) |>
    mutate(site_depth = fct_relevel(as_factor(glue('{hillslope} {depth_label}')), site_levels))
}

df_meq <- add_site_depth(df_meq)
df_ppm <- add_site_depth(df_ppm)
df_bs  <- add_site_depth(df_bs)

## ----------------------------------------- ##
#               Palettes ----
## ----------------------------------------- ##

palette_moonrise2   <- wes_palette('Moonrise2')
palette_chevalier1  <- wes_palette('Chevalier1')
my_palette <- c(palette_moonrise2, palette_chevalier1)

site_fill <- c(
  'HS 20-30cm' = '#759295',
  'HS 0-10cm'  = '#A9D9DE',
  'LS 20-30cm' = '#7E2E27',
  'LS 0-10cm'  = '#C55648'
)

## ----------------------------------------- ##
#                   Plots ----
## ----------------------------------------- ##

## ----------------------------------------- ##
#                  meq ----

df_meq.plt <- df_meq |>
  filter(!ion %in% c('total_bases', 'cec'))

# Fixed axes: all ions in one panel
plt.fixed.meq <-
  ggplot(df_meq.plt) +
  geom_col(aes(x = site_depth, y = meq, fill = ion), position = 'dodge') +
  scale_fill_manual(values = my_palette) +
  labs(
    title = 'Major Cations and Anions',
    x = 'Site ID',
    y = 'meq/100g',
    fill = 'Ion'
  )

# Free axes: one facet per ion, fill by site_depth
plt.free.meq <-
  ggplot(df_meq.plt) +
  geom_col(aes(x = site_depth, y = meq, fill = site_depth), position = 'dodge') +
  facet_wrap(~ion, scales = 'free', nrow = 1) +
  scale_fill_manual(values = site_fill) +
  labs(
    title = 'Major Cations and Anions',
    x = 'Site ID',
    y = 'meq/100g',
    fill = 'Site ID'
  ) +
  theme(axis.text.x = element_blank())

## ----------------------------------------- ##
#             Base saturation ----

plt.baseSat <-
  ggplot(df_bs) +
  geom_col(aes(x = site_depth, y = base_saturation, fill = site_depth)) +
  scale_fill_manual(values = site_fill) +
  labs(
    title = 'Base Saturation',
    x = 'Site ID',
    y = 'Base Saturation (%)',
    fill = 'Site ID'
  )

## ----------------------------------------- ##
#                  ppm ----

# Fixed axes: all elements in one panel
plt.fixed.ppm <-
  ggplot(df_ppm) +
  geom_col(aes(x = site_depth, y = ppm, fill = ion), position = 'dodge') +
  scale_fill_manual(values = my_palette) +
  labs(
    title = 'Micronutrients and Trace Elements',
    x = 'Site ID',
    y = 'ppm',
    fill = 'Element'
  )

# Free axes: one facet per element, fill by site_depth
plt.free.ppm <-
  ggplot(df_ppm) +
  geom_col(aes(x = site_depth, y = ppm, fill = site_depth), position = 'dodge') +
  facet_wrap(~ion, scales = 'free', nrow = 1) +
  scale_fill_manual(values = site_fill) +
  labs(
    title = 'Micronutrients and Trace Elements',
    x = 'Site ID',
    y = 'ppm',
    fill = 'Site ID'
  ) +
  theme(axis.text.x = element_blank())

## ----------------------------------------- ##
#                  Export ----
## ----------------------------------------- ##

# Major ions (meq)
CairoPDF(file.path(PD, 'Fig.S8A.pdf'), width = 4, height = 3.5)
print(plt.fixed.meq)
dev.off()

CairoPDF(file.path(PD, 'Fig.S8C.pdf'), width = 8, height = 3)
print(plt.free.meq)
dev.off()

# Trace elements (ppm)
CairoPDF(file.path(PD, 'Fig.S8B.pdf'), width = 4, height = 3.5)
print(plt.fixed.ppm)
dev.off()

CairoPDF(file.path(PD, 'Fig.S8D.pdf'), width = 12, height = 2)
print(plt.free.ppm)
dev.off()

# Base saturation
CairoPDF(file.path(PD, 'Fig.S9.pdf'), width = 4, height = 3.5)
print(plt.baseSat)
dev.off()

## ----------------------------------------- ##
#                 Cleanup ----
## ----------------------------------------- ##

rm(df_meq, df_ppm, df_bs, df_meq.plt,
   plt.fixed.meq, plt.free.meq, plt.fixed.ppm, plt.free.ppm, plt.baseSat,
   palette_moonrise2, palette_chevalier1, my_palette, site_fill, site_levels,
   add_site_depth, DD, PD)
