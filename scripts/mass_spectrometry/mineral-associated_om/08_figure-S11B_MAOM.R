# ------------------------------------------------------------------------------
# 08_figure-S11B_MAOM.R
#
# MAOM compound-class composition -- normalized barplot (unique formulas)
# ------------------------------------------------------------------------------
# Purpose:
#   - Import MAOM sample-average composition metrics computed from unique formulas
#   - Parse sample metadata from standardized Sample identifiers
#   - Generate a normalized (relative) compound-class barplot
#
# Input (FTICR processing pipeline):
#   - 2023_FTICR-MS_mineral-associated_om_SampleAvg_Unique.csv
#   (written by 04_unique_formula_subsetting_MAOM.R)
#
# Output:
#   - Fig.S11B_MAOM_CompClass.pdf
#
# Created:  2026-01-20
# Author:   X. Takver
# ------------------------------------------------------------------------------

## ----------------------------------------- ##
#             Package Import ----
## ----------------------------------------- ##

# Data processing + plotting
librarian::shelf(tidyverse, here, glue)

# PDF export device
librarian::shelf(Cairo)

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

# DD: FTICR outputs (unique formulas)
# PD: figure output directory for this script
DD <- here('data', 'mass_spectrometry', 'processed', 'mineral-associated_om', 'unique_formulas')
PD <- here('output', 'figures', 'mass_spectrometry', 'mineral-associated_om')

ensure_dir(PD)

## ----------------------------------------- ##
#             Set Plot Theme ----
## ----------------------------------------- ##

# Apply project theme (fonts + ggplot defaults)
set_plot_theme()

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

maom_sampleProp <- read_csv(file.path(DD, '2023_FTICR-MS_mineral-associated_om_SampleAvg_Unique.csv'))

## ----------------------------------------- ##
#             Data Preparation ----
## ----------------------------------------- ##

# Parse Sample ID structure: <HS|LS><Pit>_<June|December>_<Top|Btm>
maom_properties <- maom_sampleProp |>
  separate(Sample, into = c('Site.ID', 'Sampling.Event', 'Core.Section'), remove = FALSE) |>
  mutate(
    Site.ID = case_when(
      str_detect(Site.ID, 'HS') ~ 'HS',
      str_detect(Site.ID, 'LS') ~ 'LS'
    ),
    Core.Section = case_when(
      Core.Section == 'Top' ~ '0-10cm',
      Core.Section == 'Btm' ~ '20-30cm'
    )
  ) |>
  mutate(across(Site.ID:Core.Section, as_factor))

# Van Krevelen compound-class percent columns only. 
class_cols <- c(
  'Lipid_percent',
  'UnsatHydrocarbon_percent',
  'Protein_percent',
  'Lignin_percent',
  'Carbohydrate_percent',
  'AminoSugar_percent',
  'Tannin_percent',
  'CondHydrocarbon_percent'
)

# Convert percent columns to long form for stacked bars
maom_barplot <- maom_properties |>
  select(Site.ID:Core.Section, any_of(class_cols)) |>
  pivot_longer(
    cols = any_of(class_cols),
    names_to = 'Class',
    values_to = 'Percent'
  ) |>
  mutate(
    Class = gsub('_percent', '', Class),
    Class = glue('{Class}-like')
  ) |>
  mutate(
    Sampling.Event = fct_relevel(Sampling.Event, 'June', 'December'),
    Core.Section = fct_relevel(Core.Section, '0-10cm', '20-30cm'),
    Sample.ID = glue('{Site.ID} {Core.Section}')
  )

## ----------------------------------------- ##
#                 Plot ----
## ----------------------------------------- ##

# Palette used in your original script
palette_opaque <- c('#E0EBF2', '#A1A997', '#EFB571', '#D79B80',
                    '#B69A8C', '#C3B59A', '#DAD7B3', '#BFC6B8')

# Normalized (relative composition; sums to 1 within each Sample.ID x Sampling.Event)
df_maom_normalized <- maom_barplot |>
  group_by(Sample.ID, Sampling.Event) |>
  mutate(
    total_pct = sum(Percent, na.rm = TRUE),
    norm_pct = if_else(total_pct > 0, 100 * Percent / total_pct, NA_real_),
    relative_comp = norm_pct / 100
  ) |>
  ungroup()

plt.maom.norm <-
  ggplot(df_maom_normalized, aes(x = Sample.ID, y = relative_comp, fill = Class)) +
  geom_col(position = 'stack') +
  facet_wrap(~Sampling.Event, ncol = 1) +
  scale_fill_manual(values = palette_opaque) +
  labs(
    y = '',
    subtitle = 'MAOM',
    x = 'Site ID'
  )

## ----------------------------------------- ##
#                 Export ----
## ----------------------------------------- ##

CairoPDF(file.path(PD, 'Fig.S11B_MAOM_CompClass.pdf'), width = 4, height = 6)
print(plt.maom.norm)
dev.off()

## ----------------------------------------- ##
#              Cleanup ----
## ----------------------------------------- ##

rm(maom_sampleProp, maom_properties, maom_barplot, df_maom_normalized,
   plt.maom.norm)
