# ------------------------------------------------------------------------------
# 08_figures-S12_and_S13_MAOM.R
#
# MAOM unique formulas -- seasonal compound-class shifts and NOSC distributions
# ------------------------------------------------------------------------------
# Purpose:
#   - Delta compound-class composition (December - June) from sample-average
#     metrics (Fig. S12)
#   - NOSC kernel densities for detected (unique) formulas, PA-weighted,
#     split by hillslope (Fig. S13, HS = panel A, LS = panel B)
#
# Inputs (FTICR pipeline):
#   - unique_formulas:
#       * 2023_FTICR-MS_mineral-associated_om_Presence-Absence_Unique.csv
#       * 2023_FTICR-MS_mineral-associated_om_SampleAvg_Unique.csv
#   - replicate_aggregation (formula-level properties incl. NOSC):
#       * 2023_FTICR-MS_mineral-associated_om_ChemProperties_aggregated.csv
#   (all written by 03_replicate_aggregation_MAOM.R / 04_unique_formula_subsetting_MAOM.R)
#
# Outputs (written to PD):
#   - Fig.S12_MAOM_CompoundChange.pdf
#   - Fig.S13A_MAOM_NOSC_HS.pdf / Fig.S13B_MAOM_NOSC_LS.pdf
#
# Notes:
#   - Some MAOM strata are missing (June for HS 20-30cm / LS 0-10cm), so
#     seasonal class-change deltas can only be calculated for HS 0-10cm
#     and LS 20-30cm; rows with an undefined (NA) delta are dropped
#     rather than plotted as zero. See manuscript Results 3.3.2.
#
# Created:  2026-02-20
# Author:   X. Takver
# ------------------------------------------------------------------------------

## ----------------------------------------- ##
#             Package Import ----
## ----------------------------------------- ##

# Wrangling + plotting
librarian::shelf(tidyverse, here, glue)

# Vector-safe PDF export
librarian::shelf(Cairo, tinythemes)

## ----------------------------------------- ##
#             Source Functions ----
## ----------------------------------------- ##

# Project-wide helpers:
#   - ensure_dir(): idempotent directory creation
#   - fonts.R / plot_theme.R: consistent manuscript styling
source(here('R', 'functions', 'paths_and_directories.R'))
source(here('R', 'functions', 'fonts.R'))
source(here('R', 'functions', 'plot_theme.R'))

## ----------------------------------------- ##
#           Directory Navigation ----
## ----------------------------------------- ##

# Unique formula products
DD_unique <- here('data', 'mass_spectrometry', 'processed', 'mineral-associated_om', 'unique_formulas')

# Formula-level properties (NOSC, etc.) used to annotate PA detections
DD_props <- here('data', 'mass_spectrometry', 'processed', 'mineral-associated_om', 'replicate_aggregation')

# Figure outputs
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

fp_pa   <- file.path(DD_unique, '2023_FTICR-MS_mineral-associated_om_Presence-Absence_Unique.csv')
fp_avg  <- file.path(DD_unique, '2023_FTICR-MS_mineral-associated_om_SampleAvg_Unique.csv')
fp_prop <- file.path(DD_props,  '2023_FTICR-MS_mineral-associated_om_ChemProperties_aggregated.csv')

stopifnot(file.exists(fp_pa), file.exists(fp_avg), file.exists(fp_prop))

maom_pa_unique  <- read_csv(fp_pa)
maom_avg_unique <- read_csv(fp_avg)
maom_props      <- read_csv(fp_prop)

## ----------------------------------------- ##
#             Helper: parse IDs ----
## ----------------------------------------- ##

# Sample ID format produced by 03/04: <HS|LS><Pit>_<June|December>_<Top|Btm>
# Example: HS1_June_Top, LS1_December_Btm
parse_sample_id <- function(df, sample_col = 'Sample.ID') {
  df |>
    separate(all_of(sample_col), into = c('Site.Raw', 'Sampling.Event', 'Core.Raw'), sep = '_', remove = FALSE) |>
    mutate(
      Site.ID = case_when(
        str_detect(Site.Raw, 'HS') ~ 'HS',
        str_detect(Site.Raw, 'LS') ~ 'LS',
        TRUE ~ NA_character_
      ),
      Core.Section = case_when(
        Core.Raw == 'Top' ~ '0-10cm',
        Core.Raw == 'Btm' ~ '20-30cm',
        TRUE ~ NA_character_
      ),
      Sampling.Event = fct_relevel(as_factor(Sampling.Event), 'June', 'December'),
      Core.Section = fct_relevel(as_factor(Core.Section), '0-10cm', '20-30cm'),
      Site.ID = as_factor(Site.ID)
    ) |>
    select(-Site.Raw, -Core.Raw)
}

## ----------------------------------------- ##
#                   NOSC ----
## ----------------------------------------- ##

# Join unique detections to per-formula NOSC, then expand to long form by sample.
# Density is PA-weighted: effectively the NOSC distribution of detected formulas.
maom_nosc <- maom_pa_unique |>
  left_join(maom_props |> select(molecular_formula, NOSC), by = 'molecular_formula') |>
  pivot_longer(
    cols = -c(molecular_formula, NOSC),
    names_to = 'Sample.ID',
    values_to = 'PA'
  ) |>
  select(Sample.ID, PA, NOSC) |>
  parse_sample_id(sample_col = 'Sample.ID')

# Sample-average NOSC for vertical reference lines (computed by your pipeline)
maom_nosc_means <- maom_avg_unique |>
  rename(Sample.ID = any_of('Sample')) |>
  select(Sample.ID, NOSC) |>
  parse_sample_id(sample_col = 'Sample.ID')

## ----------------------------------------- ##
##               NOSC Plots ----

# Publication colors consistent with your hillslope palette
cols_hs <- c('0-10cm' = '#A9D9DE', '20-30cm' = '#759295')
cols_ls <- c('0-10cm' = '#C55648', '20-30cm' = '#7E2E27')

# High Storage (HS)
plt.nosc.maom.hs <-
  ggplot(maom_nosc |> filter(Site.ID == 'HS')) +
  geom_density(
    aes(x = NOSC, weight = PA, fill = Core.Section),
    color = NA,
    alpha = 0.6
  ) +
  geom_vline(
    data = maom_nosc_means |> filter(Site.ID == 'HS'),
    aes(xintercept = NOSC, color = Core.Section),
    linewidth = 0.35
  ) +
  scale_fill_manual(values = cols_hs) +
  scale_color_manual(values = cols_hs, breaks = c('0-10cm', '20-30cm')) +
  coord_cartesian(ylim = c(0, 1)) +
  facet_wrap(~Sampling.Event, ncol = 1) +
  labs(
    title = 'MAOM',
    subtitle = 'High Storage',
    y = 'Density',
    fill = 'Core Section',
    color = 'Mean'
  )

# Low Storage (LS)
plt.nosc.maom.ls <-
  ggplot(maom_nosc |> filter(Site.ID == 'LS')) +
  geom_density(
    aes(x = NOSC, weight = PA, fill = Core.Section),
    color = NA,
    alpha = 0.6
  ) +
  geom_vline(
    data = maom_nosc_means |> filter(Site.ID == 'LS'),
    aes(xintercept = NOSC, color = Core.Section),
    linewidth = 0.35
  ) +
  scale_fill_manual(values = cols_ls) +
  scale_color_manual(values = cols_ls, breaks = c('0-10cm', '20-30cm')) +
  coord_cartesian(ylim = c(0, 1)) +
  facet_wrap(~Sampling.Event, ncol = 1) +
  labs(
    title = 'MAOM',
    subtitle = 'Low Storage',
    y = 'Density',
    fill = 'Core Section',
    color = 'Mean'
  )

## ----------------------------------------- ##
#         Delta Compound Classes ----
## ----------------------------------------- ##

# Parse SampleAvg identifiers to Site/Event/Core and compute delta % by class (Dec - June).
# Missing MAOM strata are dropped (delta NA) -- see Notes above.
maom_properties <- maom_avg_unique |>
  rename(Sample.ID = any_of('Sample')) |>
  parse_sample_id(sample_col = 'Sample.ID') |>
  mutate(Sample.Label = glue('{Site.ID} {Core.Section}'))

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

maom_compound_changes <- maom_properties |>
  select(Site.ID, Sampling.Event, Core.Section, Sample.Label, any_of(class_cols)) |>
  pivot_longer(
    cols = any_of(class_cols),
    names_to = 'Class',
    values_to = 'Percent'
  ) |>
  mutate(Class = gsub('_percent', '', Class),
         # Facet order matches the manuscript figure: complex classes
         # (upper left) to simple classes (lower right), filled by
         # facet_wrap() row-by-row (top row first 4 levels, bottom row
         # last 4), not the default alphabetical order.
         Class = factor(Class, levels = c('CondHydrocarbon', 'Tannin', 'Lignin', 'UnsatHydrocarbon',
                                          'Lipid', 'AminoSugar', 'Protein', 'Carbohydrate'))) |>
  pivot_wider(names_from = Sampling.Event, values_from = Percent) |>
  mutate(
    delta_percent = December - June,
    direction = case_when(
      delta_percent > 0 ~ 'Increase',
      delta_percent < 0 ~ 'Decrease',
      delta_percent == 0 ~ 'Neutral'
    )
  ) |>
  filter(!is.na(delta_percent))

plt.maom.deltaCompound <-
  ggplot(maom_compound_changes) +
  geom_col(aes(x = Sample.Label, y = delta_percent, fill = Sample.Label)) +
  facet_wrap(~Class, nrow = 2) +
  scale_fill_manual(values = c(
    'HS 0-10cm' = '#A9D9DE',
    'HS 20-30cm' = '#759295',
    'LS 0-10cm' = '#C55648',
    'LS 20-30cm' = '#7E2E27'
  )) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = '#848482') +
  labs(
    subtitle = 'MAOM',
    x = 'Site',
    y = expression(Delta~'Composition (%)'),
    fill = 'Site'
  ) +
  theme(
    legend.position = 'none',
    strip.text.x = element_text(margin = margin(b = 2)),
    panel.spacing.x = unit(0.5, 'lines'),
    panel.spacing.y = unit(0, 'lines'),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.subtitle = element_text(margin = margin(b = 2))
  )

## ----------------------------------------- ##
#                 Export ----
## ----------------------------------------- ##

CairoPDF(file.path(PD, 'Fig.S13A_MAOM_NOSC_HS.pdf'), width = 2.75, height = 3.25)
print(plt.nosc.maom.hs)
dev.off()

CairoPDF(file.path(PD, 'Fig.S13B_MAOM_NOSC_LS.pdf'), width = 2.75, height = 3.25)
print(plt.nosc.maom.ls)
dev.off()

CairoPDF(file.path(PD, 'Fig.S12_MAOM_CompoundChange.pdf'), width = 6.5, height = 3.25)
print(plt.maom.deltaCompound)
dev.off()

## ----------------------------------------- ##
#              Cleanup ----
## ----------------------------------------- ##

rm(maom_pa_unique, maom_avg_unique, maom_props,
   maom_nosc, maom_nosc_means, maom_properties, maom_compound_changes,
   plt.nosc.maom.hs, plt.nosc.maom.ls, plt.maom.deltaCompound)
