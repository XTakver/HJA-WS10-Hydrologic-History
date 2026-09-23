# -----------------------------------------------------------------------------#
# 07_figure-5B_MAOM.R
#
# Purpose
#   Generate the MAOM NMDS ordination figure for the manuscript (Fig. 5B;
#   paired with the WEOM ordination as Fig. 5A), including significant
#   (p < 0.05) molecular-property vectors from envfit.
#
# Inputs (read from DD)
#   1) 2023_FTICR-MS_NMDS_Scores_MAOM.csv
#        - NMDS site scores (samples) with parsed metadata fields
#   2) 2023_FTICR-MS_NMDS_SampleProp_Scores_MAOM.csv
#        - envfit vector coordinates (NMDS axes), p-values, and r2 for molecular properties
#   (both written by 05_rarefied_NMDS_MAOM.R)
#
# Output (written to PD)
#   - Fig.5B_FTICR-MS_NMDS_MAOM.pdf
#
# Notes
#   - NMDS is computed upstream in the NMDS processing script; this file only plots results.
#   - Vectors are filtered to p < 0.05 and ranked by r2 (effect size); the top 15 are displayed.
#   - Vector coordinates are scaled (divided) for readability on the plot.
#
# Created:  2026-01-29
# Author:   X. Takver
# -----------------------------------------------------------------------------#


## ----------------------------------------- ##
#             Import Packages ----
## ----------------------------------------- ##

# Data processing
librarian::shelf(tidyverse, here)

# Data Visualization
librarian::shelf(ggrepel, tinythemes, Cairo)

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
#          Directory Navigation ----
## ----------------------------------------- ##

# Makes quick navigation directory commands

# DD: processed NMDS outputs
# PD: manuscript figure output directory
DD <- here('data', 'mass_spectrometry', 'processed', 'mineral-associated_om', 'NMDS')
PD <- here('output', 'figures', 'mass_spectrometry', 'mineral-associated_om')

# Ensure output directories exist (idempotent)
ensure_dir(PD)


## ----------------------------------------- ##
#             Set Plot Theme ----
## ----------------------------------------- ##

# Set global ggplot2 theme for consistent manuscript figure styling
set_plot_theme()

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

# NMDS site scores (samples) and fitted molecular-property vectors
df_NMDS <- read_csv(file.path(DD, '2023_FTICR-MS_NMDS_Scores_MAOM.csv'))
df_sampleProp <- read_csv(file.path(DD, '2023_FTICR-MS_NMDS_SampleProp_Scores_MAOM.csv'))

## ----------------------------------------- ##
#                 MAOM Plots ----
## ----------------------------------------- ##

## ----------------------------------------- ##
##             NMDS Base Plot ----

# Base ordination: samples colored by site x depth and shaped by sampling event
plt.hcl <- ggplot(df_NMDS, 
                   aes(x = NMDS1, y = NMDS2, 
                       color = Site.ID, shape = Sampling.Event)) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = '#848482', linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = '#848482', linewidth = 0.4) +
  geom_point(size = 2.25) +
  scale_color_manual(values = c('HS 0-10cm' = '#A9D9DE', 
                                'HS 20-30cm' = '#759295', 
                                'LS 0-10cm' = '#C55648', 
                                'LS 20-30cm' = '#7E2E27')) +
  scale_shape_manual(values = c('June' = 16, 'December' = 17))  +
  guides(shape = guide_legend(override.aes = list(size = c(1.75, 2))),
         color = guide_legend(override.aes = list(shape = 16, size = 1.75))) +
  labs(subtitle = 'MAOM',
       color = 'Site ID',
       shape = 'Sampling Event')

plt.hcl

## ----------------------------------------- ##
##          Add Sample Properties ----

# Filter properties to p-value < 0.05 and max properties to 15

# Display the strongest fitted molecular-property vectors:
#   - filter to significant vectors (p < 0.05)
#   - rank by r2 (effect size) and keep the top 15 for legibility
hcl_sampleProp_filter <- df_sampleProp |> 
  filter(p.value < 0.05) |> 
  mutate(vector.length = sqrt(NMDS1^2 + NMDS2^2)) |> 
  arrange(desc(r2)) |> 
  head(15)

# Add vectors and labels (scaled for readability)
plt.hcl.sampleProp <- 
  plt.hcl + 
  geom_segment(data = hcl_sampleProp_filter, 
               aes(x = 0, y = 0, xend = NMDS1 / 1.85, yend = NMDS2 / 1.85),
               arrow = arrow(length = unit(0.1, 'cm')), color =  	'#926C00', alpha = 0.7,
               inherit.aes = FALSE, lwd = 0.25) +
  geom_text_repel(data = hcl_sampleProp_filter,
                  aes(x = NMDS1 /1.2, y = NMDS2 /1.2, label = mol_properties),
                  max.overlaps = 15,
                  force = 0.75,
                  point.padding = 0.05,
                  size = 2, color = '#926C00',
                  inherit.aes = FALSE,
                  family = 'sans') +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank())

plt.hcl.sampleProp

## ----------------------------------------- ##
#               Export Plots ----
## ----------------------------------------- ##

# Export manuscript figure as PDF
CairoPDF(file.path(PD, 'Fig.5B_FTICR-MS_NMDS_MAOM'), width = 5, height = 3.25)
print(plt.hcl.sampleProp)
dev.off()
