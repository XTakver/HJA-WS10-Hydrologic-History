# ------------------------------------------------------------------------------
# WHC Subset Highlighting differences in HR -- Figure 3B
# ------------------------------------------------------------------------------
# Purpose:
#   - Read WHC data and subset to LS and HS 20 cm
#   - Restrict to a short window in mid-July 2024 used in the manuscript
#   - Export a small, publication-style PDF panel for the paper
#
# Input:
#   - data/volumetric_water_content/processed/whc.csv
#
# Output:
#   - output/figures/volumetric_water_content/Fig.3B-Detailed_HR (CairoPDF)
#
# Created:  2026-06-26
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

df_whc <- read_csv(file.path(DD, 'whc.csv'))

## ----------------------------------------- ##
#             Data Preparation ----
## ----------------------------------------- ##

whc_subset <- df_whc |> 
  filter(depth == 20,
         datetime < ymd('2024-07-20') & datetime > ymd('2024-07-16')) |> 
  mutate(hillslope = case_when(str_detect(site_id, 'HS') ~ 'High Storage',
                               str_detect(site_id, 'LS') ~ 'Low Storage'),
         .after = site_id) |> 
  group_by(hillslope, datetime) |> 
  summarize(mean_whc = mean(whc_interp))

## ----------------------------------------- ##
#               Plot Subset ----
## ----------------------------------------- ##

plt.hr_zoom <- 
  ggplot(whc_subset) +
    geom_line(aes(x = datetime, y = mean_whc, color = hillslope),
              linewidth = 0.35) +
    scale_color_manual(values = c('High Storage' = '#759295',
                                  'Low Storage'  = '#7E2E27')) +
    labs(x = '', y = 'WHC (%)')

## ----------------------------------------- ##
#                 Export ----
## ----------------------------------------- ##

CairoPDF(file.path(PD, 'Fig.3B-HR_Subset'), width = 3, height = 2.25)
print(plt.hr_zoom)
dev.off()


