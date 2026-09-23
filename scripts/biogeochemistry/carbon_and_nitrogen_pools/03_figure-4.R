# ------------------------------------------------------------ #
# 03_figure-4.R
#
# Purpose
#   Plot WEOM C, MAOM C, and microbial biomass C (MBC) response
#   ratios from the processed WS10 MONet biogeochemistry dataset.
#   Produces manuscript Figure 4.
#
# Inputs
#   2023_Carbon_and_Nitrogen_pools_responses.csv
#
# Outputs
#   plots/biogeochemistry/*.pdf
#
#
# Created: 2026-06-16
# Author: X. Takver
# ------------------------------------------------------------ #

## ----------------------------------------- ##
#             Package Import ----
## ----------------------------------------- ##

librarian::shelf(tidyverse, here, cowplot, Cairo, tinythemes, glue)

## ----------------------------------------- ##
#             Source Functions ----
## ----------------------------------------- ##

# Helper for creating output directories (idempotent)
source(here('R', 'functions', 'paths_and_directories.R'))
source(here('R', 'functions', 'fonts.R'))
source(here('R', 'functions', 'plot_theme.R'))

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

df_OCarbon <- read_csv(file.path(DD, '2023_Carbon_and_Nitrogen_pools_responses.csv'))

## ----------------------------------------- ##
#               Data Preparation ----
## ----------------------------------------- ##

# Create combined Site x Depth category for plotting
# Create a manual y axis offset based on carbon pool for better plotting.
# depth_label re-expresses core_section (top/bottom) as the depth
# ranges (0-10cm/20-30cm) shown in the published figure legend.
df_plt <- df_OCarbon |> 
  mutate(depth_label = case_when(core_section == 'top' ~ '0-10cm',
                                 core_section == 'bottom' ~ '20-30cm'),
         .after = core_section) |> 
  mutate(site_depth = factor(glue('{hillslope} {depth_label}'),
                             levels = c('HS 0-10cm',
                                        'HS 20-30cm',
                                        'LS 0-10cm',
                                        'LS 20-30cm')),
         .after = depth_label) |> 
  mutate(carbon_pool = factor(carbon_pool,
                              levels = c('MAOM C', 'WEOM C', 'Microbial Biomass C')),
         y_base = as.numeric(carbon_pool),
         y_offset = case_when(core_section == 'top' ~ 0.075,
                              core_section == 'bottom' ~ -0.075),
         y_plot = y_base + y_offset)


## ----------------------------------------- ##
#                   Plots ----
## ----------------------------------------- ##

# Error bar end caps
cap_height <- 0.075

plt.log2fc <- 
  ggplot(df_plt) +
  geom_vline(xintercept = 0, linetype = 'dashed', 
             color = '#848482', linewidth = 0.4) +
  geom_segment(aes(x = log2_ci_low,xend = log2_ci_low,
                   y = y_plot - cap_height, yend = y_plot + cap_height, 
                   color = site_depth), linewidth = 0.2) +
  geom_segment(aes(x = log2_ci_high, xend = log2_ci_high,
                   y = y_plot - cap_height, yend = y_plot + cap_height,
                   color = site_depth), linewidth = 0.2) +
  geom_segment(aes(x = log2_ci_low, xend = log2_ci_high,
                   y = y_plot, yend = y_plot,
                   color = site_depth), linewidth = 0.3) +
  geom_point(aes(x = log2_fold_change, y = y_plot, 
                 color = site_depth),
             position = position_dodge(width = 0.25, reverse = TRUE),
             size = 2) +
  scale_y_continuous(breaks = seq_along(levels(df_plt$carbon_pool)),
                     labels = levels(df_plt$carbon_pool)) +
  scale_color_manual(values = c('HS 0-10cm' = '#A9D9DE', 
                                'HS 20-30cm' = '#759295', 
                                'LS 0-10cm' = '#C55648', 
                                'LS 20-30cm' = '#7E2E27')) +
  labs(y = '',
       x = 'Log2 Fold Change',
       color = 'Site ID') +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = 'bottom',
        legend.margin = margin(t = -5),
        axis.title.x = element_text(size = 6.5)) +
  guides(color = guide_legend(label.vjust = 0.5,
                              override.aes = list(size = 1.5)))

plt.log2fc

## ----------------------------------------- ##
#               Export Plots ----
## ----------------------------------------- ##

# Export manuscript figure as PDF
CairoPDF(file.path(PD, 'Fig.4_Ocarbon_pools_Log2FC'), width = 5, height = 3)
print(plt.log2fc)
dev.off()
