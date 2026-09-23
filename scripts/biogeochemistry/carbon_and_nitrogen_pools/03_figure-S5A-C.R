# ------------------------------------------------------------ #
# 03_figure-S5.R
#
# Purpose
#   Plot WEOM C, MAOM C, and microbial biomass C (MBC) boxplots
#   from the processed WS10 MONet biogeochemistry dataset.
#   Produces manuscript Figure S5 (panels A-C).
#
# Inputs
#   2023_Carbon_and_Nitrogen_pools_processed.csv
#
# Outputs
#   output/figures/biogeochemistry/*.pdf
#
#
# Created: 2026-02-19
# Author: X. Takver
# ------------------------------------------------------------ #

## ----------------------------------------- ##
#             Package Import ----
## ----------------------------------------- ##

librarian::shelf(tidyverse, here, patchwork, cowplot, Cairo, tinythemes)

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

df_OCarbon <- read_csv(file.path(DD, '2023_Carbon_and_Nitrogen_pools_processed.csv'))

# depth_label re-expresses core_section (top/bottom) as the depth
# ranges (0-10cm/20-30cm) shown in the published figure legend
df_OCarbon <- df_OCarbon |>
  mutate(depth_label = case_when(core_section == 'top' ~ '0-10cm',
                                 core_section == 'bottom' ~ '20-30cm'),
         .after = core_section)

## ----------------------------------------- ##
#                   Plots ----
## ----------------------------------------- ##

## ----------------------------------------- ##
##                  MAOM ----

###              MAOM Carbon ----

# HS
plt.MAOM.HS.carbon <- ggplot(df_OCarbon |> filter(hillslope == 'HS'),
                             aes(x = hillslope, y = maom_c, fill = depth_label)) +
  geom_boxplot(aes(color = ifelse(sampling_month == 'June', '#515151', as.character(depth_label)),
                   alpha = ifelse(sampling_month == 'June', 1, 0)),
               linewidth = 0.25) +
  scale_fill_manual(values = c('0-10cm' = '#A9D9DE', '20-30cm' = '#759295')) +
  scale_color_manual(values = c('0-10cm' = '#A9D9DE', '20-30cm' = '#759295', '#515151' = '#515151'),
                     breaks = c('0-10cm', '20-30cm')) +
  scale_alpha_identity() +
  coord_cartesian(ylim = c(0, 3500)) +
  labs(y = '', x = '', title = 'MAOM Carbon', fill = 'Core Section') +
  guides(color = guide_legend(title = 'December'),
         fill = guide_legend(title = 'June')) +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'))

# LS
plt.MAOM.LS.carbon <- ggplot(df_OCarbon |> filter(hillslope == 'LS'),
                             aes(x = hillslope, y = maom_c, fill = depth_label)) +
  geom_boxplot(aes(color = ifelse(sampling_month == 'June', '#515151', as.character(depth_label)),
                   alpha = ifelse(sampling_month == 'June', 1, 0)),
               linewidth = 0.25) +
  scale_fill_manual(values = c('0-10cm' = '#C55648', '20-30cm' = '#7E2E27')) +
  scale_color_manual(values = c('0-10cm' = '#C55648', '20-30cm' = '#7E2E27', '#515151' = '#515151'),
                     breaks = c('0-10cm', '20-30cm')) +
  scale_alpha_identity() +
  coord_cartesian(ylim = c(0, 3500)) +
  labs(y = 'Carbon (mg/kg)', x = 'Site ID', title = '', fill = 'Core Section') +
  guides(color = guide_legend(title = 'December'),
         fill = guide_legend(title = 'June')) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), 'cm'))

###           Combine MAOM Plots ----

plt.MAOM <- (plt.MAOM.HS.carbon + plt.MAOM.LS.carbon) +
  plot_layout(guides = 'collect')

###           Export MAOM Plots ----

CairoPDF(file.path(PD, 'Fig.S5B.pdf'), width = 4, height = 3.5)
print(plt.MAOM)
dev.off()

## ----------------------------------------- ##
##                   WEOM ----

###              WEOM Carbon ----

# HS < 100 mg/kg
plt.WEOM.HS.carbon.btm <- ggplot(df_OCarbon |> filter(hillslope == 'HS'),
                                 aes(x = hillslope, y = weom_c, fill = depth_label)) +
  geom_boxplot(aes(color = ifelse(sampling_month == 'June', '#515151', as.character(depth_label)),
                   alpha = ifelse(sampling_month == 'June', 1, 0)),
               linewidth = 0.25) +
  scale_fill_manual(values = c('0-10cm' = '#A9D9DE', '20-30cm' = '#759295')) +
  scale_color_manual(values = c('0-10cm' = '#A9D9DE', '20-30cm' = '#759295', '#515151' = '#515151'),
                     breaks = c('0-10cm', '20-30cm')) +
  scale_alpha_identity() +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = '', y = '', fill = 'Core Section') +
  guides(color = guide_legend(title = 'December'),
         fill = guide_legend(title = 'June')) +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'))

# HS > 650 mg/kg
plt.WEOM.HS.carbon.top <- ggplot(df_OCarbon |> filter(hillslope == 'HS'),
                                 aes(x = hillslope, y = weom_c, fill = depth_label)) +
  geom_boxplot(aes(color = ifelse(sampling_month == 'June', '#515151', as.character(depth_label)),
                   alpha = ifelse(sampling_month == 'June', 1, 0)),
               linewidth = 0.25) +
  scale_fill_manual(values = c('0-10cm' = '#A9D9DE', '20-30cm' = '#759295')) +
  scale_color_manual(values = c('0-10cm' = '#A9D9DE', '20-30cm' = '#759295', '#515151' = '#515151'),
                     breaks = c('0-10cm', '20-30cm')) +
  scale_alpha_identity() +
  coord_cartesian(ylim = c(650, 900)) +
  labs(y = 'Carbon (mg/kg)', title = 'WEOM Carbon', fill = 'Core Section') +
  guides(color = guide_legend(title = 'December'),
         fill = guide_legend(title = 'June')) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), 'cm'))

# LS < 100 mg/kg
plt.WEOM.LS.carbon.btm <- ggplot(df_OCarbon |> filter(hillslope == 'LS'),
                                 aes(x = hillslope, y = weom_c, fill = depth_label)) +
  geom_boxplot(aes(color = ifelse(sampling_month == 'June', '#515151', as.character(depth_label)),
                   alpha = ifelse(sampling_month == 'June', 1, 0)),
               linewidth = 0.25) +
  scale_fill_manual(values = c('0-10cm' = '#C55648', '20-30cm' = '#7E2E27')) +
  scale_color_manual(values = c('0-10cm' = '#C55648', '20-30cm' = '#7E2E27', '#515151' = '#515151'),
                     breaks = c('0-10cm', '20-30cm')) +
  scale_alpha_identity() +
  coord_cartesian(ylim = c(0, 100)) +
  labs(y = 'Carbon (mg/kg)', x = '', title = '', fill = 'Core Section') +
  guides(color = guide_legend(title = 'December'),
         fill = guide_legend(title = 'June')) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), 'cm'))

# LS > 650 mg/kg
plt.WEOM.LS.carbon.top <- ggplot(df_OCarbon |> filter(hillslope == 'LS'),
                                 aes(x = hillslope, y = weom_c, fill = depth_label)) +
  geom_boxplot(aes(color = ifelse(sampling_month == 'June', '#515151', as.character(depth_label)),
                   alpha = ifelse(sampling_month == 'June', 1, 0)),
               linewidth = 0.25) +
  scale_fill_manual(values = c('0-10cm' = '#C55648', '20-30cm' = '#7E2E27')) +
  scale_color_manual(values = c('0-10cm' = '#C55648', '20-30cm' = '#7E2E27', '#515151' = '#515151'),
                     breaks = c('0-10cm', '20-30cm')) +
  scale_alpha_identity() +
  coord_cartesian(ylim = c(650, 900)) +
  labs(fill = 'Core Section', title = '') +
  guides(color = guide_legend(title = 'December'),
         fill = guide_legend(title = 'June')) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), 'cm'))

###            Combine WEOM Plots ----

plt.WEOM <-
  (plt.WEOM.HS.carbon.top + plt.WEOM.LS.carbon.top) /
  plot_spacer() /
  (plt.WEOM.HS.carbon.btm + plt.WEOM.LS.carbon.btm) +
  plot_layout(guides = 'collect', heights = c(1, 0.025, 2))

# Add axis break marker
plt.WEOM.final <-
  ggdraw(plt.WEOM) +
  draw_text('/', x = 0.1, y = 0.65, angle = 130, size = 14) +
  draw_text('/', x = 0.1, y = 0.63, angle = 130, size = 14)

###              Export WEOM Plots ----

CairoPDF(file.path(PD, 'Fig.S5A.pdf'), width = 4.5, height = 3.5)
print(plt.WEOM.final)
dev.off()

## ----------------------------------------- ##
##            Microbial Biomass ----

###        Microbial Biomass Carbon ----

# HS
plt.mbio.HS.carbon <- ggplot(df_OCarbon |> filter(hillslope == 'HS'),
                             aes(x = hillslope, y = mbc, fill = depth_label)) +
  geom_boxplot(aes(color = ifelse(sampling_month == 'June', '#515151', as.character(depth_label)),
                   alpha = ifelse(sampling_month == 'June', 1, 0)),
               linewidth = 0.25) +
  scale_fill_manual(values = c('0-10cm' = '#A9D9DE', '20-30cm' = '#759295')) +
  scale_color_manual(values = c('0-10cm' = '#A9D9DE', '20-30cm' = '#759295', '#515151' = '#515151'),
                     breaks = c('0-10cm', '20-30cm')) +
  scale_alpha_identity() +
  coord_cartesian(ylim = c(50, 1700)) +
  labs(y = 'Carbon (mg/kg)', x = '', title = 'Microbial Biomass Carbon', fill = 'Core Section') +
  guides(color = guide_legend(title = 'December'),
         fill = guide_legend(title = 'June')) +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'))

# LS
plt.mbio.LS.carbon <- ggplot(df_OCarbon |> filter(hillslope == 'LS'),
                             aes(x = hillslope, y = mbc, fill = depth_label)) +
  geom_boxplot(aes(color = ifelse(sampling_month == 'June', '#515151', as.character(depth_label)),
                   alpha = ifelse(sampling_month == 'June', 1, 0)),
               linewidth = 0.25) +
  scale_fill_manual(values = c('0-10cm' = '#C55648', '20-30cm' = '#7E2E27')) +
  scale_color_manual(values = c('0-10cm' = '#C55648', '20-30cm' = '#7E2E27', '#515151' = '#515151'),
                     breaks = c('0-10cm', '20-30cm')) +
  scale_alpha_identity() +
  coord_cartesian(ylim = c(50, 1700)) +
  labs(y = 'Carbon (mg/kg)', x = 'Site ID', title = '', fill = 'Core Section') +
  guides(color = guide_legend(title = 'December'),
         fill = guide_legend(title = 'June')) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), 'cm'))

###           Combine mBio Plots ----

plt.mbio <- (plt.mbio.HS.carbon + plt.mbio.LS.carbon) +
  plot_layout(guides = 'collect')

###           Export mBio Plots ----

CairoPDF(file.path(PD, 'Fig.S5C.pdf'), width = 4, height = 3.5)
print(plt.mbio)
dev.off()
