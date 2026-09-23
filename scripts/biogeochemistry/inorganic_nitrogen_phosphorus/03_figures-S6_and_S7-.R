# ------------------------------------------------------------ #
# 03_figures-S6_and_S7.R
#
# Purpose
#   Plot WS10 MONet extractable inorganic nitrogen (NH4-N, NO3-N),
#   bicarbonate-extractable phosphorus (P_bicarb), and derived ratios
#   (NH4:NO3 and total N:P) as publication-style paired (HS/LS) boxplots.
#   Produces manuscript Figure S6 (panels A-D) and Figure S7.
#
# Inputs
#   data/biogeochemistry/processed/2023_Inorganic_Nitrogen_and_Phosphorus_processed.csv
#
# Outputs
#   output/figures/biogeochemistry/*.pdf
#     - Fig.S6A.pdf (Ammonium)
#     - Fig.S6B.pdf (Nitrate)
#     - Fig.S6C.pdf (NH4:NO3 ratio)
#     - Fig.S6D.pdf (N:P ratio)
#     - Fig.S7.pdf  (Phosphorus)
#
# Notes
#   - Uses shared theme + font helpers (plot_theme.R, fonts.R).
#   - Visual encoding:
#       fill  = June core section colors
#       color = December core section colors (June outlines are grey)
#       alpha = hides December fill so December appears as outline-only
#
# Created:  2026-02-20
# Author:   X. Takver
# ------------------------------------------------------------ #

## ----------------------------------------- ##
#             Package Import ----
## ----------------------------------------- ##

librarian::shelf(tidyverse, here, patchwork, cowplot, Cairo, tinythemes)

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

df_NP <- read_csv(file.path(DD, '2023_Inorganic_Nitrogen_and_Phosphorus_processed.csv')) |>
  mutate(
    hillslope = as_factor(hillslope),
    sampling_month = fct_relevel(as_factor(sampling_month), 'June', 'December'),
    core_section = fct_relevel(as_factor(core_section), 'Top', 'Bottom')
  )

# depth_label re-expresses core_section (Top/Bottom) as the depth
# ranges (0-10cm/20-30cm) shown in the published figure legends
df_NP <- df_NP |>
  mutate(depth_label = case_when(core_section == 'Top' ~ '0-10cm',
                                 core_section == 'Bottom' ~ '20-30cm'),
         .after = core_section)

## ----------------------------------------- ##
#             Derived metrics ----
## ----------------------------------------- ##

# Ratios (avoids division-by-zero for NO3 = 0)
df_NP <- df_NP |>
  mutate(
    no3_adj = if_else(is.na(NO3N), NA_real_, if_else(NO3N == 0, 0.0005, NO3N)),
    nh4_no3_ratio = NH4N / no3_adj,
    total_n = NH4N + NO3N,
    np_ratio = total_n / P_bicarb
  )

## ----------------------------------------- ##
#             Plot helpers ----
## ----------------------------------------- ##

# Color Palette
cols_hs <- c('0-10cm' = '#A9D9DE', '20-30cm' = '#759295')
cols_ls <- c('0-10cm' = '#C55648', '20-30cm' = '#7E2E27')
col_june_outline <- '#515151'

# Boxplot style: June is filled; December is outline-only (colored by depth).
geom_box_bi <- function(...) {
  geom_boxplot(
    aes(
      color = ifelse(sampling_month == 'June', col_june_outline, as.character(depth_label)),
      alpha = ifelse(sampling_month == 'June', 1, 0)
    ),
    linewidth = 0.25,
    ...
  )
}

# Build paired HS/LS panels and collect guides.
paired_panels <- function(plt_hs, plt_ls) {
  (plt_hs + plt_ls) + plot_layout(guides = 'collect')
}

## ----------------------------------------- ##
#                   Plots ----
## ----------------------------------------- ##

## ----------------------------------------- ##
##                 NH4-N ----

# HS
plt.HS.NH4 <-
  ggplot(df_NP |> filter(hillslope == 'HS'),
         aes(x = hillslope, y = NH4N, fill = depth_label)) +
  geom_box_bi() +
  coord_cartesian(ylim = c(0, 65)) +
  scale_fill_manual(values = cols_hs) +
  scale_color_manual(values = c(cols_hs, col_june_outline = col_june_outline),
                     breaks = c('0-10cm', '20-30cm')) +
  scale_alpha_identity() +
  labs(
    y = expression(NH[4]^'+' ~ '(mg/kg)'),
    x = '',
    title = 'Ammonium',
    fill = 'Core Section'
  ) +
  guides(color = guide_legend(title = 'December'),
         fill  = guide_legend(title = 'June')) +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'))

# LS
plt.LS.NH4 <-
  ggplot(df_NP |> filter(hillslope == 'LS'),
         aes(x = hillslope, y = NH4N, fill = depth_label)) +
  geom_box_bi() +
  coord_cartesian(ylim = c(0, 65)) +
  scale_fill_manual(values = cols_ls) +
  scale_color_manual(values = c(cols_ls, col_june_outline = col_june_outline),
                     breaks = c('0-10cm', '20-30cm')) +
  scale_alpha_identity() +
  labs(x = 'Site ID', title = '', fill = 'Core Section') +
  guides(color = guide_legend(title = 'December'),
         fill  = guide_legend(title = 'June')) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    plot.margin  = unit(c(0, 0, 0, 0), 'cm')
  )

plt.NH4 <- paired_panels(plt.HS.NH4, plt.LS.NH4)

CairoPDF(file.path(PD, 'Fig.S6A.pdf'), width = 4, height = 3.5)
print(plt.NH4)
dev.off()

rm(plt.NH4, plt.HS.NH4, plt.LS.NH4)

## ----------------------------------------- ##
##                 NO3-N ----

# HS
plt.HS.NO3 <-
  ggplot(df_NP |> filter(hillslope == 'HS'),
         aes(x = hillslope, y = NO3N, fill = depth_label)) +
  geom_box_bi() +
  coord_cartesian(ylim = c(0, 25)) +
  scale_fill_manual(values = cols_hs) +
  scale_color_manual(values = c(cols_hs, col_june_outline = col_june_outline),
                     breaks = c('0-10cm', '20-30cm')) +
  scale_alpha_identity() +
  labs(
    y = expression(NO[3]^'-' ~ '(mg/kg)'),
    x = '',
    title = 'Nitrate',
    fill = 'Core Section'
  ) +
  guides(color = guide_legend(title = 'December'),
         fill  = guide_legend(title = 'June')) +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'))

# LS
plt.LS.NO3 <-
  ggplot(df_NP |> filter(hillslope == 'LS'),
         aes(x = hillslope, y = NO3N, fill = depth_label)) +
  geom_box_bi() +
  coord_cartesian(ylim = c(0, 25)) +
  scale_fill_manual(values = cols_ls) +
  scale_color_manual(values = c(cols_ls, col_june_outline = col_june_outline),
                     breaks = c('0-10cm', '20-30cm')) +
  scale_alpha_identity() +
  labs(x = 'Site ID', title = '', fill = 'Core Section') +
  guides(color = guide_legend(title = 'December'),
         fill  = guide_legend(title = 'June')) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    plot.margin  = unit(c(0, 0, 0, 0), 'cm')
  )

plt.NO3 <- paired_panels(plt.HS.NO3, plt.LS.NO3)

CairoPDF(file.path(PD, 'Fig.S6B.pdf'), width = 4, height = 3.5)
print(plt.NO3)
dev.off()

rm(plt.NO3, plt.HS.NO3, plt.LS.NO3)

## ----------------------------------------- ##
##              NH4:NO3 ratio ----

# HS
plt.HS.Nratio <-
  ggplot(df_NP |> filter(hillslope == 'HS'),
         aes(x = hillslope, y = nh4_no3_ratio, fill = depth_label)) +
  geom_box_bi() +
  coord_cartesian(ylim = c(0, 3750)) +
  scale_fill_manual(values = cols_hs) +
  scale_color_manual(values = c(cols_hs, col_june_outline = col_june_outline),
                     breaks = c('0-10cm', '20-30cm')) +
  scale_alpha_identity() +
  labs(
    y = expression(NH[4]^'+' ~ '/' ~ NO[3]^'-'),
    x = '',
    title = 'Ammonium–Nitrate Ratio',
    fill = 'Core Section'
  ) +
  guides(color = guide_legend(title = 'December'),
         fill  = guide_legend(title = 'June')) +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'))

# LS
plt.LS.Nratio <-
  ggplot(df_NP |> filter(hillslope == 'LS'),
         aes(x = hillslope, y = nh4_no3_ratio, fill = depth_label)) +
  geom_box_bi() +
  coord_cartesian(ylim = c(0, 3750)) +
  scale_fill_manual(values = cols_ls) +
  scale_color_manual(values = c(cols_ls, col_june_outline = col_june_outline),
                     breaks = c('0-10cm', '20-30cm')) +
  scale_alpha_identity() +
  labs(x = 'Site ID', title = '', fill = 'Core Section') +
  guides(color = guide_legend(title = 'December'),
         fill  = guide_legend(title = 'June')) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    plot.margin  = unit(c(0, 0, 0, 0), 'cm')
  )

plt.Nratio <- paired_panels(plt.HS.Nratio, plt.LS.Nratio)

CairoPDF(file.path(PD, 'Fig.S6C.pdf'), width = 4, height = 3.5)
print(plt.Nratio)
dev.off()

rm(plt.Nratio, plt.HS.Nratio, plt.LS.Nratio)

## ----------------------------------------- ##
##                 Phosphorus ----

# HS
plt.HS.P <-
  ggplot(df_NP |> filter(hillslope == 'HS'),
         aes(x = hillslope, y = P_bicarb, fill = depth_label)) +
  geom_box_bi() +
  coord_cartesian(ylim = c(0, 40)) +
  scale_fill_manual(values = cols_hs) +
  scale_color_manual(values = c(cols_hs, col_june_outline = col_june_outline),
                     breaks = c('0-10cm', '20-30cm')) +
  scale_alpha_identity() +
  labs(
    y = 'P (mg/kg)',
    x = '',
    title = 'Phosphorus',
    fill = 'Core Section'
  ) +
  guides(color = guide_legend(title = 'December'),
         fill  = guide_legend(title = 'June')) +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'))

# LS
plt.LS.P <-
  ggplot(df_NP |> filter(hillslope == 'LS'),
         aes(x = hillslope, y = P_bicarb, fill = depth_label)) +
  geom_box_bi(linewidth = 0.25) +
  coord_cartesian(ylim = c(0, 40)) +
  scale_fill_manual(values = cols_ls) +
  scale_color_manual(values = c(cols_ls, col_june_outline = col_june_outline),
                     breaks = c('0-10cm', '20-30cm')) +
  scale_alpha_identity() +
  labs(x = 'Site ID', title = '', fill = 'Core Section') +
  guides(color = guide_legend(title = 'December'),
         fill  = guide_legend(title = 'June')) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    plot.margin  = unit(c(0, 0, 0, 0), 'cm')
  )

plt.P <- paired_panels(plt.HS.P, plt.LS.P)

CairoPDF(file.path(PD, 'Fig.S7.pdf'), width = 4, height = 3.5)
print(plt.P)
dev.off()

rm(plt.P, plt.HS.P, plt.LS.P)

## ----------------------------------------- ##
##                 Total N:P ----

# HS
plt.HS.NP <-
  ggplot(df_NP |> filter(hillslope == 'HS'),
         aes(x = hillslope, y = np_ratio, fill = depth_label)) +
  geom_box_bi() +
  coord_cartesian(ylim = c(0, 25)) +
  scale_fill_manual(values = cols_hs) +
  scale_color_manual(values = c(cols_hs, col_june_outline = col_june_outline),
                     breaks = c('0-10cm', '20-30cm')) +
  scale_alpha_identity() +
  labs(
    y = 'Total N / Total P',
    x = '',
    title = 'N:P Ratio',
    fill = 'Core Section'
  ) +
  guides(color = guide_legend(title = 'December'),
         fill  = guide_legend(title = 'June')) +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'))

# LS
plt.LS.NP <-
  ggplot(df_NP |> filter(hillslope == 'LS'),
         aes(x = hillslope, y = np_ratio, fill = depth_label)) +
  geom_box_bi() +
  coord_cartesian(ylim = c(0, 25)) +
  scale_fill_manual(values = cols_ls) +
  scale_color_manual(values = c(cols_ls, col_june_outline = col_june_outline),
                     breaks = c('0-10cm', '20-30cm')) +
  scale_alpha_identity() +
  labs(x = 'Site ID', title = '', fill = 'Core Section') +
  guides(color = guide_legend(title = 'December'),
         fill  = guide_legend(title = 'June')) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    plot.margin  = unit(c(0, 0, 0, 0), 'cm')
  )

plt.NP <- paired_panels(plt.HS.NP, plt.LS.NP)

CairoPDF(file.path(PD, 'Fig.S6D.pdf'), width = 4, height = 3.5)
print(plt.NP)
dev.off()

rm(plt.NP, plt.HS.NP, plt.LS.NP)

## ----------------------------------------- ##
#                 Cleanup ----
## ----------------------------------------- ##

rm(df_NP, cols_hs, cols_ls, col_june_outline, geom_box_bi, paired_panels)
