# -----------------------------------------------------------------------------#
# Purpose
#   Core plotting logic for the soil texture triangle figures (S16C & S16D)
#   and their shared legend. Isolated into a function so it can be run in a
#   disposable subprocess (see 03_figure-S16C-D.R), since ggtern patches
#   ggplot2 internals on load in ways that cannot be reliably undone in the
#   calling R session.
#
# Inputs (read from DD)
#   - 2023_WS10-HJA_PSF_Processed-long.csv
#
# Outputs (written to PD)
#   - Fig.S16C_Texture_Triangle.pdf
#   - Fig.S16D_Texture_Triangle.pdf
#   - Fig.S16C-D_Legend.pdf
#
# Created: 01/23/2026
# Author: X. Takver
# -----------------------------------------------------------------------------#

# Writes three PDF files to PD as a side effect; has no return value.
run_triangle_figures <- function() {
  
  librarian::shelf(tidyverse, here)
  librarian::shelf(tinythemes, Cairo, ggtern, grid, cowplot)
  
  source(here('R', 'functions', 'paths_and_directories.R'))
  source(here('R', 'functions', 'fonts.R'))
  source(here('R', 'functions', 'plot_theme.R'))
  
  DD <- here('data', 'soil_texture', 'processed')
  PD <- here('output', 'figures', 'soil_properties')
  ensure_dir(PD)
  
  set_plot_theme()
  

  ## ----------------------------------------- ##
  #                 Data Import ----
  ## ----------------------------------------- ##
  
  # Read processed particle size fractions in long format
  
  df_texture <- read_csv(here(DD, '2023_WS10-HJA_PSF_Processed-long.csv')) |>
    # Ensure IDs/classes are factors for consistent ordering/labels in plots
    mutate(across(c(site_id, psize_class), as_factor)) |>
    # Standardize class names to title case and enforce plotting order
    mutate(
      psize_class = fct_recode(psize_class, Sand = 'sand', Silt = 'silt', Clay = 'clay'),
      psize_class = fct_relevel(psize_class, 'Sand', 'Silt', 'Clay'))
  
  ## ----------------------------------------- ##
  #        Prepare Texture Triangle Base ----
  ## ----------------------------------------- ##

  # Builds the reusable ternary background (USDA texture classes); independent
  # of df_texture and shared by both panels below.
  data(USDA)

  # Place USDA class labels at polygon midpoints (mean of vertices)
  USDA.LAB <- USDA |>
    group_by(Label) |>
    summarise(across(1:3, mean, .names = 'mean_{.col}'))
  
  # Adjust label angles for readability (only where needed)
  USDA.LAB$Angle = 0
  USDA.LAB$Angle[which(USDA.LAB$Label == 'Loamy Sand')] = -35
  
  # Rename to match ggtern axis aesthetics (x = Sand, y = Clay, z = Silt)
  USDA.LAB <- USDA.LAB |>
    rename(Clay = mean_Clay, Sand = mean_Sand, Silt = mean_Silt)
  
  # Base ternary plot: USDA polygons + class labels + axis arrows
  texture_triangle <-
    ggplot(data = USDA, aes(y = Clay, x = Sand, z = Silt)) +
    coord_tern(L = 'x', T = 'y', R = 'z') +
    geom_polygon(aes(fill = Label, color = Label),
                 alpha = 0.75, linewidth = 0.25, color = 'black') +
    scale_fill_manual(values = c('#B84430', '#d56f24', '#cc713b', '#d98f42',
                                 '#e6942d', '#F29A18', '#c36744', '#cd8957',
                                 '#c5a38a', '#bb8177', '#c1785c', '#BDBDBD'),
                      guide = 'none') +
    geom_text(data = USDA.LAB, aes(label = Label, angle = Angle),
              color = 'black', size = 2, family = 'sans') +
    labs(yarrow = 'Clay (%)', zarrow = 'Silt (%)', xarrow = 'Sand (%)') +
    theme_linedraw(base_size = 5.5) +
    theme_showarrows() +
    theme_clockwise() +
    theme_hidetitles() +
    theme(plot.margin = margin(l = -30, r = -30, t = -100, b = -75),
          tern.axis.arrow = element_line(linewidth = 0.75),
          tern.axis.arrow.text = element_text(size = 6.5, family = 'sans'),
          tern.axis.arrow.sep = 0.075)
  
  texture_triangle
  
  ## ----------------------------------------- ##
  #             Plot Upper 30 cm ----
  ## ----------------------------------------- ##
  
  # Convert long texture data to wide (columns: Sand/Silt/Clay), subset to upper 30 cm
  df_triangle <-
    df_texture |>
    select(site_id:mean_psf) |>
    pivot_wider(names_from = psize_class, values_from = mean_psf) |>
    filter(mid_depth_cm <= 30) |>
    mutate(core.section = case_when(mid_depth_cm <= 10 ~ '0-10cm',
                                    mid_depth_cm >= 20 ~ '20-30cm'))
  
  ## ----------------------------------------- ##
  ##             Plot on Triangle ----
  
  # Manually add small jitter to reduce overplotting, then renormalize to sum to ~100%
  # (jitter magnitudes chosen by visual inspection for the upper-30cm subset)
  set.seed(123)
  
  df_triangle_jitter <- df_triangle |>
    mutate(Sand_j = Sand + runif(n(), 0, 7.5),
           Clay_j = Clay + runif(n(), -3, 0),
           Silt_j = Silt + runif(n(), -3, 3)) |>
    rowwise() |>
    mutate(total = sum(Sand_j, Clay_j, Silt_j),
           Sand_j = Sand_j / total * 100,
           Clay_j = Clay_j / total * 100,
           Silt_j = Silt_j / total * 100) |>
    ungroup()
  
  # NOTE: ggtern intercepts the z aesthetic to generate ternary coordinates. 
  # An innocuous warning will occur as a result.
  plt.triangle <-
    texture_triangle +
    geom_point(data = df_triangle_jitter,
               aes(x = Sand_j, y = Clay_j, z = Silt_j, shape = core.section),
               # Map site IDs to fill colors (kept outside aes to avoid extra legend handling)
               fill = df_triangle_jitter$site_id |> 
                 dplyr::recode('HS1' = '#759295', 'HS2' = '#A9D9DE',
                               'LS1' = '#7E2E27', 'LS2' = '#C55648'),
      color = 'grey30',inherit.aes = FALSE, size = 2, alpha = 0.85) +
    scale_shape_manual(values = c('0-10cm' = 21, '20-30cm' = 24), guide = 'none') +
    labs(color = 'Site ID', shape = 'Depth')
  
  plt.triangle
  
  
  ###         Export upper 30 Triangle ----
  
  CairoPDF(here(PD, 'Fig.S16C_Texture_Triangle.pdf'), width = 3.25, height = 3.25)
  print(plt.triangle)
  dev.off()
  
  ## ----------------------------------------- ##
  #               Plot All Data ----
  ## ----------------------------------------- ## 
  
  # Wide format with depth bins for all samples
  df_triangle_all <-
    df_texture |>
    select(site_id:mean_psf) |>
    pivot_wider(names_from = psize_class, values_from = mean_psf) |>
    mutate(core.section = case_when(mid_depth_cm <= 10 ~ '0-10cm',
                                    mid_depth_cm > 10 & mid_depth_cm <= 30 ~ '20-30cm',
                                    mid_depth_cm > 30 & mid_depth_cm <= 60 ~ '30-60cm',
                                    mid_depth_cm > 60 & mid_depth_cm <= 100 ~ '60-100cm'))
  
  # Manually Jitter + renormalize (tuned for readability across more depth groups)
  # (larger/asymmetric ranges here vs. above to separate the additional depth classes)
  set.seed(123)
  
  df_triangle_jitter <- df_triangle_all |>
    mutate(Sand_j = Sand + runif(n(), 2, 7),
           Clay_j = Clay + runif(n(), -1, -1),
           Silt_j = Silt + runif(n(), -1, 1.25)) |>
    rowwise() |>
    mutate(total = sum(Sand_j, Clay_j, Silt_j),
           Sand_j = Sand_j / total * 100,
           Clay_j = Clay_j / total * 100,
           Silt_j = Silt_j / total * 100) |>
    ungroup()
  
  # Adjust point sizes to equalize visual weight across shapes
  df_triangle_jitter <- df_triangle_jitter |>
    mutate(size_adj = case_when(core.section == '0-10cm' ~ 1.75,
                                core.section == '20-30cm' ~ 1.75,
                                core.section == '30-60cm' ~ 2.25,
                                core.section == '60-100cm' ~ 2.25))
  
  # Plot all depths on the triangle
  plt.triangle.all <-
    texture_triangle +
    geom_point(data = df_triangle_jitter,
               aes(x = Sand_j, y = Clay_j, z = Silt_j, shape = core.section, size = size_adj),
               fill = df_triangle_jitter$site_id |>
                 dplyr::recode('HS1' = '#759295', 'HS2' = '#A9D9DE',
                               'LS1' = '#7E2E27', 'LS2' = '#C55648'),
               color = 'grey30', inherit.aes = FALSE, alpha = 0.85) +
    scale_shape_manual(values = c('0-10cm' = 21, '20-30cm' = 24, '30-60cm' = 22, '60-100cm' = 23),
                       guide = 'none') +
    scale_size_identity(guide = 'none') +
    labs(color = 'Site ID', shape = 'Depth')
  
  plt.triangle.all
  
  
  ###         Export Triangle All ----
  
  CairoPDF(here(PD, 'Fig.S16D_Texture_Triangle.pdf'), width = 3.25, height = 3.25)
  print(plt.triangle.all)
  dev.off()
  
  ## ----------------------------------------- ##
  #           Build Shared Legend ----
  ## ----------------------------------------- ##
  
  # Minimal dummy plot solely to generate a legend matching the Site ID colors
  # and Depth shapes used across both triangle panels (C & D).
  legend_data <- expand_grid(
    site_id = factor(c('HS1', 'HS2', 'LS1', 'LS2'), levels = c('HS1', 'HS2', 'LS1', 'LS2')),
    core.section = factor(c('0-10cm', '20-30cm', '30-60cm', '60-100cm'),
                          levels = c('0-10cm', '20-30cm', '30-60cm', '60-100cm')))
  
  legend_plot <- ggplot(legend_data, aes(x = 1, y = 1, fill = site_id, shape = core.section)) +
    geom_point(size = 3, color = 'grey30') +
    scale_fill_manual(name = 'Site ID',
                      values = c('HS1' = '#759295', 'HS2' = '#A9D9DE',
                                 'LS1' = '#7E2E27', 'LS2' = '#C55648'),
                      labels = c('HS1' = 'HS-1', 'HS2' = 'HS-2',
                                 'LS1' = 'LS-1', 'LS2' = 'LS-2')) +
    scale_shape_manual(name = 'Depth',
                       values = c('0-10cm' = 21, '20-30cm' = 24,
                                  '30-60cm' = 22, '60-100cm' = 23)) +
    guides(fill = guide_legend(override.aes = list(shape = 21)),
           shape = guide_legend(override.aes = list(fill = 'grey50'))) +
    theme_void() +
    theme(legend.position = 'bottom', legend.box = 'horizontal',
          legend.title = element_text(face = 'bold'))
  
  shared_legend <- cowplot::get_legend(legend_plot)
  
  ###         Export Shared Legend ----
  
  CairoPDF(here(PD, 'Fig.S16C-D_Legend.pdf'), width = 8, height = 0.75)
  grid::grid.draw(shared_legend)
  dev.off()
  

}
