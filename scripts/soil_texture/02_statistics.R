# -----------------------------------------------------------------------------#
# Purpose
#   Test for hillslope effects on soil particle-size fractions (PSF; sand,
#   silt, clay). Fits linear mixed-effects models (PSF ~ hillslope + horizon,
#   random intercept for site_id) and simple linear models (PSF ~ hillslope +
#   depth_groups) for each particle-size class, then exports ANOVA tables and
#   model coefficients as a single PDF report.
#
# Inputs (read from DD)
#   - 2023_WS10-HJA_PSF_Processed.csv
#
# Outputs (written to SD)
#   - 2023_WS10-HJA_PSF_Statistics.pdf
#       LME ANOVA/fixed effects and LM Type III ANOVA/coefficients for
#       clay, silt, and sand.
#
# Created: 01/23/2026
# Author: X. Takver
# -----------------------------------------------------------------------------#

## ----------------------------------------- ##
#             Import Packages ----
## ----------------------------------------- ##

# Data Analysis
librarian::shelf(tidyverse, here, glue, lme4, lmerTest, car,
                 broom, broom.mixed, gridExtra, grid)

## ----------------------------------------- ##
#             Source Functions ----
## ----------------------------------------- ##

# Directory helpers (ensure_dir)
source(here('R', 'functions', 'paths_and_directories.R'))

## ----------------------------------------- ##
#          Directory Navigation ----
## ----------------------------------------- ##

# Makes quick navigation directory commands

# Data Directory Variables
DD <- here('data', 'soil_texture', 'processed')

# Stats Directory Variables
SD <- here('output','statistics', 'soil_properties')

# Creates a statistical output directory
ensure_dir(SD)

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

df_psf <- read_csv(file.path(DD, '2023_WS10-HJA_PSF_Processed.csv'))

## ----------------------------------------- ##
#            Data Preparation ----
## ----------------------------------------- ##

## ----------------------------------------- ##
##                PSF ----

# Derive hillslope from the site_id prefix (HS/LS), set factor levels, and
# bin samples into 15 cm depth increments for the simple LMs below.
psf_model_data <- df_psf |> 
  mutate(
    hillslope = if_else(grepl('^HS', site_id), 'High', 'Low'),
    hillslope = factor(hillslope, levels = c('Low', 'High')),
    
    site_id   = factor(site_id),
    horizon  = factor(horizon),
    psize_class = factor(psize_class))  |> 
  
  separate(depth_range_cm, into = c('top_depth_cm', 'bottom_depth_cm'), sep = '-') |> 
  mutate(top_depth_cm = as.numeric(top_depth_cm),
         bottom_depth_cm = as.numeric(bottom_depth_cm)) |> 
         mutate(mid_depth_cm = ((bottom_depth_cm - top_depth_cm)/2) + top_depth_cm,
                .after = bottom_depth_cm) |> 
  
  mutate(depth_groups = case_when(mid_depth_cm <= 15 ~ 15,
                                  mid_depth_cm > 15 & mid_depth_cm <= 30 ~ 30,
                                  mid_depth_cm > 30 & mid_depth_cm <= 50 ~ 50,
                                  mid_depth_cm > 50 & mid_depth_cm <=75 ~ 75,
                                  mid_depth_cm > 75 & mid_depth_cm <= 100 ~ 100),
         .after = mid_depth_cm) |> 
  mutate(depth_groups = as_factor(depth_groups))

## ----------------------------------------- ##
#           Statistical Analysis ----
## ----------------------------------------- ##

# LME: fixed effects of hillslope and horizon, random intercept for site_id
# to account for repeated horizon samples within the same site.

# Separate data by PSF
clay_data <- psf_model_data |>
  filter(psize_class == 'clay')

silt_data <- psf_model_data |>
  filter(psize_class == 'silt')

sand_data <- psf_model_data |>
  filter(psize_class == 'sand')

## ----------------------------------------- ##
##                PSF LME ----

###               Clay ----
m_clay <- lmer(psf_mean ~ hillslope + horizon + (1 | site_id),
               data = clay_data)

anova(m_clay)
summary(m_clay)    

###               Silt ----
m_silt <- lmer(psf_mean ~ hillslope + horizon + (1 | site_id),
               data = silt_data)

anova(m_silt)
summary(m_silt)

###               Sand ----
m_sand <- lmer(psf_mean ~ hillslope + horizon + (1 | site_id),
               data = sand_data)

anova(m_sand)
summary(m_sand) 

## ----------------------------------------- ##
##          Simple Linear Models ----

# No hillslope x depth_groups interaction: with only 2 site_ids per hillslope,
# some hillslope-by-depth combinations have no data, producing aliased
# coefficients. Main effects only.

###             By Depth ----

# Clay
m_clay_lm <- lm(psf_mean ~ hillslope + depth_groups, data = clay_data)

Anova(m_clay_lm, type = 'III')  # tests for Hillslope and depth group
summary(m_clay_lm)

# Silt
m_silt_lm <- lm(psf_mean ~ hillslope + depth_groups, data = silt_data)

Anova(m_silt_lm, type = 'III')  # tests for Hillslope and depth group
summary(m_silt_lm)

# Sand
m_sand_lm <- lm(psf_mean ~ hillslope + depth_groups, data = sand_data)

Anova(m_sand_lm, type = 'III')  # tests for Hillslope and depth group
summary(m_sand_lm)

## ----------------------------------------- ##
#      Supplement Summary Statistics ----
## ----------------------------------------- ##

# Table S4: horizon-level summary statistics (n, mean, SD) for each
# particle-size fraction by hillslope and depth group. Descriptive only --
# not derived from a fitted model.
summarize_psf <- function(data, fraction_label) {
  data |> 
    group_by(hillslope, depth_groups) |> 
    summarise(
      pit_n = n_distinct(site_id),
      horizon_sample_n = n(),
      instrument_replicate_n = sum(psf_n_runs, na.rm = TRUE),
      mean = mean(psf_mean, na.rm = TRUE),
      sd = sd(psf_mean, na.rm = TRUE),
      .groups = 'drop') |> 
    mutate(fraction = fraction_label, .before = 1)
}

table_s4 <- bind_rows(
  summarize_psf(clay_data, 'Clay'),
  summarize_psf(silt_data, 'Silt'),
  summarize_psf(sand_data, 'Sand')) |> 
  mutate(hillslope = dplyr::recode(hillslope, 'High' = 'HS', 'Low' = 'LS')) |> 
  arrange(fraction, hillslope, depth_groups)


## ----------------------------------------- ##
#        Export Statistical Results ----
## ----------------------------------------- ##

# Helper: draw a section header on its own page area, then a table below it
add_table_page <- function(df, title) {
  grid.newpage()
  # Section title
  grid.text(title, x = 0.5, y = 0.97, gp = gpar(fontsize = 14, fontface = 'bold'))
  # Round numeric columns for readability, then draw as a table
  df_display <- df |> 
    mutate(across(where(is.numeric), ~ round(.x, 4)))
  tbl <- tableGrob(df_display, rows = NULL,
                   theme = ttheme_default(base_size = 8))
  grid.draw(tbl)
}

# Helper: add all results for one LME model
export_lme <- function(model, label) {
  add_table_page(as.data.frame(anova(model)) |> rownames_to_column('term'),
                 glue('{label}: LME ANOVA (psf_mean ~ hillslope + horizon + (1|site_id))'))
  add_table_page(broom.mixed::tidy(model, effects = 'fixed'),
                 glue('{label}: LME Fixed Effects'))
}

# Helper: add all results for one simple LM model
export_lm <- function(lm_model, label) {
  add_table_page(Anova(lm_model, type = 'III') |> as.data.frame() |> rownames_to_column('term'),
                 glue('{label}: LM Type III ANOVA (psf_mean ~ hillslope + depth_groups)'))
  add_table_page(broom::tidy(lm_model),
                 glue('{label}: LM Coefficients'))
}

# ---- Write everything to one PDF ----
pdf(file.path(SD, '2023_WS10-HJA_PSF_Statistics.pdf'), width = 8.5, height = 11)

export_lme(m_clay, 'Clay')
export_lm(m_clay_lm, 'Clay')

export_lme(m_silt, 'Silt')
export_lm(m_silt_lm, 'Silt')

export_lme(m_sand, 'Sand')
export_lm(m_sand_lm, 'Sand')

# Table S4: descriptive summary statistics (not a model result)
add_table_page(table_s4, 'Table S4: Particle-Size Fraction Summary Statistics')

dev.off()

# Write as a standalone CSV for building the formatted supplement table
write_csv(table_s4, file.path(SD, 'TableS4_PSF_Summary_Statistics.csv'))