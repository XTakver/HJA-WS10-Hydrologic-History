# -----------------------------------------------------------------------------#
# Aggregate Size Distribution (ASD): Linear Model and ANOVA
#
# Purpose: Import processed ASD geometric mean diameter summary statistics, define
# discrete depth groups, and test hillslope and depth effects using additive and
# interaction linear models (Type III ANOVA). Model outputs are exported to a
# plain-text file for reproducibility.
#
# Input:
#   - Data/aggregate_size_distribution/Processed/processed_asd_GeoMean_SummaryStats.csv
#
# Output:
#   - Statistics/Soil_Properties/ASD_Statistical_Output.txt
#
# Created: 2026-01-23
# Author: X. Takver
# -----------------------------------------------------------------------------#

## ----------------------------------------- ##
#             Import Packages ----
## ----------------------------------------- ##

# Data processing and statistical modeling packages
librarian::shelf(tidyverse, here, emmeans, car)

## ----------------------------------------- ##
#             Source Functions ----
## ----------------------------------------- ##

# Directory helpers (ensure_dir)
source(here('R', 'functions', 'paths_and_directories.R'))

## ----------------------------------------- ##
#          Directory Navigation ----
## ----------------------------------------- ##

# DD: directory containing processed ASD datasets
# SD: directory for statistical outputs written by this script
DD <- here('data', 'aggregate_size_distribution', 'processed')
SD <- here('output', 'statistics', 'soil_properties')

# Ensure output directories exist (idempotent)
ensure_dir(SD)

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

# Processed ASD geometric mean summary statistics:
#   - read_csv() imports the replicate-aggregated geometric mean diameter and
#     associated summary statistics used for downstream modeling
df_asd <- read_csv(file.path(DD, 'processed_asd_GeoMean_SummaryStats.csv'))

## ----------------------------------------- ##
#               Data Preparation ----
## ----------------------------------------- ##

# Create depth-group bins from continuous mid-depth (cm) for categorical modeling.
# Convert grouping variables to factors so they are treated as categorical predictors.
asd_model_data <- 
  df_asd |> 
  mutate(depth_groups = case_when(mid_depth <= 15 ~ 15,
                                  mid_depth > 15 & mid_depth <= 30 ~ 30,
                                  mid_depth > 30 & mid_depth <= 50 ~ 50,
                                  mid_depth > 50 & mid_depth <=75 ~ 75,
                                  mid_depth > 75 & mid_depth <= 100 ~ 100),
         .after = mid_depth) |> 
  mutate(hillslope = as_factor(hillslope),
         pit_id = as_factor(pit_id),
         site_id = as_factor(site_id),
         depth_groups = as_factor(depth_groups))

## --------------------------------------------- ##
#                 Statistics ----
## --------------------------------------------- ##

## --------------------------------------------- ##
##                Setup Export ----

# Export model outputs (including messages/warnings) to a plain-text file
out_file <- file.path(SD, 'ASD_Statistical_Output.txt')

zz <- file(out_file, open = 'wt')
sink(zz)
sink(zz, type = 'message')

cat('ASD stats export\n')
cat('Run time:', as.character(Sys.time()), '\n\n')

## --------------------------------------------- ##
##    Statistical analysis: additive model ----



# Use sum-to-zero contrasts so Type III tests are interpretable for main effects.
options(contrasts = c('contr.sum', 'contr.poly'))

# Additive model: effects of hillslope and depth group on geometric mean diameter.
cat('\n## Additive model\n')
asd_lm <- lm(final_mean_diam ~ hillslope + depth_groups, data = asd_model_data)

cat('\n### Type III ANOVA\n')
print(Anova(asd_lm, type = 'III'))

cat('\n### Model summary\n')
print(summary(asd_lm))

# Estimated marginal means for hillslope (averaged over depth_groups) and Tukey pairs.
cat('\n### EMMs (hillslope) + Tukey pairs\n')
emm_asd <- emmeans(asd_lm, ~hillslope)
print(emm_asd)
print(pairs(emm_asd, adjust = 'tukey'))

## --------------------------------------------- ##
#  Statistical analysis: interaction model ----
## --------------------------------------------- ##

# Interaction model: tests whether the hillslope effect differs by depth group.
# Depth groups are restricted to those represented in both hillslopes to avoid
# non-estimable coefficients due to empty hillslope × depth combinations.
asd_model_data <- droplevels(asd_model_data)

# Tabulate observations per Hillslope × depth cell to assess balance/sparsity.
cell_counts <- with(asd_model_data, table(hillslope, depth_groups, useNA = 'ifany'))

# Keep only depth bins represented in *both* hillslopes
keep_depths <- colnames(cell_counts)[apply(cell_counts > 0, 2, all)]

asd_balanced <-
  asd_model_data |>
  filter(depth_groups %in% keep_depths) |>
  droplevels()

# Fit the interaction model on the balanced subset.
asd_lm_int <- lm(final_mean_diam ~ hillslope * depth_groups, data = asd_balanced)

# Capture aliasing diagnostics for troubleshooting (not printed to export file).
aliased <- alias(asd_lm_int)


cat('\n## Interaction model\n')
cat('\n### Type III ANOVA\n')
print(Anova(asd_lm_int, type = 'III'))

# Hillslope EMMs within each depth group and Tukey-adjusted pairwise comparisons.
cat('\n### EMMs (hillslope | depth_groups) + Tukey pairs\n')
emm_int <- emmeans(asd_lm_int, ~ hillslope | depth_groups)
print(emm_int)
print(pairs(emm_int, adjust = 'tukey'))

## --------------------------------------------- ##
##          Close Export Environment ----

# Stop capturing and close file connection
sink(type = 'message')
sink()
close(zz)

message('Wrote stats to: ', out_file)
