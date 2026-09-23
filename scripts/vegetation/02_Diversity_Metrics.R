# ------------------------------------------------------------ #
# 02_Diversity_Metrics.R
#
# Purpose
#   Build a species community matrix and calculate alpha-diversity
#   indices (Shannon, Simpson, inverse Simpson, Fisher's alpha,
#   richness, evenness) per year x plot x hillslope for the WS10
#   long-term ground-cover survey (all plots, all years).
#
# Inputs
#   data/vegetation/processed/tidy/WS10_GroundCover_AllPlots.csv
#
# Outputs (written to OD)
#   - Diversity_DF_All.csv       (feeds 02_statistics_Vegetation_Community.R,
#                                  03_figure-S14A.R)
#   - Community_Matrix_All.csv   (feeds 02_statistics_Vegetation_Community.R,
#                                  03_figure-S14A.R)
#
# Notes
#   - Evenness = shannon / log(richness).
#
# Created:  2025-06-04
# Author:   X. Takver
# ------------------------------------------------------------ #

## ----------------------------------------- ##
#             Package Import ----
## ----------------------------------------- ##

librarian::shelf(tidyverse, here, vegan)

## ----------------------------------------- ##
#             Source Functions ----
## ----------------------------------------- ##

source(here('R', 'functions', 'paths_and_directories.R'))  # ensure_dir()

## ----------------------------------------- ##
#           Directory Creation ----
## ----------------------------------------- ##

DD <- here('data', 'vegetation', 'processed')
OD <- here('data', 'vegetation', 'processed', 'diversity')

ensure_dir(OD)

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

# All years, all plots. Groundcover
df_grndcover.all <- read_csv(file.path(DD, 'WS10_GroundCover_AllPlots_processed.csv'))

## ----------------------------------------- ##
#            Diversity Indices ----
## ----------------------------------------- ##

### Format data into a community Matrix ----

# Turn the dataframe into species count data
community_matrix <- df_grndcover.all |>
  count(year, plot, hillslope, species) |>
  pivot_wider(names_from = species,
              values_from = n,
              values_fill = 0) |>
  unite('YR_PLT_HS', year, plot, hillslope, remove = FALSE) |>
  column_to_rownames('YR_PLT_HS')

### Calculate diversity Metrics ----

# Shannon Diversity
shannon_values <- diversity(community_matrix[, -(1:3)], index = 'shannon')

# Simpson Diversity
simpson_values <- diversity(community_matrix[, -(1:3)], index = 'simpson')

# Inverse Simpson
invsimpson_values <- diversity(community_matrix[, -(1:3)], index = 'invsimpson')

# Fishers Alpha
fisher_values <- fisher.alpha(community_matrix[, -(1:3)])

# Species Richness
richness_values <- specnumber(community_matrix[, -(1:2)])

### Create a Diversity Dataframe ----

# Create df and add evenness metric
df_diversity <- community_matrix |>
  select(year, plot, hillslope) |>
  mutate(shannon = shannon_values,
         simpson = simpson_values,
         invsimpson = invsimpson_values,
         fisher = fisher_values,
         richness = richness_values) |>
  mutate(evenness = shannon/log(richness)) |>
  rownames_to_column('YR_PLT_HS')

## ----------------------------------------- ##
#                  Export ----
## ----------------------------------------- ##

write_csv(df_diversity, file.path(OD, 'Diversity_DF_All.csv'))
write_csv(community_matrix, file.path(OD, 'Community_Matrix_All.csv'))

# Cleanup workspace
rm(df_grndcover.all,
   shannon_values, simpson_values, invsimpson_values, fisher_values, richness_values)
