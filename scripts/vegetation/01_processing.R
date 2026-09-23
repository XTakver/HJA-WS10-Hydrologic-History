# ------------------------------------------------------------ #
# Purpose
#   Import WS10 ground-cover and tree survey data (1973-2016) from
#   the H.J. Andrews LTER, assign hillslope (high- vs low-storage),
#   and produce tidied subsets: all long-term plots, plots co-located
#   with our study sites (all years), and co-located plots restricted
#   to 2010-present and 1989-present.
#
# Inputs
#   TP04101_v2.csv (ground cover estimates)
#   TP04103_v2.csv (tree survey)
#
# Outputs (written to OD)
#   - WS10_GroundCover_AllPlots.csv               / WS10_TreeSurvey_AllPlots.csv
#   - WS10_GroundCover_CoLocated_AllYears.csv      / WS10_TreeSurvey_CoLocated_AllYears.csv
#   - WS10_GroundCover_CoLocated_2010-Present.csv  / WS10_TreeSurvey_CoLocated_2010-Present.csv
#   - WS10_GroundCover_CoLocated_1989-Present.csv  / WS10_TreeSurvey_CoLocated_1989-Present.csv
#
# Notes
#   - Plot 35 is excluded: it does not appear in the site metadata
#     or map on the H.J. Andrews website.
#   - Hillslope assignment (High Storage vs Low Storage) is by plot
#     number, per the watershed's hillslope layout.
#   - Ground cover species codes for non-vegetation cover classes
#     (litter, bare ground, logs, etc.) are excluded throughout.
#   - Co-located plots (11, 34, 27, 28, 26, 4, 25, 13) are the subset
#     spatially nearest our biogeochemistry/hydrology sampling sites.
#
# Created:  2025-05-07
# Author:   X. Takver
# ------------------------------------------------------------ #

## ----------------------------------------- ##
#             Package Import ----
## ----------------------------------------- ##

librarian::shelf(tidyverse, here)

## ----------------------------------------- ##
#             Source Functions ----
## ----------------------------------------- ##

source(here('R', 'functions', 'paths_and_directories.R'))  # ensure_dir()

## ----------------------------------------- ##
#             Directories ----
## ----------------------------------------- ##

DD <- here('data', 'vegetation', 'raw')
OD <- here('data', 'vegetation', 'processed')

ensure_dir(OD)

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

# Ground cover estimates
raw_GCE <- read_csv(file.path(DD, 'TP04101_v2.csv')) |>
  mutate(SAMPLEDATE = as_date(SAMPLEDATE, format = '%Y-%m-%d')) |>
  rename_with(tolower)

# Tree survey; select relevant columns
raw_Tree <- read_csv(file.path(DD, 'TP04103_v2.csv'),
                     col_types = cols(DB_NOTES = col_character())) |>
  mutate(SAMPLEDATE = as_date(SAMPLEDATE, format = '%Y-%m-%d')) |>
  select(DBCODE, WATERSHED, PLOTID:MICROPLOT, TREE_SPP:DBH, SAMPLEDATE) |>
  rename_with(tolower)

## ----------------------------------------- ##
#                   Tidy ----
## ----------------------------------------- ##

## ----------------------------------------- ##
#             All Plots ----

# Add hillslope to each data frame; plot 35 is excluded (see Notes)
tidy_tree_all <- raw_Tree |>
  filter(plot != 35) |>
  mutate(hillslope = case_when(plot %in% c(5, 7, 8:11, 17, 18, 23, 24, 27:31, 34, 36) ~ 'High Storage',
                               plot %in% c(1:4, 6, 12:16, 19:22, 25, 26, 32, 33) ~ 'Low Storage'))

# Select relevant columns; drop non-descript ground-cover classes
tidy_grndcover_all <- raw_GCE |>
  filter(plot != 35) |>
  mutate(hillslope = case_when(plot %in% c(5, 7, 8:11, 17, 18, 23, 24, 27:31, 34, 36) ~ 'High Storage',
                               plot %in% c(1:4, 6, 12:16, 19:22, 25, 26, 32, 33) ~ 'Low Storage')) |>
  select(year, plot, hillslope, species, cover) |>
  filter(!species %in% c('LICHEN', 'MOSS', 'LITTER', 'LOG', 'STUMP', 'BARE',
                         'BUTT', 'TREEB', 'TREE', 'HERB'))

# Export: all survey years
write_csv(tidy_tree_all, file.path(OD, 'WS10_TreeSurvey_AllPlots_processed.csv'))
write_csv(tidy_grndcover_all, file.path(OD, 'WS10_GroundCover_AllPlots_processed.csv'))

rm(tidy_tree_all, tidy_grndcover_all)

## ----------------------------------------- ##
#             Co-located Plots ----

# Restrict to plots co-located with our study sites (all years)
tidy_grndcover_coloc <- raw_GCE |>
  filter(plot %in% c(11, 34, 27, 28, 26, 4, 25, 13)) |>
  mutate(hillslope = case_when(plot %in% c(11, 34, 27, 28) ~ 'High Storage',
                               plot %in% c(26, 4, 25, 13) ~ 'Low Storage')) |>
  select(year, plot, hillslope, species) |>
  filter(!species %in% c('LICHEN', 'MOSS', 'LITTER', 'LOG', 'STUMP', 'BARE',
                         'BUTT', 'TREEB', 'TREE', 'HERB'))

tidy_tree_coloc <- raw_Tree |>
  filter(plot %in% c(11, 34, 27, 28, 26, 4, 25, 13)) |>
  mutate(hillslope = case_when(plot %in% c(11, 34, 27, 28) ~ 'High Storage',
                               plot %in% c(26, 4, 25, 13) ~ 'Low Storage'))

# Export: all survey years, co-located plots only
write_csv(tidy_tree_coloc, file.path(OD, 'WS10_TreeSurvey_CoLocated_AllYears_processed.csv'))
write_csv(tidy_grndcover_coloc, file.path(OD, 'WS10_GroundCover_CoLocated_AllYears_processed.csv'))

rm(tidy_tree_coloc, tidy_grndcover_coloc)

## ----------------------------------------- ##
#             Co-located, 2010-present ----

# Co-located plots, restricted to the most recent survey years
tidy_grndcover2010 <- raw_GCE |>
  filter(year >= 2010,
         plot %in% c(11, 34, 27, 28, 26, 4, 25, 13)) |>
  mutate(hillslope = case_when(plot %in% c(11, 34, 27, 28) ~ 'High Storage',
                               plot %in% c(26, 4, 25, 13) ~ 'Low Storage')) |>
  select(year, plot, hillslope, species) |>
  filter(!species %in% c('LICHEN', 'MOSS', 'LITTER', 'LOG', 'STUMP', 'BARE',
                         'BUTT', 'TREEB', 'TREE', 'HERB'))

tidy_tree2010 <- raw_Tree |>
  filter(year >= 2010,
         plot %in% c(11, 34, 27, 28, 26, 4, 25, 13)) |>
  mutate(hillslope = case_when(plot %in% c(11, 34, 27, 28) ~ 'High Storage',
                               plot %in% c(26, 4, 25, 13) ~ 'Low Storage'))

# Export: 2010-present, co-located plots only
write_csv(tidy_tree2010, file.path(OD, 'WS10_TreeSurvey_CoLocated_2010-Present_processed.csv'))
write_csv(tidy_grndcover2010, file.path(OD, 'WS10_GroundCover_CoLocated_2010-Present_processed.csv'))

rm(tidy_tree2010, tidy_grndcover2010)

## ----------------------------------------- ##
#             Co-located, 1989-present ----

# Co-located plots, restricted to 1989-present
tidy_grndcover1989 <- raw_GCE |>
  filter(year >= 1989,
         plot %in% c(11, 34, 27, 28, 26, 4, 25, 13)) |>
  mutate(hillslope = case_when(plot %in% c(11, 34, 27, 28) ~ 'High Storage',
                               plot %in% c(26, 4, 25, 13) ~ 'Low Storage')) |>
  select(year, plot, hillslope, species) |>
  filter(!species %in% c('LICHEN', 'MOSS', 'LITTER', 'LOG', 'STUMP', 'BARE',
                         'BUTT', 'TREEB', 'TREE', 'HERB'))

tidy_tree1989 <- raw_Tree |>
  filter(year >= 1989,
         plot %in% c(11, 34, 27, 28, 26, 4, 25, 13)) |>
  mutate(hillslope = case_when(plot %in% c(11, 34, 27, 28) ~ 'High Storage',
                               plot %in% c(26, 4, 25, 13) ~ 'Low Storage'))

# Export: 1989-present, co-located plots only
write_csv(tidy_tree1989, file.path(OD, 'WS10_TreeSurvey_CoLocated_1989-Present_processed.csv'))
write_csv(tidy_grndcover1989, file.path(OD, 'WS10_GroundCover_CoLocated_1989-Present_processed.csv'))

rm(tidy_tree1989, tidy_grndcover1989)
