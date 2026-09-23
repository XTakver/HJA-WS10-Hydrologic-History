# -----------------------------------------------------------------------------#
# 06_statistics_WEOM_NMDS.R
#
# Purpose
#   Run community-level statistics on WEOM FTICR-MS presence/absence data using a
#   rarefied Jaccard dissimilarity matrix:
#     1) PERMANOVA (adonis2) to test compositional differences among factors
#     2) Beta-dispersion (betadisper + permutest) to test homogeneity of multivariate
#        dispersion (a key assumption/diagnostic for PERMANOVA)
#     3) Season-stratified versions of both tests to evaluate hillslope effects
#        within June and within December separately
#
# Inputs
#   - Rarefied dissimilarity matrix (RDS):
#       WEOM_Jaccard_avgdist_rare3739_it750.rds
#   - Sample metadata / molecular-property summary (CSV):
#       2023_FTICR-MS_WEOM_Sample_Properties_MADfiltered.csv
#   (both written by 05_rarefied_NMDS_WEOM.R / 02_QAQC_WEOM.R, respectively)
#
# Outputs (written to SD)
#   - Full factorial PERMANOVA:
#       2023_FTICR-MS_WEOM_PERMANOVA_FullFactorial.csv
#       2023_FTICR-MS_WEOM_PERMANOVA_FullFactorial_ByTerm.csv
#   - Season-stratified PERMANOVA summaries:
#       WEOM_PERMANOVA_bySeason.csv
#   - Beta-dispersion permutation test summaries:
#       2023_FTICR-MS_WEOM_betadisper_permutest_summary.csv
#       2023_FTICR-MS_WEOM_betadisper_hillslope_withinSeason_summary.csv
#
# Notes
#   - The distance matrix includes sample labels; metadata are reordered to match
#     dist label order prior to any subsetting.
#   - set.seed() is applied before permutation-based tests to make p-values
#     reproducible across runs.
#   - Models reference Hillslope, Sampling.Event, and Core.Section as
#     published in 2023_FTICR-MS_WEOM_Sample_Properties_MADfiltered.csv
#     (set in 01_processing_WEOM.R).
#
# Created:  2026-01-29
# Author:   X. Takver
# -----------------------------------------------------------------------------#


## ----------------------------------------- ##
#             Import Packages ----
## ----------------------------------------- ##

# Data processing + ordination/statistical utilities
librarian::shelf(tidyverse, here, vegan)

## ----------------------------------------- ##
#             Source Functions ----
## ----------------------------------------- ##

# Directory helpers (ensure_dir)
source(here('R', 'functions', 'paths_and_directories.R'))

## ----------------------------------------- ##
#          Directory Navigation ----
## ----------------------------------------- ##

# DD: FTICR-MS processed directory containing NMDS outputs and QAQC metadata
# SD: directory for statistical outputs
DD <- here('data', 'mass_spectrometry', 'processed', 'water-extractable_om')
SD <- here('output', 'statistics', 'mass_spectrometry')

# Ensure output directories exist (idempotent)
ensure_dir(SD)

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

# Rarefied Jaccard dissimilarity matrix (avgdist output saved as RDS)
weom_distMat <- readRDS(file.path(DD, 'NMDS', 'WEOM_Jaccard_avgdist_rare3739_it750.rds'))

# Sample metadata (must include Sample + factor columns used in models)
weom_sampleProp <- read_csv(file.path(DD, 'QAQC', '2023_FTICR-MS_WEOM_Sample_Properties_MADfiltered.csv'))

## ----------------------------------------- ##
#             Data Preparation ----
## ----------------------------------------- ##

# Reorder metadata rows to match the dist label order.
# This ensures all downstream models use aligned sample order.
weom_metadata <- weom_sampleProp[match(attributes(weom_distMat)$Labels,
                                       weom_sampleProp$Sample), ]

# Confirm that all distance labels were matched to metadata rows
if (any(is.na(weom_metadata$Sample))) {
  stop('Some dist labels were not matched to metadata$Sample. Check sample naming and inputs.')
}

# Create season masks (used for within-season subsetting)
weom_june <- weom_metadata$Sampling.Event == 'June'
weom_dec  <- weom_metadata$Sampling.Event == 'December'

## ----------------------------------------- ##
#            Helper Functions ----
## ----------------------------------------- ##

# Standardize adonis2 output tables for export
tidy_adonis2 <- function(x, season, model) {
  as.data.frame(x) |>
    rownames_to_column('term') |>
    mutate(season = season,
           model = model,
           .before = 1)
}

# Run a consistent set of PERMANOVA models for a given season
run_permanova_season <- function(d, meta, season_name) {
  
  # Set seed for reproducible permutation p-values
  set.seed(55)
  pn_full <- adonis2(d ~ Hillslope + Core.Section, data = meta)
  
  # Set seed for reproducible permutation p-values
  set.seed(55)
  pn_terms <- adonis2(d ~ Hillslope + Core.Section, data = meta, by = 'terms')
  
  # Set seed for reproducible permutation p-values
  set.seed(55)
  pn_hs <- adonis2(d ~ Hillslope, data = meta)
  
  bind_rows(tidy_adonis2(pn_full,  season = season_name, model = 'Hillslope_plus_Depth'),
            tidy_adonis2(pn_terms, season = season_name, model = 'Hillslope_plus_Depth_by_terms'),
            tidy_adonis2(pn_hs,    season = season_name, model = 'Hillslope_only'))
}

# Extract a consistent summary row from permutest(betadisper())
tidy_permutest <- function(pt, grouping, season = NA_character_) {
  
  tab <- as.data.frame(pt$tab)
  tab$term <- rownames(tab)
  rownames(tab) <- NULL
  
  out <- tab |>
    filter(term %in% c('Groups', 'group', 'Group')) |>
    mutate(test = 'betadisper_permutest',
           grouping = grouping,
           season = season,
           permutations = 999)
  
  # Fallback if the table uses an unexpected row label
  if (nrow(out) == 0) {
    out <- tab[1, , drop = FALSE] |>
      mutate(test = 'betadisper_permutest',
             grouping = grouping,
             season = season,
             permutations = 999)
  }
  
  out
}

## ----------------------------------------- ##
#                PERMANOVA ----
## ----------------------------------------- ##

# Set seed for reproducible permutation p-values
set.seed(55)

# Full factorial PERMANOVA: tests main effects + interactions
weom_PN_global <-
  adonis2(weom_distMat ~ Hillslope * Sampling.Event * Core.Section,
          data = weom_metadata)

# Set seed for reproducible permutation p-values
set.seed(55)

# By-term PERMANOVA: sequential tests of each term in the model
weom_PN_term <-
  adonis2(weom_distMat ~ Hillslope * Sampling.Event * Core.Section,
          data = weom_metadata,
          by = 'terms')

# Convert to data frames for export
df_weom_PN_global <- as.data.frame(weom_PN_global)
df_weom_PN_term   <- as.data.frame(weom_PN_term)

## ----------------------------------------- ##
#            Export PERMANOVA ----
## ----------------------------------------- ##

write_csv(df_weom_PN_global, file.path(SD, '2023_FTICR-MS_WEOM_PERMANOVA_FullFactorial.csv'))
write_csv(df_weom_PN_term,   file.path(SD, '2023_FTICR-MS_WEOM_PERMANOVA_FullFactorial_ByTerm.csv'))

# Cleanup workspace
rm(weom_PN_global, weom_PN_term, df_weom_PN_global, df_weom_PN_term)

## ----------------------------------------- ##
##      WEOM PERMANOVA by season ----
## ----------------------------------------- ##

# Subset dist objects by season (dist must be subset as a matrix first)
weom_diss_matrix_full <- as.matrix(weom_distMat)

weom_diss_june <- weom_diss_matrix_full[weom_june, weom_june, drop = FALSE]
weom_diss_dec  <- weom_diss_matrix_full[weom_dec,  weom_dec,  drop = FALSE]

weom_diss_june_dist <- as.dist(weom_diss_june)
weom_diss_dec_dist  <- as.dist(weom_diss_dec)

# Split metadata by season (droplevels removes unused factor levels after subsetting)
weom_meta_june <- droplevels(weom_metadata[weom_june, , drop = FALSE])
weom_meta_dec  <- droplevels(weom_metadata[weom_dec,  , drop = FALSE])

# Run within-season PERMANOVA models (explicit calls for readability)
df_weom_permanova_by_season <- bind_rows(
  run_permanova_season(weom_diss_june_dist, weom_meta_june, 'June'),
  run_permanova_season(weom_diss_dec_dist,  weom_meta_dec,  'December'))

## ----------------------------------------- ##
#       Export Seasonal PERMANOVA ----
## ----------------------------------------- ##

write_csv(df_weom_permanova_by_season, file.path(SD, 'WEOM_PERMANOVA_bySeason.csv'))

## ----------------------------------------- ##
#               Beta Dispersion ----
## ----------------------------------------- ##

# Beta-dispersion tests whether within-group dispersions differ among factor levels.
# This is commonly used to diagnose whether PERMANOVA differences could reflect
# dispersion differences rather than centroid (composition) shifts.

# Global beta-dispersion tests (explicit calls for readability)

bd_hillslope <- betadisper(weom_distMat, group = weom_metadata$Hillslope)
set.seed(55)
pt_hillslope <- permutest(bd_hillslope, permutations = 999)

bd_season <- betadisper(weom_distMat, group = weom_metadata$Sampling.Event)
set.seed(55)
pt_season <- permutest(bd_season, permutations = 999)

bd_depth <- betadisper(weom_distMat, group = weom_metadata$Core.Section)
set.seed(55)
pt_depth <- permutest(bd_depth, permutations = 999)

bd_hs_x_season <- betadisper(weom_distMat,
                             group = interaction(weom_metadata$Hillslope, weom_metadata$Sampling.Event))
set.seed(55)
pt_hs_x_season <- permutest(bd_hs_x_season, permutations = 999)

# Combine summaries into one export table
betadisp_results <- bind_rows(
  tidy_permutest(pt_hillslope, grouping = 'Hillslope'),
  tidy_permutest(pt_season, grouping = 'Sampling.Event'),
  tidy_permutest(pt_depth, grouping = 'Core.Section'),
  tidy_permutest(pt_hs_x_season, grouping = 'Hillslope_x_Sampling.Event'))

## ----------------------------------------- ##
#            Export BetaDisp ----
## ----------------------------------------- ##

write_csv(betadisp_results, file.path(SD, '2023_FTICR-MS_WEOM_betadisper_permutest_summary.csv'))

## ----------------------------------------- ##
##      Beta dispersion within season ----
## ----------------------------------------- ##

# Test hillslope dispersion within each season separately (June vs December).
# Distances must be subset as matrices prior to conversion back to dist objects.

bd_june <- betadisper(weom_diss_june_dist, group = factor(weom_metadata$Hillslope[weom_june]))
set.seed(55)
pt_june <- permutest(bd_june, permutations = 999)

bd_dec <- betadisper(weom_diss_dec_dist, group = factor(weom_metadata$Hillslope[weom_dec]))
set.seed(55)
pt_dec <- permutest(bd_dec, permutations = 999)

weom_betadisp_by_season <- bind_rows(
  tidy_permutest(pt_june, grouping = 'Hillslope', season = 'June'),
  tidy_permutest(pt_dec, grouping = 'Hillslope', season = 'December'))

## ----------------------------------------- ##
#        Export Seasonal BetaDisp ----
## ----------------------------------------- ##

write_csv(weom_betadisp_by_season,
          file.path(SD, '2023_FTICR-MS_WEOM_betadisper_hillslope_withinSeason_summary.csv'))

# Cleanup workspace
rm(weom_diss_matrix_full, weom_diss_june, weom_diss_dec,
   weom_diss_june_dist, weom_diss_dec_dist, weom_meta_june, weom_meta_dec,
   df_weom_permanova_by_season, betadisp_results, weom_betadisp_by_season,
   bd_hillslope, bd_season, bd_depth, bd_hs_x_season,
   pt_hillslope, pt_season, pt_depth, pt_hs_x_season, bd_june, bd_dec, pt_june, pt_dec)
