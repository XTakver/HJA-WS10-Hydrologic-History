# -----------------------------------------------------------------------------#
# 06_statistics_MAOM_NMDS.R
#
# Purpose
#   Run community-level statistics on MAOM FTICR-MS presence/absence data using a
#   rarefied Jaccard dissimilarity matrix:
#     1) PERMANOVA (adonis2) to test compositional differences among factors
#     2) Beta-dispersion (betadisper + permutest) to test homogeneity of multivariate
#        dispersion (a key assumption/diagnostic for PERMANOVA)
#     3) Season-stratified versions of both tests to evaluate hillslope effects
#        within June and within December separately
#
# Inputs
#   - Rarefied dissimilarity matrix (RDS):
#       mineral-associated_om_Jaccard_avgdist_rare2195_it750.rds
#   - Sample metadata / molecular-property summary (CSV):
#       2023_FTICR-MS_mineral-associated_om_Sample_Properties_MADfiltered.csv
#   (both written by 05_rarefied_NMDS_MAOM.R / 02_QAQC_MAOM.R, respectively)
#
# Outputs (written to SD)
#   - Full factorial PERMANOVA:
#       2023_FTICR-MS_MAOM_PERMANOVA_FullFactorial.csv
#       2023_FTICR-MS_MAOM_PERMANOVA_FullFactorial_ByTerm.csv
#   - Season-stratified PERMANOVA summaries:
#       MAOM_PERMANOVA_bySeason.csv
#   - Beta-dispersion permutation test summaries:
#       2023_FTICR-MS_MAOM_betadisper_permutest_summary.csv
#       2023_FTICR-MS_MAOM_betadisper_hillslope_withinSeason_summary.csv
#
# Notes
#   - The distance matrix includes sample labels; metadata are reordered to match
#     dist label order prior to any subsetting.
#   - set.seed() is applied before permutation-based tests to make p-values
#     reproducible across runs.
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
DD <- here('data', 'mass_spectrometry', 'processed', 'mineral-associated_om')
SD <- here('output', 'statistics', 'mass_spectrometry')

# Ensure output directories exist (idempotent)
ensure_dir(SD)

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

# Rarefied Jaccard dissimilarity matrix (avgdist output saved as RDS)
hcl_distMat <- readRDS(file.path(DD, 'NMDS', 'mineral-associated_om_Jaccard_avgdist_rare2195_it750.rds'))

# Sample metadata (must include Sample + factor columns used in models)
hcl_sampleProp <- read_csv(file.path(DD, 'QAQC', '2023_FTICR-MS_mineral-associated_om_Sample_Properties_MADfiltered.csv'))

## ----------------------------------------- ##
#             Data Preparation ----
## ----------------------------------------- ##

# Reorder metadata rows to match the dist label order.
# This ensures all downstream models use aligned sample order.
hcl_metadata <- hcl_sampleProp[match(attributes(hcl_distMat)$Labels,
                                       hcl_sampleProp$Sample), ]

# Confirm that all distance labels were matched to metadata rows
if (any(is.na(hcl_metadata$Sample))) {
  stop('Some dist labels were not matched to metadata$Sample. Check sample naming and inputs.')
}

# Create season masks (used for within-season subsetting)
hcl_june <- hcl_metadata$Sampling.Event == 'June'
hcl_dec  <- hcl_metadata$Sampling.Event == 'December'

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
hcl_PN_global <-
  adonis2(hcl_distMat ~ Hillslope * Sampling.Event * Core.Section,
          data = hcl_metadata)

# Set seed for reproducible permutation p-values
set.seed(55)

# By-term PERMANOVA: sequential tests of each term in the model
hcl_PN_term <-
  adonis2(hcl_distMat ~ Hillslope * Sampling.Event * Core.Section,
          data = hcl_metadata,
          by = 'terms')

# Convert to data frames for export
df_hcl_PN_global <- as.data.frame(hcl_PN_global)
df_hcl_PN_term   <- as.data.frame(hcl_PN_term)

## ----------------------------------------- ##
#            Export PERMANOVA ----
## ----------------------------------------- ##

write_csv(df_hcl_PN_global, file.path(SD, '2023_FTICR-MS_MAOM_PERMANOVA_FullFactorial.csv'))
write_csv(df_hcl_PN_term,   file.path(SD, '2023_FTICR-MS_MAOM_PERMANOVA_FullFactorial_ByTerm.csv'))

# Cleanup workspace
rm(hcl_PN_global, hcl_PN_term, df_hcl_PN_global, df_hcl_PN_term)

## ----------------------------------------- ##
##      MAOM PERMANOVA by season ----
## ----------------------------------------- ##

# Subset dist objects by season (dist must be subset as a matrix first)
hcl_diss_matrix_full <- as.matrix(hcl_distMat)

hcl_diss_june <- hcl_diss_matrix_full[hcl_june, hcl_june, drop = FALSE]
hcl_diss_dec  <- hcl_diss_matrix_full[hcl_dec,  hcl_dec,  drop = FALSE]

hcl_diss_june_dist <- as.dist(hcl_diss_june)
hcl_diss_dec_dist  <- as.dist(hcl_diss_dec)

# Split metadata by season (droplevels removes unused factor levels after subsetting)
hcl_meta_june <- droplevels(hcl_metadata[hcl_june, , drop = FALSE])
hcl_meta_dec  <- droplevels(hcl_metadata[hcl_dec,  , drop = FALSE])

# Run within-season PERMANOVA models (explicit calls for readability)
df_hcl_permanova_by_season <- bind_rows(
  run_permanova_season(hcl_diss_june_dist, hcl_meta_june, 'June'),
  run_permanova_season(hcl_diss_dec_dist,  hcl_meta_dec,  'December'))

## ----------------------------------------- ##
#       Export Seasonal PERMANOVA ----
## ----------------------------------------- ##

write_csv(df_hcl_permanova_by_season, file.path(SD, 'MAOM_PERMANOVA_bySeason.csv'))

## ----------------------------------------- ##
#               Beta Dispersion ----
## ----------------------------------------- ##

# Beta-dispersion tests whether within-group dispersions differ among factor levels.
# This is commonly used to diagnose whether PERMANOVA differences could reflect
# dispersion differences rather than centroid (composition) shifts.

# Global beta-dispersion tests (explicit calls for readability)

bd_hillslope <- betadisper(hcl_distMat, group = hcl_metadata$Hillslope)
set.seed(55)
pt_hillslope <- permutest(bd_hillslope, permutations = 999)

bd_season <- betadisper(hcl_distMat, group = hcl_metadata$Sampling.Event)
set.seed(55)
pt_season <- permutest(bd_season, permutations = 999)

bd_depth <- betadisper(hcl_distMat, group = hcl_metadata$Core.Section)
set.seed(55)
pt_depth <- permutest(bd_depth, permutations = 999)

bd_hs_x_season <- betadisper(hcl_distMat,
                             group = interaction(hcl_metadata$Hillslope, hcl_metadata$Sampling.Event))
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

write_csv(betadisp_results, file.path(SD, '2023_FTICR-MS_MAOM_betadisper_permutest_summary.csv'))

## ----------------------------------------- ##
##      Beta dispersion within season ----
## ----------------------------------------- ##

# Test hillslope dispersion within each season separately (June vs December).
# Distances must be subset as matrices prior to conversion back to dist objects.

bd_june <- betadisper(hcl_diss_june_dist, group = factor(hcl_metadata$Hillslope[hcl_june]))
set.seed(55)
pt_june <- permutest(bd_june, permutations = 999)

bd_dec <- betadisper(hcl_diss_dec_dist, group = factor(hcl_metadata$Hillslope[hcl_dec]))
set.seed(55)
pt_dec <- permutest(bd_dec, permutations = 999)

hcl_betadisp_by_season <- bind_rows(
  tidy_permutest(pt_june, grouping = 'Hillslope', season = 'June'),
  tidy_permutest(pt_dec, grouping = 'Hillslope', season = 'December'))

## ----------------------------------------- ##
#        Export Seasonal BetaDisp ----
## ----------------------------------------- ##

write_csv(hcl_betadisp_by_season,
          file.path(SD, '2023_FTICR-MS_MAOM_betadisper_hillslope_withinSeason_summary.csv'))

# Cleanup workspace
rm(hcl_diss_matrix_full, hcl_diss_june, hcl_diss_dec,
   hcl_diss_june_dist, hcl_diss_dec_dist, hcl_meta_june, hcl_meta_dec,
   df_hcl_permanova_by_season, betadisp_results, hcl_betadisp_by_season,
   bd_hillslope, bd_season, bd_depth, bd_hs_x_season,
   pt_hillslope, pt_season, pt_depth, pt_hs_x_season, bd_june, bd_dec, pt_june, pt_dec)
