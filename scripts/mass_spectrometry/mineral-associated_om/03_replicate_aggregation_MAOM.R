# -----------------------------------------------------------------------------#
# 03_replicate_aggregation_MAOM.R
#
# Purpose
#   Filters MAOM (FTICR-MS) replicates down to the sample level:
#     - Instrument replicates (mRep) are first collapsed to technical
#       replicates (Rep) using a 2/3 consensus rule.
#     - Technical replicates (Rep) are then collapsed to the sample
#       level (Site x Sampling.Event x Core.Section) using the same
#       consensus rule.
#
# Inputs
#   2023_FTICR-MS_mineral-associated_om_PA_MADfiltered.csv
#   2023_FTICR-MS_mineral-associated_om_ChemProperties_MADfiltered.csv
#   (both written by 02_QAQC_MAOM.R)
#
# Outputs (written to OD)
#   - 2023_FTICR-MS_mineral-associated_om_aggregated.csv                
#   - 2023_FTICR-MS_mineral-associated_om_ChemProperties_aggregated.csv  
#   - 2023_FTICR-MS_mineral-associated_om_SampleAvg_aggregated.csv       
#
# Notes
#   - Consensus rule for collapsing replicates to a single presence/absence
#     call: 3 replicates -> detected in >= 2; 2 replicates -> detected in
#     both; 1 replicate -> detected in that one.
#   - Both June and December MAOM samples have 3 technical replicates
#     (Rep), each with 3 instrument replicates (mRep) -- see manuscript
#     Methods 2.4.
#
# Created:  2026-01-26
# Author:   X. Takver
# -----------------------------------------------------------------------------#

## ----------------------------------------- ##
#             Import Packages ----
## ----------------------------------------- ##

# Data Processing
librarian::shelf(tidyverse, broom, here)

## ----------------------------------------- ##
#             Source Functions ----
## ----------------------------------------- ##

# Source helper functions
source(here('R', 'functions', 'paths_and_directories.R'))
source(here('R', 'functions', 'fticr_sample_means_function.R'))

## ----------------------------------------- ##
#             Directory Setup ----
## ----------------------------------------- ##

# Base data directory and output locations
DD <- here('data', 'mass_spectrometry', 'processed', 'mineral-associated_om', 'QAQC')
OD <- file.path(dirname(DD), 'replicate_aggregation')

# Ensure output directories exist (idempotent)
ensure_dir(OD)

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

# MAOM inputs are QAQCd at the replicate level (post-MAD filtering).
hcl_pa <- read_csv(file.path(DD, '2023_FTICR-MS_mineral-associated_om_PA_MADfiltered.csv'))
hcl_chem <- read_csv(file.path(DD, '2023_FTICR-MS_mineral-associated_om_ChemProperties_MADfiltered.csv'))

## ----------------------------------------- ##
#              MAOM Aggregated ----
## ----------------------------------------- ##

## ----------------------------------------- ##
##              Prepare Data ----

# Convert wide presence/absence table to long form and parse sample metadata
# directly from the standardized identifier (expected format:
# <Fraction>_<HS|LS><Pit>_<June|December>_<Top|Btm>_R<Rep>_mR<mRep>).
hcl_long <- hcl_pa |> 
  pivot_longer(cols = -molecular_formula,
               names_to = 'Sample.ID',
               values_to = 'PA') |> 
  mutate(Site.ID = str_extract(Sample.ID, 'HS\\d+|LS\\d+'),
         Sampling.Event = str_extract(Sample.ID, 'June|December'),
         Core.Section = str_extract(Sample.ID, 'Top|Btm'),
         Rep = str_extract(Sample.ID, '(?<=_R)\\d+'),
         mRep = str_extract(Sample.ID, '(?<=_mR)\\d+'))

## ----------------------------------------- ##
##             Condense mRep ----

# Collapse instrument replicates (mRep) within each technical replicate (Rep).
# A compound is retained if detected in >= 2 instrument replicates (2/3 rule when three are present).
hcl_mRep_aggregated <- hcl_long |> 
  group_by(molecular_formula, Site.ID, Sampling.Event, Core.Section, Rep) |> 
  summarise(PA = if_else(sum(PA) >= 2, 1, 0), .groups = 'drop')

## ----------------------------------------- ##
##             Condense Rep ----

# Summarise replicate counts (diagnostic; used to confirm the expected replicate structure).
hcl_sample_summary <- hcl_mRep_aggregated |>
  group_by(molecular_formula, Site.ID, Sampling.Event, Core.Section) |>
  summarise(reps = n(), .groups = 'drop')

# Collapse technical replicates (Rep) to sample level using the same consensus rule:
#   - 3 reps: require detection in at least 2
#   - 2 reps: require detection in both
#   - 1 rep: require detection in that replicate
hcl_aggregated <- hcl_long |> 
  group_by(molecular_formula, Site.ID, Sampling.Event, Core.Section) |> 
  summarise(Total_Reps = n_distinct(Rep),
            Total_PA = sum(PA),
            PA = case_when(Total_Reps == 3 & Total_PA >= 2 ~ 1,
                           Total_Reps == 2 & Total_PA >= 2 ~ 1,
                           Total_Reps == 1 & Total_PA >= 1 ~ 1,
                           TRUE ~ 0),
            .groups = 'drop')

## ----------------------------------------- ##
##             Reconstitute Df ----

# Reconstitute a sample-level wide table and remove formulas with zero detections across all samples.
df_hcl_aggregated <- hcl_aggregated |> 
  unite('Sample.ID', Site.ID, Sampling.Event, Core.Section, sep = '_') |> 
  select(-Total_Reps, -Total_PA) |> 
  pivot_wider(names_from = Sample.ID,
              values_from = PA) |> 
  filter(rowSums(across(where(is.numeric))) > 0)

# Cleanup workspace
rm(hcl_long, hcl_mRep_aggregated, hcl_pa, hcl_sample_summary)

## ----------------------------------------- ##
##      Update Molecular Properties ----

# Subset formula-level properties to the formulas retained after replicate aggregation.
df_hcl_chem <- hcl_chem |> 
  filter(molecular_formula %in% df_hcl_aggregated$molecular_formula)

## ----------------------------------------- ##
##     Average Molecular Properties ----  

# Join sample-level detections to molecular properties to compute per-sample means and composition metrics.
df_long <- df_hcl_aggregated |>
  pivot_longer(cols = -molecular_formula,
               names_to = 'Sample',
               values_to = 'w') |>
  left_join(df_hcl_chem, by = 'molecular_formula') |>
  mutate(present = w != 0)

# Recalculate average properties based on sample-level molecular formulas.
# el_forms is limited to the P-free groups, since phosphorus is not
# present in the published MAOM data.
df_hcl_AvgProp <- 
  calc_fticr_means(df_long,
                   el_forms = c('CHO', 'CHON', 'CHOS', 'CHONS'))

## ----------------------------------------- ##
##          Export Data Frames ----  

write_csv(df_hcl_aggregated, file.path(OD, '2023_FTICR-MS_mineral-associated_om_aggregated.csv'))
write_csv(df_hcl_chem, file.path(OD, '2023_FTICR-MS_mineral-associated_om_ChemProperties_aggregated.csv'))
write_csv(df_hcl_AvgProp, file.path(OD, '2023_FTICR-MS_mineral-associated_om_SampleAvg_aggregated.csv'))

# Cleanup Workspace
rm(df_hcl_AvgProp, df_hcl_chem, df_long, df_hcl_aggregated, 
   hcl_aggregated, hcl_chem)
