# -----------------------------------------------------------------------------#
# 05_rarefied_NMDS_WEOM.R
#
# Purpose
#   Pulls in MAD-filtered WEOM presence/absence data and performs rarefied NMDS.
#   Computes a rarefied Jaccard dissimilarity matrix (avgdist) and fits NMDS ordination (vegan::metaMDS).
#   Exports NMDS site scores and envfit vector scores for downstream plotting/statistics.
#
# Inputs (read from DD)
#   - 2023_FTICR-MS_WEOM_PA_MADfiltered.csv
#   - 2023_FTICR-MS_WEOM_Sample_Properties_MADfiltered.csv
#   (both written by 02_QAQC_WEOM.R)
#
# Outputs (written to OD)
#   - 2023_FTICR-MS_NMDS_Scores_WEOM.csv
#   - 2023_FTICR-MS_NMDS_SampleProp_Scores_WEOM.csv
#   - WEOM_Jaccard_avgdist_rare3739_it750.rds
#
# Key settings
#   - Distance: Jaccard (binary)
#   - Rarefaction: sample = 3739, iterations = 750
#   - NMDS: k = 4, trymax = 100
#   - Reproducibility: set.seed(55) used prior to avgdist/metaMDS
#
# Created:  2026-01-26
# Author:   X. Takver
# -----------------------------------------------------------------------------#

## ----------------------------------------- ##
#             Import Packages ----
## ----------------------------------------- ##

# Data processing and ordination (vegan) utilities
librarian::shelf(tidyverse, vegan, here)

## ----------------------------------------- ##
#             Source Functions ----
## ----------------------------------------- ##

# Directory helpers (ensure_dir)
source(here('R', 'functions', 'paths_and_directories.R'))

## ----------------------------------------- ##
#          Directory Creation ----
## ----------------------------------------- ##

# DD: input directory containing MAD-filtered QAQC outputs
# OD: output directory for NMDS results and intermediate objects
DD <- here('data', 'mass_spectrometry', 'processed', 'water-extractable_om', 'QAQC')
OD <- file.path(dirname(DD), 'NMDS')

# Ensure output directories exist (idempotent)
ensure_dir(OD)

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

# Presence/absence matrix (formulas x samples) and per-sample molecular properties
weom_pa <- read_csv(file.path(DD, '2023_FTICR-MS_WEOM_PA_MADfiltered.csv'))
weom_sampleProp <- read_csv(file.path(DD, '2023_FTICR-MS_WEOM_Sample_Properties_MADfiltered.csv'))

## ----------------------------------------- ##
#             Compound Summary ----
## ----------------------------------------- ##

# Summarize detected compound counts per sample to set rarefaction bounds
weom_summary <- weom_pa |> 
  pivot_longer(cols = -molecular_formula,
               names_to = 'Site.ID',
               values_to = 'PA') |> 
  group_by(Site.ID) |> 
  summarise(Total.Compounds = sum(PA, na.rm = TRUE), .groups = 'drop')

min(weom_summary$Total.Compounds)
max(weom_summary$Total.Compounds)

# Minimum peaks detected: 3,739
# Maximum peaks detected: 12,376

## ----------------------------------------- ##
#             Distance Matrix ----
## ----------------------------------------- ##

# Convert to samples x formulas matrix required by avgdist/metaMDS
weom_transposed <- weom_pa |> 
  column_to_rownames('molecular_formula') |> 
  t()

# Drop samples with zero detected formulas (cannot compute distances/ordination)
empty_rows <- rowSums(weom_transposed) == 0
nmds_matrix_clean <- weom_transposed[!empty_rows, ]

# Set seed before any stochastic steps (avgdist rarefaction + NMDS starts)
set.seed(55)

# Rarefied Jaccard dissimilarities for binary data:
#   - Jaccard uses shared presences only (appropriate for PA data)
#   - avgdist rarefies to a fixed detection count (sample) across iterations
dissimilarity_matrix <- avgdist(nmds_matrix_clean, 
                                sample = 3739,
                                iterations = 750,
                                dmethod = 'jaccard',
                                binary = TRUE)

# Set seed for reproducible permutation-based metaMDS
set.seed(55)
# Screen candidate NMDS dimensions (k) by stress to choose a stable solution
stress_values <- numeric(4)  

for(i in 2:5){
  stress_values[i - 1] <- metaMDS(dissimilarity_matrix, k=i, trace = FALSE)$stress * 100
}

print(stress_values)

# Set seed for reproducible permutation-based metaMDS
set.seed(55)
# Fit final NMDS solution (k can be adjusted based on the stress screen above)
nmds_result <- metaMDS(dissimilarity_matrix, k = 4, trymax = 100)

# Inspect ordination diagnostics
print(nmds_result)

# Ordination stability check:
#   compare two runs with different random starts via Procrustes
set.seed(1); n1 <- metaMDS(dissimilarity_matrix, k = 4, trymax = 100)
set.seed(2); n2 <- metaMDS(dissimilarity_matrix, k = 4, trymax = 100)
pro <- procrustes(scores(n1), scores(n2))
pro_sum <- summary(pro)

# Extract NMDS site scores and parse sample metadata directly from the
# standardized identifier (expected format:
# <Fraction>_<HS|LS><Pit>_<June|December>_<Top|Btm>_R<Rep>_mR<mRep>).
weom_nmds_scores <- as.data.frame(scores(nmds_result)) |>
  rownames_to_column('Sample') |>
  mutate(Site = str_extract(Sample, 'HS\\d+|LS\\d+'),
         Sampling.Event = str_extract(Sample, 'June|December'),
         Section = str_extract(Sample, 'Top|Btm'),
         Rep = str_extract(Sample, '(?<=_R)\\d+'),
         mRep = str_extract(Sample, '(?<=_mR)\\d+')) |>
  mutate(across(c(Site, Sampling.Event, Section, Rep, mRep), as_factor)) |>
  mutate(Sample.ID = as_factor(str_c(Site, Section, sep = '_')),
         .after = mRep) |>
  mutate(Site = case_when(str_detect(Site, 'HS') ~ 'HS',
                          str_detect(Site, 'LS') ~ 'LS'),
         Section = case_when(Section == 'Top' ~ '0-10cm',
                             Section == 'Btm' ~ '20-30cm'),
         Site.ID = as_factor(paste(Site, Section, sep = ' ')),
         .after = Sample.ID) |>
  mutate(Site.ID = fct_relevel(Site.ID,
                               'HS 0-10cm',
                               'HS 20-30cm',
                               'LS 0-10cm',
                               'LS 20-30cm'))

## ----------------------------------------- ##
##            Calculate EnvFit ----
## ----------------------------------------- ##

# Fit molecular property vectors to the NMDS ordination using permutation tests.
# Align sample properties to the NMDS site-score order and drop constant variables (envfit requires variation).

site_order <- rownames(scores(nmds_result, display = 'sites'))

nmds_weom_sampleProp <- weom_sampleProp |>
  filter(Sample %in% site_order) |>
  mutate(Sample = factor(Sample, levels = site_order)) |>
  arrange(Sample) |>
  select(mass:CondHydrocarbon_percent) |>
  select(where(~ !is.numeric(.) | any(. != 0, na.rm = TRUE)))

# Set seed for reproducible permutation-based p-values.
set.seed(55)

fit_mol_prop <- envfit(nmds_result, nmds_weom_sampleProp, permutations = 999, choices = 1:4)

# Extract fitted vectors and statistics for plotting and reporting.
weom_sampleProp_nmdsScores <- as.data.frame(scores(fit_mol_prop, 'vectors')) |>
  rownames_to_column('mol_properties') |>
  mutate(mol_properties = str_replace(mol_properties, '_percent', ' %')) |>
  mutate(p.value = format(fit_mol_prop$vector$pvals, digits = 10, nsmall = 10),
         r2 = fit_mol_prop$vectors$r)

## ----------------------------------------- ##
##               Export Data ----
## ----------------------------------------- ##

# NMDS site scores (samples) and envfit vectors for downstream plotting and statistics
write_csv(weom_nmds_scores, file.path(OD, '2023_FTICR-MS_NMDS_Scores_WEOM.csv'))
write_csv(weom_sampleProp_nmdsScores, file.path(OD, '2023_FTICR-MS_NMDS_SampleProp_Scores_WEOM.csv'))

# Save rarefied dissimilarity matrix as an intermediate object for reproducibility
saveRDS(dissimilarity_matrix, file.path(OD, 'WEOM_Jaccard_avgdist_rare3739_it750.rds'),
        compress = 'xz')
