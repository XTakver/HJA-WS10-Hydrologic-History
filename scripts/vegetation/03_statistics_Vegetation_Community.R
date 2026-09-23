# ------------------------------------------------------------ #
# 03_statistics_Vegetation_Community.R
#
# Purpose
#   PERMANOVA and diversity modeling for the WS10 long-term
#   vegetation community (H.J. Andrews LTER, 1973-2016): tests
#   hillslope, year, and hillslope x year effects on community
#   composition (all plots and co-located plots), models alpha-
#   diversity metrics by year x hillslope, and identifies
#   indicator species by hillslope x decade.
#
# Inputs
#   data/vegetation/processed/diversity/Community_Matrix_All.csv
#   data/vegetation/processed/diversity/Diversity_DF_All.csv
#
# Outputs (written to SD)
#   - Veg_Permanova_noInteraction.csv
#   - Veg_Permanova_Interaction.csv
#   - Veg_Permanova_byYear.csv             (feeds 03_figure-S14B.R)
#   - coLoc_Permanova_byYear.csv           (feeds 03_figure-S14B.R)
#   - PSME_removed_Permanova_byYear.csv
#   - Indicator_Species_HillslopeEra.csv
#   - {metric}_Anova.txt / {metric}_lmer_summary.txt
#   - {metric}_emmeans_pairwise_by_year.csv
#     (metric = shannon, simpson, invsimpson, fisher, evenness, richness)
#
# Notes
#   - Diversity models: lmer(metric ~ year * hillslope + (1 | plot)),
#     Type III ANOVA, and emmeans pairwise contrasts by year (see
#     Table S3).
#   - Indicator-species analysis (multipatt) groups samples by
#     hillslope x decade ('80s/'90s/'00s/'10s').
#   - The PSME/ACCI/ANDE3 sensitivity check removes the dominant
#     species from the community matrix and re-runs the per-year
#     PERMANOVA, matching the manuscript's "sensitivity analyses
#     with dominant species removed."
#   - Community_Matrix_All.csv and Diversity_DF_All.csv are produced
#     by 01_processing_Vegetation_Diversity.R, not by this script.
#
# Created:  2025-11-23
# Author:   X. Takver
# ------------------------------------------------------------ #

## ----------------------------------------- ##
#             Package Import ----
## ----------------------------------------- ##

# Data Processing
librarian::shelf(tidyverse, broom, vegan, lme4, lmerTest, car, emmeans,
                 here, glue, indicspecies)

## ----------------------------------------- ##
#             Source Functions ----
## ----------------------------------------- ##

source(here('R', 'functions', 'paths_and_directories.R'))  # ensure_dir()

## ----------------------------------------- ##
#          Directory Navigation ----
## ----------------------------------------- ##

# Data Directory Variables
DD <- here('data', 'vegetation', 'processed', 'diversity')

# Stats Directory Variables
SD <- here('output', 'statistics', 'vegetation')

ensure_dir(SD)

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

comm_matrix <- read_csv(file.path(DD, 'Community_Matrix_All.csv'))
diversity_data <- read_csv(file.path(DD, 'Diversity_DF_All.csv'))

## ----------------------------------------- ##
#             Data Preparation ----
## ----------------------------------------- ##

## ----------------------------------------- ##
##             Dissimilarity ----

# Create a distance matrix
dist_matrix <- vegdist(comm_matrix[, -(1:3)], method = 'bray')

# Create an environmental metadata df
metadata <- comm_matrix |>
  select(year, plot, hillslope) |>
  unite('YR_PLT_HS', year, plot, hillslope, remove = FALSE)

## ----------------------------------------- ##
#           Statistical Analysis ----
## ----------------------------------------- ##

##               PERMANOVA ----

# non-interaction
perma_fixed <- adonis2(dist_matrix ~ year + hillslope,
                       data   = metadata,
                       strata = metadata$plot,
                       by = 'terms')   # respects repeated measures per plot

# Interaction
perma_interaction <- adonis2(dist_matrix ~ year * hillslope,
                             data   = metadata,
                             strata = metadata$plot,
                             by = 'terms')   # respects repeated measures per plot

# Pairwise by-year comparison
years <- sort(unique(metadata$year))

per_year_results <- map_df(years, ~ {
  yr <- .x

  # subset metadata and community to this year
  idx      <- metadata$year == yr
  meta_yr  <- metadata[idx, , drop = FALSE]
  comm_yr  <- comm_matrix[idx, -(1:3), drop = FALSE]     # species matrix for that year

  # distance matrix for this year only
  dist_yr <- vegdist(comm_yr, method = 'bray')

  # PERMANOVA: hillslope effect within this year
  ad <- adonis2(dist_yr ~ hillslope,
                data = meta_yr,
                permutations = 999,
                by = 'terms')

  hillslope_row <- ad['hillslope', ]

  tibble(
    year    = yr,
    R2      = hillslope_row$R2,
    F       = hillslope_row$F,
    p_value = hillslope_row$`Pr(>F)`
  )
})

# optional multiple-testing correction
per_year_results <- per_year_results |>
  mutate(p_adj = p.adjust(p_value, method = 'BH'))

per_year_results

# Save Permanova results
write_csv(perma_fixed, file.path(SD, 'Veg_Permanova_noInteraction.csv'))
write_csv(perma_interaction, file.path(SD, 'Veg_Permanova_Interaction.csv'))
write_csv(per_year_results, file.path(SD, 'Veg_Permanova_byYear.csv'))

##               Diversity ----

# Shannon
m_shannon <- lmer(shannon ~ year * hillslope + (1 | plot), data = diversity_data)
Anova(m_shannon, type = 'III')
summary(m_shannon)

emm_shan <- emmeans(
  m_shannon,
  ~ hillslope | year,
  at = list(year = sort(unique(diversity_data$year))))

summary(emm_shan)
pairs(emm_shan, by = 'year')

writeLines(capture.output(Anova(m_shannon, type = 'III')), file.path(SD, 'shannon_Anova.txt'))
writeLines(capture.output(summary(m_shannon)), file.path(SD, 'shannon_lmer_summary.txt'))
write_csv(as.data.frame(pairs(emm_shan, by = 'year')),
          file.path(SD, 'shannon_emmeans_pairwise_by_year.csv'))

# Simpson
m_simpson<- lmer(simpson ~ year * hillslope + (1 | plot), data = diversity_data)
Anova(m_simpson, type = 'III')
summary(m_simpson)

emm_simp <- emmeans(
  m_simpson,
  ~ hillslope | year,
  at = list(year = sort(unique(diversity_data$year))))

pairs(emm_simp, by = 'year')

writeLines(capture.output(Anova(m_simpson, type = 'III')), file.path(SD, 'simpson_Anova.txt'))
writeLines(capture.output(summary(m_simpson)), file.path(SD, 'simpson_lmer_summary.txt'))
write_csv(as.data.frame(pairs(emm_simp, by = 'year')),
          file.path(SD, 'simpson_emmeans_pairwise_by_year.csv'))

# Inv Simpson
m_invsimpson<- lmer(invsimpson ~ year * hillslope + (1 | plot), data = diversity_data)
Anova(m_invsimpson, type = 'III')
summary(m_invsimpson)

emm_invsimp <- emmeans(
  m_invsimpson,
  ~ hillslope | year,
  at = list(year = sort(unique(diversity_data$year))))

pairs(emm_invsimp, by = 'year')

writeLines(capture.output(Anova(m_invsimpson, type = 'III')), file.path(SD, 'invsimpson_Anova.txt'))
writeLines(capture.output(summary(m_invsimpson)), file.path(SD, 'invsimpson_lmer_summary.txt'))
write_csv(as.data.frame(pairs(emm_invsimp, by = 'year')),
          file.path(SD, 'invsimpson_emmeans_pairwise_by_year.csv'))

# Fisher
m_fisher<- lmer(fisher ~ year * hillslope + (1 | plot), data = diversity_data)
Anova(m_fisher, type = 'III')
summary(m_fisher)

emm_fisher <- emmeans(
  m_fisher,
  ~ hillslope | year,
  at = list(year = sort(unique(diversity_data$year))))

pairs(emm_fisher, by = 'year')

writeLines(capture.output(Anova(m_fisher, type = 'III')), file.path(SD, 'fisher_Anova.txt'))
writeLines(capture.output(summary(m_fisher)), file.path(SD, 'fisher_lmer_summary.txt'))
write_csv(as.data.frame(pairs(emm_fisher, by = 'year')),
          file.path(SD, 'fisher_emmeans_pairwise_by_year.csv'))

# Evenness
m_evenness<- lmer(evenness ~ year * hillslope + (1 | plot), data = diversity_data)
Anova(m_evenness, type = 'III')
summary(m_evenness)

emm_evenness <- emmeans(
  m_evenness,
  ~ hillslope | year,
  at = list(year = sort(unique(diversity_data$year))))

pairs(emm_evenness, by = 'year')

writeLines(capture.output(Anova(m_evenness, type = 'III')), file.path(SD, 'evenness_Anova.txt'))
writeLines(capture.output(summary(m_evenness)), file.path(SD, 'evenness_lmer_summary.txt'))
write_csv(as.data.frame(pairs(emm_evenness, by = 'year')),
          file.path(SD, 'evenness_emmeans_pairwise_by_year.csv'))

# Richness
m_richness<- lmer(richness ~ year * hillslope + (1 | plot), data = diversity_data)
Anova(m_richness, type = 'III')
summary(m_richness)

emm_richness <- emmeans(
  m_richness,
  ~ hillslope | year,
  at = list(year = sort(unique(diversity_data$year))))

pairs(emm_richness, by = 'year')

writeLines(capture.output(Anova(m_richness, type = 'III')), file.path(SD, 'richness_Anova.txt'))
writeLines(capture.output(summary(m_richness)), file.path(SD, 'richness_lmer_summary.txt'))
write_csv(as.data.frame(pairs(emm_richness, by = 'year')),
          file.path(SD, 'richness_emmeans_pairwise_by_year.csv'))

## ----------------------------------------- ##
#           Additional Analysis ----
## ----------------------------------------- ##

## ----------------------------------------- ##
##           Indicator Species ----

# Prepare data
comm_only <- comm_matrix |>
  select(-year, -plot, -hillslope)

ind_metadata <- metadata |>
  mutate(era = cut(year,
                   c(1979, 1990, 2000, 2010, 2020),
                   c('80s', '90s', '00s', '10s')))

# Group hillslope and era
ind_metadata <- ind_metadata |>
  mutate(hs_era = interaction(hillslope, era, drop = TRUE))

group_factor <- ind_metadata$hs_era

ind <- multipatt(comm_only, group_factor, func = 'IndVal.g', control = how(nperm = 999))
summary(ind)

write_csv(as_tibble(ind$sign, rownames = 'species'),
          file.path(SD, 'Indicator_Species_HillslopeEra.csv'))

## ----------------------------------------- ##
##          Co-Located Convergence ----

# Filter to co located sites
coLoc_matrix <- comm_matrix |>
  filter(plot %in% c(11, 34, 27, 28, 26, 4, 25, 13))

coLoc_metadata <- metadata |>
  filter(plot %in% c(11, 34, 27, 28, 26, 4, 25, 13))

# Pairwise comparisons
# Pairwise by-year comparison
years <- sort(unique(coLoc_metadata$year))

per_year_results <- map_df(years, ~ {
  yr <- .x

  # subset metadata and community to this year
  idx      <- coLoc_metadata$year == yr
  meta_yr  <- coLoc_metadata[idx, , drop = FALSE]
  comm_yr  <- coLoc_matrix[idx, -(1:3), drop = FALSE]     # species matrix for that year

  # distance matrix for this year only
  dist_yr <- vegdist(comm_yr, method = 'bray')

  # PERMANOVA: hillslope effect within this year
  ad <- adonis2(dist_yr ~ hillslope,
                data = meta_yr,
                permutations = 999,
                by = 'terms')

  hillslope_row <- ad['hillslope', ]

  tibble(
    year    = yr,
    R2      = hillslope_row$R2,
    F       = hillslope_row$F,
    p_value = hillslope_row$`Pr(>F)`
  )
})

# optional multiple-testing correction
per_year_results <- per_year_results |>
  mutate(p_adj = p.adjust(p_value, method = 'BH'))

per_year_results

# Export Data
write_csv(per_year_results, file.path(SD, 'coLoc_Permanova_byYear.csv'))

## ----------------------------------------- ##
##              PSME Driver ----

# Remove PSME (and co-dominant ACCI, ANDE3) to test sensitivity to
# the most abundant species
PSME_test_matrix <- comm_matrix |>
  select(-c(PSME, ACCI, ANDE3))

# Pairwise by-year comparison
years <- sort(unique(metadata$year))

per_year_results <- map_df(years, ~ {
  yr <- .x

  # subset metadata and community to this year
  idx      <- metadata$year == yr
  meta_yr  <- metadata[idx, , drop = FALSE]
  comm_yr  <- PSME_test_matrix[idx, -(1:3), drop = FALSE]     # species matrix for that year

  # distance matrix for this year only
  dist_yr <- vegdist(comm_yr, method = 'bray')

  # PERMANOVA: hillslope effect within this year
  ad <- adonis2(dist_yr ~ hillslope,
                data = meta_yr,
                permutations = 999,
                by = 'terms')

  hillslope_row <- ad['hillslope', ]

  tibble(
    year    = yr,
    R2      = hillslope_row$R2,
    F       = hillslope_row$F,
    p_value = hillslope_row$`Pr(>F)`
  )
})

# optional multiple-testing correction
per_year_results <- per_year_results |>
  mutate(p_adj = p.adjust(p_value, method = 'BH'))

per_year_results

write_csv(per_year_results, file.path(SD, 'PSME_removed_Permanova_byYear.csv'))
