# -----------------------------------------------------------------------------#
# 06_statistics_WEOM_NOSC.R
#
# Purpose
#   Test seasonal (June vs. December) shifts in the nominal oxidation state of
#   carbon (NOSC) for unique WEOM molecular formulas, within each hillslope x
#   core-section stratum:
#     1) Two-sample Kolmogorov-Smirnov test on the NOSC distributions, with a
#        Monte Carlo p-value (NOSC is rounded, so many formulas share identical
#        values, which violates the exact KS test's no-ties assumption).
#     2) A permutation test on the mean NOSC difference (December - June), with
#        a bootstrapped 95% CI for that difference.
#
# Inputs
#   2023_FTICR-MS_WEOM_Presence-Absence_Unique.csv
#     (written by 04_unique_formula_subsetting_WEOM.R)
#   2023_FTICR-MS_WEOM_ChemProperties_aggregated.csv
#     (written by 03_replicate_aggregation_WEOM.R; NOSC is a per-formula
#     property and is unaffected by the uniqueness filtering in 04, so the
#     aggregated -- not unique-filtered -- properties table is the correct source)
#
# Outputs (written to SD)
#   - NOSC_KS_wMonteCarlo_WEOM.csv
#   - NOSC_mean-permutation_WEOM.csv
#
# Notes
#   - MAOM NOSC is not included here: the manuscript reports MAOM NOSC
#     descriptively only ("Because MAOM replication was limited, we
#     interpret these patterns descriptively"), with no formal KS test
#     or bootstrap. This script covers WEOM only, matching Methods 2.5.3
#     and the Results text, which report the KS/bootstrap tests for
#     WEOM specifically.
#   - B = 9999 permutations/bootstraps throughout; set.seed(55) for
#     reproducibility of the mean-permutation/bootstrap block. 
#
# Created:  2025-11-22
# Author:   X. Takver
# -----------------------------------------------------------------------------#

## ----------------------------------------- ##
#             Import Packages ----
## ----------------------------------------- ##

# Data Analysis
librarian::shelf(tidyverse, here, glue)

## ----------------------------------------- ##
#             Source Functions ----
## ----------------------------------------- ##

# Directory helpers (ensure_dir)
source(here('R', 'functions', 'paths_and_directories.R'))

## ----------------------------------------- ##
#          Directory Creation ----
## ----------------------------------------- ##

# DD_unique: unique-formula WEOM presence/absence
# DD_props:  formula-level properties (incl. NOSC), pre-uniqueness filtering
DD_unique <- here('data', 'mass_spectrometry', 'processed', 'water-extractable_om', 'unique_formulas')
DD_props  <- here('data', 'mass_spectrometry', 'processed', 'water-extractable_om', 'replicate_aggregation')

# Stats Directory Variables
SD <- here('output', 'statistics', 'mass_spectrometry')

ensure_dir(SD)

## ----------------------------------------- ##
#                 Data Import ----
## ----------------------------------------- ##

## ----------------------------------------- ##
##                WEOM Data ----

weom_pa <- read_csv(file.path(DD_unique, '2023_FTICR-MS_WEOM_Presence-Absence_Unique.csv'))
weom_compounds <- read_csv(file.path(DD_props, '2023_FTICR-MS_WEOM_ChemProperties_aggregated.csv'))

## ----------------------------------------- ##
#               NOSC Calcs ----
## ----------------------------------------- ##

# Join detections to per-formula NOSC, expand to long form, and parse sample
# metadata directly from the standardized identifier (format produced by
# 04_unique_formula_subsetting_WEOM.R: <HS|LS><Pit>_<June|December>_<Top|Btm>).
weom_nosc <- weom_pa |> 
  left_join(weom_compounds |> select(molecular_formula, NOSC), by = 'molecular_formula') |> 
  pivot_longer(cols = -c(molecular_formula, NOSC),
               names_to = 'Sample.ID',
               values_to = 'PA') |> 
  separate(Sample.ID, into = c('Site.ID', 'Sampling.Event', 'Core.Section'), sep = '_') |> 
  mutate(Site.ID = case_when(str_detect(Site.ID, 'HS') ~ 'HS',
                             str_detect(Site.ID, 'LS') ~ 'LS'),
         Core.Section = case_when(Core.Section == 'Top' ~ '0-10cm',
                                  Core.Section == 'Btm' ~ '20-30cm'),
         Sampling.Event = fct_relevel(Sampling.Event, 'June', 'December'))

## ----------------------------------------- ##
#           Statistical Analysis ----
## ----------------------------------------- ##

## ----------------------------------------- ##
##         KS (Distribution Test) ----

# Because NOSC is rounded, many values share the same NOSC. 
# This produces a warning in the KS test below so they are simulated with monte-carlo.

ks_results_weom <- weom_nosc |>
  filter(PA == 1) |>
  group_by(Site.ID, Core.Section) |>
  group_modify(~ {
    nosc_june <- .x$NOSC[.x$Sampling.Event == "June"]
    nosc_dec  <- .x$NOSC[.x$Sampling.Event == "December"]
    
    ks <- ks.test(
      nosc_june,
      nosc_dec,
      exact = FALSE,             # don't try exact (ties violate it)
      simulate.p.value = TRUE,   # Monte Carlo p-value
      B = 9999                   # number of simulations
    )
    
    tibble(
      D       = as.numeric(ks$statistic),
      p_value = ks$p.value,
      n_june  = length(nosc_june),
      n_dec   = length(nosc_dec)
    )
  }) |>
  ungroup()
  
### Export Monte Carlo ----

write_csv(ks_results_weom, file.path(SD, 'NOSC_KS_wMonteCarlo_WEOM.csv'))

## ----------------------------------------- ##
##        Testing Means (Bootstrap) ----

set.seed(55)   # for reproducibility
B <- 9999       # number of permutations / bootstraps

nosc_change_weom <- weom_nosc |>
  filter(PA == 1) |>
  group_by(Site.ID, Core.Section) |>
  group_modify(~ {
    dat <- .x
    
    # Split by season
    nosc_june <- dat$NOSC[dat$Sampling.Event == "June"]
    nosc_dec  <- dat$NOSC[dat$Sampling.Event == "December"]
    
    # Observed difference in means (Dec - June)
    obs_diff <- mean(nosc_dec, na.rm = TRUE) - mean(nosc_june, na.rm = TRUE)
    
    ## -------- Permutation test on mean -------- ##
    all_vals <- c(nosc_june, nosc_dec)
    n_june   <- length(nosc_june)
    
    perm_diffs <- replicate(B, {
      idx      <- sample(seq_along(all_vals))      # shuffle labels
      june_new <- all_vals[idx[1:n_june]]
      dec_new  <- all_vals[idx[(n_june + 1):length(all_vals)]]
      mean(dec_new) - mean(june_new)
    })
    
    perm_p <- mean(abs(perm_diffs) >= abs(obs_diff))  # two-sided p
    
    ## -------- Bootstrap CI for mean difference -------- ##
    boot_diffs <- replicate(B, {
      june_boot <- sample(nosc_june, length(nosc_june), replace = TRUE)
      dec_boot  <- sample(nosc_dec,  length(nosc_dec),  replace = TRUE)
      mean(dec_boot) - mean(june_boot)
    })
    
    tibble(
      n_june       = length(nosc_june),
      n_dec        = length(nosc_dec),
      obs_diff     = obs_diff,
      perm_pvalue  = perm_p,
      boot_mean    = mean(boot_diffs),
      boot_ci_low  = quantile(boot_diffs, 0.025, na.rm = TRUE),
      boot_ci_high = quantile(boot_diffs, 0.975, na.rm = TRUE)
    )
  }) |>
  ungroup()

### Export Bootstraped Mean ----

write_csv(nosc_change_weom, file.path(SD, 'NOSC_mean-permutation_WEOM.csv'))
