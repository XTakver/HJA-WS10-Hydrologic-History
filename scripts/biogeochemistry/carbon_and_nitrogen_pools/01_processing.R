# ------------------------------------------------------------ #
#
# Purpose
#   Process WEOM C, MAOM C, and microbial biomass C (MBC) for
#   WS10 MONet biogeochemistry outputs.
#
# Inputs
#   HJA_WS10_2023_Soil_Carbon_and_Nitrogen_Pools.csv
#
# Outputs
#   data/Biogeochemistry/processed/2023_Carbon_and_Nitrogen_pools_responses.csv
#   data/Biogeochemistry/processed/2023_Carbon_and_Nitrogen_pools_processed.csv
#
# Created: 2026-02-19
# Author: X. Takver
# ------------------------------------------------------------ #

## ----------------------------------------- ##
#             Package Import ----
## ----------------------------------------- ##

librarian::shelf(tidyverse, here)

## ----------------------------------------- ##
#             Source Functions ----
## ----------------------------------------- ##

# Directory helpers (ensure_dir)
source(here('R', 'functions', 'paths_and_directories.R'))

## ----------------------------------------- ##
#        Directory helpers + paths ----
## ----------------------------------------- ##

DD <- here('data', 'biogeochemistry', 'raw')
OD <- here('data', 'biogeochemistry', 'processed')

ensure_dir(OD)

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

df_OCarbon <- read_csv(file.path(DD, 'HJA_WS10_2023_Soil_Carbon_and_Nitrogen_Pools.csv'))

## ----------------------------------------- ##
#            Factor Conversion ----
## ----------------------------------------- ##

# Sets grouping/plotting order for downstream models and figures
df_OCarbon <- df_OCarbon |>
  mutate(hillslope = factor(hillslope, levels = c('HS', 'LS')),
         sampling_month = factor(sampling_month, levels = c('June', 'December')),
         core_section = factor(core_section, levels = c('top', 'bottom')))

## ----------------------------------------- ##
#          Calculate log Change ----
## ----------------------------------------- ##

df_OCarbon_Long <- 
  df_OCarbon |> 
  select(doe_proposal_id:maom_c, weom_c, mbc) |> 
  pivot_longer(cols = c(maom_c, weom_c, mbc),
               names_to = 'carbon_pool',
               values_to = 'mg_per_kg') |> 
  mutate(carbon_pool = dplyr::recode(carbon_pool,
                              maom_c = 'MAOM C',
                              weom_c = 'WEOM C',
                              mbc = 'Microbial Biomass C'))

# Bootstrap Confidence intervals
#
# Computes the December-vs-June response ratio (and its log2 fold
# change) for one hillslope x core_section x carbon_pool group, with
# a percentile bootstrap CI on the log2 fold change. Applied per
# group via group_modify() below rather than run once on the whole
# dataset, since each stratum needs its own June/December samples.
set.seed(55)

boot_log2fc <- function(dat, n_boot = 5000){
  june <- dat |> 
    filter(sampling_month == 'June') |> 
    pull(mg_per_kg) |> 
    na.omit()
  
  dec <- dat |> 
    filter(sampling_month == 'December') |> 
    pull(mg_per_kg) |> 
    na.omit()
  
  # Bail out to NAs (still reporting n and raw means) when a group
  # can't support a stable log ratio: too few reps to resample, or a
  # non-positive mean makes log()/log2() undefined
  if(length(june) < 2 | length(dec) < 2 | any(june <= 0) | any(dec <= 0)){
    return(tibble(June = mean(june, na.rm = TRUE),
                  December = mean(dec, na.rm = TRUE),
                  response_ratio = NA_real_,
                  log_response_ratio = NA_real_,
                  log2_fold_change = NA_real_,
                  percent_change = NA_real_,
                  log2_ci_low = NA_real_,
                  log2_ci_high = NA_real_,
                  log_ci_low = NA_real_,
                  log_ci_high = NA_real_,
                  n_June = length(june),
                  n_December = length(dec)
    )
    )
  }
  
  # Resample June and December independently, with replacement, at
  # their original sample sizes; recompute log2(December/June) on
  # each resample to build the null-free sampling distribution of
  # the fold change
  boot_log2 <- replicate(n_boot,
                         {
                           june_mean <- mean(sample(june, size = length(june), replace = TRUE))
                           dec_mean  <- mean(sample(dec,  size = length(dec),  replace = TRUE))
                           
                           log2(dec_mean / june_mean)
                         }
  )
  
  tibble(June = mean(june),
         December = mean(dec),
         response_ratio = December / June,
         log_response_ratio = log(December / June),
         log2_fold_change = log2(December / June),
         percent_change = ((December - June) / June) * 100,
         # 95% percentile CI: 2.5th/97.5th percentiles of the
         # bootstrap distribution, then rescaled from log2 to
         # natural-log units (log2_ci * log(2)) for log_response_ratio
         log2_ci_low = quantile(boot_log2, 0.025, na.rm = TRUE),
         log2_ci_high = quantile(boot_log2, 0.975, na.rm = TRUE),
         log_ci_low = log2_ci_low * log(2),
         log_ci_high = log2_ci_high * log(2),
         n_June = length(june),
         n_December = length(dec)
  )
}

# Applies boot_log2fc() within each hillslope x core_section x
# carbon_pool stratum, so June/December are always resampled from
# the same group rather than pooled across groups
OCarbon_response <- df_OCarbon_Long |> 
  group_by(hillslope, core_section, carbon_pool) |> 
  group_modify(~ boot_log2fc(.x, n_boot = 5000)) |> 
  ungroup() |> 
  mutate(site_depth = paste(hillslope, core_section),
         carbon_pool = factor(carbon_pool,
                              levels = c('WEOM C', 'MAOM C', 'Microbial Biomass C')))


## ----------------------------------------- ##
#                  Export ----
## ----------------------------------------- ##

write_csv(df_OCarbon, file.path(OD, '2023_Carbon_and_Nitrogen_pools_processed.csv'))
write_csv(OCarbon_response, file.path(OD, '2023_Carbon_and_Nitrogen_pools_responses.csv'))
