# ------------------------------------------------------------------------------
# Hydraulic redistribution summary plots -- Figure 3C and Figure S3
# ------------------------------------------------------------------------------
# Purpose:
#   - Read HR calculations, bin time as "days post-rain," summarize recharge (HR)
#   - Plot cumulative HR vs time bins as bubble plots (all sites, and hillslope-aggregated)
#
# Input:
#   - data/volumetric_water_content/processed/HR_Calculations_June-Aug.csv
#
# Outputs:
#   - Fig.S3-Cumulative-HR_All-sites (all sites, by site_id)
#   - Fig.3C-Cumulative-HR_Aggregated (hillslope-aggregated)
#
# Created:  2026-02-10
# Author:   X. Takver
# ------------------------------------------------------------------------------

## ----------------------------------------- ##
#             Package Import ----
## ----------------------------------------- ##

# Data Analysis
librarian::shelf(tidyverse, pracma, quantmod, zoo, mgcv, here)

# Data visualization
librarian::shelf(Cairo, tinythemes)

## ----------------------------------------- ##
##            Source functions ----
## ----------------------------------------- ##

# Project-wide helpers:
#   - ensure_dir(): idempotent directory creation
#   - fonts.R / plot_theme.R: consistent manuscript styling
source(here('R', 'functions', 'paths_and_directories.R')) 
source(here('R', 'functions', 'fonts.R'))                  
source(here('R', 'functions', 'plot_theme.R'))        

## ----------------------------------------- ##
#             Set Plot Theme ----
## ----------------------------------------- ##

# Set global ggplot2 theme for consistent manuscript figure styling
set_plot_theme()

## ----------------------------------------- ##
#           Directory Creation ----
## ----------------------------------------- ##

# DD: processed VWC/HR outputs
# PD: figure output directory for this script
DD <- here('data', 'volumetric_water_content', 'processed')
PD <- here('output', 'figures', 'volumetric_water_content')

# Ensure output directory exists
ensure_dir(PD)

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

# HR calculations include at least:
#   - datetime (time)
#   - site_id (sensor/site ID)
#   - recharge (HR / water input metric used here)
df_hr <- read_csv(file.path(DD, 'HR_Calculations_June-Aug.csv'))

## ----------------------------------------- ##
#              Data Processing ----
## ----------------------------------------- ##

# Create a time axis relative to the start of the period (days since first timestamp).
# Offset by +2 (script-specific convention) and round to whole days for binning.
df_hr <- df_hr |> 
  arrange(datetime) |>
  mutate(days_post_rain = 2 + time_length(datetime - first(datetime), "days"),
         days_round = round(days_post_rain))

# Bin days into 6-day intervals (0-54 by 6), with a final bin ending at 62.
# Labels represent the upper end (or nominal center) used for plotting.
breaks_vec <- c(seq(0, 54, by = 6), 62)        
labels_vec <- c(seq(6, 54, by = 6), 60)  

df_binned <- df_hr |>
  mutate(
    days_bin = cut(
      days_round,
      breaks = breaks_vec,
      include_lowest = TRUE,
      labels = labels_vec
    )
  )

## ----------------------------------------- ##
##    Calculate Average (aggregated)  ----


# Summarize HR by site_id within each time bin:
#   - cum_HR: cumulative sum of recharge over time (within site_id)
#   - average_HR / sd_HR: bin-level recharge summary (used for bubble size)
#   - cum_HR: max cumulative value reached within each bin (for plotting trajectory)
HR_weekly <- df_binned |> 
  group_by(site_id) |> 
  mutate(cum_HR = cumsum(replace_na(recharge, 0))) |> 
  ungroup() |> 
  group_by(days_bin, site_id) |> 
  summarise(n = n(), 
            average_HR = mean(recharge, na.rm = TRUE),
            sd_HR = sd(recharge, na.rm = TRUE),
            cum_HR = max(cum_HR, na.rm = TRUE),
            .groups = 'drop') 

# Bubble plot: x = time bin, y = cumulative HR; bubble size = |avg HR| in bin
plt.aggregated.bubble <- 
  HR_weekly |>
  mutate(HR_mag = abs(average_HR)) |>  
  ggplot(aes(x = days_bin,
             y = cum_HR,
             color = site_id,
             group = site_id)) + 
  geom_point(aes(size = HR_mag), alpha = 0.5) +
  scale_size_area(max_size = 6, name = "Avg HR") +
  scale_color_manual(values = c('HS1' = '#759295',
                                'HS2' = '#A9D9DE',
                                'LS1' = '#7E2E27',
                                'LS2' = '#C55648')) +
  labs(x = "Days post-Rain",
       y = "Cumulative hydraulic redistribution",
       color = "Site ID")

plt.aggregated.bubble

###             Export Plots ----

# Export site_id-level panel (Fig. S3)
CairoPDF(file.path(PD, 'Fig.S3-Cumulative-HR_All-sites'), width = 6.5, height = 3)
print(plt.aggregated.bubble)
dev.off()

## ----------------------------------------- ##
##   Average (aggregated Hillslope)  ----
## ----------------------------------------- ##

# Repeat summarization but aggregate site_ids to hillslope groups (HS vs LS):
#   - cum_HR is computed within hillslope and then averaged by bin (draft choice)
HR_weekly <- df_binned |> 
  mutate(hillslope = case_when(str_detect(site_id, 'HS') ~ 'High Storage',
                               str_detect(site_id, 'LS') ~ 'Low Storage')) |> 
  group_by(hillslope) |> 
  mutate(cum_HR = cumsum(replace_na(recharge, 0))) |> 
  ungroup() |> 
  group_by(days_bin, hillslope) |> 
  summarise(n = n(), 
            average_HR = mean(recharge, na.rm = TRUE),
            sd_HR = sd(recharge, na.rm = TRUE),
            cum_HR = mean(cum_HR, na.rm = TRUE),
            .groups = 'drop') 

# Bubble plot: hillslope-level aggregation
plt.aggregatedHS.bubble <- 
  HR_weekly |>
  mutate(HR_mag = abs(average_HR)) |>  
  ggplot(aes(x = days_bin,
             y = cum_HR,
             color = hillslope,
             group = hillslope)) + 
  geom_point(aes(size = HR_mag), alpha = 0.5) +
  scale_size_area(max_size = 6, name = "Avg HR") +
  scale_color_manual(values = c('High Storage' = '#759295',
                                'Low Storage' = '#7E2E27')) +
  labs(x = "Days post-Rain",
       y = "Cumulative hydraulic redistribution",
       color = "Hillslope")

plt.aggregatedHS.bubble

###             Export Plots ----

# Export hillslope-aggregated panel (Fig. 3C)
CairoPDF(file.path(PD, 'Fig.3C-Cumulative-HR_Aggregated'), width = 6.5, height = 3)
print(plt.aggregatedHS.bubble)
dev.off()
