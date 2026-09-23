# ------------------------------------------------------------------------------
# Water holding content time series (Upper 30 cm) -- Figure 3A

# Purpose:
#   - Read smoothed VWC time series (site_id already matches manuscript convention)
#   - Normalize to each sensor's max VWC to compute WHC (%)
#   - Linearly interpolate missing WHC (for plotting/aggregation)
#   - Aggregate to hillslope mean WHC by depth and timestamp
#   - Plot June-October 2024 upper 30 cm mean WHC with:
#       * linetype = number of active sensors
#       * shaded 50-70% WHC reference band
#
# Input:
#   - data/volumetric_water_content/processed/smoothed_vwc.csv
#
# Output:
#   - output/figures/volumetric_water_content/Fig.3A-WHC (CairoPDF)
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

# Project-wide helpers (consistent with manuscript scripts):
#   - ensure_dir(): create output directory if needed
#   - fonts.R: font registration (if used by plot theme)
#   - set_plot_theme(): global ggplot styling
source(here('R', 'functions', 'paths_and_directories.R')) 
source(here('R', 'functions', 'fonts.R'))                  
source(here('R', 'functions', 'plot_theme.R'))             

## ----------------------------------------- ##
#           Directory Creation ----
## ----------------------------------------- ##

# DD = directory containing processed VWC data
# PD = directory where figure outputs from this script will be saved
DD <- here('data', 'volumetric_water_content', 'processed')
PD <- here('output', 'figures', 'volumetric_water_content')

# Ensure output directory exists (idempotent)
ensure_dir(PD)

## ----------------------------------------- ##
#             Set Plot Theme ----
## ----------------------------------------- ##

# Set global ggplot2 theme for consistent manuscript figure styling
set_plot_theme()

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

# Read in smoothed VWC time series (site_id already matches manuscript labels).
df_vwc <- read_csv(file.path(DD, 'smoothed_vwc.csv'))

# ----------------------------------------- ##
#             Data Processing ----
## ----------------------------------------- ##

# Step 1: For each sensor (site_id) and depth, compute the maximum observed VWC.
# This max value becomes the normalization reference for WHC (% of maximum).
df_max <- df_vwc |> 
  group_by(site_id, depth) |> 
  summarize(max_vwc = max(Smooth_VWC_Avg, na.rm = TRUE))


# Step 2: Join max values back onto the full time series and compute WHC.
# WHC here is expressed as percent of each site_id x depth maximum VWC.
df_whc <- df_vwc |> 
  left_join(df_max) |> 
  mutate(whc = (Smooth_VWC_Avg/max_vwc * 100))

# Interpolate WHC across time to fill gaps for plotting/aggregation.
# - zoo::na.approx performs linear interpolation.
# - datetime is converted to numeric to provide a proper x-axis for interpolation.
# - na.rm = FALSE keeps leading/trailing NAs as NA (no extrapolation).
# - maxgap = Inf allows interpolation across arbitrarily large gaps (draft choice).
df_whc_interp <- df_whc |>
  arrange(site_id, depth, datetime) |>
  group_by(site_id, depth) |>
  mutate(whc_interp = na.approx(whc, x = as.numeric(datetime),
                                na.rm = FALSE, maxgap = Inf))

## ----------------------------------------- ##
#          Gap Filled Mean Plots ----
## ----------------------------------------- ##

###          Mean WHC for plots ----

# Add a hillslope label derived from site_id (HS vs LS).
# Then aggregate across sensors within each hillslope x depth x time:
#   - mean_whc: mean of interpolated WHC
#   - sd_whc / max_whc / min_whc: simple descriptive stats for context/QC
#   - n_sensors: number of non-missing *original* WHC values at that timestamp
#                (used as a "data completeness" indicator)

whc_mean <- df_whc_interp |> 
  mutate(hillslope = case_when(str_detect(site_id, 'HS') ~ 'High Storage',
                               str_detect(site_id, 'LS') ~ 'Low Storage'),
         .after = site_id) |> 
  group_by(hillslope, depth, datetime) |> 
  summarise(mean_whc = mean(whc_interp),
            sd_whc = sd(whc_interp),
            max_whc = max(whc_interp),
            min_whc = min(whc_interp),
            n_sensors = sum(!is.na(whc)),
            .groups   = 'drop')

# Subset for plotting: focus on upper 30 cm
whc_upper30 <- whc_mean |> 
  filter(depth <= 30) 

# QC check: identify timestamps where neither sensor contributed data (n_sensors == 0)
# (These times are later excluded from plotting so lines don't bridge missing periods.)
whc_check <- whc_upper30 |> 
  filter(n_sensors == 0)

# Define time bounds for full plotting period (used for shaded rectangle extent)
x_rng <- range(whc_upper30$datetime, na.rm = TRUE)


# Drop periods with 0 sensors so the plot doesn't connect across fully missing intervals.
# Create a segment ID (seg_id) that increments when the number of active sensors changes.
# This seg_id is later used to break line connections when data completeness shifts (1->2 or 2->1)
whc_upper30_plot <- whc_upper30 |>
  filter(n_sensors > 0) |>                                
  arrange(hillslope, depth, datetime) |>
  group_by(hillslope, depth) |>
  mutate(
    seg_id = cumsum(
      n_sensors != dplyr::lag(n_sensors, default = first(n_sensors))
    )
  ) |>
  ungroup()

###           Upper 30 Plots ----

# Plot mean WHC time series by depth (facets) and hillslope (color).
#   - Shaded band highlights an "optimal condition" window (50-70% WHC)
#   - linetype encodes data completeness: dashed = 1 sensor, solid = 2 sensors
#   - group includes seg_id to avoid connecting across completeness-change boundaries
plt.mean.upper <- 
  ggplot(whc_upper30_plot) +
  annotate('rect', 
           xmin = x_rng[1], xmax = x_rng[2], 
           ymin = 50, ymax = 70,
           fill = '#a3b88c', alpha = 0.45) +
  geom_line(aes(x = datetime, y = mean_whc,
                color = hillslope, linetype = factor(n_sensors), 
                group = interaction(hillslope, seg_id)),  
            linewidth = 0.35) +
  facet_wrap(~ depth, ncol = 1) +
  scale_color_manual(values = c('High Storage' = '#759295',
                                'Low Storage'  = '#7E2E27')) +
  scale_linetype_manual(values = c(`1` = 'dashed', `2` = 'solid'),
                        name   = 'Data Completeness') +
  labs(subtitle = 'Water Holding Content: June-October 2024',
       x = 'Date',
       y = 'Water Holding Content (%)',
       color = 'Site ID') +
  theme(panel.spacing   = unit(0.1, 'lines'),
        legend.position = 'bottom',
        legend.margin   = margin(t = -12.5))

plt.mean.upper

###             Export Plots ----

# Export as PDF using Cairo (reliable text rendering)
CairoPDF(file.path(PD, 'Fig.3A-WHC'), width = 6.5, height = 3)
print(plt.mean.upper)
dev.off()
