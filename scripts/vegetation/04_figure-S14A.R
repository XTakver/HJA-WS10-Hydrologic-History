# ------------------------------------------------------------ #
# 04_figure-S14A.R
#
# Purpose
#   NMDS ordination of species abundance across all 36 long-term
#   WS10 vegetation microplots (1983-2016), with the four microplots
#   co-located with our study sites outlined, a sampling-year color
#   gradient by hillslope, and fitted species vectors. Produces
#   manuscript Figure S14, panel A.
#
# Inputs
#   data/vegetation/processed/WS10_GroundCover_AllPlots.csv
#   data/vegetation/processed/WS10_TreeSurvey_AllPlots.csv
#   data/vegetation/processed/diversity/Community_Matrix_All.csv
#   data/vegetation/processed/diversity/Diversity_DF_All.csv
#
# Outputs
#   output/figures/vegetation/Fig.S14A.pdf
#
# Notes
#   - NMDS fit at k = 3 dimensions (Bray-Curtis dissimilarity);
#     stress values for k = 2-5 are checked first.
#   - Fitted species vectors are filtered to p < 0.005, then the 12
#     longest vectors are kept. The Fig. S14 caption describes this
#     as "p < 0.05" with no fixed count -- this discrepancy has not
#     yet been resolved and should be checked before finalizing.
#   - Community_Matrix_All.csv and Diversity_DF_All.csv are produced
#     by 01_processing_Vegetation_Diversity.R, not by this script.
#
# Created:  2025-07-17
# Author:   X. Takver
# ------------------------------------------------------------ #

## ----------------------------------------- ##
#             Package Import ----
## ----------------------------------------- ##

librarian::shelf(tidyverse, glue, here, vegan, ggrepel, scales, Cairo)

## ----------------------------------------- ##
#             Source Functions ----
## ----------------------------------------- ##

source(here('R', 'functions', 'paths_and_directories.R'))  # ensure_dir()
source(here('R', 'functions', 'fonts.R'))                  # e.g., register_fonts()
source(here('R', 'functions', 'plot_theme.R'))              # e.g., set_plot_theme()

## ----------------------------------------- ##
#           Directory Creation ----
## ----------------------------------------- ##

DD <- here('data', 'vegetation', 'processed')
DivD <- here('data', 'vegetation', 'processed', 'diversity')

PD <- here('output', 'figures', 'vegetation')

ensure_dir(PD)


## ----------------------------------------- ##
#             Set Plot Theme ----
## ----------------------------------------- ##

# Set global ggplot2 theme for consistent manuscript figure styling
set_plot_theme()

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

# Import Surface ground cover estimates
tidy_grndcover <- read_csv(file.path(DD, 'WS10_GroundCover_AllPlots_processed.csv'))

# Import tree data and select relevant columns
tidy_tree <- read_csv(file.path(DD, 'WS10_TreeSurvey_AllPlots_processed.csv'))

community_matrix <- read_csv(file.path(DivD, 'Community_Matrix_All.csv')) |>
  unite(col = 'YR_PLT_HS', year, plot, hillslope, sep = '_',
        remove = FALSE) |>
  column_to_rownames('YR_PLT_HS')

diversity_df <- read_csv(file.path(DivD, 'Diversity_DF_All.csv'))

## ----------------------------------------- ##
#       All Years Groundcover NMDS ----
## ----------------------------------------- ##

## ----------------------------------------- ##
#             Species Summary ----

# Gets general summary statistics (min and max) for rarefaction bounds

community_summary <- community_matrix |>
  group_by(year, plot, hillslope) |>
  summarize(Total.Individuals = rowSums(across(-1:-3)))

min(community_summary$Total.Individuals)
max(community_summary$Total.Individuals)

# Minimum peaks detected: 42
# Maximum peaks detected: 480

## ----------------------------------------- ##
##      All Years NMDS Calculations ----

### Prepare data ----

# Create a distance matrix
dist_matrix <- vegdist(community_matrix[, -(1:3)], method = 'bray')


# Create an environmental metadata df
metadata <- tidy_grndcover |>
  distinct(year, plot, hillslope) |>
  unite('YR_PLT_HS', year, plot, hillslope, remove = FALSE) |>
  column_to_rownames('YR_PLT_HS')

### Calculate NMDS ----
# Set a seed for reproducibility
set.seed(42)

# Try different dimensions to determine lowest stress
stress_values <- numeric(4)

for(i in 2:5){
  stress_values[i - 1] <- metaMDS(dist_matrix, distance = 'bray', k=i, trace = FALSE)$stress * 100
}

print(stress_values)

# Run NMDS
nmds_result <- metaMDS(dist_matrix, k = 3, trymax = 100)

# Check the result
print(nmds_result)

# Assign nmds scores to df
nmds_scores <- as.data.frame(scores(nmds_result)) |>
  rownames_to_column('YR_PLT_HS') |>
  full_join(diversity_df)

### Create Species Variable Dataframe ----
# Add Species drivers

df_species <- as.data.frame(community_matrix[, -(1:3)])

fit_species <- envfit(nmds_result, df_species, permutations = 999)

# Too many significant species. Chose just the 15 longest vector lengths
species_scores <- as.data.frame(scores(fit_species, 'vectors')) |>
  mutate(pvalue = fit_species$vectors$pvals) |>
  rownames_to_column('species') |>
  filter(pvalue < 0.005) |>
  mutate(vector_length = sqrt(NMDS1^2 + NMDS2^2)) |>
  arrange(desc(vector_length)) |>
  head(12)

# ### Export data frames
#
# # Save NMDS Results
# write_csv(nmds_scores, file.path(DivD, 'NMDS_Scores_All_Years.csv'))
#
# # Save Species Variable Results
# write_csv(species_scores, file.path(DivD, 'Species_Scores_All_Years.csv'))

## ----------------------------------------- ##
##                Plot NMDS ----

###         Color=Site & No Ann ----

# Combine hillslope and year
nmds_scores_Grad <- nmds_scores |>
  mutate(hillslope_year = paste(hillslope, year, sep = '_')) |>
  # Normalize year within the data
  mutate(year_norm = rescale(year)) |>
  # generate custom colors
  mutate(color = case_when(
    hillslope == 'High Storage' ~ seq_gradient_pal('#d3e8ea', '#759295')(year_norm),
    hillslope == 'Low Storage'  ~ seq_gradient_pal('#f2d6cf', '#7E2E27')(year_norm)
  ))

# Create a co-located DF
nmds_scores_coLoc <-
  nmds_scores |>
  separate(YR_PLT_HS, into = c('Year', 'Plot', 'Hillslope1', 'Hillslope2'), remove = FALSE) |>
  select(-Hillslope1, -Hillslope2) |>
  filter(Plot %in% c(11, 34, 27, 28, 26, 4, 25, 13))

# Create plot with co-located layer
# Add co-located outlines
plt.nmds.coLoc <-
  ggplot(nmds_scores_coLoc) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = '#848482', linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = '#848482', linewidth = 0.4) +
  geom_point(aes(x = NMDS1, y = NMDS2, shape = hillslope),
             fill = NA, color = 'grey30', size = 1.4, stroke = 1) +
  scale_shape_manual(values = c('High Storage' = 21, 'Low Storage' = 24))

plt.nmds.coLoc

# Plot with Gradient Legend and Colored Hillslope
plt.nmds.combined <-
  plt.nmds.coLoc +
  geom_point(data = nmds_scores_Grad,
             aes(x = NMDS1, y = NMDS2, color = I(color), fill = I(color), shape = hillslope),
             size = 1.25) +
  coord_cartesian(xlim = c(-1.25, 1.15)) +
  scale_color_identity() +

  # Dummy layer for year legend using fill
  geom_point(aes(x = NMDS1, y = NMDS2, fill = year), alpha = 0) +
  scale_fill_gradient(low = 'grey90', high = 'grey30', name = 'Sampling Year') +
  guides(shape = guide_legend(override.aes = list(size = c(1.35, 1.25),
                              color = c('#759295', '#7E2E27'), fill = c('#759295', '#7E2E27'))),
         color = guide_legend(override.aes = list(shape = 21, size = 1.25))) +

  labs(subtitle = 'All Vegetation Plots: 1983–2016',
       shape = 'Hillslope')

plt.nmds.combined

### Add Annotations ----

plt.nmds.annotation <-
  plt.nmds.combined +
  annotate('segment', x = -1, xend = -1, y = -0.6, yend = 0.6,
           arrow = arrow(type = 'closed',
                         length = unit(0.15, 'cm'),
                         ends = 'last'),
           linewidth = 0.4, linejoin = 'mitre',
           color = '#7E2E27') +
  annotate('segment', x = 0.6, xend = -0.5, y = -0.9, yend = -0.9,
           arrow = arrow(type = 'closed',
                         length = unit(0.15, 'cm')),
           linewidth = 0.4, linejoin = 'mitre',
           color = '#759295') +
  annotate('text', x = -1, y = 0.75, label = 'Low Storage \nSuccessional Trajectory',
           size = 2.5, hjust = 0.5,
           color = '#7E2E27') +
  annotate('text', x = -0.9, y = -0.9, label = 'High Storage \nSuccessional Trajectory',
           size = 2.5, hjust = 0.5,
           color = '#759295')

plt.nmds.annotation

### Species Variables ----
# Add Species Variables

plt.nmds.EnvVar <-
  plt.nmds.annotation +
  geom_segment(data = species_scores,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.1, 'cm')), color = '#926C00', alpha = 0.7,
               inherit.aes = FALSE,
               lwd = 0.25,
               lineend = 'round', linejoin = 'mitre') +
  geom_text_repel(data = species_scores,
                  aes(x = NMDS1*1.25, y = NMDS2*1.25, label = species),
                  max.overlaps = 25,
                  force = 0.1,
                  point.padding = 0.0,
                  size = 2, hjust = 0.5,
                  color = '#926C00',
                  inherit.aes = FALSE)

plt.nmds.EnvVar

### Final Touches ----
# Final theme alterations

plt.nmds.final <-
  plt.nmds.EnvVar +
  theme(plot.subtitle = element_text(margin = margin(b = 2)),
        axis.title.y = element_text(margin = margin(r = -3)),
        axis.title.x = element_text(margin = margin(t = -2)),
        legend.margin = margin(l = -5))

### Save Annotated Plots ----

CairoPDF(file.path(PD, 'Fig.S14A.pdf'), width = 5, height = 3.5)
print(plt.nmds.final)
dev.off()

# Cleanup workspace
rm(plt.nmds.annotation, plt.nmds.combined, plt.nmds.EnvVar)
