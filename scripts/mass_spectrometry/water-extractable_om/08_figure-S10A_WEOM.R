# ------------------------------------------------------------------------------
# 08_figure-S10A_WEOM.R
#
# WEOM condensed heatmap (FTICR-MS) -- post-QAQC + 2/3 consensus
# ------------------------------------------------------------------------------
# Purpose:
#   - Read replicate-aggregated WEOM presence/absence (post-MAD filtering + 2/3 rule)
#   - Join per-formula compound class annotations
#   - Build a ComplexHeatmap presence/absence heatmap with class splits
#   - Export a publication-ready PDF for manuscript figures/supplement
#
# Inputs (from FTICR processing pipeline):
#   - 2023_FTICR-MS_WEOM_aggregated.csv
#   - 2023_FTICR-MS_WEOM_ChemProperties_aggregated.csv
#   (both written by 03_replicate_aggregation_WEOM.R)
#
# Output:
#   - Fig.S10A_WEOM_Heatmap.pdf
#
# Created:  2025-08-20
# Author:   X. Takver
# ------------------------------------------------------------------------------

## ----------------------------------------- ##
#             Import Packages ----
## ----------------------------------------- ##

# Data processing + distance calculation
librarian::shelf(tidyverse, vegan, here, glue)

# Heatmap rendering + PDF export
librarian::shelf(Cairo, ComplexHeatmap, circlize, grid)

## ----------------------------------------- ##
##            Source functions ----
## ----------------------------------------- ##

# Project-wide helpers:
#   - ensure_dir(): idempotent directory creation
#   - fonts.R / plot_theme.R: consistent manuscript styling (fonts used in gpar below)
source(here('R', 'functions', 'paths_and_directories.R'))
source(here('R', 'functions', 'fonts.R'))
source(here('R', 'functions', 'plot_theme.R'))

## ----------------------------------------- ##
#           Directory Creation ----
## ----------------------------------------- ##

# DD: replicate-aggregated (sample-level) FTICR outputs
# PD: figure output directory for this script
DD <- here('data', 'mass_spectrometry', 'processed', 'water-extractable_om', 'replicate_aggregation')
PD <- here('output', 'figures', 'mass_spectrometry', 'water-extractable_om')

ensure_dir(PD)

## ----------------------------------------- ##
#             Set Plot Theme ----
## ----------------------------------------- ##

# Apply project-wide styling (ggplot theme is not required for ComplexHeatmap,
set_plot_theme()

## ----------------------------------------- ##
#               Data Import ----
## ----------------------------------------- ##

# Presence/absence after:
#   - replicate outlier screening (MAD)
#   - instrument + technical replicate consensus (2/3 rule)
weom_pa <- read_csv(file.path(DD, '2023_FTICR-MS_WEOM_aggregated.csv'))

# Per-formula properties table (includes compound Class used for column splits)
weom_chemProp <- read_csv(file.path(DD, '2023_FTICR-MS_WEOM_ChemProperties_aggregated.csv'))

## ----------------------------------------- ##
#             Data Preparation ----
## ----------------------------------------- ##

# Join only the class annotation onto the PA table (heatmap is detection-based)
weom_combined <- weom_pa |>
  left_join(weom_chemProp |> select(molecular_formula, class),
            by = 'molecular_formula')

# Long format (one row per formula x sample) for relabeling / restructuring.
# Sample ID format produced upstream: <HS|LS><Pit>_<June|December>_<Top|Btm>
weom_long <- weom_combined |>
  pivot_longer(cols = -c(molecular_formula, class),
               names_to = 'Sample.ID',
               values_to = 'PA') |>
  separate(Sample.ID, into = c('Site.ID', 'Sampling.Event', 'Core.Section')) |>
  mutate(
    Site.ID = case_when(
      str_detect(Site.ID, 'HS') ~ 'HS',
      str_detect(Site.ID, 'LS') ~ 'LS'
    ),
    Sampling.Event = case_when(
      Sampling.Event == 'June' ~ 'June',
      Sampling.Event == 'December' ~ 'Dec.'
    ),
    Core.Section = case_when(
      Core.Section == 'Top' ~ '0-10cm',
      Core.Section == 'Btm' ~ '20-30cm'
    )
  )

# Wide heatmap input: columns are sample labels, rows are molecular formulas
df_cHM <- weom_long |>
  mutate(Sample.ID = glue('{Site.ID}, {Sampling.Event} \n{Core.Section}')) |>
  select(-Site.ID, -Core.Section, -Sampling.Event) |>
  pivot_wider(names_from = 'Sample.ID',
              values_from = 'PA') |>
  arrange(class, molecular_formula) |>
  ungroup()

## ----------------------------------------- ##
#           WEOM Condensed Heatmap ----
## ----------------------------------------- ##

# Annotation dataframe: formula -> compound class (used for column splits + top bar)
dfcHM_ann <- df_cHM |>
  select(molecular_formula, class) |>
  column_to_rownames('molecular_formula')

# ComplexHeatmap expects a matrix; transpose so:
#   - rows = samples (for clustering / ordered display)
#   - columns = formulas (split/annotated by compound class)
mat_cHM <- df_cHM |>
  select(-class) |>
  column_to_rownames('molecular_formula') |>
  as.matrix() |>
  t() |>
  replace_na(0)

## ----------------------------------------- ##
#      Plot Aesthetics + Annotations ----
## ----------------------------------------- ##

# Compound colors (compound class bar + column splits)
compound_colors <- c(
  'Carbohydrate'      = '#798E87',
  'Protein'           = '#C27D38',
  'Lipid'             = '#CCC591',
  'Amino Sugar'       = '#29211F',
  'Unsat Hydrocarbon' = '#446455',
  'Cond Hydrocarbon'  = '#FDD262',
  'Tannin'            = '#D3DDDC',
  'Lignin'            = '#DD8D29',
  'Other'             = '#C7B19C'
)

# Column annotation (compound class). Legend suppressed for compact panel layouts.
column_ann <- columnAnnotation(
  ` ` = dfcHM_ann$class,
  col = list(` ` = compound_colors),
  show_legend = FALSE,
  simple_anno_size = unit(0.4, 'cm')
)

# Heatmap colors (binary detection)
HM_colors <- colorRamp2(c(0, 1), c('black', 'red'))

# Split columns by compound class to group related formulas
column_split_factor <- dfcHM_ann$class

# Left annotation spacer (visual buffer between row labels and heatmap body)
left_annot <- rowAnnotation(
  spacer = anno_block(gp = gpar(fill = 'white', col = NA)),
  gap = unit(2, 'mm'),
  border = FALSE,
  width = unit(0.5, 'cm'),
  height = unit(0.5, 'cm')
)

## ----------------------------------------- ##
#             Sample Clustering ----
## ----------------------------------------- ##

# Cluster samples by Jaccard dissimilarity in detection space
row_dist <- vegdist(mat_cHM, 'jaccard')
row_clust <- hclust(row_dist, method = 'complete')
row_dend <- as.dendrogram(row_clust)

# Enforce a manuscript-consistent sample order in the dendrogram display
desired_order <- c(
  'LS, June \n0-10cm',
  'LS, June \n20-30cm',
  'HS, June \n0-10cm',
  'HS, June \n20-30cm',
  'LS, Dec. \n0-10cm',
  'LS, Dec. \n20-30cm',
  'HS, Dec. \n0-10cm',
  'HS, Dec. \n20-30cm'
)

row_dend_ordered <- reorder(
  row_dend,
  order(match(labels(row_dend), desired_order))
)

## ----------------------------------------- ##
#                 Plot ----
## ----------------------------------------- ##

plt.heatmap <-
  Heatmap(
    mat_cHM,
    cluster_rows = row_dend_ordered,
    cluster_columns = FALSE,
    
    # Annotations + layout
    top_annotation = column_ann,
    left_annotation = left_annot,
    column_split = column_split_factor,
    column_gap = unit(1, 'mm'),
    row_split = 8,
    
    # Labels / legend
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_gp = gpar(family = 'Goldman Sans Condensed', fontsize = 5.5),
    show_heatmap_legend = FALSE,
    
    # Dendrogram styling
    row_dend_reorder = FALSE,
    row_dend_gp = gpar(lwd = 1),
    row_dend_width = unit(1, 'cm'),
    
    # Panel styling
    column_title = '   ',
    column_title_gp = gpar(family = 'Goldman Sans Condensed', fontsize = 8.5),
    row_title = NULL,
    border = TRUE,
    col = HM_colors,
    
    # Output sizing (tuned for manuscript panel dimensions)
    heatmap_width = unit(4.25, 'in'),
    heatmap_height = unit(2.25, 'in')
  )

## ----------------------------------------- ##
#                 Export ----
## ----------------------------------------- ##

# Export as vector PDF (rasterize heatmap body for file size / rendering stability)
CairoPDF(file.path(PD, 'Fig.S10A_WEOM_Heatmap.pdf'), width = 4.5, height = 2.5)
draw(
  plt.heatmap,
  use_raster = TRUE,
  raster_quality = 2
)
dev.off()
