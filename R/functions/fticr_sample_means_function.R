# -----------------------------------------------------------------------------#
# fticr_sample_means_function.R
#
# Purpose
#   Helper functions for summarizing FTICR-MS peak-level data to
#   sample-level molecular property means and detection-based
#   composition metrics (percent by elemental group and by van
#   Krevelen compound class).
#
# Functions
#   wmean(x, w)
#     Weighted mean, NA-safe: excludes pairs where x or w is NA or
#     w <= 0; returns NA_real_ if nothing remains.
#
#   calc_fticr_means(df_long, el_forms, class_levels, prop_cols)
#     Given a long-format table (one row per formula x sample, with a
#     detection flag `present` and a weight `w`), returns one row per
#     sample with:
#       - weighted means of prop_cols (w-weighted; w is either
#         presence/absence or intensity, set by the caller)
#       - percent of detected formulas in each of el_forms (elemental
#         group, e.g. CHO/CHON/CHOS)
#       - percent of detected formulas in each of class_levels (van
#         Krevelen compound class, e.g. Lignin/Protein/Lipid)
#     Percent metrics are always detection-based (over `present`),
#     regardless of what `w` represents.
#
# Created:  2026-01-25
# Author:   X. Takver
# -----------------------------------------------------------------------------#

## ----------------------------------------- ##
#          Weighted Mean Helper ----
## ----------------------------------------- ##

# Weighted mean with NA-safe behavior.
# If all weights are 0 (or all values are NA), returns NA_real_.
wmean <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

## ----------------------------------------- ##
#      Sample-Level Means + Composition ----
## ----------------------------------------- ##

# Sample-level means and detection-based composition metrics.
# - 'present' defines detections for percent metrics (presence/absence denominator).
# - 'w' defines weights for molecular property means (PA- or intensity-weighted).
calc_fticr_means <- function(df_long,
                             el_forms,
                             class_levels = c('Lipid', 'Unsat Hydrocarbon', 'Protein', 'Lignin',
                                              'Carbohydrate', 'Amino Sugar', 'Tannin', 'Cond Hydrocarbon'),
                             prop_cols = c('mass', 'C', 'H', 'O', 'N', 'S', 'ai_mod', 'DBE', 'NOSC', 'OC', 'HC')) {
  
  # Core property means + denominator for detection-based percentages
  out <- df_long |>
    dplyr::group_by(Sample) |>
    dplyr::summarise(
      dplyr::across(
      dplyr::all_of(prop_cols),
        ~ wmean(.x, w),
        .names = '{.col}'),
      n_detected = sum(present, na.rm = TRUE),
      .groups = 'drop')
  
  # Elemental diversity (detection-based)
  el_df <- df_long |>
    dplyr::group_by(Sample) |>
    dplyr::summarise(
      n_detected = sum(present, na.rm = TRUE),
      .groups = 'drop')
  
  for (ef in el_forms) {
    col_nm <- paste0(ef, '_percent')
    tmp <- df_long |>
      dplyr::group_by(Sample) |>
      dplyr::summarise(
        n = sum(present & el_form == ef, na.rm = TRUE),
        .groups = 'drop') |>
      dplyr::left_join(el_df, by = 'Sample') |>
      dplyr::mutate(!!col_nm := 100 * n / n_detected) |>
      dplyr::select(Sample, dplyr::all_of(col_nm))
    
    out <- out |>
      dplyr::left_join(tmp, by = 'Sample')
  }
  
  # Compound class diversity (detection-based)
  class_df <- df_long |>
    dplyr::group_by(Sample) |>
    dplyr::summarise(
      n_detected = sum(present, na.rm = TRUE),
      .groups = 'drop')
  
  for (cl in class_levels) {
    col_nm <- paste0(gsub(' ', '', cl), '_percent')
    tmp <- df_long |>
      dplyr::group_by(Sample) |>
      dplyr::summarise(
        n = sum(present & class == cl, na.rm = TRUE),
        .groups = 'drop') |>
      dplyr::left_join(class_df, by = 'Sample') |>
      dplyr::mutate(!!col_nm := 100 * n / n_detected) |>
      dplyr::select(Sample, dplyr::all_of(col_nm))
    
    out <- out |>
      dplyr::left_join(tmp, by = 'Sample')
  }
  
  # Drop denominator (kept only as an internal helper)
  out |>
    dplyr::select(-n_detected)
}
