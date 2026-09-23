# -----------------------------------------------------------------------------#
# fticr_qaqc_outlier_functions.R
#
# Purpose
#   Helper functions for FTICR-MS replicate-level QA/QC: summarizing
#   presence/absence data to per-sample peak counts with parsed sample
#   metadata, and flagging/removing outlier samples via a median +/-
#   k*MAD rule within Site x Sampling.Event x Core.Section strata.
#
# Functions
#   make_peak_counts(df_PA)
#     Collapses a formula x sample presence/absence table to one row
#     per sample (peak_no = total detections), then parses Site
#     (High/Low Storage), Sampling.Event (June/December), Core.Section
#     (Top/Bottom), Rep, and mRep from the standardized sample
#     identifier.
#
#   mad_filter_samples(df_peaks, k = 3)
#     Flags samples whose peak_no falls outside median +/- k*MAD of
#     their Site x Sampling.Event x Core.Section group, and returns
#     the flagged table, the retained (kept) samples, and the
#     group-level summary statistics.
#
# Notes
#   - Sample identifiers are expected in the published format
#     <Fraction>_<HS|LS><Pit>_<June|December>_<Top|Btm>_R<Rep>_mR<mRep>
#     (e.g. WEOM_HS1_December_Btm_R2_mR3).
# Created:  2026-01-25
# Author:   X. Takver
# -----------------------------------------------------------------------------#

## ----------------------------------------- ##
#            Peak Count Summary ----
## ----------------------------------------- ##

# Summarise presence/absence data into per-sample peak counts with parsed metadata
# (used to standardise outlier checks across fractions and years).
make_peak_counts <- function(df_PA) {
  df_PA |>
    # Convert presence/absence matrix to per-sample peak counts
    tibble::column_to_rownames('molecular_formula') |>
    dplyr::summarize(dplyr::across(dplyr::everything(), ~ sum(.x))) |>
    tidyr::pivot_longer(cols = dplyr::everything(),
                        names_to = 'SampleID',
                        values_to = 'peak_no') |>
    # Parse sample metadata directly from the standardized identifier
    # (expected format: <Fraction>_<HS|LS><Pit>_<June|December>_<Top|Btm>_R<Rep>_mR<mRep>)
    dplyr::mutate(
      Site = dplyr::case_when(stringr::str_detect(SampleID, 'HS') ~ 'High Storage',
                              stringr::str_detect(SampleID, 'LS') ~ 'Low Storage',
                              TRUE ~ NA_character_),
      Sampling.Event = stringr::str_extract(SampleID, 'June|December'),
      Rep  = as.integer(stringr::str_extract(SampleID, '(?<=_R)\\d+')),
      mRep = as.integer(stringr::str_extract(SampleID, '(?<=_mR)\\d+')),
      Core.Section = dplyr::case_when(stringr::str_detect(SampleID, 'Top') ~ 'Top',
                                      stringr::str_detect(SampleID, 'Btm') ~ 'Bottom',
                                      TRUE ~ NA_character_)
    ) |>
    dplyr::mutate(Core.Section = forcats::as_factor(Core.Section),
                  Site = forcats::as_factor(Site),
                  Sampling.Event = forcats::as_factor(Sampling.Event))
}

## ----------------------------------------- ##
#              MAD Outlier Filter ----
## ----------------------------------------- ##

# Flag and remove outlier samples using median ± k*MAD within grouping strata
# (used to keep QAQC rules consistent while allowing k to be tuned if needed).
mad_filter_samples <- function(df_peaks, k = 3) {
  # Summarise peak counts within Site x Sampling.Event x Core.Section groups
  df_peaks_summary <- df_peaks |>
    dplyr::group_by(Site, Sampling.Event, Core.Section) |>
    dplyr::summarise(avg_peak_no = mean(peak_no),
                     median_peak_no = median(peak_no),
                     SD1 = sd(peak_no),
                     IQR = IQR(peak_no),
                     MAD = mad(peak_no),
                     .groups = 'drop') |>
    dplyr::mutate(MADk = k * MAD)
  
  # Flag samples outside median ± k*MAD
  df_peaks_flagged <- df_peaks |>
    dplyr::left_join(df_peaks_summary, by = c('Site', 'Sampling.Event', 'Core.Section')) |>
    dplyr::mutate(FLAG.MAD = dplyr::case_when(
      peak_no > median_peak_no + MADk | peak_no < median_peak_no - MADk ~ 'FLAG',
      TRUE ~ 'NA'
    ))
  
  # Keep only samples not flagged by MAD
  df_kept <- df_peaks_flagged |>
    dplyr::filter(FLAG.MAD != 'FLAG')
  
  # Return flagged peaks table, retained samples, and group-level summary stats
  list(df_peaks_flagged = df_peaks_flagged,
       df_kept = df_kept,
       summary = df_peaks_summary)
}
