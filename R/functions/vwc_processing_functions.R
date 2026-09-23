# ------------------------------------------------------------------------------
# vwc_processing_functions.R
#
# Purpose
#   Helper functions for processing WS10 Sentek volumetric water content
#   (VWC) data: gap-limited interpolation + loess smoothing per site x
#   depth, peak/valley (drydown/recharge) identification, and hydraulic
#   redistribution (HR) calculations from the smoothed series.
#
# Functions
#   interp_and_smooth(df, pit_id, start_date, end_date, Depth)
#     Subsets df to one site_id x depth x date range, gap-limited
#     linear interpolation (na.approx, maxgap = 5) followed by loess
#     smoothing (span = 0.005). Returns a list(data, plot).
#
#   pad_with_na(x, target_length)
#     Left-pads a (possibly list-like) numeric vector with NA to a
#     target length; used to align slider::slide_period_* output back
#     to the original row count.
#
#   PV_ID(df, threshold, minpeakdistance)
#     Identifies peaks/valleys in Smooth_VWC_Avg using both
#     pracma::findpeaks() (threshold/distance-based) and
#     quantmod::findPeaks() (derivative-based), for more robust
#     detection than either method alone.
#
#   combined_PV_Deriv(df)
#     Combines the two peak/valley detection methods from PV_ID() into
#     a single PeakValley series, resolving adjacent same-type flags.
#
#   PV_excess_occurence(df, window, start_time, end_time, step)
#     Flags and removes peaks/valleys that occur too close together in
#     time (>= 2 within +/- window hours), which likely reflect noise
#     rather than genuine drydown/recharge events.
#
#   calc_HR(df)
#     Converts smoothed VWC to mm H2O and computes hydraulic
#     redistribution as the peak-to-valley and valley-to-peak change
#     in stored water.
#
#   calc_vwc_hr(pit_id, depth_cm, start_time, end_time, df_vwc, ...)
#     Wrapper that runs the full peak/valley + HR pipeline for one
#     site_id x depth x time window and returns the combined table
#     (smoothed VWC, peak/valley flags, drawdown/recharge, and the
#     predicted-no-HR series).
#
# Created:  2025-06-09
# Author:   X. Takver
# ------------------------------------------------------------------------------

## ----------------------------------------- ##
#           Processing Functions ----
## ----------------------------------------- ##

## ----------------------------------------- ##
##        Interpolate and Smooth ----

# This function interpolates missing data, where applicable and uses a loess smoothing function on data.
# Outputs a subset df with smoothed and interpolated VWC (Avg) values and a plot of smoothed v original data
interp_and_smooth <- function(df, pit_id, start_date, end_date, Depth) {
  
  # Subset and keep full time range
  df_subset <- df %>%
    filter(site_id == pit_id,
           depth == Depth,
           datetime >= start_date & datetime <= end_date) %>%
    arrange(datetime) %>%
    mutate(
      # Interpolate short gaps only (max 5 time steps), longer gaps stay NA
      VW_Interp = na.approx(vwc, maxgap = 5, na.rm = FALSE),
      datetime_num = as.numeric(datetime))
  
  # Only keep rows with interpolated data for loess model
  df_loess_fit <- df_subset %>%
    filter(!is.na(VW_Interp), !is.na(datetime_num))
  
  # Fit loess
  loess_model <- loess(VW_Interp ~ datetime_num, data = df_loess_fit, span = 0.005)
  
  # Predict only where interpolated data exists
  df_subset <- df_subset %>%
    mutate(
      Smooth_VWC_Avg = ifelse(
        !is.na(VW_Interp),
        predict(loess_model, newdata = data.frame(datetime_num = datetime_num)),
        NA)) %>%
    select(datetime, site_id, depth, vwc, VW_Interp, Smooth_VWC_Avg)
  
  # Plot result
  plot <- ggplot(df_subset, aes(x = datetime)) +
    geom_line(aes(y = vwc), color = "black", alpha = 0.4) +
    geom_line(aes(y = Smooth_VWC_Avg), color = "red") +
    labs(title = paste(pit_id, Depth, "Smoothed vs Raw"), y = "Water Content", x = "Date")
  
  return(list(data = df_subset, plot = plot))
}

## ----------------------------------------- ##
##            Pad with NA ----

pad_with_na <- function(x, target_length) {
  # slider::slide_* often returns a list-like object; flatten it
  if (is.list(x)) {
    x <- unlist(x, use.names = FALSE)
  }
  
  x <- as.numeric(x)
  
  len <- length(x)
  
  if (len < target_length) {
    return(c(rep(NA_real_, target_length - len), x))
  } else if (len > target_length) {
    # just in case it's too long, trim to target_length
    return(x[seq_len(target_length)])
  } else {
    return(x)
  }
}

## ----------------------------------------- ##
##      Peak/Valley Identification----

# Function uses both pracma and quantmod findpeaks to isolate peaks and valleys
# Using two methods makes for more robust peak and valley identification

PV_ID <- function(df, threshold, minpeakdistance) {
  # Keep a stable index for mapping back
  df <- df %>% 
    dplyr::mutate(index = dplyr::row_number())
  
  # Indices of non-NA smoothed values
  valid_idx <- which(!is.na(df$Smooth_VWC_Avg))
  
  # If too few points to detect anything, just return with 'none'
  if (length(valid_idx) < 3) {
    return(
      df %>%
        dplyr::mutate(
          value = NA_real_,
          type  = "none",
          Max   = NA_real_,
          Min   = NA_real_))
  }
  
  # Work on the non-NA series only
  s <- df$Smooth_VWC_Avg[valid_idx]
  
  ## ---- pracma peaks/valleys on non-NA data ----
  peaks   <- pracma::findpeaks(s,  threshold = threshold, minpeakdistance = minpeakdistance)
  valleys <- pracma::findpeaks(-s, threshold = threshold, minpeakdistance = minpeakdistance)
  
  # Helper to safely build tibbles even when no peaks/valleys found
  make_extremes <- function(mat, type_label, sign_flip = 1) {
    if (is.null(mat) || nrow(mat) == 0) {
      tibble::tibble(
        value = numeric(0),
        index = integer(0),
        type  = character(0)
      )
    } else {
      tibble::tibble(
        value = sign_flip * mat[, 1],
        # map local index back to original row index
        index = valid_idx[mat[, 2]],
        type  = type_label
      )
    }
  }
  
  df_peaks   <- make_extremes(peaks,   "peak",   sign_flip =  1)
  df_valleys <- make_extremes(valleys, "valley", sign_flip = -1)
  
  all_extremes <- dplyr::bind_rows(df_peaks, df_valleys)
  
  ## ---- join pracma peaks/valleys back to full df (including NAs) ----
  df <- df %>%
    dplyr::left_join(all_extremes, by = "index") %>%
    dplyr::mutate(
      type = dplyr::if_else(is.na(type), "none", type)
    )
  
  ## ---- quantmod peaks/valleys (derivative-based) on non-NA data ----
  max_idx_local <- quantmod::findPeaks(s)
  min_idx_local <- quantmod::findPeaks(-s)
  
  # Map back to global row indices
  max_idx <- valid_idx[max_idx_local]
  min_idx <- valid_idx[min_idx_local]
  
  df <- df %>%
    dplyr::mutate(
      Max = dplyr::if_else(index %in% max_idx, Smooth_VWC_Avg, NA_real_),
      Min = dplyr::if_else(index %in% min_idx, Smooth_VWC_Avg, NA_real_)
    )
  
  return(df)
}

## ----------------------------------------- ##
##     Create a Peak/Valley DataFrame ----

# Function to create two dfs (Peak_Valley and Derivative) and bind them together. 

combined_PV_Deriv <- function(df){
  
  # Create a pracma peak-valley df and flag regions where 2+ peaks or 2+ valleys occur in a row
  df_PV <- df %>% 
    filter(type == 'peak' | type == 'valley') %>% 
    mutate(flag = case_when(type == lag(type) ~ 'FLAG',
                            type != lag(type) ~ 'OK'))
  
  # Create a quantmod peak-valley df
  df_Deriv <-  df %>% 
    filter(!is.na(Max) | !is.na(Min)) %>% 
    mutate(flag = NA)
  
  # combine dataframes
  df_MaxMin <- rbind(df_Deriv, df_PV) %>% 
    arrange(as_datetime(datetime)) %>% 
    mutate(PeakValley = case_when(lead(flag) == 'FLAG' & lead(type) == 'valley' ~ 'peak',
                                  lead(flag) == 'FLAG' & lead(type) == 'peak' ~ 'valley',
                                  TRUE ~ type)) %>% 
    filter(PeakValley == 'peak' | PeakValley == 'valley') %>% 
    arrange(as_datetime(datetime)) 
  
  return(df_MaxMin)
}

## ----------------------------------------- ##
##            Peak/Valley Check ----

# Function to check data and see if there are too many peaks and valleys too close together
# Final cleanup of PV data - removes peaks and valleys when there are too many. 
# Sets Peak = Valley when there are too many

PV_excess_occurence <- function(df,
                                window,
                                start_time = NULL,
                                end_time   = NULL,
                                step       = "10 min") {
  
  df <- df %>%
    dplyr::mutate(datetime = lubridate::as_datetime(datetime))
  
  # if no explicit start/end, use range of the data
  if (is.null(start_time)) {
    start_time <- min(df$datetime, na.rm = TRUE)
  }
  if (is.null(end_time)) {
    end_time <- max(df$datetime, na.rm = TRUE)
  }
  
  # regular 10-min grid
  blank_ts <- seq.POSIXt(start_time, end_time, by = step)
  blank_df <- tibble::tibble(datetime = blank_ts)
  
  # ensure complete time grid
  df <- blank_df %>%
    dplyr::left_join(df, by = "datetime")
  
  n <- nrow(df)
  
  # counts of peaks / valleys in +/- `window` hours around each time
  peak_counts <- slider::slide_period_dbl(
    .x       = df$PeakValley,
    .i       = df$datetime,
    .period  = "second",
    .before  = window * 60 * 60,
    .after   = window * 60 * 60,
    .complete = TRUE,
    .f       = ~ sum(.x == "peak",   na.rm = TRUE)
  )
  
  valley_counts <- slider::slide_period_dbl(
    .x       = df$PeakValley,
    .i       = df$datetime,
    .period  = "second",
    .before  = window * 60 * 60,
    .after   = window * 60 * 60,
    .complete = TRUE,
    .f       = ~ sum(.x == "valley", na.rm = TRUE)
  )
  
  df <- df %>%
    dplyr::mutate(
      PeakCount   = pad_with_na(peak_counts,   n),
      ValleyCount = pad_with_na(valley_counts, n)
    )
  
  df_MaxMin <- df %>%
    dplyr::filter(PeakValley %in% c("peak", "valley")) %>%
    dplyr::mutate(
      Flag_PV = dplyr::case_when(
        PeakCount >= 2 | ValleyCount >= 2 ~ "FLAG",
        TRUE                              ~ "OK"
      ),
      value = dplyr::case_when(
        !is.na(value) ~ value,
        !is.na(Max)   ~ Max,
        !is.na(Min)   ~ Min,
        TRUE          ~ NA_real_
      )
    ) %>%
    dplyr::filter(!(PeakCount >= 2 & ValleyCount >= 2)) %>%
    dplyr::select(datetime:type, PeakValley)
  
  return(df_MaxMin)
}

## ----------------------------------------- ##
##            HR Calculations ----

calc_HR <- function(df){
  df_MaxMin <- df %>% 
    mutate(
      mm_h2o  = Smooth_VWC_Avg * 10,
      delta_mm = mm_h2o - lag(mm_h2o),
      flag = case_when(
        PeakValley == "peak"   & sign(delta_mm) == -1 ~ "FLAG",
        PeakValley == "valley" & sign(delta_mm) ==  1 ~ "FLAG",
        TRUE                                        ~ NA_character_  # <-- key change
      )
    ) %>%
    # Look for delta_mm where peaks are negative; adjust mm_h2o and recompute delta_mm
    mutate(
      mm_h2o = case_when(
        PeakValley == "peak" & sign(delta_mm) == -1 ~ lag(mm_h2o),
        TRUE                                        ~ mm_h2o
      ),
      delta_mm = mm_h2o - lag(mm_h2o)
    )
  
  return(df_MaxMin)
}

## ----------------------------------------- ##
##    Process files and bind together ----

calc_vwc_hr <- function(pit_id,
                        depth_cm,
                        start_time,
                        end_time,
                        df_vwc,
                        threshold       = 0.0003,
                        minpeakdistance = 100,
                        window          = 6,
                        smooth_col      = "Smooth_VWC_Avg") {
  
  # 1. Filter smoothed VWC for this pit/depth/window
  df_VWC <- df_vwc %>%
    dplyr::filter(
      site_id == pit_id,
      depth == depth_cm,
      datetime >= start_time,
      datetime <= end_time
    ) %>%
    dplyr::mutate(
      site_id = forcats::as_factor(site_id),
      depth = forcats::as_factor(depth)
    )
  
  if (nrow(df_VWC) == 0) {
    return(tibble::tibble())
  }
  
  # 2. Peak/valley + HR pipeline on this filtered series
  df_pv_full <- df_VWC %>%
    PV_ID(threshold = threshold, minpeakdistance = minpeakdistance) %>%
    combined_PV_Deriv() %>%
    PV_excess_occurence(window = window,
                        start_time = start_time,
                        end_time   = end_time) %>%
    calc_HR()
  
  # Keep only the HR-related bits we actually need to join back
  df_pv <- df_pv_full %>%
    dplyr::select(
      datetime,
      PeakValley,
      mm_h2o,
      delta_mm
      # add 'value' or other HR-specific cols here if needed
    )
  
  # 3. Create full 10-min grid for this window
  blank_ts <- seq.POSIXt(start_time, end_time, by = "10 min")
  blank_df <- tibble::tibble(datetime = blank_ts)
  
  # 4. Spread HR info over full time grid (no VWC columns here!)
  df_hr_full <- blank_df %>%
    dplyr::left_join(df_pv, by = "datetime") %>%
    dplyr::mutate(
      site_id = forcats::as_factor(pit_id),
      depth = forcats::as_factor(depth_cm)
    )
  
  # 5. Join HR info back to VWC; VWC columns only come from df_VWC
  df_join <- df_VWC %>%
    dplyr::left_join(
      df_hr_full,
      by = c("datetime", "site_id", "depth")
    ) %>%
    dplyr::mutate(
      flag = dplyr::case_when(
        PeakValley == "peak"   ~ 0L,
        PeakValley == "valley" ~ 1L,
        TRUE                   ~ NA_integer_
      )
    )
  
  # Sanity: make sure smoothed column is present
  if (!smooth_col %in% names(df_join)) {
    stop("Column '", smooth_col,
         "' not found in df_join. Names are: ",
         paste(names(df_join), collapse = ", "))
  }
  
  # 6. Drawdown / recharge / predicted-no-HR calcs
  df_out <- df_join %>%
    dplyr::mutate(
      drawdown = dplyr::if_else(PeakValley == "valley", delta_mm, NA_real_),
      recharge = dplyr::if_else(PeakValley == "peak",   delta_mm, NA_real_),
      total_recharge = cumsum(tidyr::replace_na(recharge, 0)),
      total_drawdown = cumsum(tidyr::replace_na(drawdown, 0)),
      mm_filled      = mm_h2o,
      perc_recharge  = total_recharge / 10
    ) %>%
    tidyr::fill(mm_filled, .direction = "down") %>%
    tidyr::fill(flag,      .direction = "down") %>%
    dplyr::mutate(
      vwc_adjusted       = .data[[smooth_col]] - perc_recharge,
      predicted_vwc_noHR = dplyr::case_when(
        flag == 1L ~ (mm_filled / 10) - perc_recharge,
        TRUE       ~ vwc_adjusted
      )
    )
  
  return(df_out)
}
