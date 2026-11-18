#' Function to decompose monthly atmospheric time series and prepare data for 
#' carbon budgeting.
#' 
#' @author Stuart K. Grange
#' 
#' @param df A data frame or tibble of monthly atmospheric observations 
#' containing columns `id`, `date`, `variable` (e.g., `"co2"`, `"o2"`), and
#' `value`.
#' 
#' @param df_emissions A data frame or tibble of annual global carbon emissions 
#' (`year`, `carbon_emissions_fossil_fuels`, `carbon_emissions_land_use_change`) 
#' to join by year.
#' 
#' @param co2_uncertainty Numeric value representing the 1-sigma uncertainty 
#' standard deviation to apply to \eqn{\text{CO}_2} observations via random normal 
#' sampling. Default is `0`.
#' 
#' @param o2_uncertainty Numeric value representing the 1-sigma uncertainty 
#' standard deviation to apply to \eqn{\text{O}_2} observations. Default is `0`.
#' 
#' @param span Numeric smoother span parameter. 
#' 
#' @param n_min Integer specifying the minimum number of valid monthly 
#' observations required within a year; annual means derived from fewer 
#' observations are set to `NA`.
#' 
#' @param invalidate_first_year Logical. If `TRUE`, sets annual observations for 
#' the initial starting year of each site (`id`) to `NA` to avoid edge artefacts 
#' from trend smoothing. 
#' 
#' @return Tibble. 
#' 
#' @seealso [apply_o2_co2_partitioning()], 
#' [apply_o2_co2_partitioning_with_table()]
#' 
#' @export
decompose_time_series_and_prepare_for_budget <- function(
    df, df_emissions, co2_uncertainty = 0, o2_uncertainty = 0, span = 1.5,
    n_min = 10L, invalidate_first_year = TRUE) {
  
  # Check for monthly resolution time series
  stopifnot(detect_date_interval(df$date, text_return = TRUE) == "month")
  
  # Sample the values based on the given uncertainties, warning suppression is
  # for missing values
  df_sampled <- df |> 
    mutate(
      sd = dplyr::case_when(
        variable == "co2" ~ co2_uncertainty,
        variable == "o2" ~ o2_uncertainty
      )
    ) |> 
    mutate(value = suppressWarnings(rnorm(dplyr::n(), mean = value, sd)))
  
  # Decompose the monthly time series, smooth the trend component
  df_trends <- df_sampled |> 
    dplyr::nest_by(id, 
                   variable,
                   .key = "observations") |> 
    mutate(
      trends = list(
        decompose_with_stlplus_extra(
          observations, window = 35, na_preserve = FALSE, span = span
        )
      )
    ) |> 
    reframe(trends)
  
  # Calculate annual means from smooth trend component and join emission
  # data for budget calculations
  df_trends_annual <- df_trends |> 
    select(id,
           variable,
           date,
           value = trend_smooth) |> 
    threadr::calculate_date_summaries(id,
                                      variable,
                                      interval = "year") |> 
    mutate(value = if_else(n <= !!n_min, NA_real_, value)) |> 
    select(-n) |> 
    tidyr::pivot_wider(names_from = variable) |> 
    mutate(year = lubridate::year(date)) |>
    left_join(df_emissions, by = join_by(year)) |>
    relocate(year,
             .after = date_end)
  
  # Invalidate the first non-missing year
  if (invalidate_first_year) {
    
    # Filter and get starting non-missing year for each id
    df_first_years <- df_trends_annual |> 
      filter(!is.na(co2) | !is.na(o2)) |> 
      summarise(year_min = min(year, na.rm = TRUE),
                .by = id)
    
    # Join first years and invalidate these years
    df_trends_annual <- df_trends_annual |> 
      left_join(df_first_years, by = join_by(id)) |> 
      mutate(
        across(
          c(co2, o2), \(x) if_else(year <= year_min, NA_real_, x)
        )
      ) |> 
      select(-year_min)
    
  }
  
  return(df_trends_annual)
  
}


# Pulled from threadr
detect_date_interval <- function(date, skip = 1L, n = 100L, text_return = FALSE) {
  
  # Check
  if (!lubridate::is.POSIXt(date)) {
    cli::cli_abort("`date` must be a POSIXt date.")
  }
  
  # A catch for vectors with fewer elements than skip
  if (length(date) <= skip) skip <- 0L
  
  # Get vectors
  # Skip if needed
  date <- date[skip:length(date)]
  
  # Filter
  date <- head(date, n)
  date_lag <- dplyr::lag(date)
  
  # Calculate difference
  seconds <- difftime(date, date_lag, units = "secs")
  seconds <- as.integer(seconds)
  seconds <- seconds[!is.na(seconds)]
  
  # Most common occurrence
  seconds <- mode_average(seconds, na.rm = TRUE)
  
  if (!text_return) {
    return(seconds)
  } else {
    
    # Default
    period <- "unknown"
    
    # Missing-ness test, when length is one
    if (length(seconds) == 1L && is.na(seconds)) {
      return(period)
    }
    
    # Known periods
    period <- dplyr::case_when(
      all(seconds == 1L) ~ "second",
      all(seconds == 60L) ~ "minute",
      all(seconds == 300L) ~ "five_minute",
      all(seconds == 600L) ~ "ten_minute",
      all(seconds == 900L) ~ "fifteen_minute",
      all(seconds == 1800L) ~ "half_hour",
      all(seconds == 3600L) ~ "hour",
      all(seconds == 86400L) ~ "day",
      all(seconds %in% c(2419200L, 2678400L, 2592000L, 2505600L)) ~ "month",
      .default = "unknown"
    )
    
    # Return text
    return(period)
    
  }
  
}


decompose_with_stlplus <- function(df, window = 35, na_preserve = TRUE) {
  
  if (!"date" %in% names(df)) {
    cli::cli_abort("`date` must be present in input.")
  }
  
  if (!lubridate::is.POSIXct(df$date)) {
    cli::cli_abort("`date` must be a parsed date (POSIXct).")
  }
  
  if (!threadr::detect_date_interval(df$date, text_return = TRUE) == "month") {
    cli::cli_abort("Time series must be at monthly resolution.")
  }
  
  if (!"value" %in% names(df)) {
    cli::cli_abort("Variable to be decomposed must be named `value`.")
  }
  
  if (nrow(df) <= 3) {
    cli::cli_warn("Too few observations to decompose...")
    return(tibble())
  } else {
    
    df <- tryCatch({
      # Warning suppression is for when there are missing values, rather common
      suppressWarnings(
        df |> 
          data_frame_to_time_series() |> 
          stlplus::stlplus(s.window = window, s.degree = 0, robust = TRUE) |> 
          stl_to_data_frame(na_preserve = na_preserve) 
      )
    }, error = function(e) {
      cli::cli_warn("Decomposition failed...")
      tibble()
    })
    
  }
  
  return(df)
  
}


decompose_with_stlplus_extra <- function(df, window = 35, na_preserve = FALSE, 
                                         span = 0.75, as_long = FALSE) {
  
  # Decompose the time series and return many variables
  df_decompose <- df |> 
    decompose_with_stlplus(window = window, na_preserve = na_preserve) |> 
    mutate(missing = is.na(raw), 
           trend_and_seasonal = trend + seasonal) |> 
    relocate(missing,
             .before = raw)
  
  # Train model on trend component, not trend + remainder
  model <- stats::loess(
    trend ~ as.numeric(date), data = df_decompose, span = span
  )
  
  # Add prediction to input, no non-missing observations in prediction
  df_decompose <- mutate(
    df_decompose, trend_smooth = unname(stats::predict(model, df))
  )
  
  # Make long data
  if (as_long) {
    df_decompose <- df_decompose |> 
      tidyr::pivot_longer(-date:-missing, names_to = "value_type") |> 
      arrange(value_type,
              date)
  }
  
  return(df_decompose)
  
}


data_frame_to_time_series <- function(df, variable = "value", period = "month") {
  
  # Get start date
  date_start <- min(df$date)
  
  if (period == "month") {
    frequency <- 12L
  } else if (period == "day") {
    frequency <- 365L
  } else {
    cli::cli_abort("`frequency` not supported.")
  }
  
  # Create the time series object
  ts <- ts(
    data = pull(df, !!variable),
    start = c(lubridate::year(date_start), lubridate::month(date_start)),
    frequency = frequency,
    names = "value"
  )
  
  return(ts)
  
}


stl_to_data_frame <- function(stl, date_round = TRUE, na_preserve = FALSE) {
  
  if (inherits(stl, "stlplus")) {
    
    # Get dates and transform into POSIXct
    date <- stl$time
    date_offset <- (date - floor(date)) * 31536000
    date <- lubridate::ymd(floor(date), truncated = 2, tz = "UTC") + date_offset
    if (date_round) date <- lubridate::round_date(date, "month")
    
    # To tibble
    df <- stl$data |> 
      mutate(date = !!date,
             trend_and_remainder = trend + remainder) |> 
      as_tibble() |> 
      relocate(date) |> 
      select(-sub.labels,
             -weights)
    
    if (na_preserve) {
      df <- dplyr::mutate(
        df, across(-c(1, 2), \(x) if_else(is.na(raw), NA_real_, x))
      )
    }
    
  } else {
    cli::cli_abort("Input not recognised.")
  }
  
  return(df)
  
}


mode_average <- function(x, na.rm = FALSE) {
  
  # Drop missing elements
  if (na.rm) x <- na.omit(x)
  
  # Compute the mode of the vector
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
  
}
