# Project: Tern SEM
# Author: Jessie Golding, Erica Mills, Gemini
# Date: 10/28/2026

#################################### Intro #####################################

# Function name: calc_winter_ocean_avg
# Description:  Calculates winter averages (Jan–Mar) for large-scale ocean metrics.

################################## Function #####################################
calc_winter_ocean_avg <- function(df, date_col = "date", metric_cols = c("enso", "pdo", "noi")) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (!requireNamespace("lubridate", quietly = TRUE)) stop("Package 'lubridate' is required.")
  
  df %>%
    dplyr::mutate(
      parsed_date = lubridate::parse_date_time(.data[[date_col]], orders = c("mdy", "ymd", "dmy", "Ymd")),
      year  = lubridate::year(parsed_date),
      month = lubridate::month(parsed_date)
    ) %>%
    dplyr::filter(month %in% c(1, 2, 3)) %>%
    dplyr::group_by(year) %>%
    dplyr::summarize(dplyr::across(dplyr::all_of(metric_cols), ~ mean(.x, na.rm = TRUE), .names = "winter_{.col}"), .groups = "drop")
} 