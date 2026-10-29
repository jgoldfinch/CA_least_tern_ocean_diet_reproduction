# Project: Tern SEM
# Author: Jessie Golding, Erica Mills, Gemini
# Date: 10/28/2026

#################################### Intro #####################################

# Function name: clean_SEM_yearly_data
# Description:  Scales continuous variables globally across the full dataset 
#               (or by site if scale_by_group = TRUE), appends _z scores, and 
#               assigns data frames to the Global Environment.

################################## Function #####################################
clean_SEM_yearly_data <- function(d, scale_by_group = FALSE) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (!("site" %in% colnames(d))) stop("The input data frame must contain a 'site' column.")
  
  target_cols <- c("fledge_rate", "diet_prop_C_enriched", "d15N_avg", 
                   "pCUI_seasonal", "HCI_seasonal", "pCUI_lagged", "HCI_lagged")
  
  if (scale_by_group) {
    scaled_data <- d %>%
      dplyr::group_by(site) %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(target_cols), ~ as.numeric(scale(.)), .names = "{.col}_z")) %>%
      dplyr::ungroup()
  } else {
    scaled_data <- d %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(target_cols), ~ as.numeric(scale(.)), .names = "{.col}_z"))
  }
  
  assign("SEM_data_yearly_clean", scaled_data, envir = .GlobalEnv)
  
  unique_sites <- sort(unique(scaled_data$site))
  for (s in unique_sites) {
    assign(paste0("SEM_data_yearly_", s), dplyr::filter(scaled_data, site == s), envir = .GlobalEnv)
  }
}