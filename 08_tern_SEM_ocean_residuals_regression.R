# Project: Tern SEM
# Author: Jessie Golding, Erica Mills, Gemini
# Date: 10/28/2026

#################################### Intro #####################################

# Function name: regress_residuals_ocean
# Description:  Regresses model residuals against ocean climate metrics.

################################## Function #####################################
regress_residuals_ocean <- function(sem_resids, ocean_df, fit_object, source_data, metric_col, resid_var = "fledge_rate_z", year_col = "year") {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (!requireNamespace("lmtest", quietly = TRUE)) stop("Package 'lmtest' is required.")
  
  case_indices <- lavaan::lavInspect(fit_object, "case.idx")
  grp_labels   <- lavaan::lavInspect(fit_object, "group.label")
  if (!is.list(case_indices)) { case_indices <- list(case_indices); grp_labels <- "Group_1" }
  
  summary_list <- list()
  for (i in seq_along(sem_resids)) {
    site <- names(sem_resids)[i]
    df_site_resids <- sem_resids[[i]]
    grp_idx <- match(site, grp_labels)
    if (is.na(grp_idx)) grp_idx <- i
    
    site_source <- source_data[case_indices[[grp_idx]], ]
    df_site_resids$year_extracted <- site_source[[year_col]]
    
    df_combined <- df_site_resids %>%
      dplyr::select(year = year_extracted, residual = dplyr::all_of(resid_var)) %>%
      dplyr::inner_join(ocean_df, by = "year") %>%
      dplyr::filter(!is.na(residual) & !is.na(.data[[metric_col]]))
    
    if (nrow(df_combined) < 4) next
    
    fit_lm   <- stats::lm(stats::as.formula(paste("residual ~", metric_col)), data = df_combined)
    lm_sum   <- summary(fit_lm)
    dw_test  <- lmtest::dwtest(fit_lm)
    coef_mat <- lm_sum$coefficients
    
    summary_list[[site]] <- data.frame(
      site = site, predictor = metric_col, n_obs = nrow(df_combined),
      estimate = coef_mat[2, "Estimate"], std_error = coef_mat[2, "Std. Error"],
      p_value = coef_mat[2, "Pr(>|t|)"], r_squared = lm_sum$r.squared,
      dw_statistic = as.numeric(dw_test$statistic), dw_p_value = dw_test$p.value,
      stringsAsFactors = FALSE
    )
  }
  return(do.call(rbind, summary_list))
}