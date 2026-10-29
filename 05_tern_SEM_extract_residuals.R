# Project: Tern SEM
# Author: Jessie Golding, Erica Mills, Gemini
# Date: 10/28/2026

#################################### Intro #####################################

# Function name: extract_sem_residuals
# Description:  Extracts casewise residuals from fitted lavaan objects.

################################## Function #####################################
extract_sem_residuals <- function(fit_object = NULL, target_var = "all", type = "casewise") {
  if (!requireNamespace("lavaan", quietly = TRUE)) stop("Package 'lavaan' is required.")
  if (is.null(fit_object)) fit_object <- get("fit_final", envir = .GlobalEnv)
  
  obs_data   <- lavaan::lavInspect(fit_object, "data")
  grp_labels <- lavaan::lavInspect(fit_object, "group.label")
  if (!is.list(obs_data)) { obs_data <- list(obs_data); grp_labels <- "Group_1" }
  
  pe <- lavaan::parameterEstimates(fit_object)
  
  resids <- lapply(seq_along(obs_data), function(g) {
    data_g <- as.data.frame(obs_data[[g]])
    pe_g   <- dplyr::filter(pe, group == g | is.na(group))
    res_df <- data_g
    endog_vars <- unique(pe_g$lhs[pe_g$op == "~"])
    
    for (v in endog_vars) {
      intercept_row <- pe_g[pe_g$lhs == v & pe_g$op == "~1", ]
      intercept     <- if (nrow(intercept_row) > 0) intercept_row$est else 0
      paths         <- pe_g[pe_g$lhs == v & pe_g$op == "~", ]
      hat_y         <- rep(intercept, nrow(data_g))
      for (i in seq_len(nrow(paths))) {
        pred_name  <- paths$rhs[i]
        slope_coef <- paths$est[i]
        if (pred_name %in% colnames(data_g)) hat_y <- hat_y + (slope_coef * data_g[[pred_name]])
      }
      res_df[[v]] <- data_g[[v]] - hat_y
    }
    res_df$site_label <- grp_labels[g]
    return(res_df)
  })
  names(resids) <- grp_labels
  
  if (target_var != "all") {
    resids <- lapply(resids, function(df) df[, intersect(c(target_var, "year", "site_label"), colnames(df)), drop = FALSE])
  }
  return(resids)
}