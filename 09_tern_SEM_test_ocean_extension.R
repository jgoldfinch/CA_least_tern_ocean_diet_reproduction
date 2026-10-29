# Project: Tern SEM
# Author: Jessie Golding, Erica Mills, Gemini
# Date: 10/28/2026

#################################### Intro #####################################

# Function name: test_ocean_sem_extension
# Description:  Tests SEMExtensions adding direct paths for winter PDO and NOI.

################################## Function #####################################
test_ocean_sem_extension <- function(d, ocean_df, group_var = "site", year_col = "year") {
  if (!requireNamespace("lavaan", quietly = TRUE)) stop("Package 'lavaan' is required.")
  
  d_merged <- dplyr::inner_join(d, ocean_df, by = year_col)
  
  syntax_base <- '
    fledge_rate_z ~ diet_prop_C_enriched_z + d15N_avg_z
    diet_prop_C_enriched_z ~ pCUI_seasonal_z + HCI_lagged_z
    d15N_avg_z ~ HCI_seasonal_z + pCUI_lagged_z
  '
  syntax_pdo  <- paste(syntax_base, "\nfledge_rate_z ~ winter_pdo")
  syntax_noi  <- paste(syntax_base, "\nfledge_rate_z ~ winter_noi")
  syntax_both <- paste(syntax_base, "\nfledge_rate_z ~ winter_pdo + winter_noi")
  
  fit_base <- lavaan::sem(syntax_base, data = d_merged, group = group_var, meanstructure = TRUE, fixed.x = FALSE)
  fit_pdo  <- lavaan::sem(syntax_pdo,  data = d_merged, group = group_var, meanstructure = TRUE, fixed.x = FALSE)
  fit_noi  <- lavaan::sem(syntax_noi,  data = d_merged, group = group_var, meanstructure = TRUE, fixed.x = FALSE)
  fit_both <- lavaan::sem(syntax_both, data = d_merged, group = group_var, meanstructure = TRUE, fixed.x = FALSE)
  
  extract_metrics <- function(fit_obj, name) {
    m <- lavaan::fitMeasures(fit_obj, c("aic", "bic", "cfi", "rmsea", "chisq", "df", "pvalue"))
    data.frame(Model = name, AIC = round(m["aic"], 2), BIC = round(m["bic"], 2), CFI = round(m["cfi"], 3), RMSEA = round(m["rmsea"], 3), stringsAsFactors = FALSE)
  }
  
  comp_df <- rbind(
    extract_metrics(fit_base, "0. Base Model"),
    extract_metrics(fit_pdo,  "1. Base + Winter PDO"),
    extract_metrics(fit_noi,  "2. Base + Winter NOI"),
    extract_metrics(fit_both, "3. Base + PDO + NOI")
  )
  
  return(invisible(list(fit_base = fit_base, fit_pdo = fit_pdo, fit_noi = fit_noi, fit_both = fit_both, fit_table = comp_df)))
}