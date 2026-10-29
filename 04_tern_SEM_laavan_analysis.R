# Project: Tern SEM
# Author: Jessie Golding, Erica Mills, Gemini
# Date: 10/28/2026

#################################### Intro #####################################

# Function name: run_stepwise_multigroup_sem
# Description:  Evaluates all nested model constraints step-by-step (free, single 
#               constrained pathways, 5 constrained pathways, fully constrained), 
#               generates LRT statistics for the manuscript supplement, and returns 
#               fitted lavaan model objects.

################################## Function #####################################
run_stepwise_multigroup_sem <- function(d, group_var = "site") {
  if (!requireNamespace("lavaan", quietly = TRUE)) stop("Package 'lavaan' is required.")
  
  d[[group_var]] <- factor(d[[group_var]])
  
  # Base Model Syntax
  m_free <- '
    fledge_rate_z ~ diet_prop_C_enriched_z + d15N_avg_z
    diet_prop_C_enriched_z ~ pCUI_seasonal_z + HCI_lagged_z
    d15N_avg_z ~ HCI_seasonal_z + pCUI_lagged_z
  '
  fit_free <- lavaan::sem(m_free, data = d, group = group_var, meanstructure = TRUE)
  
  # Single path equality models
  fit_diet_equal     <- lavaan::sem('fledge_rate_z ~ c(b,b)*diet_prop_C_enriched_z + d15N_avg_z\n diet_prop_C_enriched_z ~ pCUI_seasonal_z + HCI_lagged_z\n d15N_avg_z ~ HCI_seasonal_z + pCUI_lagged_z', data = d, group = group_var)
  fit_d15N_equal     <- lavaan::sem('fledge_rate_z ~ diet_prop_C_enriched_z + c(b,b)*d15N_avg_z\n diet_prop_C_enriched_z ~ pCUI_seasonal_z + HCI_lagged_z\n d15N_avg_z ~ HCI_seasonal_z + pCUI_lagged_z', data = d, group = group_var)
  fit_pCUI_diet_equal<- lavaan::sem('fledge_rate_z ~ diet_prop_C_enriched_z + d15N_avg_z\n diet_prop_C_enriched_z ~ c(b,b)*pCUI_seasonal_z + HCI_lagged_z\n d15N_avg_z ~ HCI_seasonal_z + pCUI_lagged_z', data = d, group = group_var)
  fit_HCI_diet_equal <- lavaan::sem('fledge_rate_z ~ diet_prop_C_enriched_z + d15N_avg_z\n diet_prop_C_enriched_z ~ pCUI_seasonal_z + c(b,b)*HCI_lagged_z\n d15N_avg_z ~ HCI_seasonal_z + pCUI_lagged_z', data = d, group = group_var)
  fit_HCI_d15N_equal <- lavaan::sem('fledge_rate_z ~ diet_prop_C_enriched_z + d15N_avg_z\n diet_prop_C_enriched_z ~ pCUI_seasonal_z + HCI_lagged_z\n d15N_avg_z ~ c(b,b)*HCI_seasonal_z + pCUI_lagged_z', data = d, group = group_var)
  fit_pCUI_d15N_equal<- lavaan::sem('fledge_rate_z ~ diet_prop_C_enriched_z + d15N_avg_z\n diet_prop_C_enriched_z ~ pCUI_seasonal_z + HCI_lagged_z\n d15N_avg_z ~ HCI_seasonal_z + c(b,b)*pCUI_lagged_z', data = d, group = group_var)
  
  # Final model (5 constrained, 1 free)
  m_final <- '
    fledge_rate_z ~ diet_prop_C_enriched_z + c(d15N, d15N)*d15N_avg_z
    diet_prop_C_enriched_z ~ c(pCUI_d, pCUI_d)*pCUI_seasonal_z + c(HCI_d, HCI_d)*HCI_lagged_z
    d15N_avg_z ~ c(HCI_n, HCI_n)*HCI_seasonal_z + c(pCUI_n, pCUI_n)*pCUI_lagged_z
  '
  fit_final <- lavaan::sem(m_final, data = d, group = group_var)
  
  # Fully constrained model
  m_all_equal <- '
    fledge_rate_z ~ c(diet, diet)*diet_prop_C_enriched_z + c(d15N, d15N)*d15N_avg_z
    diet_prop_C_enriched_z ~ c(pCUI_d, pCUI_d)*pCUI_seasonal_z + c(HCI_d, HCI_d)*HCI_lagged_z
    d15N_avg_z ~ c(HCI_n, HCI_n)*HCI_seasonal_z + c(pCUI_n, pCUI_n)*pCUI_lagged_z
  '
  fit_all_equal <- lavaan::sem(m_all_equal, data = d, group = group_var)
  
  models <- list(
    "Fully unconstrained"                      = fit_free,
    "Diet -> Fledge constrained"               = fit_diet_equal,
    "d15N -> Fledge constrained"               = fit_d15N_equal,
    "pCUI -> Diet constrained"                 = fit_pCUI_diet_equal,
    "HCI -> Diet constrained"                  = fit_HCI_diet_equal,
    "HCI -> d15N constrained"                  = fit_HCI_d15N_equal,
    "pCUI -> d15N constrained"                 = fit_pCUI_d15N_equal,
    "Final: 5 constrained, Diet -> Fledge free"= fit_final,
    "Fully constrained"                        = fit_all_equal
  )
  
  comp_df <- data.frame(
    Model      = names(models),
    AIC        = sapply(models, function(x) lavaan::fitMeasures(x, "aic")),
    BIC        = sapply(models, function(x) lavaan::fitMeasures(x, "bic")),
    Chi_Square = sapply(models, function(x) lavaan::fitMeasures(x, "chisq")),
    df         = sapply(models, function(x) lavaan::fitMeasures(x, "df"))
  )
  
  comp_df$Delta_AIC <- round(comp_df$AIC - comp_df$AIC[1], 2)
  comp_df$Delta_BIC <- round(comp_df$BIC - comp_df$BIC[1], 2)
  
  get_lrt <- function(mod) {
    t <- lavaan::anova(fit_free, mod)
    c(Delta_ChiSq = t$`Chisq diff`[2], Delta_df = t$`Df diff`[2], P_val = t$`Pr(>Chisq)`[2])
  }
  
  lrt_res <- rbind(c(NA, NA, NA), t(sapply(models[-1], get_lrt)))
  comp_df <- cbind(comp_df, lrt_res)
  
  assign("fit_free", fit_free, envir = .GlobalEnv)
  assign("fit_final", fit_final, envir = .GlobalEnv)
  
  return(invisible(list(models = models, comparison_table = comp_df)))
}