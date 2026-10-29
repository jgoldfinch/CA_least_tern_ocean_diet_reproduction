# Project: Tern SEM
# Author: Jessie Golding, Erica Mills, Gemini
# Date: 10/28/2026

# Code to run analyses for the manuscript
# Mills, Golding, and Lewison (2026/2027)
# Climate-mediated diet shifts shape reproductive outcomes for an endangered seabird
# Ecology and Evolution

# Primary File
# Run all data loading, cleaning, analysis, and output creation from here

################################################################################
## 00. Setup

# Change filepath in line below with the appropriate file path for your file system
# where data files and R scripts are saved
setwd("filepath")

# Load functions
source("./01_tern_SEM_load_packages.R")
source("./02_tern_SEM_clean_data.R")
source("./03_tern_SEM_simmr_diet.R")
source("./04_tern_SEM_laavan_analysis.R")
source("./05_tern_SEM_extract_residuals.R")
source("./06_tern_SEM_plot_residuals.R")
source("./07_tern_SEM_ocean_summary.R")
source("./08_tern_SEM_ocean_residuals_regression.R")
source("./09_tern_SEM_test_ocean_extension.R")

# NOTE: You need JAGS to run this code. It is free to download at 
# https://sourceforge.net/projects/mcmc-jags/

################################################################################
## 01. Run SIMMR Isotope Mixing Models & Generate Diet Proportions

# Load raw isotopic data
ALL_EGGS_STACKED     <- read.csv("ALL_EGGS_STACKED.csv")
simmr_sources_FINAL  <- read.csv("SIMMR_sources_FINAL.csv")
TEFs_Buckner         <- read.csv("TEFs_Buckner_10.22.26.csv")

# Run SIMMR estimation, combining, and export carbon diet proportions
carbon_proportions_table <- run_simmr_diet_pipeline(
  all_eggs_df   = ALL_EGGS_STACKED,
  sources_df    = simmr_sources_FINAL,
  tefs_df       = TEFs_Buckner
)

# Creates "carbon_proportions_table" which contains carbon enriched diet percentage data
# for both sites

################################################################################
## 02. Format and Clean SEM Data (Global Dataset Scaling)

# Load raw productivity (fledge rate), isotope, and ocean metric data for 2 sites
SEM_data_yearly <- read.csv("SEM_data_yearly.csv")

# Run cleaning and standardization function
clean_SEM_yearly_data(SEM_data_yearly, scale_by_group = FALSE)

# Creates "SEM_data_yearly_clean" file for analysis

################################################################################
## 03. Step-by-Step Multigroup SEM Model Comparisons

# Run stepwise multigroup SEM analysis
sem_stepwise_results <- run_stepwise_multigroup_sem(d = SEM_data_yearly_clean, group_var = "site")

# View summary table (for manuscript supplement)
print(sem_stepwise_results$comparison_table)
write.csv(sem_stepwise_results$comparison_table, "./sem_model_comparisons_supplement.csv", row.names = FALSE)

# Creates sem_stepwise_results, a list with model results ("model"), and 
# the result comparison table ("comparison_table"), which contains info in supplement
# Table S4 in the manuscript.

################################################################################
## 04. Residual Extraction & Diagnostics

resids_fledge <- extract_sem_residuals(fit_object = fit_final, target_var = "fledge_rate_z", type = "casewise")
resids_c      <- extract_sem_residuals(fit_object = fit_final, target_var = "diet_prop_C_enriched_z", type = "casewise")
resids_n      <- extract_sem_residuals(fit_object = fit_final, target_var = "d15N_avg_z", type = "casewise")

custom_names <- c(
  "fledge_rate_z"          = "Fledge Rate",
  "diet_prop_C_enriched_z" = "% C-Enriched Diet",
  "d15N_avg_z"             = "D15N"
)

plots_fledge <- plot_sem_residuals(resids_fledge, var_labels = custom_names)
plots_c      <- plot_sem_residuals(resids_c, var_labels = custom_names)
plots_n      <- plot_sem_residuals(resids_n, var_labels = custom_names)

print(plots_fledge$raw_plot)
print(plots_fledge$acf_plot)

# Creates residual lists (resids_fledge, resids_c, resids_n) and lists with 
# plots of residuals (plots_fledge, plots_c, plots_n) that contain raw and summarized plots
# Summarized plots are in Supplement Figure S4

################################################################################
## 05. Ocean Metrics & Residual Regressions

tern_large_scale_ocean <- read.csv("tern_large_scale_ocean.csv")
winter_ocean_df <- calc_winter_ocean_avg(
  df          = tern_large_scale_ocean, 
  date_col    = "date", 
  metric_cols = c("enso", "pdo", "noi")
)

ocean_vars <- c("winter_enso", "winter_pdo", "winter_noi")

reg_c <- do.call(rbind, lapply(ocean_vars, function(m) {
  regress_residuals_ocean(resids_c, winter_ocean_df, fit_final, SEM_data_yearly_clean, m, "diet_prop_C_enriched_z")
}))
reg_c$variable <- "C"

reg_n <- do.call(rbind, lapply(ocean_vars, function(m) {
  regress_residuals_ocean(resids_n, winter_ocean_df, fit_final, SEM_data_yearly_clean, m, "d15N_avg_z")
}))
reg_n$variable <- "d15N"

reg_fledge <- do.call(rbind, lapply(ocean_vars, function(m) {
  regress_residuals_ocean(resids_fledge, winter_ocean_df, fit_final, SEM_data_yearly_clean, m, "fledge_rate_z")
}))
reg_fledge$variable <- "fledge_r"

all_results <- rbind(reg_c, reg_n, reg_fledge)
write.csv(all_results, "./all_results.csv", row.names = FALSE)

# Creates regression models for residuals for SEM model residuals relative to 
# larger scale ocean variables. Information in manuscript Supplement Table S5.

################################################################################
## 06. Ocean Extension SEM Tests

ocean_extension_results <- test_ocean_sem_extension(
  d         = SEM_data_yearly_clean,
  ocean_df  = winter_ocean_df,
  group_var = "site",
  year_col  = "year"
)
print(ocean_extension_results$fit_table)