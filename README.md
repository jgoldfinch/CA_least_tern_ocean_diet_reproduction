<img width="960" height="768" alt="image" src="https://github.com/user-attachments/assets/96a49d23-d79f-4ad3-86c2-23a4e5d046e2" />

**California least tern (*Sternula antillarum browni*)**

*Photo credit: U.S. Fish and Wildlife Service, National Conservation Training Center Pubilcations and Training Material*

# California least tern (*Sternula antillarum browni*) ocean, diet, and reproduction analysis
A repository to recreate the analysis in the Mills et al. (2026, In Revision) publication examining the influence of ocean conditions on California least tern diet and reproductive success.

## Programs needed to run:
* R
* JAGS

## Code files needed to run:
* 00_tern_SEM_primary.R
* 01_tern_SEM_load_packages.R
* 02_tern_SEM_clean_data.R
* 03_tern_SEM_simmr_diet.R
* 04_tern_SEM_laavan_analysis.R
* 05_tern_SEM_extract_residuals.R
* 06_tern_SEM_plot_residuals.R
* 07_tern_SEM_ocean_summary.R"
* 8_tern_SEM_ocean_residuals_regression.R
* 09_tern_SEM_test_ocean_extension.R

## Data files needed to run:
* SEM_data_yearly.csv
* ALL_EGGS_STACKED.csv
* SIMMR_sources_FINAL.csv
* TEFs_Buckner_10.22.26.csv
* tern_large_scale_ocean.csv
* All data and code files hosted on https://doi.org/10.5281/zenodo.20435101


## Complete the following steps in this order:

1) Open the *00_tern_SEM_primary.R* file (all R code will be run from here and line numbers below refer to this file).

2) Set your working directory where you have all code files and data files saved (L18).

3) Load functions (L20-29).

4) Run SIMMR isotope mixing models and generate diet proportions (L35-51).

5) Format and clean SEM data (L53-L62).
   
6) Run and compare multiple iterations of multigroup SEM models (L64-L75).

7) Extract residuals and diagnostics (L78-99).

8) Examine large-scale ocean metrics relationship to residuals (L102-143).
