# Project: Tern SEM
# Author: Jessie Golding, Erica Mills, Gemini
# Date: 10/28/2026

#################################### Intro #####################################

# Function name: load_tern_sem_packages
# Description:  Checks for, installs (if missing), and loads required R packages.

################################## Function #####################################
load_tern_sem_packages <- function() {
  packages <- c("tidyverse", "piecewiseSEM", "ggplot2", "lme4", "lavaan", "simmr", 
                "lmtest", "lubridate", "tibble")
  
  lapply(packages, function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  })
}

load_tern_sem_packages()