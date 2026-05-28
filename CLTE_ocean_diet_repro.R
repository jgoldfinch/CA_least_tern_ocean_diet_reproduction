# R script for California least tern 2026 ocean, diet, reproduction analysis
# Author: Erica Mills
# Date: 05/20/2026

################################ Load Packages #################################

library(readr)
library(tidyr)
library(dplyr)
library(piecewiseSEM)

############################ Set Working Directory #############################

setwd() # Insert file path that contains "SEM_data_yearly.csv" in ()

################################# Load Data ####################################

SEM_data_yearly <- read.csv("SEM_data_yearly.csv")
SEM_data_yearly_MCBCP <- subset(SEM_data_yearly, site == "MCBCP")
SEM_data_yearly_NBC <- subset(SEM_data_yearly, site == "NBC")

########################### Scale all variables ################################

# Pendleton
SEM_data_yearly_MCBCP <- SEM_data_yearly_MCBCP  %>%
  mutate(across(
    c(
      fledge_rate,
      diet_prop_C_enriched,
      d15N_avg,
      pCUI_seasonal,
      HCI_seasonal,
      pCUI_lagged,
      HCI_lagged
    ),
    ~ as.numeric(scale(.)),
    .names = "{.col}_z"
  ))

# Coronado
SEM_data_yearly_NBC <- SEM_data_yearly_NBC  %>%
  mutate(across(
    c(
      fledge_rate,
      diet_prop_C_enriched,
      d15N_avg,
      pCUI_seasonal,
      HCI_seasonal,
      pCUI_lagged,
      HCI_lagged
    ),
    ~ as.numeric(scale(.)),
    .names = "{.col}_z"
  ))

############################## Build Models ####################################

## Pendleton model with breeding year and prior year ocean metrics 
yearly_MCBCP_both_timescales_pSEM <- psem(
  glm(
    fledge_rate_z ~ diet_prop_C_enriched_z + d15N_avg_z,
    data = SEM_data_yearly_MCBCP,
    na.action = na.exclude
  ),
  glm(
    diet_prop_C_enriched_z ~ pCUI_seasonal_z + HCI_lagged_z,
    data = SEM_data_yearly_MCBCP,
    na.action = na.exclude
  ),
  glm(
    d15N_avg_z ~ HCI_seasonal_z + pCUI_lagged_z,
    data = SEM_data_yearly_MCBCP,
    na.action = na.exclude
  )
)

summary(yearly_MCBCP_both_timescales_pSEM) 

# d-sep test shows significant missing pathway between d15N and pCUI_lagged
yearly_MCBCP_both_timescales_pSEM2 <- psem(
  glm(
    fledge_rate_z ~ diet_prop_C_enriched_z + d15N_avg_z,
    data = SEM_data_yearly_MCBCP,
    na.action = na.exclude
  ),
  glm(
    diet_prop_C_enriched_z ~ pCUI_seasonal_z + HCI_lagged_z,
    data = SEM_data_yearly_MCBCP,
    na.action = na.exclude
  ),
  glm(
    d15N_avg_z ~ HCI_seasonal_z + pCUI_seasonal_z + pCUI_lagged_z,
    data = SEM_data_yearly_MCBCP,
    na.action = na.exclude
  )
)

summary(yearly_MCBCP_both_timescales_pSEM2)

## Coronado model with breeding year and prior year ocean metrics
yearly_NBC_both_timescales_pSEM <- psem(
  glm(
    fledge_rate_z ~ diet_prop_C_enriched_z + d15N_avg_z,
    data = SEM_data_yearly_NBC,
    na.action = na.exclude
  ),
  glm(
    diet_prop_C_enriched_z ~ pCUI_seasonal_z + HCI_lagged_z,
    data = SEM_data_yearly_NBC,
    na.action = na.exclude
  ),
  glm(
    d15N_avg_z ~ HCI_seasonal_z + pCUI_lagged_z,
    data = SEM_data_yearly_NBC,
    na.action = na.exclude
  )
)

summary(yearly_NBC_both_timescales_pSEM)

#d-sep test shows significant missing pathway between d15N and C-enriched diet 
yearly_NBC_both_timescales_pSEM2 <- psem(
  glm(
    fledge_rate_z ~ diet_prop_C_enriched_z + d15N_avg_z,
    data = SEM_data_yearly_NBC,
    na.action = na.exclude
  ),
  glm(
    diet_prop_C_enriched_z ~ pCUI_seasonal_z + HCI_lagged_z,
    data = SEM_data_yearly_NBC,
    na.action = na.exclude
  ),
  glm(
    d15N_avg_z ~ diet_prop_C_enriched_z + HCI_seasonal_z + pCUI_lagged_z,
    data = SEM_data_yearly_NBC,
    na.action = na.exclude
  )
)

summary(yearly_NBC_both_timescales_pSEM2)

# d-sep test shows significant (relaxed significance assumption in manuscript
# of p <= 0.1) missing pathway  between d15N and HCI lagged
yearly_NBC_both_timescales_pSEM3 <- psem(
  glm(
    fledge_rate_z ~ diet_prop_C_enriched_z + d15N_avg_z,
    data = SEM_data_yearly_NBC,
    na.action = na.exclude
  ),
  glm(
    diet_prop_C_enriched_z ~ pCUI_seasonal_z + HCI_lagged_z,
    data = SEM_data_yearly_NBC,
    na.action = na.exclude
  ),
  glm(
    d15N_avg_z ~ diet_prop_C_enriched_z + HCI_seasonal_z + pCUI_lagged_z + HCI_lagged_z,
    data = SEM_data_yearly_NBC,
    na.action = na.exclude
  )
)

summary(yearly_NBC_both_timescales_pSEM3)
