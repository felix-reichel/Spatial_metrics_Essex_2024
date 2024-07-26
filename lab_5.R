# 4 practicals ----

rm(list = ls())
library(haven)
library(MASS)
library(spatialreg)
library(spdep)
library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## import homicide_South_1990.dta
## y: hrate
## X: ln_population, ln_pdensity, gini

data <- read_dta("data/homicide_South_1990.dta")
hist(data$hrate)

weights_matrix <- mat2listw(as.matrix(read_csv("data/W_South_1990.csv", col_select = -v1)), style = "W")

## 3.1 estimate SAR model via
### 3.1.1 naive OLS

data$W_hrate <- lag.listw(weights_matrix, data$hrate)
sar_ols_model <- lm(hrate ~ W_hrate + ln_population + ln_pdensity + gini, data)
summary(sar_ols_model)

### 3.1.2 maximum likelihood (ML)

sar_ml_model <- lagsarlm(hrate ~ ln_population + ln_pdensity + gini, data = data, listw = weights_matrix)
summary(sar_ml_model)

### 3.1.3 two-stage least squares with simple and robust SEs respectively — would point estimates become different?

sar_2sls_model <- stsls(hrate ~ ln_population + ln_pdensity + gini, data, weights_matrix)
summary(sar_2sls_model)

sar_2sls_robust_model <- stsls(hrate ~ ln_population + ln_pdensity + gini, data, weights_matrix, robust = TRUE)
summary(sar_2sls_robust_model)

## 3.2 estimate SEM model

sem_model <- errorsarlm(hrate ~ ln_population + ln_pdensity + gini, data, weights_matrix)
summary(sem_model)

## 3.3 estimate SLX (lagging all RHS variables) model with OLS using lm()

data <- data %>% mutate(across(c(ln_population, ln_pdensity, gini), ~ lag.listw(weights_matrix, .x), .names = "W_{.col}"))
slx_model <- lm(hrate ~ ln_population + ln_pdensity + gini + W_ln_population + W_ln_pdensity + W_gini, data)
summary(slx_model)

## 3.4 estimate the following models
### 3.4.1 SAC

sac_model <- sacsarlm(hrate ~ ln_population + ln_pdensity + gini, data, weights_matrix)
summary(sac_model)

### 3.4.2 SDM

sdm_model <- lagsarlm(hrate ~ ln_population + ln_pdensity + gini, data, weights_matrix, Durbin = TRUE)
summary(sdm_model)

### 3.4.3 SDEM

sdem_model <- errorsarlm(hrate ~ ln_population + ln_pdensity + gini, data, weights_matrix, Durbin = TRUE)
summary(sdem_model)

### 3.4.4 GNS

gns_model <- sacsarlm(hrate ~ ln_population + ln_pdensity + gini, data, weights_matrix, Durbin = TRUE)
summary(gns_model)

## 3.5 compare different models (using SAR-ML for the SAR model)

compare_models <- function(model1, model2) {
  list(
    LR_Test = LR.Sarlm(model1, model2),
    AIC = AIC(model1, model2)
  )
}

compare_models(gns_model, sac_model)
compare_models(gns_model, sdem_model)
compare_models(gns_model, sdm_model)
compare_models(gns_model, sar_ml_model)
compare_models(gns_model, slx_model)
compare_models(gns_model, sem_model)

compare_models(sac_model, sar_ml_model)
compare_models(sac_model, sem_model)

compare_models(sdem_model, slx_model)
compare_models(sdem_model, sem_model)

compare_models(sdm_model, slx_model)
compare_models(sdm_model, sar_ml_model)

## 3.6 draw 1,000 rho (based upon the SAR-ML estimate) from the distribution it's 
#assumed to follow under the null hypothesis and plot the histogram

set.seed(2022207)
rho_values <- mvrnorm(1000, coef(sar_ml_model), vcov(sar_ml_model))[,"rho"]
hist(rho_values, main = "Rho Parametric Simulation")

