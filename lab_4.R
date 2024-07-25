# 3 practicals ----

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

library(cshapes)
library(countrycode)
library(haven)
library(spdep)
library(sf)
library(tidyverse)

## 3.1 ----
### load Brexit_chippy.RData (Breixt_chippy — data, sf_GBR — simple features, 
#W_GB_constituency — spatial weights matrix)

load("data/Brexit_chippy.RData")


### 3.1.1 use the code given below to drop "islands"
#### islands <- rownames(W_GB_constituency[rowSums(W_GB_constituency) == 0,])
#### Brexit_chippy <- filter(Brexit_chippy, !(pcon15cd %in% islands))
#### W_GB_constituency <- W_GB_constituency[!(rownames(W_GB_constituency) %in% islands), !(colnames(W_GB_constituency) %in% islands)]
#### sf_GBR <- filter(sf_GBR, !(pcon15cd %in% islands))

islands <- rownames(W_GB_constituency[rowSums(W_GB_constituency) == 0,])
Brexit_chippy <- filter(Brexit_chippy, !(pcon15cd %in% islands))
W_GB_constituency <- W_GB_constituency[!(rownames(W_GB_constituency) %in% islands), !(colnames(W_GB_constituency) %in% islands)]
sf_GBR <- filter(sf_GBR, !(pcon15cd %in% islands))

### 3.1.2 create a weights list, row-standardized
weights_list <- mat2listw(W_GB_constituency, style = "W")

### 3.1.3 calculate Moran's I (analytical), Geary's C, and Getis-Ord G for the following two variables: leave and chippy
### 3.1.4 perform Monte Carlo-based Moran's I test for these two variables with 10,000 simulations and visually present the statistical inference
### 3.1.5 create two Moran scatterplots for these two variables

variables <- c("leave", "chippy")

for (var in variables) {
  moran.test(Brexit_chippy[[var]], weights_list)
  geary.test(Brexit_chippy[[var]], weights_list)
  globalG.test(Brexit_chippy[[var]], weights_list)
  plot(moran.mc(Brexit_chippy[[var]], weights_list, nsim = 10000))
  moran.plot(Brexit_chippy[[var]], weights_list)
}

### 3.1.6 create spatial lags for these two variables and replicate the two Moran scatterplots 
#(fitted line not needed) only using the base R plot function
Brexit_chippy <- Brexit_chippy %>%
  mutate(across(all_of(variables), ~ lag.listw(weights_list, .), .names = "W{.col}"))

plot(Brexit_chippy$leave, Brexit_chippy$Wleave, xlab = "leave", ylab = "Wleave", main = "Spatial Autocorrelation of Leave Votes")
plot(Brexit_chippy$chippy, Brexit_chippy$Wchippy, xlab = "chippy", ylab = "Wchippy", main = "Spatial Autocorrelation of Chippies")


### 3.1.7 **perform local Moran's I test (right-tailed) for leave, access attributes of the 
#resultant object, visualize the four types of local spatial autocorrelation (corresponding to four quadrants on a Moran scatterplot) 
#on the map (using sf_GBR)

quadr_localI <- localmoran(Brexit_chippy$leave, weights_list, alternative = "greater") %>%
  attr("quadr")

bind_cols(sf_GBR, quadr_localI) %>%
  ggplot() + geom_sf(aes(fill = mean), alpha = 0.5, size = 0.25) +
  scale_fill_manual(values = c("red", "pink", "white", "cyan", "steelblue"), name = "LISA") +
  theme_void() + 
  theme(text = element_text(family = "Georgia"))


## 3.2 ----

rm(list = ls())
dev.off()

### import (read_dta from the haven package) senate_1840_votes.dta (data) as well as senate_1840_floor.dta (W)

senate_data <- read_dta("data/senate_1840_votes.dta")
senate_weights <- as.matrix(read_dta("data/senate_1840_floor.dta"))

### 3.2.1 use senate_1840_votes.dta to regress (OLS) nom1 on population, dem_2party_inv, mfct_pct, and V61
#### nom1: a senator's score on the first dimension of the scaled NOMINATE measure
#### population: a state's population in the 1840 census
#### dem_2party_inv: Republican share of the two-party vote
#### mfct_pct: percentage of the population employed in manufacturing and trade sectors
#### V61: a senator's total years of senatorial service

print(table(rowSums(senate_weights)))  # 30 senators had two colleagues, 12 had one
print(isSymmetric(senate_weights))


### 3.2.2 perform Moran's I test on the regression residuals with a row-standardized W from senate_1840_floor

model <- lm(nom1 ~ population + dem_2party_inv + mfct_pct + V61, data = senate_data)
lm.morantest(model, mat2listw(senate_weights, style = "W"))


### 3.2.3 perform robust Lagrange Multiplier tests for spatially lagged error and spatially lagged 
#dependent variable respectively with the same W

lm.LMtests(model, mat2listw(senate_weights, style = "W"), test = c("RLMerr", "RLMlag"))


