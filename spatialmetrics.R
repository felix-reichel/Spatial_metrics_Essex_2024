# Author: Felix Reichel
library(sf)
library(spdep)
library(haven)
library(spatialreg)
library(fourPNO)

# Info: input data has to be placed at the same level as source file.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Q. 1
us_states <- st_read("cb_2017_us_state_500k/cb_2017_us_state_500k.shp")
# 56 states

# Since there are 56 states in the data 
# and there are only 48 states part of the lower 48, define them manuelly.
# retrieved from: https://worldpopulationreview.com/state-rankings/lower-48-states
lower_48 <- c(
  "WA", "MT", "ND", "MN", "WI", "MI", "NY", "VT", "NH", "ME",
  "RI", "MA", "CT", "NJ", "DE", "OH", "IL", "IN", "MD", "PA",
  "WV", "VA", "SC", "GA", "TN", "KY", "MO", "CO", "NV", "OR",
  "ID", "IA", "NE", "SD", "WY", "NC", "FL", "AL", "MS", "CA",
  "UT", "AZ", "NM", "OK", "LA", "TX", "KS", "AR"
)

contiguous_states <- us_states[us_states$STUSPS %in% lower_48, ] # filter for lower 48
centroids <- st_centroid(contiguous_states)  # compute centroids

dist_matrix <- st_distance(centroids)
dist_matrix <- as.matrix(dist_matrix) / 1000  # from meters to kilometers

inv_dist_matrix <- matrix(0, nrow = nrow(dist_matrix), ncol = ncol(dist_matrix))
for (i in 1:nrow(dist_matrix)) {
  for (j in 1:ncol(dist_matrix)) {
    if (i != j) {
      inv_dist_matrix[i, j] <- 1 / dist_matrix[i, j]
    }
  }
}

states <- contiguous_states$STUSPS
alabama_idx <- which(states == "AL")
westvirginia_idx <- which(states == "WV")
inv_dist_al_wv <- inv_dist_matrix[alabama_idx, westvirginia_idx]
inv_dist_al_wv
inv_dist_matrix[westvirginia_idx, alabama_idx] 


?st_distance
1/inv_dist_al_wv
# Notice: Using the Haversine package would result in a distance (between centroids) of 
# 859.0499 kilometres instead of 859.288 kilometres which would be more accurate,
# but this package 'sf' only computes the Euclidean and the dist_fun parameter is deprecated.

# Thus:
inv_dist_al_wv
# The inverse distance between Alabama and West Virginia is: 0.001163754.




# Q. 2
oxford_data <- read_dta("oxford_data.dta")  # social policy variables at US state level (contiguous lower 48).
oxford_w <- read_dta("oxford_w.dta")        # contains a row-standardized queen contiguity matrix for the contiguous lower 48 US states.

oxford_w <- oxford_w[, -1]  # remove index column

oxford_data$ben95_z <- scale(oxford_data$ben95)   # uses Z-score norm

listw <- mat2listw(as.matrix(oxford_w), style = "W") # row-standardized matrix 

moran_test <- moran.test(oxford_data$ben95_z, listw)

print(moran_test)

# Moran I test under randomisation

# data:  oxford_data$ben95_z  
# weights: listw    

# Moran I statistic standard deviate = 5.7757, p-value = 3.831e-09
# alternative hypothesis: greater
# sample estimates:
#  Moran I statistic       Expectation          Variance 
#        0.540686607      -0.021276596       0.009466741 

moran_plot <- moran.plot(as.vector(oxford_data$ben95_z), listw, main = "Moran's I Plot for AFDC Benefits")




# Q. 3

## Q. 3.1.
lm_tests <- lm.LMtests(lm(ben95 ~ rskpovpc + wage95 + instcoad + ipcfold + teitrend + match, data = oxford_data), listw)

print(lm_tests)

# Rao's score (a.k.a Lagrange multiplier) diagnostics for spatial dependence

# data:  
# model: lm(formula = ben95 ~ rskpovpc + wage95 + instcoad + ipcfold + teitrend + match,
# data = oxford_data)
# test weights: listw

# RSerr = 5.8447, df = 1, p-value = 0.01562

# p-value < 0.05 => Reject null => OLS residuals are spatial autocorrelated.


ols_model <- lm(ben95 ~ rskpovpc + wage95 + instcoad + ipcfold + teitrend + match, 
                data = oxford_data)

lm_tests <- lm.LMtests(ols_model, listw, test="all")

lm_tests
# RSerr = 5.8447, df = 1, p-value = 0.01562 => SAR/SAC
# RSlag = 11.606, df = 1, p-value = 0.0006574 => SAR/SAC

# Robust Tests
# adjRSerr = 0.71624, df = 1, p-value = 0.3974 => Insignificant => not SAC
# adjRSlag = 6.4774, df = 1, p-value = 0.01093 => Significant => SAR, 
#  because the robust LM-Lag test stays significant and would be sufficient.

# SARMA = 12.322, df = 2, p-value = 0.00211



## Q. 3.2.
hist(oxford_data$ben95)

# Fit a SAC
sac_model <- sacsarlm(ben95 ~ rskpovpc + wage95 + instcoad + ipcfold + teitrend + match, data = oxford_data, listw)
summary(sac_model)


effect_ideology <- impacts(sac_model, listw = listw, R = 1000) # 1000 monte carlo sims
effect_ideology

# Impact measures (sac, exact):
#               Direct     Indirect        Total
# rskpovpc   5.13409572    8.8007929    13.9348887
# wage95     0.05935024    0.1017373    0.1610876
# instcoad   1.50169572    2.5741852    4.0758809
# ipcfold    627.04353101  1074.8689878 1701.9125188
# teitrend   1.79986967    3.0853107    4.8851803
# match     -6.30137076   -10.8017190  -17.1030897



## Q.  3.3.
W <- oxford_w / rowSums(oxford_w) # row standardizes weights

N <- nrow(W)

I <- diag(as.matrix(W)) # 48 zeros (neighbours contig.)

rho <- sac_model$rho

M <- solve(I - rho * W) # calc spatial multiplier


# instcoad effects
direct_effect <- sum(diag(M * coef(sac_model)["instcoad"])) / N
# 1.824513

indirect_effect <- mean(rowSums(M * coef(sac_model)["instcoad"]) - diag(M * coef(sac_model)["instcoad"]))
# -3.637528

total_effect <- mean(rowSums(M * coef(sac_model)["instcoad"]))
# -1.813015



# Neighboring Impact Effects
illinois_row <- which(oxford_data$statenm == "IL")

# Simulate impacts
nSims <- 1000
set.seed(08142024)
coef <- coef(sac_model)
vcov <- vcov(sac_model)
coefs <- rmvnorm(n = nSims, coef, vcov)

neighbor_indices <- which(oxford_w[illinois_row, ] != 0)

neighbor_states <- oxford_data$statenm[neighbor_indices]
# 5 neighbor_states: WI = Wisconsin, IA = Iowa, KY = Kentucky, IN = Indiana, MO = Missouri, Michigan -> Lake -> doesn't count


eff_from_IL <- matrix(
  data = NA,
  nrow = length(neighbor_indices),
  ncol = nSims,
  dimnames = list(neighbor_states, NULL)
)

for (i in 1:nSims) {
  eff_from_IL[, i] <- 
    coefs[i, 6] * solve(I - coefs[i, 1] * W)[neighbor_indices, illinois_row]     # rho is at idx 1, instcoad is at idx 6.
}


eff_summary <- apply(eff_from_IL, 
                     MARGIN = 1, 
                     quantile, 
                     c(0.05, 0.50, 0.95)) %>% t(.)

print(eff_summary)




# Q. 4
# MLE estimation
sar_model_mle <- lagsarlm(ben95 ~ rskpovpc + wage95 + instcoad + ipcfold + teitrend + match, data = oxford_data, listw)
summary(sar_model_mle)

# S2SLS estimation
sar_model_s2sls <- stsls(ben95 ~ rskpovpc + wage95 + instcoad + ipcfold + teitrend + match, data = oxford_data, listw)
summary(sar_model_s2sls)



