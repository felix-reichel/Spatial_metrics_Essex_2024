# practical ----
# 1. import the shapefile for Westminster Parliamentary Constituencies, merge it with the Brexit data, and visualize the Leave votes;
# 2. create a contiguity W based on the imported data with useful row and column names;
# 3. check the distribution of the number of neighbors and figure out which constituencies(s) are contiguous with Colchester;
# 4. check the invertibility of W by two means without using solve();
# 5. have a new W without "islands," row-standardized, and check its symmetry;
# 6. simulate a rho that ~ uniform(0.05, 0.95) and is rounded to 3 decimal points;
# 7. have a spatial multiplier based on rho and the new W and calculate its trace.

library(sf)
library(spdep)
library(tidyverse)
library(Matrix)
library(rstudioapi)
library(viridis)

setwd(dirname(getActiveDocumentContext()$path))

# Load spatial and Brexit data
load("data/sf_GBR.RData")  # Boundary data
brexit_data <- read_csv("data/Brexit.csv")
spatial_data <- left_join(sf_GBR, brexit_data, by = c("pcon15cd" = "PCON11CD"))

# 1. Visualize Leave votes
ggplot(spatial_data) + 
  geom_sf(aes(fill = leave), color = "black", size = 0.01) +
  scale_fill_viridis_c(option = "C", name = "Leave votes in Brexit") +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", size = 10),
    legend.position = "right"
  ) +
  coord_sf()

# 2. Create spatial weight matrix (W)
spatial_weights <- nb2mat(poly2nb(spatial_data), style = "B", zero.policy = TRUE)
dimnames(spatial_weights) <- list(spatial_data$pcon15nm, spatial_data$pcon15nm)

# 3. Check distribution of number of neighbors and find Colchester's neighbors
hist(rowSums(spatial_weights), xlab = "Number of neighbours", col = "whitesmoke")
spatial_weights["Colchester", spatial_weights["Colchester", ] != 0]

# 4. Check invertibility of spatial weight matrix without using solve()
det(spatial_weights)
rankMatrix(spatial_weights)

# 5. Remove islands, row-standardize spatial weight matrix, and check its symmetry
spatial_weights_no_islands <- spatial_weights[rowSums(spatial_weights) != 0, colSums(spatial_weights) != 0]
spatial_weights_standardized <- spatial_weights_no_islands / rowSums(spatial_weights_no_islands)
isSymmetric(spatial_weights_standardized)

# 6. Simulate rho from uniform(0.05, 0.95) rounded to 3 decimal points
set.seed(202307)
rho <- round(runif(1, 0.05, 0.95), 3)

# 7. Calculate spatial multiplier matrix and its trace
spatial_multiplier <- solve(diag(nrow(spatial_weights_standardized)) - rho * spatial_weights_standardized)
sum(diag(spatial_multiplier))

