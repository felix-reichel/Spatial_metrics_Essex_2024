# Load required libraries
# install.packages("spdep")
# install.packages("spatialreg")
library(spdep)
library(spatialreg)

# Parameters
b <- c(5, 10)  # True coefficient vector (beta)
rho <- 0.7     # Spatial autocorrelation parameter (rho), 0.7 for indicating strong correlation
n <- 900       # Number of observations

# Create an artificial DGP
# Generate artificial covariate matrix X with an intercept and a uniformly distributed covariate
X <- matrix(c(rep(1, n), runif(n, -10, 10)), n, 2)

# Create an artificial spatial weight matrix W for a rectangular grid
W <- cell2nb(sqrt(n), sqrt(n))

# Invert the matrix I - rho * W
IW <- invIrM(W, rho)

# Number of bootstrap simulations
r <- 420

# Matrix to store bootstrap estimates of beta
b_hat <- matrix(NA, r, 2)

# Set seed for reproducibility
set.seed(420)

# Bootstrap to estimate beta coefficients
for (i in 1:r) {
  e <- rnorm(n, mean = 0, sd = 1)  # Generate random error terms from N(0, 1)
  y_sim <- IW %*% X %*% b + IW %*% e  # Simulate dependent variable with spatial autocorrelation
  ols <- lm(y_sim ~ X[, 2] - 1)  # Ordinary Least Squares regression (no intercept)
  b_hat[i, ] <- coef(ols)  # Store the estimated coefficients
}

# Plot the sampling distribution of the second coefficient (beta_2)
plot(density(b_hat[, 2]), xlab = expression(widehat(beta[2])), main = "Sampling Distribution of beta_2")

# Summary statistics for bootstrap estimates
summary(b_hat)

# Note: The confidence intervals of the bootstrap estimates do not include the true beta
#       indicating that positive spatial correlation biases the estimates of beta.

