library(ggplot2)
library(spdep)
library(spatialreg)
library(Matrix)

set.seed(42)
n <- 100  # Number of neighborhoods

coords <- matrix(runif(n * 2, 0, 10), ncol = 2)

cluster_coords <- rbind(
  cbind(runif(20, 2, 3), runif(20, 2, 3)),
  cbind(runif(20, 6, 7), runif(20, 6, 7)),
  cbind(runif(20, 8, 9), runif(20, 2, 3))
)

coords <- rbind(coords[-(1:60), ], cluster_coords)

# Create spatial weight matrix
nb <- dnearneigh(coords, 0, 2)  # Neighbors within distance of 2 units
W <- nb2listw(nb, style = "W")  # Spatial weights list

# True coefficients
beta <- c(100000, 2, 3)

median_income <- rnorm(n, mean = 40000, sd = 10000)
population_density <- rnorm(n, mean = 3000, sd = 500)

cluster_indices <- (n - 59):n
median_income[cluster_indices] <- rnorm(60, mean = 50000, sd = 12000)
population_density[cluster_indices] <- rnorm(60, mean = 12000, sd = 300)

X <- cbind(1, median_income, population_density)

results <- data.frame(rho = numeric(), model = character(),
                      term = character(), estimate = numeric(), 
                      conf.low = numeric(), conf.high = numeric())

rho_values <- seq(0, 0.9, by = 0.1)

for (rho in rho_values) {
  I <- diag(n)
  W_matrix <- as(as_dgRMatrix_listw(W), "CsparseMatrix")
  e <- rnorm(n, mean = 0, sd = 10000)
  epsilon <- solve(I - rho * W_matrix) %*% e
  
  housing_prices <- as.vector(X %*% beta + epsilon)
  housing_prices[cluster_indices] <- housing_prices[cluster_indices] + 150000
  
  data <- data.frame(housing_prices, median_income, population_density)

  sar_model <- lagsarlm(housing_prices ~ median_income + population_density, 
                        data = data, listw = W, method = "Matrix", tol.solve = 1e-22)
  
  sar_coefs <- coef(summary(sar_model))
  sar_ci <- confint(sar_model)
  
  for (term in c("median_income", "population_density")) {
    results <- rbind(results, data.frame(rho = rho, model = "SAR",
                                         term = term, 
                                         estimate = sar_coefs[term, "Estimate"],
                                         conf.low = sar_ci[term, 1], 
                                         conf.high = sar_ci[term, 2]))
  }
  
  ols_model <- lm(housing_prices ~ median_income + population_density, data = data)
  
  ols_coefs <- coef(summary(ols_model))
  ols_ci <- confint(ols_model)
  
  for (term in c("median_income", "population_density")) {
    results <- rbind(results, data.frame(rho = rho, model = "OLS",
                                         term = term, 
                                         estimate = ols_coefs[term, "Estimate"],
                                         conf.low = ols_ci[term, 1], 
                                         conf.high = ols_ci[term, 2]))
  }
}

ggplot(results, aes(x = factor(rho), y = estimate, color = model)) +
  geom_line(aes(group = interaction(term, model)), size = 1) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1) +
  facet_wrap(~ term, scales = "free_y") +
  labs(title = "Coefficient Estimates and Confidence Intervals for Different Values of rho",
       x = "rho", y = "Coefficient Estimate") +
  theme_minimal() +
  theme(legend.position = "bottom")

cluster_labels <- rep("Other", n)
cluster_labels[cluster_indices] <- rep(c("Cluster 1", "Cluster 2", "Cluster 3"), each = 20)

coords_df <- data.frame(coords, cluster = cluster_labels)
colnames(coords_df) <- c("x", "y", "cluster")

ggplot(coords_df, aes(x = x, y = y, color = cluster)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("Other" = "grey", "Cluster 1" = "red", "Cluster 2" = "blue", "Cluster 3" = "green")) +
  labs(title = "Cluster Map with Neighborhoods",
       x = "X Coordinate", y = "Y Coordinate") +
  theme_minimal() +
  theme(legend.position = "bottom")
