##############
# TUTORIAL 1 #
##############

# 1. draw two random variables X1 and X2 from a multivariate normal distribution (seed = 84956) that:
#### N = 10,000
#### mean X1 = -2,
#### mean X2 = 2,
#### variance X1 = variance X2 = 1,
#### covariance between X1 and X2 = 0.3;

# 2. draw random error (epsilon) from a normal distribution with mean = 0 and variance = 1, with the same seed value;

# 3. generate y = -5 * X1 + 4 * X2 + epsilon;

# 4. apply OLS to estimate that model (including intercept) with matrix algebra and 
#    present both point and uncertainty estimate of regression coefficients.


set.seed(84956)
N <- 10000
mean <- c(-2, 2)
cov_matrix <- matrix(c(1, 0.3, 0.3, 1), ncol = 2)

library(MASS)
X <- mvrnorm(N, mu = mean, Sigma = cov_matrix)
X1 <- X[, 1]
X2 <- X[, 2]
eps <- rnorm(N, mean = 0, sd = 1)
y <- -5 * X1 + 4 * X2 + eps
model <- lm(y ~ X1 + X2)
summary(model)




# essex summer school: spatial econometrics
# solution 1
# Muzhou Zhang & Martin Steinwand
# 22 July 2024

library(MASS)

N <- 10^4
mu <- c(-2, 2)
Sigma <- matrix(c(1, 0.3, 0.3, 1), 2, 2, byrow = FALSE)

set.seed(84956)
X <- mvrnorm(N, mu, Sigma)
epsilon <- rnorm(N)

beta <- matrix(c(-5, 4), nrow = 2)
y <- X %*% beta + epsilon

X <- cbind(1, X) # include intercept

beta_hat <- solve(t(X) %*% X) %*% t(X) %*% y

e <- y - X %*% beta_hat

beta_vcov <- c(t(e) %*% e / (N - ncol(X))) * (solve(t(X) %*% X))

est <- cbind(beta_hat, sqrt(diag(beta_vcov)))

dimnames(est) <- list(c("Intercept", "X1", "X2"), c("beta", "se"))

est

summary(lm(y ~ -1 + X))

rm(list = ls())

