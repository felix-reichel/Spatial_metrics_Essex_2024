library(maxLik)


sim_reg_data <- function(n_obs, betas, sgma) {
  X <- matrix(
    rnorm( n_obs * length(betas) ), 
    nrow = n_obs
  )
  eps <- rnorm(n_obs, sd = sgma)
  y <- X %*% betas + eps
  return(list(X = X, y = y))
}

logLikFunc <- function(pars, data) {
  n <- nrow(data)
  K_regressors <- ncol(data) - 1
  b_hats <- pars[ 1:K_regressors ]
  sgma <- pars[ K_regressors + 1]
  
  y <- data[, 1]
  X <- data[, 2:(K_regressors + 1)]
  
  err_resid <- y - X %*% b_hats

  # logLikFunc for SLR is given by:
  logLikFunc <- -n * log(sgma) - 0.5 * (t(err_resid) %*% err_resid) / sgma^2
  return(logLikFunc)
}

# MLE for SLR using the maxLik package
regression_mle <- function(X, y, method = "NR") { # defaults to Newton-Raphson method for maxLik
  data <- cbind(y, X)
  
  startVals <- c(
    lm(y ~ X - 1)$coefficients, 
    1)
  
  result <- maxLik(logLik = logLikFunc, start = startVals, method = method, data = data)
  return(result)
}

compare_OLS_MLE_estimates <- function(n_sim, 
                              n, 
                              betas, 
                              sgma, 
                              methods = c("BFGS", "BFGSR","NM", "SANN", "CG", "NR")) {
  results <- list()
  lm_times <- numeric(n_sim)
  mle_times <- matrix(0, n_sim, length(methods))
  colnames(mle_times) <- methods
  
  for (i in 1:n_sim) {
    
    data <- sim_reg_data(n, betas, sgma)
    X <- data$X
    y <- data$y
    
    lm_start <- Sys.time()
    lm_fit <- lm(y ~ X - 1)
    lm_end <- Sys.time()
    lm_times[i] <- as.numeric(difftime(lm_end, lm_start, units = "secs"))
    
    lm_coef <- coef(lm_fit)
    lm_sigma <- summary(lm_fit)$sigma
    
    for (method in methods) {

      mle_start <- Sys.time()
      mle_fit <- regression_mle(X, y, method = method)
      mle_end <- Sys.time()
      mle_times[i, method] <- as.numeric(difftime(mle_end, mle_start, units = "secs"))

      mle_coef <- coef(mle_fit)[1:ncol(X)]
      mle_sigma <- coef(mle_fit)[ncol(X) + 1]
      
      if (is.null(results[[method]])) {
        results[[method]] <- list(lm_coef = matrix(0, n_sim, length(betas)),
                                  mle_coef = matrix(0, n_sim, length(betas)),
                                  lm_sigma = numeric(n_sim),
                                  mle_sigma = numeric(n_sim))
      }
      
      results[[method]]$lm_coef[i, ] <- lm_coef
      results[[method]]$mle_coef[i, ] <- mle_coef
      results[[method]]$lm_sigma[i] <- lm_sigma
      results[[method]]$mle_sigma[i] <- mle_sigma
    }
  }
  

  lm_time_mean <- mean(lm_times)
  lm_time_sd <- sd(lm_times)
  
  mle_time_means <- colMeans(mle_times)
  mle_time_sds <- apply(mle_times, 2, sd)
  
  return(list(results = results,
       lm_time_mean = lm_time_mean,
       lm_time_sd = lm_time_sd,
       mle_time_means = mle_time_means,
       mle_time_sds = mle_time_sds))
}

# Run Simulations
n_sim <- 10
n <- 10000
betas <- c(5, 2, 4)
sgma <- 1.75

comparison <- compare_OLS_MLE_estimates(n_sim, n, betas, sgma)
print(comparison)
