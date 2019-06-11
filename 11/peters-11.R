# Name: Christian Peters

library(MASS)

# No. 30)
# =======

n <- 20000

X <- cbind(1, 1:50)
X_inv <- ginv(X)
beta_real <- c(0, 1)
alpha = 0.1

results <- sapply(seq(-0.8, 0.8, 0.2), function(p) {
  results_p <- replicate(n, {
    
    # simulate the error vector
    e <- numeric(50)
    for (i in 2:length(e)) {
      e[i] <- p * e[i-1] + rnorm(1)
    }
    
    # generate y
    y <- X %*% beta_real + e
    
    # estimate coefficients using OLS
    beta_estimate <- X_inv %*% y
    
    # do the t test
    c <- c(0, 1)
    residuals <- y - X %*% beta_estimate
    var_estimate <- sum(residuals**2) / (50 - 2)
    test_stat <- (beta_estimate[2] - 1) / sqrt(var_estimate * t(c) %*% X_inv %*% t(X_inv) %*% c)
    test_res <- 0
    if (abs(test_stat) > qt(1 - alpha/2, df = 48)) {
      test_res <- 1
    }
    
    return(c(test_res = test_res, beta_0 = beta_estimate[1], beta_1 = beta_estimate[2]))
  })
  
  # calculate the simulation statistics
  return(c(rejection_percent = mean(results_p['test_res', ]), 
           mean_b0 = mean(results_p['beta_0', ]), 
           mean_b1 = mean(results_p['beta_1', ])))
})

colnames(results) <- seq(-0.8, 0.8, 0.2)
print(results)