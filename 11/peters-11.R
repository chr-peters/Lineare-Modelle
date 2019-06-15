# Name: Christian Peters

library(MASS)

set.seed(123)

# No. 30)
# =======

n <- 20000

X <- cbind(1, 1:50)
X_inv <- ginv(X)
beta_real <- c(0, 1)
alpha = 0.05

# simulate for each rho (p)
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

# a)
# As we can see by looking at the results, the OLS estimator for beta is still
# unbiased, because the mean values for beta_0 and beta_1 are still 0 and 1.

# b)
# The rejection percentages of H0 can be seen in the console output.
# Only when rho = 0, which means that the errors are independent, the test
# seems to reject H0 roughly 5% of the time. As soon as rho gets bigger, 
# the amount of rejection grows as well. When rho = 0.8, H0 is rejected about 50% of the time.
# Contrary, when rho is -0.8, H0 is never rejected.