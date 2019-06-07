# Name: Christian Peters

library(MASS)

# No. 26)
# ======

cobbdouglas <- read.csv2('cobbdouglas.csv')

# a)

# get the OLS estimate of beta
X <- cbind(1, log(cobbdouglas$xL), log(cobbdouglas$xK))
X_inv <- ginv(X)
y <- cobbdouglas$y

beta <- X_inv %*% y

# estimate the error variance
residuals <- y - X %*% beta
var_error <- sum(residuals**2) / (length(y) - 3)

# estimate the variances of the beta_i
vars_beta <- sapply(1:3, function(i) {
  C <- numeric(3)
  C[i] <- 1
  return(var_error * t(C) %*% X_inv %*% t(X_inv) %*% C)
})

# estimate the 90% confidence bounds
confidence_beta <- t(sapply(1:3, function(i) {
  qnorm(c(0.05, 0.95), mean = beta[i], sd = sqrt(vars_beta[i]))
}))
colnames(confidence_beta) <- c('lower_bound', 'upper_bound')
rownames(confidence_beta) <- c('beta_0', 'beta_1', 'beta_2')
print(confidence_beta)

# d)

# Function to carry out the test. Returns 1 if it rejects the zero
# hypothesis, otherwise 0.
test_beta <- function(z) {
  # create C matrix
  C <- matrix(0, nrow = 2, ncol = 3)
  C[1, 2] <- 1
  C[2, 3] <- 1
  
  # calculate test statistic using the confidence ellipsoid for beta
  test_val <- t(C %*% beta - z) %*% ginv(C %*% X_inv %*% t(X_inv) %*% t(C)) %*% (C %*% beta - z)
  
  # now check if z is inside of the confidence ellipsoid
  rg_cb <- 2 # rank of C_B is 2
  ne <- length(y) - 3 # three coefficients were estimated
  if (test_val <= rg_cb * var_error * qf(0.9, rg_cb, ne)) {
    # no rejection of the zero hypothesis
    return(0)
  }
  # rejection of the zero hypothesis
  return(1)
}

# carry out the tests
(test_beta(c(0.3, 1)))
(test_beta(c(0.35, 0.5)))

# as we can see, in both cases the zero hypothesis is rejected