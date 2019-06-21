# Name: Christian Peters

library(MASS)
library(car)

# No. 31)
# =======

cobbdouglas <- read.csv2('cobbdouglas.csv')

# a)

# find the coefficients
X_cobbdouglas <- cbind(1, log(cobbdouglas$xL), log(cobbdouglas$xK))
y_cobbdouglas <- log(cobbdouglas$y)
X_inv_cobbdouglas <- ginv(X_cobbdouglas)
beta_cobbdouglas <- X_inv_cobbdouglas %*% y_cobbdouglas
residuals_cobbdouglas <- y_cobbdouglas - X_cobbdouglas %*% beta_cobbdouglas
ne_cobbdouglas <- length(residuals_cobbdouglas) - 3
var_cobbdouglas <- sum(residuals_cobbdouglas**2) / ne_cobbdouglas

# find the 90% bonferroni condfidence intervals
C <- list(matrix(c(0, 1, 0), nrow = 1), matrix(c(0, 0, 1), nrow = 1))
sd_bonferroni <- sapply(1:2, function(i) {
  sqrt(var_cobbdouglas * C[[i]] %*% X_inv_cobbdouglas %*%
         t(X_inv_cobbdouglas) %*% t(C[[i]]))
})
intervals_bonferroni <- sapply(1:2, function(i) {
  lower <- C[[i]] %*% beta_cobbdouglas - sd_bonferroni[i] * qt(1-0.1/4, ne_cobbdouglas)
  upper <- C[[i]] %*% beta_cobbdouglas + sd_bonferroni[i] * qt(1-0.1/4, ne_cobbdouglas)
  return(c(lower = lower, upper = upper))
})
intervals_bonferroni <- t(intervals_bonferroni)
rownames(intervals_bonferroni) <- c('beta1', 'beta2')
print('Bonferroni Intervals:')
print(intervals_bonferroni)

# find the 90% scheffe confidence intervals, Rg(C) = 2
intervals_scheffe <- sapply(1:2, function(i) {
  lower <- C[[i]] %*% beta_cobbdouglas - sqrt(var_cobbdouglas * C[[i]] %*% X_inv_cobbdouglas %*%
                                                t(X_inv_cobbdouglas) %*% t(C[[i]]) * 2 * 
                                                qf(0.9, 2, ne_cobbdouglas))
  upper <- C[[i]] %*% beta_cobbdouglas + sqrt(var_cobbdouglas * C[[i]] %*% X_inv_cobbdouglas %*%
                                                t(X_inv_cobbdouglas) %*% t(C[[i]]) * 2 * 
                                                qf(0.9, 2, ne_cobbdouglas))
  return(c(lower = lower, upper = upper))
})
intervals_scheffe <- t(intervals_scheffe)
rownames(intervals_scheffe) <- c('beta1', 'beta2')
cat('\n')
print('Scheffe Intervals:')
print(intervals_scheffe)

# b)

# draw the confidence ellipse using lm
model_cobbdouglas <- lm(log(y) ~ log(xL) + log(xK), data = cobbdouglas)
confidenceEllipse(model_cobbdouglas, xlab = 'beta1', ylab = 'beta2', grid = FALSE)

# draw the bonferroni intervals (red)
abline(v = intervals_bonferroni['beta1', 'lower'], col = 'red', lwd = 2)
abline(v = intervals_bonferroni['beta1', 'upper'], col = 'red', lwd = 2)
abline(h = intervals_bonferroni['beta2', 'lower'], col = 'red', lwd = 2)
abline(h = intervals_bonferroni['beta2', 'upper'], col = 'red', lwd = 2)

# draw the scheffe intervals (green)
abline(v = intervals_scheffe['beta1', 'lower'], col = 'green', lwd = 2)
abline(v = intervals_scheffe['beta1', 'upper'], col = 'green', lwd = 2)
abline(h = intervals_scheffe['beta2', 'lower'], col = 'green', lwd = 2)
abline(h = intervals_scheffe['beta2', 'upper'], col = 'green', lwd = 2)

# draw the single confidence intervals (yellow)
abline(v = 0.5576, col = 'yellow', lwd = 2)
abline(v = 1.0569, col = 'yellow', lwd = 2)
abline(h = 0.1237, col = 'yellow', lwd = 2)
abline(h = 0.3424, col = 'yellow', lwd = 2)

# c)

# Because we are dealing with a multiple test problem, we have to look at the
# confidence intervals which take this circumstance into account, namely
# Bonferroni and Scheffe.
# The first hypothesis (beta1 = 0.55) lies whithin the bounds of both confidence intervals
# which means that we can't reject it at the 10% niveau.