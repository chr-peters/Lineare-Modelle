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
intervals_bonferroni <- sapply(1:2, function(i) {
  sd_est <- sqrt(var_cobbdouglas * C[[i]] %*% X_inv_cobbdouglas %*%
                   t(X_inv_cobbdouglas) %*% t(C[[i]]))
  lower <- C[[i]] %*% beta_cobbdouglas - sd_est * qt(1-0.1/4, ne_cobbdouglas)
  upper <- C[[i]] %*% beta_cobbdouglas + sd_est * qt(1-0.1/4, ne_cobbdouglas)
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
confidenceEllipse(model_cobbdouglas, xlab = 'beta1', ylab = 'beta2', main = 'Confidence Ellipse and Intervals for a = 10%', grid = FALSE, levels = 0.9)

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

# As we can see, the Scheffe confidence intervals completely contain the
# confidence ellipse. The intervals by Bonferroni are a bit smaller and leave
# some parts of the ellipse out. The single confidence intervals (as drawn in yellow) which
# don't factor in the multiple testing scenario are clearly the smallest.

# c)

# draw the new point given by the hypotheses
points(0.55, 0.1, pch = 16)

# Because we are dealing with a multiple test problem, we have to look at the
# confidence intervals which take this circumstance into account, namely
# Bonferroni and Scheffe (the red and the green one, not the yellow one).
# As we can see in the plot, Bonferroni would reject H_0 for beta2 (very closely)
# and not for beta1. Scheffe doesn't reject H_0 for both coefficients.

# No. 32)
# =======

# a)

# calculate beta manually
y <- c(29, 32, 19, 20, 27, 24)
X <- cbind(1, c(1, 1, 0, 0, 0, 0), c(0, 0, 1, 1, 0, 0), c(0, 0, 0, 0, 1, 1))
B <- matrix(c(0, 1, 1, 1), nrow = 1)
X_B <- X %*% (diag(4) - ginv(B) %*% B)
X_B_inv <- ginv(X_B)
beta <- X_B_inv %*% y

cat('\n')
print('Model coefficients:')
print(beta)

# now use lm
data <- data.frame(y = y, a = as.factor(c(1, 1, 2, 2, 3, 3)))
model <- lm(y ~ a, data = data, contrasts = list(a = 'contr.sum'))

cat('\n')
print('lm coefficients')
print(model$coefficients)

# We can see that for whatever reason, lm doesn't use a3.

# b)
confidenceEllipse(model, main = 'Confidence Ellipse for a1 and a2 coefficients', levels = 0.95)

# No. 33)
# =======

cuckoo <- read.table('kuckuck.txt', header = TRUE)

# a)

# Fit a balanced ANOVA model:
# y = mu + a_i + eij
# y := length of egg
# mu := mean of egg length
# a_i := influence of i'th species on egg length
# eij := Error
# Side condition: sum(n_i * a_i) = 0

# create the design matrix and the response vector
n_WP <- sum(!is.na(cuckoo$WP))
n_BP <- sum(!is.na(cuckoo$BP))
n_RK <- sum(!is.na(cuckoo$RK))
n_ZK <- sum(!is.na(cuckoo$ZK))
n <- n_WP + n_BP + n_RK + n_ZK
X_cuckoo <- cbind(1, c(rep(1, n_WP), numeric(n - n_WP)), c(numeric(n_WP), rep(1, n_BP), numeric(n_RK + n_ZK)),
           c(numeric(n_WP + n_BP), rep(1, n_RK), numeric(n_ZK)), c(numeric(n - n_ZK), rep(1, n_ZK)))
colnames(X_cuckoo) <- c('mu', 'a_WP', 'a_BP', 'a_RK', 'a_ZK')           
y_cuckoo <- c(cuckoo$WP[!is.na(cuckoo$WP)], cuckoo$BP[!is.na(cuckoo$BP)], cuckoo$RK[!is.na(cuckoo$RK)], cuckoo$ZK[!is.na(cuckoo$ZK)])

# Find the coefficients for a balanced model
B <- matrix(c(0, n_WP, n_BP, n_RK, n_ZK), nrow = 1)
X_B_cuckoo <- X_cuckoo %*% (diag(5) - ginv(B) %*% B)
X_B_inv_cuckoo <- ginv(X_B_cuckoo)
beta_cuckoo <- X_B_inv_cuckoo %*% y_cuckoo
rownames(beta_cuckoo) <- c('mu', 'a_WP', 'a_BP', 'a_RK', 'a_ZK')
cat('\n')
print('Model coefficients:')
print(beta_cuckoo)

# b)

# Use an F-Test to test if the effects of bird species are different.
# Precondition: The parametric function that is used is LE-estimable (which is the case here).

# Hypothesis: C * beta = 0
C_cuckoo <- matrix(c(0, 1, -1, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1), nrow = 3, byrow = TRUE)

# compute the test statistic
residuals_cuckoo <- y_cuckoo - X_cuckoo %*% beta_cuckoo
ne_cuckoo <- length(residuals_cuckoo) - 5
var_cuckoo <- sum(residuals_cuckoo**2) / ne_cuckoo
# rank(CB) = 3
teststat <- t(C_cuckoo %*% beta_cuckoo) %*% ginv(C_cuckoo %*% X_B_inv_cuckoo %*% t(X_B_inv_cuckoo) %*% t(C_cuckoo)) %*%
  (C_cuckoo %*% beta_cuckoo) / (3 * var_cuckoo)
# 13.32955

# compute the critical value
critical <- qf(0.9, 3, ne_cuckoo)
# 2.14859

# We can see that the test statistic is greater than the critical value, which means
# that there is significant evidence that the lengths of the cuckoo eggs differ depending
# on the species of bird which breeds them.
# We reject H0 at the 10% niveau.

# c)

# iterate over all possible combinations of bird species and perform a bonferroni
# test for each
invisible(sapply(combn(c('a_WP', 'a_BP', 'a_RK', 'a_ZK'), 2, simplify = FALSE), function(pair) {
  cat('\nTest ', pair[1], ' vs. ', pair[2], ':\n')
  
  # create hypothesis matrix
  pos_1 <- match(pair[1], rownames(beta_cuckoo))
  pos_2 <- match(pair[2], rownames(beta_cuckoo))
  C <- matrix(0, ncol = 5)
  C[1, pos_1] = 1
  C[1, pos_2] = -1
  
  # calculate test statistic
  teststat <- abs(C %*% beta_cuckoo)

  # estimate the standard deviation
  sd_est <- sqrt(var_cuckoo * C %*% X_B_inv_cuckoo %*% t(X_B_inv_cuckoo) %*% t(C))

  # estimate the critical value:
  critical <- qt(1 - 0.1 / (2 * 6), ne_cuckoo) * sd_est
  
  cat('Value of teststat: ', teststat, ', Critical value: ', critical)
  cat('\nResult: ', ifelse(teststat > critical, 'Rejection (H1)', 'No rejection (H0)'), '\n')
}))
