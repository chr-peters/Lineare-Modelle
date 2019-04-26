# Name: Christian Peters

library(mvtnorm)

set.seed(1234)

# No. 9)
# ======

# b)

# define density functions
f_u <- function(u1, u2) {
  dmvnorm(x=c(u1, u2), mean = c(0, 0), sigma = diag(2))
}
f_v <- function(v1, v2) {
  dmvnorm(x=c(v1, v2), mean = c(1, 1), sigma = matrix(c(2, 1, 1, 1), nrow = 2))
}
f_y <- function(y1, y2) {
  dmvnorm(x=c(y1, y2), mean = c(0, 1), sigma = matrix(c(2, -1, -1, 1), nrow = 2))
}

# create contour plots
x <- seq(-1, 1, 0.1)
y <- seq(-1, 1, 0.1)
z <- sapply(y, function(u2) sapply(x, function (u1) f_u(u1, u2)))
contour(x=x, y=y, z=z, xlab = 'u1', ylab = 'u2', main = 'Density function of u')

x <- seq(-2, 4, 0.1)
y <- seq(-1.5, 3.5, 0.1)
z <- sapply(y, function(v2) sapply(x, function (v1) f_v(v1, v2)))
contour(x=x, y=y, z=z, xlab = 'v1', ylab = 'v2', main = 'Density function of v')

x <- seq(-3, 3, 0.1)
y <- seq(-1.5, 3.5, 0.1)
z <- sapply(y, function(y2) sapply(x, function (y1) f_y(y1, y2)))
contour(x=x, y=y, z=z, xlab = 'y1', ylab = 'y2', main = 'Density function of y')

# correlation matrices
corr_u <- matrix(c(1, 0, 0, 1), nrow = 2)
print(corr_u)
corr_v <- matrix(c(1, 1/sqrt(2), 1/sqrt(2), 1), nrow=2)
print(corr_v)
corr_y <- matrix(c(1, -1/sqrt(2), -1/sqrt(2), 1), nrow=2)
print(corr_y)

# Description:
# u is centered around the origin with equal variance in each direction. The
# components u1 and u2 are uncorrelated.
# v is centered around (1, 1) with a variance of 2 in the direction of v1 and a
# variance of 1 in the direction of v2. v1 and v2 have a correlation coefficient
# of ~0.71.
# The center of y is located at (0, 1) and it has a variance of 1 in the direction
# of y1 and a variance of 1 in the direction of y2. The components y1 and y2
# are negatively correlated with a correlation coefficient of ~-0.71.

# c)

# define the functions to generate random numbers using the cholesky decomposition
# of the covariance matrix (2.35)
generate_u <- function() {
  diag(2) %*% rnorm(2) + c(0, 0)
}
generate_v <- function() {
  matrix(c(1, 0, 1, 1), nrow = 2)  %*% rnorm(2) + c(1, 1)
}
generate_y <- function() {
  matrix(c(1, 0, -1, 1), nrow = 2) %*% rnorm(2) + c(0, 1)
}

# generate plots
samples_u <- replicate(10000, generate_u(), simplify = TRUE)
plot(samples_u[1,], samples_u[2,], xlab = 'u1', ylab = 'u2', main = '10000 samples of u')

samples_v <- replicate(10000, generate_v(), simplify = TRUE)
plot(samples_v[1,], samples_v[2,], xlab = 'v1', ylab = 'v2', main = '10000 samples of v')

samples_y <- replicate(10000, generate_y(), simplify = TRUE)
plot(samples_y[1,], samples_y[2,], xlab = 'y1', ylab = 'y2', main = '10000 samples of y')

# No. 11)
# =======

# a)

generate_chisq <- function(df) {
  sum(replicate(df, rnorm(1)) ** 2)
}

samples_chisq_2 <- replicate(1000, generate_chisq(2))
plot(ecdf(samples_chisq_2), main = 'Distribution of 1000 chi squared samples with df=2')
curve(pnorm(x), col='red', add = TRUE)
legend('bottomright', inset=0.1, legend=c('chi squared', 'standard normal'),
       col = c('black', 'red'), lty=c(1, 1))

samples_chisq_4 <- replicate(1000, generate_chisq(4))
plot(ecdf(samples_chisq_4), main = 'Distribution of 1000 chi squared samples with df=4')
curve(pnorm(x), col='red', add = TRUE)
legend('bottomright', inset=0.1, legend=c('chi squared', 'standard normal'),
       col = c('black', 'red'), lty=c(1, 1))

samples_chisq_6 <- replicate(1000, generate_chisq(6))
plot(ecdf(samples_chisq_6), main = 'Distribution of 1000 chi squared samples with df=6')
curve(pnorm(x), col='red', add = TRUE)
legend('bottomright', inset=0.1, legend=c('chi squared', 'standard normal'),
       col = c('black', 'red'), lty=c(1, 1))

samples_chisq_8 <- replicate(1000, generate_chisq(8))
plot(ecdf(samples_chisq_8), main = 'Distribution of 1000 chi squared samples with df=8')
curve(pnorm(x), col='red', add = TRUE)
legend('bottomright', inset=0.1, legend=c('chi squared', 'standard normal'),
       col = c('black', 'red'), lty=c(1, 1))

# Comparison:
# As we can see, the chi squared distribution doesn't converge to a standard
# normal distribution with increasing degrees of freedom.

# b)

generate_t <- function(df) {
  rnorm(1) / sqrt(generate_chisq(df) / df)
}

samples_t_2 <- replicate(1000, generate_t(2))
plot(ecdf(samples_t_2), main = 'Distribution of 1000 t samples with df=2')
curve(pnorm(x), col='red', add = TRUE)
legend('topleft', inset=0.05, legend=c('t', 'standard normal'),
       col = c('black', 'red'), lty=c(1, 1))

samples_t_4 <- replicate(1000, generate_t(4))
plot(ecdf(samples_t_4), main = 'Distribution of 1000 t samples with df=4')
curve(pnorm(x), col='red', add = TRUE)
legend('topleft', inset=0.05, legend=c('t', 'standard normal'),
       col = c('black', 'red'), lty=c(1, 1))

samples_t_6 <- replicate(1000, generate_t(6))
plot(ecdf(samples_t_6), main = 'Distribution of 1000 t samples with df=6')
curve(pnorm(x), col='red', add = TRUE)
legend('topleft', inset=0.05, legend=c('t', 'standard normal'),
       col = c('black', 'red'), lty=c(1, 1))

samples_t_8 <- replicate(1000, generate_t(8))
plot(ecdf(samples_t_8), main = 'Distribution of 1000 t samples with df=8')
curve(pnorm(x), col='red', add = TRUE)
legend('topleft', inset=0.05, legend=c('t', 'standard normal'),
       col = c('black', 'red'), lty=c(1, 1))

# Comparison:
# In this case we can see that the t distribution converges to a standard normal
# distribution with increasing degrees of freedom.
