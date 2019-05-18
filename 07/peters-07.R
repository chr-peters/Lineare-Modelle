# Name: Christian Peters

# No. 19)
# =======

data <- read.csv2('cars.csv')

# a)

# create design matrix and response vector
x <- cbind(1, data$Hubraum, data$Leistung, data$Gewicht, data$Verbrauch)
y <- data$Preis

# get the least squares estimate (gauss-markov-estimate)
beta <- solve(t(x) %*% x, t(x) %*% y)
# -13614.6216356691, 526.935980973853, 145.011255933093, 19141.205906062, 739.892449453102

# Interpretation:
# Intercept: no useful way of interpretation
# beta_1: 526.94 euros per liter of engine displacement
# beta_2: 145.01 euros per horsepower
# beta_3: 19141.21 euros per ton of weight
# beta_4: 739.89 euros per liter of fuel consumption

# b)

# estimate the variance of the error term
residuals <- y - x %*% beta
n <- length(y)
r <- 5 # rank of design matrix
var_estimate <- 1 / (n - r) * sum(residuals**2)

# use this to estimate the covariance matrix
cov_estimate <- var_estimate * solve(t(x) %*% x)
cor_estimate <- diag(1/sqrt(diag(cov_estimate))) %*% cov_estimate %*% diag(1/sqrt(diag(cov_estimate)))
print(cor_estimate)

# c)

# estimate the new beta vector, the residuals and the variance
beta_ridge <- solve(t(x) %*% x + 0.1 * diag(5), t(x) %*% y)
residuals_ridge <- y - x %*% beta_ridge
