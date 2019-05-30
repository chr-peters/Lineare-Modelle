# Name: Christian Peters

# No. 25)
# ======

data <- read.table('mietspiegel2015.txt', header = TRUE)
data <- data[data$bj >= 1970 & data$bj <= 1980, ]

# a)

# Model:
# nm = beta0 + beta1 * wfl + e

# get the least squares estimate using the normal equations
X <- cbind(1, data$wfl)
y <- data$nm
beta <- solve(t(X) %*% X, t(X) %*% y)
# beta0 = 138.751308, beta1 = 8.075865

lm(data$nm ~ data$wfl)
# Same results.

# Interpretation:
# Intercept (beta0): No useful way of interpreting the rent for flats without living space.
# wfl (beta1): One square meter of living space costs around 8.08 euro a month.

# b)
residuals <- y - X %*% beta
plot(data$wfl, residuals, xlab = 'wfl', main = 'Residual Plot')

# I don't think that the assumption of equal variance among the errors is justified
# because the deviations from the mean of the residuals clearly increase with
# increasing wfl value.

sse <- sum(residuals**2)
sst <- sum((y - mean(y))**2)
mse <- sse / length(y)
# Mean squared error: 22725.49

r_squared <- 1 - sse / sst
# 0.5984767

# c)

# Estimate the coefficients for the new model using least squares
X <- cbind(1, log(data$wfl))
y <- data$nmqm
beta <- solve(t(X) %*% X, t(X) %*% y)
# beta0 = 23.307780, beta1 = -3.094016

plot(log(data$wfl), data$nmqm, xlab = 'wfl', ylab = 'nmqm', main = 'Square Meter Costs by Living Space')
abline(a = beta[1], b = beta[2])

# d)

# add the new variable and estimate the new coefficients using least squares
X_location <- cbind(X, data$wohngut | data$wohnbest)
beta_location <- solve(t(X_location) %*% X_location, t(X_location) %*% y)
# beta0 = 23.103903, beta1 = -3.172943, beta2 = 1.723611

# e)

# calculate the adjusted coefficients of determination (r_squared)
residuals <- y - X %*% beta
residuals_location <- y - X_location %*% beta_location
sse <- sum(residuals**2)
sse_location <- sum(residuals_location**2)
sst <- sum((y-mean(y))**2)
r_squared <- 1 - sse / sst
r_squared_location <- 1 - sse_location / sst

n <- length(y)
r_squared_adj <- 1 - (n-1) / (n - 2) * (1 - r_squared)
# 0.2142187

r_squared_location_adj <- 1 - (n-1) / (n - 3) * (1 - r_squared_location)
# 0.3173157

# As we can see, the adjusted r_squared of the new model containing the location
# is around 10 percentage points higher than that of the old model.
# So it's safe to say that adding the new variable lead to a significant improvement
# of the model in (c).