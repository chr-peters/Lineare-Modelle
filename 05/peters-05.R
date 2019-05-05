# Name: Christian Peters

# No. 12)
# =======

data <- read.csv('campus.csv')

# a)

# Model (Analysis of variance):
# y = mu + a_i + b_j + e_ij; i = 1, 2; j = 1, 2
#
# y    := Mean time per km in seconds
# a_i  := Fixed effect of distance on the mean (1=2.5km, 2=10km)
# b_j  := Fixed effect of gender on the mean (1=male, 2=female)
# e_ij := Random error term

# Side condition to account for differential effects:
# 79 * a_1 + 231 * a_2 = 0
# 226 * b_1 + 84 * b_2 = 0
# where 79 is the number of data points with distance = 2500,
# 231 is the number of data points with distance = 10000,
# 226 is the number of data points with gender = m and
# 84 is the number of data points with gender = w

# create the design matrix
x <- cbind(1, ifelse(data$Distanz==2500, 1, 0), ifelse(data$Distanz==10000, 1, 0),
           ifelse(data$Geschlecht=="m", 1, 0), ifelse(data$Geschlecht=="w", 1, 0))
# side condition
x <- rbind(x, c(0, 79, 231, 0, 0), c(0, 0, 0, 226, 84))

# create target vector with side condition
y <- c(data$Sekunden, 0, 0)

# estimate parameters using least squares:
estimates <- solve(t(x) %*% x, t(x) %*% y)
# mu = 329.634194
# a_1 = 32.177786
# a_2 = -11.004524
# b_1 = -9.618911
# b_2 = 25.879451

# Interpretation:
# We can see that the overall mean of running times is roughly 329.6 seconds
# per km. Interestingly, people who run the 2.5km distance seem to be about
# 32 seconds slower than the average and people who run the longer distance
# are 11 seconds faster. This could be the case because people who are better
# in shape tend to seek the challenge of the longer run while running novices
# prefer the shorter distance. But that's speculation of course.
# Regarding the gender influence, we can see that men are about 9.6 seconds
# faster than the average and women are 25.9 seconds slower.
# This means that on average, men are roughly 35 seconds faster than women.

# Prediction for Stacey:
target_stacey <- c(1, 1, 0, 0, 1) %*% estimates * 2.5
# 969.2286

# According to the model, Stacey would have to finish the 2.5km within
# about 969 seconds.

# b)

# get the sorted times of the 30 fastest men and women
times_men <- sort(data$Sekunden[data$Geschlecht=='m'])[1:30]
times_women <- sort(data$Sekunden[data$Geschlecht=='w'])[1:30]

# Model:
# time_woman = intercept + beta * time_man + error

# create the design matrix
x <- cbind(1, times_men)

# get parameter estimates
params <- solve(t(x) %*% x, t(x) %*% times_women)
# intercept = -18.417027
# beta = 1.326073

# Here we can see that women are on average 1.33 times slower than men.
# The interesting part about this model is, that the intercept value is negative.
# I was expecting this value to be positive which would support the thesis that
# women are slower on average.

# c)

# define new loss function
loss_function <- function(beta) {
  sum(abs(times_women - beta[1] - beta[2] * times_men))
}

# get parameter vector
beta <- optim(c(0, 0), loss_function)$par
# -25.091352   1.354165

# Here we can see that the '1.3 rule' still holds true, but the intercept
# has become even more negative.

# Let's compare the models by creating a plot:
plot(times_men, times_women, main = 'Comparison of running times between men and women')
curve(params[1] + params[2] * x, add = TRUE, col = 'blue')
curve(beta[1] + beta[2] * x, add = TRUE, col = 'red')
legend('topleft', inset = 0.05, legend = c('OLS', 'LAD'), col = c('blue', 'red'), lty = c(1, 1))

# In the plot we can see that the models are almost the same.

# d)

# manipulate data
fake_data <- data
fake_data[1, 'Sekunden'] <- 100

times_men_fake <- sort(fake_data$Sekunden[fake_data$Geschlecht=='m'])[1:30]
times_women_fake <- sort(fake_data$Sekunden[fake_data$Geschlecht=='w'])[1:30]

# OLS
x <- cbind(1, times_men_fake)
beta_ols <- solve(t(x) %*% x, t(x) %*% times_women_fake)
# -171.888769  1.957291

# LAD
loss_function_fake <- function(beta) {
  sum(abs(times_women_fake - beta[1] - beta[2] * times_men_fake))
}
beta_lad <- optim(c(0, 0), loss_function_fake)$par
# -49.644752   1.453525

# create plot
plot(times_men_fake, times_women_fake, main = 'Running times with measurement error')
curve(beta_ols[1] + beta_ols[2] * x, add = TRUE, col = 'blue')
curve(beta_lad[1] + beta_lad[2] * x, add = TRUE, col = 'red')
legend('topleft', inset = 0.05, legend = c('OLS', 'LAD'), col = c('blue', 'red'), lty = c(1, 1))

# As we can see in the plot, the LAD estimates are more robust than the OLS estimates regarding outliers.