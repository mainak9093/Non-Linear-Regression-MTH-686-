
# Load data
data <- read.table("set-49.dat")
head(data)

# data
t <- as.numeric(unlist(data[1])) # replace with your time values
t <- t/100
y1 <- as.numeric(unlist(data[2])) # replace with your observed y(t) values

# Plot the data
plot(data$V1, data$V2, type = "l", main = "Time Series Data", xlab = "Time", ylab = "Value")

# Residual Plot after fitting a simple model
fit <- lm(V2 ~ V1, data = data) # Example linear model
residuals <- resid(fit)
plot(data$V1, residuals, main = "Residual Plot", xlab = "Time", ylab = "Residuals")
abline(h = 0, col = "red")

# Calculate mean and standard deviation to estimate SNR
signal <- mean(y1, na.rm = TRUE)
noise <- sd(y1, na.rm = TRUE)
SNR <- signal / noise
cat("Signal-to-Noise Ratio:", SNR)

# ACF plot of residuals to check for white noise
acf(residuals, main = "ACF of Residuals")


model_function <- function(t, alpha0, alpha1, alpha2, beta1, beta2) {
  a <- alpha0 + alpha1 * exp(beta1 * t) + alpha2 * exp(beta2 * t)
  return(a)
}

###############################################


# Create delayed values for t = 1 to (length(y) - 2)
Y1 <- y1[1:(length(y1) - 2)]
Y2 <- y1[2:(length(y1) - 1)]
Y3 <- y1[3:length(y1)]

# Formulate the linear system
A <- cbind(Y2, Y1)  # Matrix of delayed terms
B <- Y3             # Response vector

# Solve for coefficients (c1, c2) in the equation Y3 ??? c1 * Y2 + c2 * Y1
coef <- lm(B ~ A - 1)$coefficients  # -1 removes the intercept
c1 <- coef[1]
c2 <- coef[2]

# Find roots of the characteristic polynomial
roots <- polyroot(c(1, -c1, -c2))

beta_estimates <- log(roots)
beta1 <- Re(beta_estimates[1])  # Take real part
beta2 <- Re(beta_estimates[2])  # Take real part

# Construct the matrix with estimated beta values
X <- cbind(1, exp(beta1 * t), exp(beta2 * t))

# Solve for alpha values
alpha_estimates <- lm(y1 ~ X - 1)$coefficients
alpha0 <- alpha_estimates[1]
alpha1 <- alpha_estimates[2]
alpha2 <- alpha_estimates[3]

cat("Estimated parameters:\n")
cat("alpha0:", alpha0, "\n")
cat("alpha1:", alpha1, "\n")
cat("alpha2:", alpha2, "\n")
cat("beta1:", beta1, "\n")
cat("beta2:", beta2, "\n")



y2 <- model_function(t,alpha0,alpha1,alpha2,beta1,beta2)

plot(t, y1, type = "l", col = "blue", main = "Observed vs Fitted Data", xlab = "Time", ylab = "Values")
lines(t, y2, col = "red", lwd = 2)
legend("topright", legend = c("Observed Data", "Prony's Model Fit"), col = c("blue", "red"), lty = 1, lwd = 2)

rmse <- sqrt(mean((y1 - y2)^2, na.rm = TRUE))
rss <- sum(y1 - y2)^2                  # Residual sum of squares
tss <- sum((y1 - mean(y1))^2)                         # Total sum of squares
r_squared <- 1 - (rss / tss)

# Display model accuracy metrics
cat("Root Mean Squared Error (RMSE):", rmse, "\n")
cat("R-squared:", r_squared, "\n")



#########################################################################
# Step 1: Estimate the intercept (mean of the data as a starting guess)
intercept <- mean(data$V2)

# Step 2: Adjust data by subtracting the intercept
adjusted_data <- data$V2-intercept 

# Create the Toeplitz matrix
n <- length(adjusted_data)
p <- 2  # Number of exponential terms
T_matrix <- matrix(0, nrow = n - p, ncol = p + 1)

for (i in 1:(n - p)) {
  T_matrix[i, ] <- c(1, adjusted_data[i:(i + p - 1)])
}

# Response vector
y <- adjusted_data[(p + 1):n]

# Estimate parameters using least squares
params <- solve(t(T_matrix) %*% T_matrix, t(T_matrix) %*% y)

# Ensure params length corresponds to the number of parameters you expect
if (length(params) < 3) {
  stop("Not enough parameters estimated.")
}



# Define the characteristic polynomial based on your coefficients
# Ensure you are correctly referencing params based on the model structure
char_poly <- c(1, -c(params[2], params[3]))  # Make sure these correspond to beta coefficients
beta_values <- log(polyroot(char_poly))

# Output beta values
cat("Beta Values:\n")
cat("Beta1:", Re(beta_values[1]), "\n")
cat("Beta2:", Re(beta_values[2]), "\n")

# Construct the matrix with estimated beta values
X <- cbind(1, exp(Re(beta_values[1]) * t), Re(beta_values[1]) * t)

# Solve for alpha values
alpha_estimates <- lm(y1 ~ X - 1)$coefficients
alpha0 <- alpha_estimates[1]
alpha1 <- alpha_estimates[2]
alpha2 <- alpha_estimates[3]

cat("Estimated Parameters:\n")
cat("Alpha0:", alpha0, "\n")
cat("Alpha1:", alpha1, "\n")
cat("Alpha2:", alpha2, "\n")

# Reconstruct the model using the estimated parameters
reconstructed_y <- alpha0 + alpha1 * exp(Re(beta_values[1]) * (1:length(data$V2))) + alpha2 * exp(Re(beta_values[2]) * (1:length(data$V2)))

# Plot original vs. reconstructed data
plot(data$V2, type = 'l', col = 'blue', lty = 1, main = 'Original vs Reconstructed Data', ylab = 'Value', xlab = 'Time')
lines(reconstructed_y, col = 'red', lty = 2)
legend('topright', legend = c('Original Data', 'Reconstructed Model'), col = c('blue', 'red'), lty = 1:2)


#################################################
# method 2

smoothed_y <- stats::filter(y, rep(1/5, 5), sides = 2)  # Example: 5-point moving average
plot(t, y, type = "l", col = "gray", main = "Original vs. Smoothed Data")
lines(t, smoothed_y, col = "blue")

log_smoothed_y <- log(smoothed_y)  # Log transformation of smoothed data

# Estimate slopes manually between selected points (for example, points 10-20 and 20-30)
beta1_est <- coef(lm(log_smoothed_y[10:20] ~ t[10:20]))[2]
beta2_est <- coef(lm(log_smoothed_y[20:30] ~ t[20:30]))[2]

X <- cbind(1, exp(beta1_est * t), exp(beta2_est * t))

alpha_est <- lm(y ~ X - 1)$coefficients  # -1 to exclude intercept as alpha_0 is already included
alpha0_est <- alpha_est[1]
alpha1_est <- alpha_est[2]
alpha2_est <- alpha_est[3]

# Load the minpack.lm package
#install.packages("minpack.lm")  # Uncomment if not installed
library(minpack.lm)

# Refine parameters with nlsLM using initial guesses from Obson's method
fit <- nlsLM(y ~ alpha0 + alpha1 * exp(beta1 * t) + alpha2 * exp(beta2 * t),
             start = list(alpha0 = alpha0_est, alpha1 = alpha1_est, alpha2 = alpha2_est,
                          beta1 = beta1_est, beta2 = beta2_est),
             data = data.frame(t = t, y = y))
summary(fit)
alpha0
alpha_est


# Example data
original_data <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

# Using a simple moving average with a window size of 3
library(zoo)
smoothed_data <- rollmean(original_data, k = 3, fill = NA, align = "center")

# Check the length of the original and smoothed data
length(original_data)  # Original length
length(smoothed_data)   # Smoothed length

