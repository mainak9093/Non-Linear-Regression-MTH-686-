setwd("D:/Sem 5/MTH686/MTH 686 Non Linear Reg/Project")  
library(readr)
getwd()
list.files()
file.exists("set-49.dat")
# then:
# data <- readLines("data-42.txt")
 data <- read.table("set-49.dat")
 head(data)
 
 # Extract time and observed values
 t <- as.numeric(unlist(data[1]))  # Replace with time values
 t <- t/100
 y1 <- as.numeric(unlist(data[2])) # Replace with observed y(t) values
 
 
 # Plot the time series data
 plot(data$V1, data$V2, type = "l", main = "Time Series Data", xlab = "Time", ylab = "Value")
 
 # Fit a simple linear model for residual analysis
 fit <- lm(V2 ~ V1, data = data) 
 residuals <- resid(fit)
 
 # Residual plot
 plot(data$V1, residuals, main = "Residual Plot", xlab = "Time", ylab = "Residuals")
 abline(h = 0, col = "red")
 
 # Calculate Signal-to-Noise Ratio (SNR)
 signal <- mean(y1, na.rm = TRUE)
 noise <- sd(y1, na.rm = TRUE)
 SNR <- signal / noise
 cat("Signal-to-Noise Ratio (SNR):", SNR, "\n")
 
 
 ####### Model 1 #########
 
 # Define the model function for Prony's method
 model1_function <- function(t, alpha0, alpha1, alpha2, beta1, beta2) {
   return(alpha0 + alpha1 * exp(beta1 * t) + alpha2 * exp(beta2 * t))
 }
 
 # Constructing the delayed values for Prony's method
 Y1 <- y1[1:(length(y1) - 2)]
 Y2 <- y1[2:(length(y1) - 1)]
 Y3 <- y1[3:length(y1)]
 
 # Formulate the linear system for Prony's method
 A <- cbind(Y2, Y1)  
 B <- Y3             
 
 # Solve for coefficients (c1, c2) in the characteristic polynomial equation
 coef <- lm(B ~ A - 1)$coefficients  
 c1 <- coef[1]
 c2 <- coef[2]
 
 # Find roots of the characteristic polynomial to estimate beta parameters
 roots <- polyroot(c(1, -c1, -c2))
 beta_estimates <- log(roots)
 
 # Extract real parts of beta estimates (assuming real roots for simplicity)
 beta1_m1 <- Re(beta_estimates[1])  
 beta2_m1 <- Re(beta_estimates[2])
 
 # Construct the X matrix with the estimated beta values
 X <- cbind(1, exp(beta1_m1 * t), exp(beta2_m1 * t))
 
 # Solve for alpha parameters using least squares
 alpha_estimates <- lm(y1 ~ X -1 )$coefficients
 alpha0_m1 <- alpha_estimates[1]
 alpha1_m1 <- alpha_estimates[2]
 alpha2_m1 <- alpha_estimates[3]
 
 # Display estimated parameters
 cat("Estimated parameters:\n")
 cat("alpha0:", alpha0_m1, "\n")
 cat("alpha1:", alpha1_m1, "\n")
 cat("alpha2:", alpha2_m1, "\n")
 cat("beta1:", beta1_m1, "\n")
 cat("beta2:", beta2_m1, "\n")
 
 # Model predictions using the estimated parameters
 y2 <- model1_function(t, alpha0_m1, alpha1_m1, alpha2_m1, beta1_m1, beta2_m1) 
 
 # Plot observed vs. fitted data
 plot(t, y1, type = "l", col = "blue", main = "Observed vs. Fitted Data", xlab = "Time", ylab = "Values")
 lines(t, y2, col = "red", lwd = 2)
 legend("topright", legend = c("Observed Data", "Prony's Model Fit"), col = c("blue", "red"), lty = 1, lwd = 2)
 
 # Calculate model accuracy metrics
 rmse_1 <- sqrt(mean((y1 - y2)^2, na.rm = TRUE))
 rss_1 <- sum((y1 - y2)^2)                  
 tss_1 <- sum((y1 - mean(y1))^2)  
 r_squared_1 <- 1 - (rss_1 / tss_1)
 adjusted_r_squared_1 <- 1 - ((1 - r_squared_1) * (n - 1) / (n - 3))
 
 
 
 # Display model accuracy metrics
 cat("Model Accuracy Metrics:\n")
 cat("Root Mean Squared Error (RMSE):", rmse_1, "\n")
 cat("R-squared:", r_squared_1, "\n")
 cat("Adjusted R-squared:", adjusted_r_squared_1, "\n")
 
 ####### Model 3 #########
 
 
 # Fit a 4th-degree polynomial regression model
 model3 <- lm(y1 ~ poly(t, 4, raw = TRUE))
 
 # Display summary of the model
 summary(model3)
 
 # Extract coefficients (OLSE for ??0, ??1, ??2, ??3, and ??4)
 beta_estimates <- coef(model3)
 cat("Estimated coefficients:\n")
 cat("??0:", beta_estimates[1], "\n")
 cat("??1:", beta_estimates[2], "\n")
 cat("??2:", beta_estimates[3], "\n")
 cat("??3:", beta_estimates[4], "\n")
 cat("??4:", beta_estimates[5], "\n")
 
 # Calculate fitted values and residuals
 fitted_values_3 <- fitted(model3)
 residuals_3 <- resid(model3)
 
 # Plot Observed vs. Fitted Data
 plot(t, y1, type = "l", col = "blue", main = "Observed vs. Fitted Data", xlab = "Time", ylab = "Values")
 lines(t, fitted_values_3, col = "red", lwd = 2)
 legend("topright", legend = c("Observed Data", "Fitted Polynomial Model"), col = c("blue", "red"), lty = 1, lwd = 2)
 
 # Residual Plot to check for patterns
 plot(t, residuals_3, main = "Residuals vs. Time", xlab = "Time", ylab = "Residuals")
 abline(h = 0, col = "red")
 
 
 # Calculate model accuracy metrics
 rmse_3 <- sqrt(mean(residuals_3^2, na.rm = TRUE))
 rss_3 <- sum(residuals_3^2)                  
 tss_3 <- sum((y1 - mean(y1))^2)  
 r_squared_3 <- 1 - (rss_3 / tss_3)
 adjusted_r_squared_3 <- 1 - ((1 - r_squared_3) * (n - 1) / (n - 5))
 
 # Display model accuracy metrics
 cat("Model Accuracy Metrics:\n")
 cat("Root Mean Squared Error (RMSE):", rmse_3, "\n")
 cat("R-squared:", r_squared_3, "\n")
 cat("Adjusted R-squared:", adjusted_r_squared_3, "\n")
 
 ####### Model 2 #########
 
 
 # Define grid ranges for beta0 and beta1, excluding zero
 beta0_range <- seq(-10, -0.1, by = 0.1)  # Values from -1 to -0.1
 beta0_range <- c(beta0_range, seq(0.1, 10, by = 0.1))  # Values from 0.1 to 1
 beta1_range <- seq(-5, -0.05, by = 0.05)  # Values from -0.5 to -0.05
 beta1_range <- c(beta1_range, seq(0.05, 5, by = 0.05))  # Values from 0.05 to 0.5
 
 # Initialize a data frame to store beta values and corresponding RSS
 results <- data.frame(beta0_2 = numeric(0), beta1_2 = numeric(0), rss_2 = numeric(0))
 
 # Track the best RSS and parameter values
 best_rss_2 <- Inf
 best_params_2 <- c()
 
 # Perform grid search
 for (beta0_2 in beta0_range) {
   for (beta1_2 in beta1_range) {
     # Calculate transformed y(t) for each combination of beta0 and beta1
     y_trans <- y1 * (beta0_2 + beta1_2 * t)
     
     # Apply linear regression to estimate alpha0 and alpha1
     fit <- lm(y_trans ~ t)
     alpha0_2 <- coef(fit)[1]
     alpha1_2 <- coef(fit)[2]
     
     # Predict y(t) using the estimated alphas and current beta values
     y_pred <- (alpha0_2 + alpha1_2 * t) / (beta0_2 + beta1_2 * t)
     
     # Calculate RSS
     rss_2 <- sum((y1 - y_pred)^2, na.rm = TRUE)
     
     # Store results in the data frame
     results <- rbind(results, data.frame(beta0_2 = beta0_2, beta1_2 = beta1_2, rss_2 = rss_2))
     
     # Update if this is the best set of parameters found
     if (rss_2 < best_rss_2) {
       best_rss_2 <- rss_2
       best_params <- c(alpha0_2, alpha1_2, beta0_2, beta1_2)
       best_y_pred <- y_pred  # Save best predictions
     }
   }
 }
 
 
 # Set initial parameter estimates from the best grid search results
 alpha0 <- best_params[1]
 alpha1 <- best_params[2]
 beta0 <- best_params[3]
 beta1 <- best_params[4]
 
 # Define convergence criteria and maximum iterations
 tolerance <- 1e-6
 max_iter <- 100
 converged <- FALSE
 
 # Track RSS in each iteration for convergence
 previous_rss <- best_rss_2
 
 # Gauss-Newton Iteration
 for (iter in 1:max_iter) {
   # Compute predicted values using current parameter estimates
   y_pred <- (alpha0 + alpha1 * t) / (beta0 + beta1 * t)
   
   # Calculate residuals
   residuals <- y1 - y_pred
   
   # Construct the Jacobian matrix (partial derivatives of y_pred with respect to parameters)
   J <- cbind(
     1 / (beta0 + beta1 * t),                # ???y_pred/???alpha0
     t / (beta0 + beta1 * t),                # ???y_pred/???alpha1
     -(alpha0 + alpha1 * t) / (beta0 + beta1 * t)^2,  # ???y_pred/???beta0
     -t * (alpha0 + alpha1 * t) / (beta0 + beta1 * t)^2  # ???y_pred/???beta1
   )
   
   # Calculate the step adjustment using least squares solution
   lambda <- 1e-5  # Small regularization factor
   delta <- solve(t(J) %*% J + lambda * diag(ncol(J))) %*% t(J) %*% residuals
   
   # Update parameter estimates
   alpha0 <- alpha0 + delta[1]
   alpha1 <- alpha1 + delta[2]
   beta0 <- beta0 + delta[3]
   beta1 <- beta1 + delta[4]
   
   # Compute new RSS
   y_pred <- (alpha0 + alpha1 * t) / (beta0 + beta1 * t)
   rss_2 <- sum((y1 - y_pred)^2)
   
   # Check for convergence
   if (abs(previous_rss - rss_2) < tolerance) {
     converged <- TRUE
     break
   }
   
   # Update previous RSS for next iteration
   previous_rss <- rss_2
 }
 
 # Display final parameter estimates
 if (converged) {
   cat("Gauss-Newton converged after", iter, "iterations.\n")
 } else {
   cat("Gauss-Newton did not converge within", max_iter, "iterations.\n")
 }
 
 cat("Final parameter estimates:\n")
 cat("alpha0:", alpha0, "\n")
 cat("alpha1:", alpha1, "\n")
 cat("beta0:", beta0, "\n")
 cat("beta1:", beta1, "\n")
 cat("RSS for final model:", rss_2, "\n")
 
 
 # Model accuracy metrics
 tss_2 <- sum((y1 - mean(y1))^2)
 r_squared_2 <- 1 - (rss_2 / tss_2)
 adjusted_r_squared_2<- 1 - ((1 - r_squared_2) * (length(y1) - 1) / (length(y1) - 3))  # Adjusted R^2 for 2 predictors
 
 cat("RSS for best model:", rss_2, "\n")
 cat("R-squared:", r_squared_2, "\n")
 cat("Adjusted R-squared:", adjusted_r_squared_2, "\n")
 
 # Plot observed vs. fitted values
 plot(t, y1, col = "blue", type = "l", main = "Observed vs. Fitted Data", xlab = "Time", ylab = "Values")
 lines(t, y_pred, col = "red", lwd = 2)
 legend("topright", legend = c("Observed", "Fitted"), col = c("blue", "red"), lty = 1, lwd = 2)
 
 # Residual plot
 residuals <- y1 - y_pred
 plot(t, residuals, main = "Residual Plot", xlab = "Time", ylab = "Residuals")
 abline(h = 0, col = "red")
 
 ############ 
 # Assuming your model y(t) = (alpha_0 + alpha_1 * t) / (beta_0 + beta_1 * t) + e(t)
 # Let's define the Jacobian matrix, Fisher Information, and Confidence Interval calculation steps
 
 # Define the model function
 model <- function(t, alpha0, alpha1, beta0, beta1) {
   return((alpha0 + alpha1 * t) / (beta0 + beta1 * t))
 }
 
 # Define the Jacobian of the model (partial derivatives with respect to alpha0, alpha1, beta0, beta1)
 jacobian <- function(t, alpha0, alpha1, beta0, beta1) {
   d_alpha0 <- 1 / (beta0 + beta1 * t)  # ???y/???alpha0
   d_alpha1 <- t / (beta0 + beta1 * t)  # ???y/???alpha1
   d_beta0 <- -(alpha0 + alpha1 * t) / ((beta0 + beta1 * t)^2)  # ???y/???beta0
   d_beta1 <- -(alpha0 + alpha1 * t) * t / ((beta0 + beta1 * t)^2)  # ???y/???beta1
   
   return(c(d_alpha0, d_alpha1, d_beta0, d_beta1))  # Return as a row vector
 }
 
 # Fisher Information matrix calculation (assuming error variance ??^2 is known or estimated)
 fisher_information <- function(t, alpha0, alpha1, beta0, beta1, sigma2) {
   n <- length(t)  # number of data points
   J <- matrix(0, nrow = n, ncol = 4)  # Initialize the Jacobian matrix
   
   for (i in 1:n) {
     J[i, ] <- jacobian(t[i], alpha0, alpha1, beta0, beta1)  # Store Jacobian row for each t
   }
   
   # Calculate Fisher Information Matrix (sum of squared Jacobian)
   fisher_matrix <- matrix(0, nrow = 4, ncol = 4)  # Initialize Fisher Information Matrix
   
   for (i in 1:n) {
     fisher_matrix <- fisher_matrix + outer(J[i, ], J[i, ])  # Add outer product of Jacobian row
   }
   
   # Multiply by 1/n to get the Fisher Information
   fisher_matrix <- fisher_matrix / n
   
   # Assuming we have variance estimate (sigma2) available, scale by sigma^2
   fisher_matrix <- fisher_matrix / sigma2
   
   return(c(fisher_matrix,J))
 }
 
 # Load the MASS package for pseudo-inverse
 library(MASS)
 
 
 # Confidence Intervals calculation (assuming the Fisher Information matrix is available)
 confidence_intervals <- function(fisher_matrix, parameter_estimates, alpha = 0.05) {
   # Inverse Fisher Information Matrix (for variance)
   fisher_inv <- ginv(fisher_matrix)
   
   # Standard errors (diagonal of the inverse Fisher Information)
   se <- sqrt(diag(fisher_inv))
   
   # Z value for 95% confidence
   z_value <- qnorm(1 - alpha / 2)
   
   # Confidence intervals for each parameter
   ci_lower <- parameter_estimates - z_value * se
   ci_upper <- parameter_estimates + z_value * se
   
   return(list(lower = ci_lower, upper = ci_upper))
 }
 
 # Example values (you should replace with your actual model parameters and data)
 alpha0 <- alpha0_2
 alpha1 <- alpha1_2
 beta0 <- beta0_2
 beta1 <- beta1_2
 sigma2 <- rss_2/(60-2)  # Assume sigma^2 is known or estimated from residuals
 
 # Fisher Information Matrix for the model parameters
 fisher_matrix <- fisher_information(t, alpha0, alpha1, beta0, beta1, sigma2)[1]
 
 # Assume parameter estimates from your model fitting procedure (use the results from your fitting method)
 parameter_estimates <- c(alpha0, alpha1, beta0, beta1)
 
 # Calculate Confidence Intervals
 ci <- confidence_intervals(fisher_matrix, parameter_estimates)
 
 # Print the results
 cat("Confidence Intervals:\n")
 print(ci)
 print(fisher_matrix)
 ##############################
 
 # Define the model function
 model <- function(t, alpha0, alpha1, beta0, beta1) {
   return((alpha0 + alpha1 * t) / (beta0 + beta1 * t))
 }
 
 # Define the Jacobian of the model (partial derivatives with respect to alpha0, alpha1, beta0, beta1)
 jacobian <- function(t, alpha0, alpha1, beta0, beta1) {
   d_alpha0 <- 1 / (beta0 + beta1 * t)  # ???y/???alpha0
   d_alpha1 <- t / (beta0 + beta1 * t)  # ???y/???alpha1
   d_beta0 <- -(alpha0 + alpha1 * t) / ((beta0 + beta1 * t)^2)  # ???y/???beta0
   d_beta1 <- -(alpha0 + alpha1 * t) * t / ((beta0 + beta1 * t)^2)  # ???y/???beta1
   
   return(c(d_alpha0, d_alpha1, d_beta0, d_beta1))  # Return as a row vector
 }
 
 # Fisher Information matrix calculation (assuming error variance ??^2 is known or estimated)
 fisher_information <- function(t, alpha0, alpha1, beta0, beta1, sigma2) {
   n <- length(t)  # number of data points
   J <- matrix(0, nrow = n, ncol = 4)  # Initialize the Jacobian matrix
   
   for (i in 1:n) {
     J[i, ] <- jacobian(t[i], alpha0, alpha1, beta0, beta1)  # Store Jacobian row for each t
   }
   
   # Calculate Fisher Information Matrix (sum of squared Jacobian)
   fisher_matrix <- matrix(0, nrow = 4, ncol = 4)  # Initialize Fisher Information Matrix
   
   for (i in 1:n) {
     fisher_matrix <- fisher_matrix + outer(J[i, ], J[i, ])  # Add outer product of Jacobian row
   }
   
   # Multiply by 1/n to get the Fisher Information
   fisher_matrix <- fisher_matrix / n
   
   # Assuming we have variance estimate (sigma2) available, scale by sigma^2
   fisher_matrix <- fisher_matrix / sigma2
   
   return(fisher_matrix)  # Return the Fisher Information Matrix
 }
 
 # Load the MASS package for pseudo-inverse
 library(MASS)
 
 # Confidence Intervals calculation (assuming the Fisher Information matrix is available)
 confidence_intervals <- function(fisher_matrix, parameter_estimates, alpha = 0.05) {
   # Inverse Fisher Information Matrix (for variance)
   fisher_inv <- ginv(fisher_matrix)
   
   # Standard errors (diagonal of the inverse Fisher Information)
   se <- sqrt(diag(fisher_inv))
   
   # Z value for 95% confidence
   z_value <- qnorm(1 - alpha / 2)
   
   # Confidence intervals for each parameter
   ci_lower <- parameter_estimates - z_value * se
   ci_upper <- parameter_estimates + z_value * se
   
   return(list(lower = ci_lower, upper = ci_upper))
 }
 
 # Example values (replace with actual model parameters and data)
 alpha0 <- alpha0_2
 alpha1 <- alpha1_2
 beta0 <- beta0_2
 beta1 <- beta1_2
 sigma2 <- rss_2/(60-2)   
 
 # Calculate Fisher Information Matrix for the model parameters
 fisher_matrix <- fisher_information(t, alpha0, alpha1, beta0, beta1, sigma2)
 
 # Assume parameter estimates from your model fitting procedure
 parameter_estimates <- c(alpha0, alpha1, beta0, beta1)
 
 # Calculate Confidence Intervals
 ci <- confidence_intervals(fisher_matrix, parameter_estimates)
 
 # Print the results
 cat("Fisher Information Matrix:\n")
 print(fisher_matrix)
 
 cat("\nConfidence Intervals:\n")
 print(ci)
 
 ##############################
 
 