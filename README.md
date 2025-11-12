MTH686: Non-Linear Regression Project

This repository contains the R solution for the MTH686 Non-Linear Regression project at IIT Kanpur. The project involves fitting and analyzing a given time-series dataset (set-49.dat) against three different candidate models.

Models Analyzed

Model 1 (Exponential): $y(t)=\alpha_{0}+\alpha_{1}e^{\beta_{1}t}+\alpha_{2}e^{\beta_{2}t}+\epsilon(t)$

Model 2 (Rational): $y(t)=\frac{\alpha_{0}+\alpha_{1}t}{\beta_{0}+\beta_{1}t}+\epsilon(t)$

Model 3 (Polynomial): $y(t)=\beta_{0}+\beta_{1}t+\beta_{2}t^{2}+\beta_{3}t^{3}+\beta_{4}t^{4}+\epsilon(t)$

Core Tasks

Least Squares Estimation (LSE): Implements robust Non-Linear Least Squares (NLS) using R's optim function (L-BFGS-B method) to find the parameters for the non-linear models (Model 1 and 2).

Numerical Stability: Uses the original time variable (t_orig) and stable methods (NLS, poly(t, 4)) to avoid numerical overflow and ill-conditioned matrices.

Statistical Analysis:

Calculates the Fisher Information Matrix (FIM) and its inverse to determine parameter standard errors.

Finds the 95% Confidence Intervals (CIs) for all parameters based on asymptotic normality.

Estimates the error variance ($\sigma^2$) for each model.

Residual Diagnostics:

Generates Residuals vs. Time plots, Q-Q plots, and Histograms.

Performs the Kolmogorov-Smirnov (KS) test to check for normality of residuals.

Model Comparison: Compares all three models using their RSS, R-Squared, and Adjusted R-Squared values to determine the "best" fitted model.
