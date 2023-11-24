#LOADING THE DATA AFTER INSTALLING NECESSARY PACKAGES AND LIBRARY

library(readxl)
#LOADING THE DATA AND CHANGING THE COLUMN NAMES
data <- read.csv("gene_data.csv", header = TRUE)
data
colnames(data) <- c("Time", "x1", "x2", "x3", "x4", "x5")
data

#SET X AND Y
x <- data[,c("x1", "x3", "x4", "x5")]
x <- as.matrix(x)
y <- data$x2
t <- data$Time
####################################
#TASK1: PRELIMINARY DATA ANALYSIS
####################################


#TASK 1.1: TIME SERIES PLOTS OF EACH GENE AGAINST TIME

par(mfrow = c(1, 3))  # Set up a 1x3 grid for time series plots

for (gene in c("x1", "x2", "x3", "x4", "x5")) {
  ts_gene <- ts(data[[gene]], start = c(-1, 1), frequency = 1)
  plot(ts_gene, main = paste("Time Series Analysis of", gene))
}


##################################################
#TASK 1.2: DISTRIBUTION FOR EACH GENE(TIME-SERIES)

# Setting up a 2x3 grid for distribution plots
par(mfrow = c(2, 3))

# Loop through each gene and create distribution plots
for (gene in c("x1", "x2", "x3", "x4", "x5")) {
  # Create a time-series object
  ts_gene <- ts(data[[gene]], start = c(-1, 1), frequency = 1)
  
  # Plot distribution of the gene over time
  hist(ts_gene, main = paste("Distribution of", gene, "(Time Series)"),
       xlab = gene, col = "skyblue", border = "black")
}

#########################################################################
#TASK 1.3:CORRELATION AND SCATTER PLOTS BETWEEN DIFFERENT COMBINATION OF TWO GENES

pairs(data[, c("x1", "x2", "x3", "x4", "x5")], 
      main = "Scatterplot Matrix", pch = 16, col = "skyblue")

# Calculate pairwise correlation coefficients
cor_matrix <- cor(data[, c("x1", "x2", "x3", "x4", "x5")])

# Display correlation matrix
print("Correlation Matrix:")
print(cor_matrix)

###########################################################################
#*******************END OF TASK 1*****************************************#


#TASK 2: Task 2: Regression – modelling the relationship between gene expression
################################################################################

#EXTRACTING RELEVANT COLUMNS
x1 <- data$x1
x3 <- data$x3
x4 <- data$x4
x5 <- data$x5
y <- data$x2  # AS X2 WILL BE THE OUTPUT VARIABLE

#DEFINING THE GIVEN MODELS

model1 <- lm(y ~ x4 + I(x3^2), data = data)
model2 <- lm(y ~ x4 + I(x3^2) + x5, data = data)
model3 <- lm(y ~ x3 + x4 + I(x5^3), data = data)
model4 <- lm(y ~ x4 + I(x3^2) + I(x5^3), data = data)
model5 <- lm(y ~ x4 + I(x1^2) + I(x3^2), data = data)

models <- list(model1, model2, model3, model4, model5)

#####################################################

#TASK 2.1: ESTIMATING MODEL PARAMETERS FOR EACH CANDIDATE MODEL USING LEAST SQUARE

theta_hat <- lapply(models, coef)
for (i in seq_along(models)) {
  cat("Model", i, "Parameters:", theta_hat[[i]], "\n")
}
#######################################################

#TASK 2.2: COMPUTE MODEL RESIDUAL (ERROR) SUM OF SQUARED ERRORS(RSS) FOR EVERY CANDIDATE MODEL

# Get the residuals for model 1
residuals_model1 <- residuals(model1)
# Compute RSS for model 1
RSS_model1 <- sum(residuals_model1^2)
RSS_model1


# Get the residuals for model 2
residuals_model2 <- residuals(model2)
# Compute RSS for model 2
RSS_model2 <- sum(residuals_model2^2)
RSS_model2



# Get the residuals for model 3
residuals_model3 <- residuals(model3)
# Compute RSS for model 3
RSS_model3 <- sum(residuals_model3^2)
RSS_model3

# Get the residuals for model 4
residuals_model4 <- residuals(model4)
# Compute RSS for model 4
RSS_model4 <- sum(residuals_model4^2)
RSS_model4

# Get the residuals for model 5
residuals_model5 <- residuals(model5)
# Compute RSS for model 5
RSS_model5 <- sum(residuals_model5^2)
RSS_model5

######################################################################

#TASK 2.3:COMPUTE THE LOG-LIKELIHOOD FUNCTION FOR EVERY CANDIDATE MODEL
##FIRST LETS DEFINE A FUNCTION FOR LOG-LIKEHOOD CALCULATION
compute_log_likelihood <- function(model, data) {
  #The residuals model
  residuals_model <- residuals(model)
  RSS_model <- sum(residuals_model^2)
  n_obs <- nrow(data)
  
  # Estimate of the residual standard deviation
  sigma_hat <- sqrt(RSS_model / (n_obs - length(coef(model))))
  
  # Compute log-likelihood
  log_likelihood <- -n_obs/2 * log(2 * pi) - n_obs/2 * log(sigma_hat^2) - 1/(2 * sigma_hat^2) * RSS_model
  
  return(log_likelihood)
}

#LOG-LIKELIHOOD FOR MODEL1
log_likelihood_model1 <- compute_log_likelihood(model1, data)
print(log_likelihood_model1)

#LOG-LIKELIHOOD FOR MODEL2
log_likelihood_model2 <- compute_log_likelihood(model2, data)
print(log_likelihood_model2)

#LOG-LIKELIHOOD FOR MODEL3
log_likelihood_model3 <- compute_log_likelihood(model3, data)
print(log_likelihood_model3)

#LOG-LIKELIHOOD FOR MODEL4
log_likelihood_model4 <- compute_log_likelihood(model4, data)
print(log_likelihood_model4)

#LOG-LIKELIHOOD FOR MODEL5
log_likelihood_model5 <- compute_log_likelihood(model5, data)
print(log_likelihood_model5)

####################################################################
#TASK 2.4: COMPUTING THE AKAIKE INFORMATION CRITERION (AIC) AND BAYESIAN INFORMATION CRITERION(BIC)

#WE ALREADY HAVE THE LOG-LIKELIHOOD FUCNTION FOR EACH PARAMETER SO NOW LETS ADJUST THE NO OF PARAMETERS FOR EACH MODEL

k_m1 <- 3
k_m2 <- 4
k_m3 <- 3
k_m4 <- 4
k_m5 <- 4

#NOW LETS CALCULATE AIC AND BIC FOR EACH CANDIDATE
AIC_m1 <- 2 * k_m1 - 2 * log_likelihood_model1
AIC_m2 <- 2 * k_m2 - 2 * log_likelihood_model2
AIC_m3 <- 2 * k_m3 - 2 * log_likelihood_model3
AIC_m4 <- 2 * k_m4 - 2 * log_likelihood_model4
AIC_m5 <- 2 * k_m5 - 2 * log_likelihood_model5
# Print AIC values
cat("AIC (Model 1):", AIC_m1, "\n")
cat("AIC (Model 2):", AIC_m2, "\n")
cat("AIC (Model 3):", AIC_m3, "\n")
cat("AIC (Model 4):", AIC_m4, "\n")
cat("AIC (Model 5):", AIC_m5, "\n")



BIC_m1 <- k_m1 * log(n) - 2 * log_likelihood_model1
BIC_m2 <- k_m2 * log(n) - 2 * log_likelihood_model2
BIC_m3 <- k_m3 * log(n) - 2 * log_likelihood_model3
BIC_m4 <- k_m4 * log(n) - 2 * log_likelihood_model4
BIC_m5 <- k_m5 * log(n) - 2 * log_likelihood_model5
#Print BIC values
cat("BIC (Model 1):", BIC_m1, "\n")
cat("BIC (Model 2):", BIC_m2, "\n")
cat("BIC (Model 3):", BIC_m3, "\n")
cat("BIC (Model 4):", BIC_m4, "\n")
cat("BIC (Model 5):", BIC_m5, "\n")

################################################################
#TASK 2.5: CHECKING THE DISTRIBUTION OF MODEL PREDICTION ERROR



# Function to create Q-Q plot for model residuals
qq_plot <- function(model, data) {
  # Get residuals
  residuals <- residuals(model)
  
  # Create Q-Q plot
  qqnorm(residuals, main = paste("Q-Q Plot for Model Residuals:", deparse(substitute(model))))
  qqline(residuals)
}

# Create Q-Q plots for each model
qq_plot(model1, data)
qq_plot(model2, data)
qq_plot(model3, data)
qq_plot(model4, data)
qq_plot(model5, data)


##############################################################
#TASK 2.7:

# Load necessary libraries
library(ggplot2)

# Model 4 is the selected as  best model
best_model <- lm(x2 ~ x4 + I(x3^2) + I(x5^3), data = train_data)  

# Step 1: Split the dataset into training and testing sets (70% for training, 30% for testing)
set.seed(123)  # for reproducibility
train <- sample(1:nrow(data), 0.7 * nrow(data))
train_data <- data[train, ]
test_data <- data[-train, ]

# Step 2: Fit the best model on the training data
best_model <- lm(x2 ~ x4 + I(x3^2) + I(x5^3), data = train_data) #y is now x2

# Step 3: Predict on the testing data
predictions <- predict(best_model, newdata = test_data, interval = "prediction", level = 0.95)

# Step 4: Plot model predictions with error bars
ggplot() +
  geom_point(data = test_data, aes(x = Time, y = x2), color = "blue", size = 2, shape = 16) +
  geom_line(data = test_data, aes(x = Time, y = predictions[, 1]), color = "red", size = 1) +
  geom_ribbon(data = data.frame(Time = test_data$Time, fit = predictions[, 1], lwr = predictions[, 2], upr = predictions[, 3]),
              aes(x = Time, y = fit, ymin = lwr, ymax = upr), alpha = 0.3, fill = "red") +
  labs(title = "Model Predictions on Testing Data with 95% Confidence Intervals",
       x = "Time",
       y = "Gene Expression") +
  theme_minimal()

#******************END OF TASK 2*******************************************************##

#TASK 3: DESCRIPTIVE/ INFERENTIAL STATISTICS
############################################################

#TASK: 3.1

# Compute sample mean and variance for each gene
mean_values <- apply(data[, c("x1", "x2", "x3", "x4", "x5")], 2, mean)
var_values <- apply(data[, c("x1", "x2", "x3", "x4", "x5")], 2, var)

# Number of observations
n <- nrow(data)

# Confidence levels
confidence_levels <- c(0.90, 0.95, 0.99)

# Degrees of freedom for t-distribution (mean)
df_mean <- n - 1

# Confidence intervals for mean
mean_conf_intervals <- sapply(confidence_levels, function(alpha) {
  margin_of_error <- qt((1 - alpha) / 2, df_mean) * sqrt(var_values / n)
  lower <- mean_values - margin_of_error
  upper <- mean_values + margin_of_error
  return(data.frame(lower, upper))
}, simplify = FALSE)

# Degrees of freedom for chi-square distribution (variance)
df_var <- n - 1

# Confidence intervals for variance
var_conf_intervals <- sapply(confidence_levels, function(alpha) {
  lower <- (n - 1) * var_values / qchisq((1 + alpha) / 2, df_var)
  upper <- (n - 1) * var_values / qchisq((1 - alpha) / 2, df_var)
  return(data.frame(lower, upper))
}, simplify = FALSE)

# Print results
for (i in seq_along(confidence_levels)) {
  cat(paste("Gene", 1:5, "Mean  ", confidence_levels[i] * 100, "% CI: ", 
            sprintf("[%0.3f, %0.3f]\n", mean_conf_intervals[[i]]$lower, mean_conf_intervals[[i]]$upper)))
  cat(paste("Gene", 1:5, "Var   ", confidence_levels[i] * 100, "% CI: ", 
            sprintf("[%0.3f, %0.3f]\n", var_conf_intervals[[i]]$lower, var_conf_intervals[[i]]$upper)))
  cat("\n")
}
#########################################

#TASK 3.2

install.packages("moments")
library(moments)


# Compute skewness for each gene
skewness_values <- apply(data[, c("x1", "x2", "x3", "x4", "x5")], 2, skewness)

# Plot histogram with normal distribution overlay
par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))

for (i in 1:5) {
  # Plot histogram
  hist(data[, i], main = paste("Gene", i), xlab = "Value", col = "lightblue", ylim = c(0, 30), freq = FALSE)
  
  # Overlay normal distribution
  mu <- mean(data[, i])
  sigma <- sd(data[, i])
  x <- seq(min(data[, i]), max(data[, i]), length = 100)
  y <- dnorm(x, mean = mu, sd = sigma)
  lines(x, y * diff(hist(data[, i], plot = FALSE)$counts[1:2]) * nrow(data), col = "red", lwd = 2)
  
  # Mark mode
  mode <- x[which.max(y)]
  points(mode, dnorm(mode, mean = mu, sd = sigma) * diff(hist(data[, i], plot = FALSE)$counts[1:2]) * nrow(data), col = "green", pch = 16)
}

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

#***********************END OF TASK 3***************************#######