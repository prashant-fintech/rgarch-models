# Load required libraries
library(rugarch)
library(fGarch)

# Configuration parameters
input_file <- "reliance.csv"
price_column <- 2
output_prefix <- "reliance_vol_estimates_EGARCH"

# Function to load and prepare data
load_data <- function(file_path, col_index) {
  price <- read.csv(file_path, header = TRUE)
  price_series <- price[, col_index]
  n <- length(price_series)
  returns <- diff(price_series) / lag(price_series)
  returns <- returns[1:(n - 1)]
  return(returns)
}

# Function to fit EGARCH model and generate results
fit_egarch_model <- function(returns, dist_model, output_file) {
  spec <- ugarchspec(
    variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
    distribution.model = dist_model
  )
  
  garch_fit <- ugarchfit(spec = spec, data = returns)
  print(garch_fit)
  
  sigma <- garch_fit@fit$sigma
  write.csv(sigma, file = output_file)
  
  forecast <- ugarchforecast(garch_fit, n.ahead = 1)
  forecast_sigma <- forecast@forecast$sigmaFor[[1]]
  print(forecast_sigma)
  
  return(garch_fit)
}

# Function to calculate and print quantiles
calculate_quantiles <- function(garch_fit, dist_type) {
  resid <- residuals(garch_fit, standardize = TRUE)
  
  if (dist_type == "std") {
    fit_params <- stdFit(resid)
    nu <- fit_params$par[3]
    print(qstd(0.05, mean = 0, sd = 1, nu = nu))
    print(qstd(0.01, mean = 0, sd = 1, nu = nu))
    print(qstd(0.005, mean = 0, sd = 1, nu = nu))
  } else if (dist_type == "ged") {
    fit_params <- gedFit(resid)
    nu <- fit_params$par[3]
    print(qged(0.05, mean = 0, sd = 1, nu = nu))
    print(qged(0.01, mean = 0, sd = 1, nu = nu))
    print(qged(0.005, mean = 0, sd = 1, nu = nu))
  } else if (dist_type == "snorm") {
    fit_params <- snormFit(resid)
    xi <- fit_params$par[3]
    print(qsnorm(0.05, mean = 0, sd = 1, xi = xi))
    print(qsnorm(0.01, mean = 0, sd = 1, xi = xi))
    print(qsnorm(0.005, mean = 0, sd = 1, xi = xi))
    print(qsnorm(0.95, mean = 0, sd = 1, xi = xi))
    print(qsnorm(0.99, mean = 0, sd = 1, xi = xi))
    print(qsnorm(0.995, mean = 0, sd = 1, xi = xi))
  } else if (dist_type == "sstd") {
    fit_params <- sstdFit(resid)
    nu <- fit_params$estimate[3]
    xi <- fit_params$estimate[4]
    print(qsstd(0.05, mean = 0, sd = 1, nu = nu, xi = xi))
    print(qsstd(0.01, mean = 0, sd = 1, nu = nu, xi = xi))
    print(qsstd(0.005, mean = 0, sd = 1, nu = nu, xi = xi))
    print(qsstd(0.95, mean = 0, sd = 1, nu = nu, xi = xi))
    print(qsstd(0.99, mean = 0, sd = 1, nu = nu, xi = xi))
    print(qsstd(0.995, mean = 0, sd = 1, nu = nu, xi = xi))
  } else if (dist_type == "sged") {
    fit_params <- sgedFit(resid)
    nu <- fit_params$par[3]
    xi <- fit_params$par[4]
    print(qsged(0.05, mean = 0, sd = 1, nu = nu, xi = xi))
    print(qsged(0.01, mean = 0, sd = 1, nu = nu, xi = xi))
    print(qsged(0.005, mean = 0, sd = 1, nu = nu, xi = xi))
    print(qsged(0.95, mean = 0, sd = 1, nu = nu, xi = xi))
    print(qsged(0.99, mean = 0, sd = 1, nu = nu, xi = xi))
    print(qsged(0.995, mean = 0, sd = 1, nu = nu, xi = xi))
  }
}

# Load data once
returns <- load_data(input_file, price_column)

# EGARCH model with Normal distribution
cat("\n=== EGARCH with Normal distribution ===\n")
garch_fit_norm <- fit_egarch_model(returns, "norm", paste0(output_prefix, ".csv"))

# EGARCH model with Student t distribution
cat("\n=== EGARCH with Student t distribution ===\n")
garch_fit_std <- fit_egarch_model(returns, "std", paste0(output_prefix, "_t.csv"))
cat("\nStudent t quantiles:\n")
calculate_quantiles(garch_fit_std, "std")

# EGARCH model with GED distribution
cat("\n=== EGARCH with GED distribution ===\n")
garch_fit_ged <- fit_egarch_model(returns, "ged", paste0(output_prefix, "_ged.csv"))
cat("\nGED quantiles:\n")
calculate_quantiles(garch_fit_ged, "ged")

# EGARCH model with Skewed Normal distribution
cat("\n=== EGARCH with Skewed Normal distribution ===\n")
garch_fit_snorm <- fit_egarch_model(returns, "snorm", paste0(output_prefix, "_snorm.csv"))
cat("\nSkewed Normal quantiles:\n")
calculate_quantiles(garch_fit_snorm, "snorm")

# EGARCH model with Skewed Student t distribution
cat("\n=== EGARCH with Skewed Student t distribution ===\n")
garch_fit_sstd <- fit_egarch_model(returns, "sstd", paste0(output_prefix, "_st.csv"))
cat("\nSkewed Student t quantiles:\n")
calculate_quantiles(garch_fit_sstd, "sstd")

# EGARCH model with Skewed GED distribution
cat("\n=== EGARCH with Skewed GED distribution ===\n")
garch_fit_sged <- fit_egarch_model(returns, "sged", paste0(output_prefix, "_sged.csv"))
cat("\nSkewed GED quantiles:\n")
calculate_quantiles(garch_fit_sged, "sged")
