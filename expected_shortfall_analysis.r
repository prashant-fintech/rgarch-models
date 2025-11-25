# Expected Shortfall (ES) Analysis for GARCH Models
# Also known as Conditional Value-at-Risk (CVaR)

library(rugarch)
library(fGarch)
library(openxlsx)
library(ggplot2)

# Configuration
input_file <- "c:/projects/rgarch-models/reliance.csv"
price_column <- 2
confidence_levels <- c(0.90, 0.95, 0.99, 0.995)
portfolio_value <- 1000000

# Load data
load_data <- function(file_path, col_index) {
  price <- read.csv(file_path, header = TRUE)
  price_series <- price[, col_index]
  n <- length(price_series)
  returns <- diff(price_series) / lag(price_series)
  returns <- returns[1:(n - 1)]
  return(returns)
}

# Fit GARCH model
fit_garch_model <- function(returns, dist_model, model_type = "sGARCH") {
  spec <- ugarchspec(
    variance.model = list(model = model_type, garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
    distribution.model = dist_model
  )
  garch_fit <- ugarchfit(spec = spec, data = returns)
  return(garch_fit)
}

# Calculate VaR and ES
calculate_risk_measures <- function(garch_fit, dist_type, confidence_level) {
  sigma <- garch_fit@fit$sigma
  resid <- residuals(garch_fit, standardize = TRUE)
  alpha <- 1 - confidence_level
  
  if (dist_type == "norm") {
    var_quantile <- qnorm(alpha, mean = 0, sd = 1)
    es_quantile <- -dnorm(qnorm(alpha)) / alpha
    
  } else if (dist_type == "std") {
    fit_params <- stdFit(resid)
    nu <- fit_params$par[3]
    var_quantile <- qstd(alpha, mean = 0, sd = 1, nu = nu)
    es_quantile <- -dstd(var_quantile, mean = 0, sd = 1, nu = nu) / alpha * 
                   (nu + var_quantile^2) / (nu - 1)
    
  } else if (dist_type == "ged") {
    fit_params <- gedFit(resid)
    nu <- fit_params$par[3]
    var_quantile <- qged(alpha, mean = 0, sd = 1, nu = nu)
    # Numerical integration for ES
    es_quantile <- integrate(function(p) qged(p, mean = 0, sd = 1, nu = nu), 
                            lower = 0, upper = alpha)$value / alpha
    
  } else if (dist_type == "snorm") {
    fit_params <- snormFit(resid)
    xi <- fit_params$par[3]
    var_quantile <- qsnorm(alpha, mean = 0, sd = 1, xi = xi)
    es_quantile <- integrate(function(p) qsnorm(p, mean = 0, sd = 1, xi = xi),
                            lower = 0, upper = alpha)$value / alpha
    
  } else if (dist_type == "sstd") {
    fit_params <- sstdFit(resid)
    nu <- fit_params$estimate[3]
    xi <- fit_params$estimate[4]
    var_quantile <- qsstd(alpha, mean = 0, sd = 1, nu = nu, xi = xi)
    es_quantile <- integrate(function(p) qsstd(p, mean = 0, sd = 1, nu = nu, xi = xi),
                            lower = 0, upper = alpha)$value / alpha
    
  } else if (dist_type == "sged") {
    fit_params <- sgedFit(resid)
    nu <- fit_params$par[3]
    xi <- fit_params$par[4]
    var_quantile <- qsged(alpha, mean = 0, sd = 1, nu = nu, xi = xi)
    es_quantile <- integrate(function(p) qsged(p, mean = 0, sd = 1, nu = nu, xi = xi),
                            lower = 0, upper = alpha)$value / alpha
  }
  
  var <- -sigma * var_quantile
  es <- -sigma * es_quantile
  
  return(data.frame(
    VaR = var,
    ES = es,
    ES_VaR_Ratio = es / var
  ))
}

# Load data
cat("Loading data...\n")
returns <- load_data(input_file, price_column)

# Models and distributions
models <- c("sGARCH", "gjrGARCH", "eGARCH")
model_names <- c("GARCH", "GJR-GARCH", "EGARCH")
dist_types <- c("norm", "std", "ged", "snorm", "sstd", "sged")
dist_names <- c("Normal", "Student-t", "GED", "Skewed Normal", "Skewed Student-t", "Skewed GED")

# Initialize results
all_es_results <- data.frame()

cat("\n=== Calculating Expected Shortfall for all models ===\n")

# Calculate for all combinations
for (m in seq_along(models)) {
  for (d in seq_along(dist_types)) {
    cat(sprintf("\nProcessing %s with %s distribution...\n", model_names[m], dist_names[d]))
    
    tryCatch({
      garch_fit <- fit_garch_model(returns, dist_types[d], models[m])
      
      for (conf in confidence_levels) {
        risk_measures <- calculate_risk_measures(garch_fit, dist_types[d], conf)
        
        # Calculate actual losses
        actual_losses <- -returns * portfolio_value
        var_monetary <- risk_measures$VaR * portfolio_value
        es_monetary <- risk_measures$ES * portfolio_value
        
        # Identify breaches
        breaches <- actual_losses > var_monetary
        actual_tail_loss <- ifelse(sum(breaches) > 0, mean(actual_losses[breaches]), 0)
        
        # ES metrics
        avg_var <- mean(var_monetary)
        avg_es <- mean(es_monetary)
        es_var_ratio <- mean(risk_measures$ES_VaR_Ratio)
        
        # ES accuracy
        es_forecast_error <- actual_tail_loss - avg_es
        es_forecast_error_pct <- ifelse(avg_es > 0, (es_forecast_error / avg_es) * 100, 0)
        
        all_es_results <- rbind(all_es_results, data.frame(
          Model = model_names[m],
          Distribution = dist_names[d],
          Confidence_Level = conf,
          Avg_VaR = avg_var,
          Avg_ES = avg_es,
          ES_VaR_Ratio = es_var_ratio,
          Actual_Tail_Loss = actual_tail_loss,
          ES_Forecast_Error = es_forecast_error,
          ES_Forecast_Error_Pct = es_forecast_error_pct,
          Num_Breaches = sum(breaches),
          stringsAsFactors = FALSE
        ))
      }
    }, error = function(e) {
      cat(sprintf("Error: %s\n", e$message))
    })
  }
}

# Create Excel workbook
cat("\nCreating Excel output...\n")
wb <- createWorkbook()

# Summary sheet
addWorksheet(wb, "ES Analysis")
writeData(wb, "ES Analysis", all_es_results)

# Format header
headerStyle <- createStyle(
  fontColour = "#FFFFFF",
  fgFill = "#4F81BD",
  halign = "center",
  valign = "center",
  textDecoration = "bold",
  border = "TopBottomLeftRight"
)
addStyle(wb, "ES Analysis", headerStyle, rows = 1, cols = 1:ncol(all_es_results), gridExpand = TRUE)
setColWidths(wb, "ES Analysis", cols = 1:ncol(all_es_results), widths = "auto")

# Best ES models
addWorksheet(wb, "Best ES Models")
best_es <- data.frame()
for (conf in confidence_levels) {
  conf_data <- all_es_results[all_es_results$Confidence_Level == conf, ]
  # Rank by absolute forecast error
  conf_data <- conf_data[order(abs(conf_data$ES_Forecast_Error)), ]
  best_es <- rbind(best_es, conf_data[1, ])
}
writeData(wb, "Best ES Models", best_es)
addStyle(wb, "Best ES Models", headerStyle, rows = 1, cols = 1:ncol(best_es), gridExpand = TRUE)
setColWidths(wb, "Best ES Models", cols = 1:ncol(best_es), widths = "auto")

# ES/VaR Ratio Analysis
addWorksheet(wb, "ES-VaR Ratio")
ratio_summary <- aggregate(ES_VaR_Ratio ~ Model + Distribution, data = all_es_results, FUN = mean)
ratio_summary <- ratio_summary[order(-ratio_summary$ES_VaR_Ratio), ]
writeData(wb, "ES-VaR Ratio", ratio_summary)
addStyle(wb, "ES-VaR Ratio", headerStyle, rows = 1, cols = 1:ncol(ratio_summary), gridExpand = TRUE)
setColWidths(wb, "ES-VaR Ratio", cols = 1:ncol(ratio_summary), widths = "auto")

# Save workbook
output_file <- "c:/projects/rgarch-models/Expected_Shortfall_Analysis.xlsx"
saveWorkbook(wb, output_file, overwrite = TRUE)

cat(sprintf("\n\nResults saved to: %s\n", output_file))

# Print summary
cat("\n=== EXPECTED SHORTFALL SUMMARY ===\n")
cat("\nBest ES Models by Confidence Level:\n")
print(best_es[, c("Model", "Distribution", "Confidence_Level", "Avg_ES", 
                  "Actual_Tail_Loss", "ES_Forecast_Error_Pct")])

cat("\n=== KEY INSIGHTS ===\n")
cat(sprintf("1. Average ES/VaR Ratio across all models: %.3f\n", mean(all_es_results$ES_VaR_Ratio)))
cat(sprintf("2. Most accurate ES model (lowest error): %s - %s\n", 
            best_es$Model[1], best_es$Distribution[1]))
cat(sprintf("3. ES typically %.1f%% higher than VaR\n", 
            (mean(all_es_results$ES_VaR_Ratio) - 1) * 100))

cat("\nNote: ES provides a more complete risk measure than VaR by capturing tail risk.\n")
cat("ES Coverage Ratio = 1 indicates perfect tail loss prediction.\n")
