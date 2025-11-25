# Load required libraries
library(rugarch)
library(fGarch)
library(openxlsx)

# Configuration parameters
# Use absolute path or relative path
input_file <- "c:/projects/rgarch-models/reliance.csv"
price_column <- 2
confidence_levels <- c(0.95, 0.99, 0.995)
portfolio_value <- 1000000  # Initial portfolio value

# Print current working directory for debugging
cat(sprintf("Current working directory: %s\n", getwd()))

# Function to load and prepare data
load_data <- function(file_path, col_index) {
  price <- read.csv(file_path, header = TRUE)
  price_series <- price[, col_index]
  n <- length(price_series)
  returns <- diff(price_series) / lag(price_series)
  returns <- returns[1:(n - 1)]
  return(returns)
}

# Function to fit GARCH model
fit_garch_model <- function(returns, dist_model, model_type = "sGARCH") {
  spec <- ugarchspec(
    variance.model = list(model = model_type, garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
    distribution.model = dist_model
  )
  
  garch_fit <- ugarchfit(spec = spec, data = returns)
  return(garch_fit)
}

# Function to calculate VaR and ES for different distributions
calculate_var_es <- function(garch_fit, dist_type, confidence_level) {
  sigma <- garch_fit@fit$sigma
  resid <- residuals(garch_fit, standardize = TRUE)
  alpha <- 1 - confidence_level
  
  # Helper function to calculate ES using numerical integration
  calculate_es_numeric <- function(alpha, pdf_func, quantile_val, params) {
    # ES = E[X | X < VaR] for the distribution
    integrand <- function(p) {
      do.call(pdf_func, c(list(p), params))
    }
    
    if (alpha > 0) {
      es_quantile <- integrate(integrand, lower = 0, upper = alpha)$value / alpha
    } else {
      es_quantile <- quantile_val
    }
    return(es_quantile)
  }
  
  if (dist_type == "norm") {
    quantile_val <- qnorm(alpha, mean = 0, sd = 1)
    # Analytical ES for normal distribution
    es_quantile <- -dnorm(qnorm(alpha)) / alpha
    
  } else if (dist_type == "std") {
    fit_params <- stdFit(resid)
    nu <- fit_params$par[3]
    quantile_val <- qstd(alpha, mean = 0, sd = 1, nu = nu)
    # Analytical ES for Student-t
    es_quantile <- -dstd(quantile_val, mean = 0, sd = 1, nu = nu) / alpha * 
                   (nu + quantile_val^2) / (nu - 1)
    
  } else if (dist_type == "ged") {
    fit_params <- gedFit(resid)
    nu <- fit_params$par[3]
    quantile_val <- qged(alpha, mean = 0, sd = 1, nu = nu)
    es_quantile <- calculate_es_numeric(alpha, qged, quantile_val, 
                                       list(mean = 0, sd = 1, nu = nu))
    
  } else if (dist_type == "snorm") {
    fit_params <- snormFit(resid)
    xi <- fit_params$par[3]
    quantile_val <- qsnorm(alpha, mean = 0, sd = 1, xi = xi)
    es_quantile <- calculate_es_numeric(alpha, qsnorm, quantile_val,
                                       list(mean = 0, sd = 1, xi = xi))
    
  } else if (dist_type == "sstd") {
    fit_params <- sstdFit(resid)
    nu <- fit_params$estimate[3]
    xi <- fit_params$estimate[4]
    quantile_val <- qsstd(alpha, mean = 0, sd = 1, nu = nu, xi = xi)
    es_quantile <- calculate_es_numeric(alpha, qsstd, quantile_val,
                                       list(mean = 0, sd = 1, nu = nu, xi = xi))
    
  } else if (dist_type == "sged") {
    fit_params <- sgedFit(resid)
    nu <- fit_params$par[3]
    xi <- fit_params$par[4]
    quantile_val <- qsged(alpha, mean = 0, sd = 1, nu = nu, xi = xi)
    es_quantile <- calculate_es_numeric(alpha, qsged, quantile_val,
                                       list(mean = 0, sd = 1, nu = nu, xi = xi))
  }
  
  var <- -sigma * quantile_val
  es <- -sigma * es_quantile
  
  return(list(var = var, es = es))
}

# Function to calculate backtesting metrics
calculate_backtesting_metrics <- function(returns, var, es, confidence_level, portfolio_val) {
  # Calculate actual losses
  actual_losses <- -returns * portfolio_val
  
  # Calculate VaR and ES in monetary terms
  var_monetary <- var * portfolio_val
  es_monetary <- es * portfolio_val
  
  # Identify breaches (when actual loss exceeds VaR)
  breaches <- actual_losses > var_monetary
  num_breaches <- sum(breaches)
  
  # Calculate breach rate
  breach_rate <- num_breaches / length(returns)
  
  # Expected breach rate
  expected_breach_rate <- 1 - confidence_level
  
  # Calculate extra capital (average VaR)
  extra_capital <- mean(var_monetary)
  
  # Calculate average ES
  avg_es <- mean(es_monetary)
  
  # Calculate extra losses (losses beyond VaR when breaches occur)
  extra_losses_per_breach <- actual_losses[breaches] - var_monetary[breaches]
  total_extra_losses <- sum(extra_losses_per_breach)
  avg_extra_loss <- ifelse(num_breaches > 0, mean(extra_losses_per_breach), 0)
  
  # ES backtesting - calculate actual average loss in tail
  actual_es <- ifelse(num_breaches > 0, mean(actual_losses[breaches]), 0)
  
  # ES coverage ratio - how well ES predicts tail losses
  es_coverage <- ifelse(avg_es > 0, actual_es / avg_es, 0)
  
  return(list(
    num_breaches = num_breaches,
    breach_rate = breach_rate,
    expected_breach_rate = expected_breach_rate,
    extra_capital = extra_capital,
    avg_es = avg_es,
    total_extra_losses = total_extra_losses,
    avg_extra_loss = avg_extra_loss,
    actual_es = actual_es,
    es_coverage = es_coverage
  ))
}

# Function to perform backtesting for a model
backtest_model <- function(returns, garch_fit, dist_type, model_name, portfolio_val) {
  results_list <- list()
  
  for (conf_level in confidence_levels) {
    risk_measures <- calculate_var_es(garch_fit, dist_type, conf_level)
    metrics <- calculate_backtesting_metrics(returns, risk_measures$var, risk_measures$es, 
                                             conf_level, portfolio_val)
    
    results_list[[as.character(conf_level)]] <- data.frame(
      Model = model_name,
      Distribution = dist_type,
      Confidence_Level = conf_level,
      Num_Breaches = metrics$num_breaches,
      Breach_Rate = metrics$breach_rate,
      Expected_Breach_Rate = metrics$expected_breach_rate,
      Extra_Capital = metrics$extra_capital,
      Avg_ES = metrics$avg_es,
      Total_Extra_Losses = metrics$total_extra_losses,
      Avg_Extra_Loss = metrics$avg_extra_loss,
      Actual_ES = metrics$actual_es,
      ES_Coverage_Ratio = metrics$es_coverage,
      stringsAsFactors = FALSE
    )
  }
  
  return(do.call(rbind, results_list))
}

# Load data
cat(sprintf("Loading data from %s...\n", input_file))
if (!file.exists(input_file)) {
  stop(sprintf("Error: File '%s' not found!", input_file))
}
returns <- load_data(input_file, price_column)
cat(sprintf("Data loaded successfully. Number of returns: %d\n", length(returns)))

# Initialize results dataframe
all_results <- data.frame()

# Distribution types
dist_types <- c("norm", "std", "ged", "snorm", "sstd", "sged")
dist_names <- c("Normal", "Student-t", "GED", "Skewed Normal", "Skewed Student-t", "Skewed GED")

cat("\n=== GARCH Model Backtesting ===\n")

# GARCH models
for (i in seq_along(dist_types)) {
  cat(sprintf("\nProcessing GARCH with %s distribution...\n", dist_names[i]))
  garch_fit <- fit_garch_model(returns, dist_types[i], "sGARCH")
  results <- backtest_model(returns, garch_fit, dist_types[i], "GARCH", portfolio_value)
  all_results <- rbind(all_results, results)
}

cat("\n=== GJR-GARCH Model Backtesting ===\n")

# GJR-GARCH models
for (i in seq_along(dist_types)) {
  cat(sprintf("\nProcessing GJR-GARCH with %s distribution...\n", dist_names[i]))
  garch_fit <- fit_garch_model(returns, dist_types[i], "gjrGARCH")
  results <- backtest_model(returns, garch_fit, dist_types[i], "GJR-GARCH", portfolio_value)
  all_results <- rbind(all_results, results)
}

cat("\n=== EGARCH Model Backtesting ===\n")

# EGARCH models
for (i in seq_along(dist_types)) {
  cat(sprintf("\nProcessing EGARCH with %s distribution...\n", dist_names[i]))
  garch_fit <- fit_garch_model(returns, dist_types[i], "eGARCH")
  results <- backtest_model(returns, garch_fit, dist_types[i], "EGARCH", portfolio_value)
  all_results <- rbind(all_results, results)
}

# Create Excel workbook
wb <- createWorkbook()

# Add summary sheet
addWorksheet(wb, "Summary")
writeData(wb, "Summary", all_results)

# Format the summary sheet
headerStyle <- createStyle(
  fontColour = "#FFFFFF",
  fgFill = "#4F81BD",
  halign = "center",
  valign = "center",
  textDecoration = "bold",
  border = "TopBottomLeftRight"
)

addStyle(wb, "Summary", headerStyle, rows = 1, cols = 1:ncol(all_results), gridExpand = TRUE)

# Add conditional formatting for breach rates
conditionalFormatting(wb, "Summary", 
                     cols = 5, 
                     rows = 2:(nrow(all_results) + 1),
                     rule = ">0.05", 
                     style = createStyle(bgFill = "#FFC7CE", fontColour = "#9C0006"))

conditionalFormatting(wb, "Summary", 
                     cols = 5, 
                     rows = 2:(nrow(all_results) + 1),
                     rule = "<=0.05", 
                     style = createStyle(bgFill = "#C6EFCE", fontColour = "#006100"))

# Set column widths
setColWidths(wb, "Summary", cols = 1:ncol(all_results), widths = "auto")

# Add sheets for each model
for (model in c("GARCH", "GJR-GARCH", "EGARCH")) {
  model_data <- all_results[all_results$Model == model, ]
  addWorksheet(wb, model)
  writeData(wb, model, model_data)
  addStyle(wb, model, headerStyle, rows = 1, cols = 1:ncol(model_data), gridExpand = TRUE)
  setColWidths(wb, model, cols = 1:ncol(model_data), widths = "auto")
}

# Add model ranking analysis
addWorksheet(wb, "Model Ranking")

# Calculate performance scores for each model
all_results$Breach_Deviation <- abs(all_results$Breach_Rate - all_results$Expected_Breach_Rate)
all_results$Efficiency_Ratio <- all_results$Total_Extra_Losses / all_results$Extra_Capital

# Rank models by confidence level
ranking_results <- data.frame()

for (conf in confidence_levels) {
  subset_data <- all_results[all_results$Confidence_Level == conf, ]
  
  # Rank by breach deviation (lower is better)
  subset_data$Breach_Rank <- rank(subset_data$Breach_Deviation)
  
  # Rank by extra capital (lower is better - less capital required)
  subset_data$Capital_Rank <- rank(subset_data$Extra_Capital)
  
  # Rank by total extra losses (lower is better)
  subset_data$Loss_Rank <- rank(subset_data$Total_Extra_Losses)
  
  # Combined score (lower is better)
  subset_data$Overall_Score <- subset_data$Breach_Rank + subset_data$Capital_Rank + subset_data$Loss_Rank
  
  # Sort by overall score
  subset_data <- subset_data[order(subset_data$Overall_Score), ]
  subset_data$Overall_Rank <- 1:nrow(subset_data)
  
  ranking_results <- rbind(ranking_results, subset_data)
}

# Write ranking results
ranking_display <- ranking_results[, c("Overall_Rank", "Model", "Distribution", "Confidence_Level", 
                                       "Num_Breaches", "Breach_Rate", "Breach_Deviation",
                                       "Extra_Capital", "Total_Extra_Losses", "Overall_Score")]
writeData(wb, "Model Ranking", ranking_display)
addStyle(wb, "Model Ranking", headerStyle, rows = 1, cols = 1:ncol(ranking_display), gridExpand = TRUE)
setColWidths(wb, "Model Ranking", cols = 1:ncol(ranking_display), widths = "auto")

# Highlight top 3 models for each confidence level
for (conf in confidence_levels) {
  conf_rows <- which(ranking_display$Confidence_Level == conf & ranking_display$Overall_Rank <= 3)
  if (length(conf_rows) > 0) {
    addStyle(wb, "Model Ranking", 
             createStyle(bgFill = "#D9EAD3", fontColour = "#274E13"),
             rows = conf_rows + 1, cols = 1:ncol(ranking_display), gridExpand = TRUE)
  }
}

# Create Best Models Summary
best_models_summary <- data.frame()
for (conf in confidence_levels) {
  subset_data <- ranking_results[ranking_results$Confidence_Level == conf, ]
  subset_data <- subset_data[order(subset_data$Overall_Score), ]
  best_model <- subset_data[1, ]
  
  best_models_summary <- rbind(best_models_summary, data.frame(
    Confidence_Level = conf,
    Best_Model = paste(best_model$Model, "-", best_model$Distribution),
    Num_Breaches = best_model$Num_Breaches,
    Breach_Rate = best_model$Breach_Rate,
    Expected_Breach_Rate = best_model$Expected_Breach_Rate,
    Breach_Deviation = best_model$Breach_Deviation,
    Extra_Capital = best_model$Extra_Capital,
    Avg_ES = best_model$Avg_ES,
    Total_Extra_Losses = best_model$Total_Extra_Losses,
    Avg_Extra_Loss = best_model$Avg_Extra_Loss,
    Actual_ES = best_model$Actual_ES,
    ES_Coverage_Ratio = best_model$ES_Coverage_Ratio
  ))
}

addWorksheet(wb, "Best Models")
writeData(wb, "Best Models", best_models_summary)
addStyle(wb, "Best Models", headerStyle, rows = 1, cols = 1:ncol(best_models_summary), gridExpand = TRUE)
setColWidths(wb, "Best Models", cols = 1:ncol(best_models_summary), widths = "auto")

# Save workbook
output_file <- "c:/projects/rgarch-models/GARCH_Backtesting_Results.xlsx"
saveWorkbook(wb, output_file, overwrite = TRUE)

cat(sprintf("\n\nBacktesting results saved to: %s\n", output_file))
cat("\n=== BEST MODELS BY CONFIDENCE LEVEL ===\n")
print(best_models_summary)

cat("\n=== TOP 5 MODELS OVERALL (95% Confidence) ===\n")
top5_95 <- ranking_results[ranking_results$Confidence_Level == 0.95, ][1:5, 
          c("Overall_Rank", "Model", "Distribution", "Breach_Deviation", "Extra_Capital", "Total_Extra_Losses")]
print(top5_95)

cat("\n=== TOP 5 MODELS OVERALL (99% Confidence) ===\n")
top5_99 <- ranking_results[ranking_results$Confidence_Level == 0.99, ][1:5, 
          c("Overall_Rank", "Model", "Distribution", "Breach_Deviation", "Extra_Capital", "Total_Extra_Losses")]
print(top5_99)

cat("\n=== TOP 5 MODELS OVERALL (99.5% Confidence) ===\n")
top5_995 <- ranking_results[ranking_results$Confidence_Level == 0.995, ][1:5, 
           c("Overall_Rank", "Model", "Distribution", "Breach_Deviation", "Extra_Capital", "Total_Extra_Losses")]
print(top5_995)

