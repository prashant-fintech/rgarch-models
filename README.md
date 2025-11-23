# GARCH Models in R

This repository contains implementations of various GARCH (Generalized Autoregressive Conditional Heteroskedasticity) models for volatility modeling and forecasting.

## Models Included

- **GARCH** - Standard GARCH model
- **GJR-GARCH** - Glosten-Jagannathan-Runkle GARCH model (captures asymmetric volatility)
- **EGARCH** - Exponential GARCH model (captures leverage effects)

## Features

Each model is implemented with multiple error distribution assumptions:
- Normal distribution
- Student t distribution
- Generalized Error Distribution (GED)
- Skewed Normal distribution
- Skewed Student t distribution
- Skewed GED distribution

## Requirements

```r
install.packages("rugarch")
install.packages("fGarch")
```

## Usage

1. Place your price data CSV file in the project directory
2. Configure the parameters in the script:
   - `input_file`: Your CSV filename
   - `price_column`: Column index containing price data
   - `output_prefix`: Prefix for output files
3. Run the script

## Output

- Volatility estimates saved as CSV files for each distribution
- One-step ahead volatility forecasts
- Quantile calculations for risk management

## Data Format

Input CSV file should contain:
- Header row
- Price data in the specified column (default: column 2)

## License

MIT
