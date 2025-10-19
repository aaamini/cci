# Conditional Copula Imputation – VaR Backtest

This repository runs the multi-asset Value-at-Risk experiment built around conditional copula based weekend imputation. The primary entry point is `run_var_test.R`, which executes the full workflow across the equity/crypto pairs configured in the script.

## Repository layout
- `run_var_test.R` – main driver that prepares data, runs rolling VaR backtests, and writes results to `output/`.
- `R/` – helper modules split by responsibility (data pipeline, copula imputation, regression imputation, registry wiring).
- `output/` – generated artefacts (RDS files with results, optional LaTeX tables). Safe to delete/regenerate.

## Requirements
1. R 4.2+ (tested with 4.3).
2. Packages: `GAS`, `dplyr`, `tibble`, `VineCopula`, `gamlss`, `gamlss.dist`, `quantmod`, `xts`, `lubridate`, `zoo`, `MASS`.  
   Optional: `mcprogress` (pretty logging) and `knitr` (LaTeX tables).  
   Install with e.g.
   ```r
   install.packages(c(
     "GAS", "dplyr", "tibble", "VineCopula",
     "gamlss", "gamlss.dist", "quantmod", "xts",
     "lubridate", "zoo", "MASS", "mcprogress", "knitr"
   ))
   ```
3. Access to Yahoo Finance via `quantmod::getSymbols()` (requires internet connectivity when the script runs).

## Running the experiment
```r
Rscript run_var_test.R
```
The script automatically sources all helpers from `R/`, downloads the required price data, fits the selected imputation methods, and evaluates the rolling GAS VaR backtest. Outputs are written to `output/`, notably:
- `final_results_df.RDS` – full run-by-run metrics.
- `copulas_<EQUITY>_<CRYPTO>.RDS` – fitted copulas (first run, per pair).
- `var_tables.tex` – optional LaTeX tables if `knitr` is available.

## Customising the run
Edit the configuration block near the top of `scripts/run_var_test.R`:
- `cfg$eq_syms`, `cfg$cr_syms` – symbols to analyse.
- `cfg$intv` and `cfg$tws` – backtest window and training length.
- `cfg$methods` – enabled imputation methods (see registry in `R/imputation_registry.R`).
- `cfg$N_runs` – Monte Carlo repetitions.
- `cfg$impute_scope` – `"train"` (default) to impute training only, or `"both"` to impute both train/test segments.

For quick smoke tests, reduce the number of symbols and runs and shorten the time window before scaling back up.

## Notes
- Intermediate data (price downloads, copulas, final metrics) can be safely removed by deleting the `output/` directory.

