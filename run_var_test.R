# VaR backtest driver for the conditional copula imputation project.
# Run this script from the project root directory.

# --- 1. PROJECT ROOT & MODULE LOADING ---------------------------------------
project_root <- getwd()
if (!file.exists(file.path(project_root, "R", "utils_data.R"))) {
  stop("Please run this script from the project root directory.")
}
source(file.path(project_root, "R", "load_all.R"), chdir = TRUE)

# --- 1.5. OUTPUT DIRECTORY SETUP ---------------------------------------------
output_dir <- file.path(project_root, "output")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# --- 2. DEPENDENCY CHECKS ----------------------------------------------------
ensure_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing required packages: ",
      paste(missing, collapse = ", "),
      ". Please install them before running the script."
    )
  }
}

ensure_packages(c(
  "GAS", "dplyr", "tibble", "VineCopula", "gamlss", "gamlss.dist",
  "parallel", "quantmod", "xts", "lubridate", "zoo", "MASS"
))
library(GAS)
library(dplyr)
library(tibble)
library(VineCopula)
library(gamlss)
library(parallel)
library(mcprogress)
library(xts)

cfg <- list(
  # eq_syms = c("^IXIC", "GC=F", "SPY", "XLK", "XLV", "XLE", "XLF"),
  eq_syms = c("SPY"),
  # eq_syms = c("GC=F"),
  # cr_syms = c("BTC-USD", "ETH-USD"),
  cr_syms = c("BTC-USD"),
  # cr_syms = c("ETH-USD"),
  intv = as.Date(c("2020-02-28", "2025-10-01")),
  tws = 365 * 3,
  refit_frequency = 60,
  N_runs = 1,
  portfolio_type = "EW",
  methods =  c("None", 
    "CCI_Det", "CCI_Stoch", 
    "CCI_4D_Det", "CCI_4D_Stoch", 
    "CCI_VAR_Det", "CCI_VAR_Stoch", 
    "Reg_Det_W", "Reg_Det"),
  alpha = 0.05,
  n_impute = 1000,
  impute_scope = "train" # "train" or "both" --- default is "train"; when selecting "both", GAS is rolled over the imputed test data, but VaR forecasts are only evaluated on Tuesday, Wednesday, Thursday, and Friday which have real realized returns
)

num_cores <- detectCores() - 1
cat("--- Starting parallel processing on", num_cores, "cores ---\n")
n_methods <- length(cfg$methods)

# --- 3. HELPER FUNCTIONS ---
backtest_var <- function(y_test, VaR, alpha = 0.05) {
  pinball_loss <- (y_test - VaR) * (alpha - (y_test < VaR))
  BT <- BacktestVaR(y_test, VaR, alpha)
  data.frame(
    pinball_loss = mean(pinball_loss),
    hit_rate = mean(as.numeric(y_test) <= VaR),
    n = NROW(y_test),
    LRuc = as.numeric(BT$LRuc[1]), LRuc_p = as.numeric(BT$LRuc[2]),
    LRcc = as.numeric(BT$LRcc[1]), LRcc_p = as.numeric(BT$LRcc[2]),
    DQ = as.numeric(BT$DQ[1]), DQ_p = as.numeric(BT$DQ[2])
  )
}

# Helper function (from your code)
tbl_to_xts <- function(tbl) {
  xts(tbl %>% select(-date, -wday), order.by = tbl$date)
}

# Function to homogenize training data
homogenize_imputed_data <- function(imputed_data) {
  imputed_data %>%
    xts_to_tibble_with_wday() %>%
    impute_prices_from_returns() %>%
    mutate(r_eq = log(p_eq / dplyr::lag(p_eq))) %>%
    tbl_to_xts()
}

# --- 4. SYMBOL GRID ---
symbol_grid <- expand.grid(
  eq_sym = cfg$eq_syms,
  cr_sym = cfg$cr_syms,
  stringsAsFactors = FALSE
)

all_experiment_results <- vector("list", nrow(symbol_grid))

# --- 5. MAIN LOOP OVER SYMBOL COMBINATIONS ---
# for (grid_idx in seq_len(nrow(symbol_grid))) {
all_experiment_results <- mclapply(seq_len(nrow(symbol_grid)), function(grid_idx) {
  eq_sym <- symbol_grid$eq_sym[grid_idx]
  cr_sym <- symbol_grid$cr_sym[grid_idx]
  cat("\n====================================================\n")
  cat("====== Starting experiments for:", eq_sym, "x", cr_sym, "======\n")
  cat("====================================================\n")

  master_data <- get_prices_and_returns(
    eq_sym,
    cr_sym,
    from = cfg$intv[1],
    to = cfg$intv[2]
  )

  run_results <- mclapply(1:cfg$N_runs, function(run_i) {
  # run_results <- lapply(1:cfg$N_runs, function(run_i) {
    cat("\n----------------------------------------------------\n")
    cat("Run", run_i, "/", cfg$N_runs, "for", eq_sym, "x", cr_sym, "\n")
    cat("----------------------------------------------------\n\n")

    # Collect fitted CCI copulas across refits (first run only) as a tibble with list-column
    copula_records_run <- if (run_i == 1) tibble::tibble(
      eq_sym = character(),
      cr_sym = character(),
      refit_date = as.Date(character()),
      train_start = as.Date(character()),
      train_end = as.Date(character()),
      family = character(),
      par = numeric(),
      par2 = numeric(),
      tau = numeric(),
      copula = list()
    ) else NULL

    test_start <- cfg$intv[1] + cfg$tws
    test_end <- cfg$intv[2]
    refit_dates <- seq.Date(from = test_start, to = test_end, by = paste(cfg$refit_frequency, "days"))

    forecast_dates <- seq(test_start, test_end, by = "day")
    all_VaR_forecasts <- xts::xts(
      matrix(NA_real_, nrow = length(forecast_dates), ncol = length(cfg$methods)),
      order.by = as.Date(forecast_dates)
    )
    colnames(all_VaR_forecasts) <- cfg$methods

    spec <- UniGASSpec(
      Dist = "std",
      ScalingType = "Identity",
      GASPar = list(location = FALSE, scale = TRUE, shape = FALSE)
    )

    for (refit_date in refit_dates) {
      refit_date <- as.Date(refit_date)
      train_end <- refit_date - 1
      train_start <- train_end - cfg$tws + 1
      test_chunk_start <- refit_date
      test_chunk_end <- min(refit_date + cfg$refit_frequency - 1, test_end)

      cat_parallel(
        "Run ", run_i,
        " || Symbols: ", eq_sym, "/", cr_sym,
        " || Refit: ", as.character(train_start), " to ", as.character(train_end), "\n"
      )

      train_data <- window(master_data, start = train_start, end = train_end)
      test_chunk_data <- window(master_data, start = test_chunk_start, end = test_chunk_end)

      # Small fit cache to avoid redundant fits across methods within the same refit window
      fit_cache <- new.env(parent = emptyenv())
      get_fit_key <- function(method) {
        if (grepl("^CCI_4D", method)) return("CCI_4D")
        if (grepl("^CCI_VAR", method)) return("CCI_VAR")
        if (grepl("^CCI_", method)) return("CCI")
        if (grepl("^Reg_", method)) return("Reg")
        return(method)
      }
      get_fit_obj <- function(method) {
        key <- get_fit_key(method)
        if (exists(key, envir = fit_cache, inherits = FALSE)) return(get(key, envir = fit_cache))
        fo <- fit_imputation(method, train_data, cfg)
        assign(key, fo, envir = fit_cache)
        fo
      }

      # Flag to ensure we record the CCI copula once per refit window
      cci_copula_recorded <- FALSE

      for (method in cfg$methods) {
        fit_obj <- get_fit_obj(method) # This is where we fit the imputation model

        # ---- The whole logic here is for recording the fitted copula 
        # Record fitted copula for 2D CCI methods (once per refit window) on first run only
        if (run_i == 1 && !cci_copula_recorded && grepl("^CCI_", method) && identical(get_fit_key(method), "CCI")) {
          cop <- NULL
          if (!is.null(fit_obj) && !is.null(fit_obj$copula)) {
            cop <- fit_obj$copula
          }
          if (!is.null(cop)) {
            copula_records_run <- dplyr::bind_rows(
              copula_records_run,
              tibble::tibble(
                eq_sym = eq_sym,
                cr_sym = cr_sym,
                refit_date = refit_date,
                train_start = train_start,
                train_end = train_end,
                family = VineCopula::BiCopName(cop$family, short = FALSE),
                par = as.numeric(cop$par),
                par2 = as.numeric(cop$par2),
                tau = as.numeric(cop$tau),
                copula = list(cop)
              )
            )
            cci_copula_recorded <- TRUE
          }
        }
        # ---- End of the logic for recording the fitted copula ----

        # Apply the imputation method to the training data
        train_imputed <- apply_imputation(method, fit_obj, train_data, cfg)

        train_homogenized <- homogenize_imputed_data(train_imputed)

        y_train <- na.omit(get_portfolio_returns(train_homogenized, cfg$portfolio_type))
        fit <- UniGASFit(spec, y_train)
        
        test_used <- if (identical(cfg$impute_scope, "both")) homogenize_imputed_data(apply_imputation(method, fit_obj, test_chunk_data, cfg)) else test_chunk_data
        y_test_chunk <- na.omit(get_portfolio_returns(test_used, cfg$portfolio_type))
      
        frc <- UniGASFor(fit, Roll = TRUE, out = y_test_chunk)
        VaR_chunk <- quantile(frc, cfg$alpha)
        all_VaR_forecasts[zoo::index(y_test_chunk), method] <- VaR_chunk
      }
    }

    # Save collected CCI copulas only once (on first run) to a sidecar RDS
    if (run_i == 1) {
      safe_eq <- gsub("[^A-Za-z0-9]+", "_", eq_sym)
      safe_cr <- gsub("[^A-Za-z0-9]+", "_", cr_sym)
      copula_outfile <- file.path(project_root, "output", sprintf("copulas_%s_%s.RDS", safe_eq, safe_cr))
      saveRDS(copula_records_run, file = copula_outfile)
    }

    actual_returns <- get_portfolio_returns(window(master_data, start = test_start, end = test_end))
    colnames(actual_returns) <- "realized"
    comparison_data <- na.omit(merge(all_VaR_forecasts, actual_returns))
    
    if (identical(cfg$impute_scope, "both")) {
      # Only keep Tuesday, Wednesday, Thursday, and Friday which have real realized returns
      is_tue_fri <- weekdays(zoo::index(comparison_data)) %in% c("Tuesday", "Wednesday", "Thursday", "Friday")
      comparison_data <- comparison_data[is_tue_fri]
    }

    run_results_list <- vector("list", length(cfg$methods))
    for (method in cfg$methods) {
      res <- backtest_var(comparison_data$realized, comparison_data[, method], cfg$alpha)
      run_results_list[[method]] <- cbind(
        data.frame(
          eq_sym = eq_sym,
          cr_sym = cr_sym,
          run = run_i,
          Method = method,
          stringsAsFactors = FALSE
        ),
        res
      )
    }

    dplyr::bind_rows(run_results_list)
  }, mc.cores = num_cores)
  # })
  # all_experiment_results[[grid_idx]] <- dplyr::bind_rows(run_results)
  return(dplyr::bind_rows(run_results))
}, mc.cores = 1)

# --- 6. FINAL AGGREGATION ---
final_results_df <- dplyr::bind_rows(all_experiment_results)
saveRDS(final_results_df, file = file.path(project_root, "output", "final_results_df.RDS"))


summary_stats <- final_results_df %>%
  group_by(eq_sym, cr_sym, Method) %>%
  summarise(
    across(
      .cols = where(is.numeric),
      .fns = list(mean = mean, sd = sd),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  select(
    eq_sym,
    cr_sym,
    Method,
    n_mean,
    pinball_loss_mean, pinball_loss_sd,
    hit_rate_mean, hit_rate_sd,
    LRuc_mean, LRuc_sd, LRuc_p_mean, LRuc_p_sd,
    LRcc_mean, LRcc_sd, LRcc_p_mean, LRcc_p_sd,
    DQ_mean, DQ_sd, DQ_p_mean, DQ_p_sd
  )

cat("\n\n--- Final Aggregated Backtesting Results across all symbol pairs ---\n")
print(as.data.frame(summary_stats))

# --- 7. FORMATTED REPORTING TABLE ---
# method_display_order <- c("None", "Simple", "CCI", "CCI_Det", "Regression_Det", "Regression_Stochastic")
method_display_order <- cfg$methods

mean_results <- final_results_df %>%
  group_by(eq_sym, cr_sym, Method) %>%
  summarise(
    hit_rate = mean(hit_rate),
    LRuc = mean(LRuc),
    LRuc_p = mean(LRuc_p),
    LRcc = mean(LRcc),
    LRcc_p = mean(LRcc_p),
    DQ = mean(DQ),
    DQ_p = mean(DQ_p),
    .groups = "drop"
  ) %>%
  mutate(Method = factor(Method, levels = method_display_order)) %>%
  arrange(cr_sym, eq_sym, Method)

format_stat_with_p <- function(stat, pval) {
  if (is.na(stat) || is.na(pval)) {
    return(NA_character_)
  }
  sprintf("%-7s (%s)", signif(stat, 3), sprintf("%.2f", round(pval, 2)))
}

format_hit_rate <- function(x) {
  if (is.na(x)) {
    return(NA_character_)
  }
  as.character(signif(x * 100, 3))
}

formatted_results <- mean_results %>%
  rowwise() %>%
  mutate(
    HitRate = format_hit_rate(hit_rate),
    LRuc_fmt = format_stat_with_p(LRuc, LRuc_p),
    LRcc_fmt = format_stat_with_p(LRcc, LRcc_p),
    DQ_fmt = format_stat_with_p(DQ, DQ_p)
  ) %>%
  select(
    cr_sym,
    Asset = eq_sym,
    Method,
    `Hit Rate (%)` = HitRate,
    LRuc = LRuc_fmt,
    LRcc = LRcc_fmt,
    DQ = DQ_fmt
  )

cat("\n\n--- Formatted reporting tables ---\n")
for (crypto in unique(formatted_results$cr_sym)) {
  cat("\n", crypto, "results:\n", sep = "")
  crypto_table <- formatted_results %>%
    filter(cr_sym == crypto) %>%
    select(-cr_sym)
  print(as.data.frame(crypto_table), row.names = FALSE)
}

if (requireNamespace("knitr", quietly = TRUE)) {
  split_tables <- split(formatted_results, formatted_results$cr_sym)
  latex_tables <- lapply(
    split_tables,
    function(df) {
      knitr::kable(
        df %>% select(-cr_sym),
        format = "latex",
        booktabs = TRUE,
        align = "llrrrr",
        linesep = c(rep("", n_methods - 1), "\\addlinespace"),
        caption = paste0("VaR backtesting results for ", unique(df$cr_sym), " portfolio.")
      )
    }
  )
  names(latex_tables) <- names(split_tables)
  latex_strings <- lapply(
    latex_tables,
    function(tbl) paste(tbl, collapse = "\n")
  )
  latex_output_path <- file.path(project_root, "output", "var_tables.tex")
  latex_body <- paste(latex_strings, collapse = "\n\n")
  latex_document <- paste0(
    "\\documentclass{article}\n",
    "\\usepackage[margin=1in]{geometry}\n",
    "\\usepackage{booktabs}\n",
    "\\usepackage{multirow}\n",
    "\\usepackage{caption}\n",
    "\\begin{document}\n",
    latex_body,
    "\n\\end{document}\n"
  )
  cat(latex_document, file = latex_output_path)
  cat("\nSaved LaTeX tables to", latex_output_path, "\n")
} else {
  latex_tables <- NULL
  cat("\n(knitr not available: skipping LaTeX table generation)\n")
}

