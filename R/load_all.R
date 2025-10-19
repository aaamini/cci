# Source all helper modules for the VaR backtest.

load_cci_modules <- function(root_dir = NULL) {
  candidates <- character()
  if (!is.null(root_dir)) {
    candidates <- c(candidates, normalizePath(root_dir, winslash = "/", mustWork = FALSE))
  }

  caller <- sys.frame(1)$ofile
  if (!is.null(caller)) {
    caller_path <- normalizePath(caller, winslash = "/", mustWork = FALSE)
    caller_dir <- normalizePath(dirname(caller_path), winslash = "/", mustWork = FALSE)
    candidates <- c(candidates, caller_dir, file.path(caller_dir, ".."))
  }

  candidates <- c(
    candidates,
    normalizePath(getwd(), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE)
  )

  candidates <- unique(candidates)
  resolved_root <- NULL
  for (cand in candidates) {
    if (is.na(cand)) next
    if (!dir.exists(cand)) next
    if (file.exists(file.path(cand, "R", "utils_data.R"))) {
      resolved_root <- normalizePath(cand, winslash = "/", mustWork = TRUE)
      break
    }
  }

  if (is.null(resolved_root)) {
    stop("Unable to locate project root. Provide 'root_dir' when calling load_cci_modules().")
  }

  module_dir <- file.path(resolved_root, "R")
  modules <- c(
    "utils_data.R",
    "data_pipeline.R",
    "cci_core.R",
    "regression_imputation.R",
    "cci_gaussian.R",
    "cci_var_gauss.R",
    "imputation_registry.R"
  )

  for (module in modules) {
    source(file.path(module_dir, module), chdir = TRUE)
  }

  invisible(TRUE)
}

load_cci_modules()
