# Registry for imputation strategies used by the VaR experiment.

imputation_registry <- list(
  "None" = list(
    fit = function(train_data, cfg) NULL,
    impute = function(fit_obj, data, cfg) data
  ),
  "Reg_Det" = list(
    fit = function(train_data, cfg) fit_t_regression(train_data),
    impute = function(fit_obj, data, cfg) {
      if (is.null(fit_obj)) return(data)
      impute_with_t_model(fit_obj, data, stochastic = FALSE, weekend_constraint = FALSE)
    }
  ),
  "Reg_Stoch" = list(
    fit = function(train_data, cfg) fit_t_regression(train_data),
    impute = function(fit_obj, data, cfg) {
      if (is.null(fit_obj)) return(data)
      impute_with_t_model(fit_obj, data, stochastic = TRUE, weekend_constraint = FALSE)
    }
  ),
  "Reg_Det_W" = list(
    fit = function(train_data, cfg) fit_t_regression(train_data),
    impute = function(fit_obj, data, cfg) {
      if (is.null(fit_obj)) return(data)
      impute_with_t_model(fit_obj, data, stochastic = FALSE, weekend_constraint = TRUE)
    }
  ),
  "Reg_Stoch_W" = list(
    fit = function(train_data, cfg) fit_t_regression(train_data),
    impute = function(fit_obj, data, cfg) {
      if (is.null(fit_obj)) return(data)
      impute_with_t_model(fit_obj, data, stochastic = TRUE, weekend_constraint = TRUE)
    }
  ),
  "CCI_Det" = list(
    fit = function(train_data, cfg) cci_fit(train_data),
    impute = function(fit_obj, data, cfg) cci_impute(fit_obj, data, m = cfg$n_impute)
  ),
  "CCI_Stoch" = list(
    fit = function(train_data, cfg) cci_fit(train_data),
    impute = function(fit_obj, data, cfg) cci_impute(fit_obj, data, m = 1)
  ),
  "CCI_4D_Det" = list(
    fit = function(train_data, cfg) cci_fit_joint4_gauss(train_data),
    impute = function(fit_obj, data, cfg) cci_impute_joint4_gauss(fit_obj, data, m = cfg$n_impute, use_rW = "var")
  ),
  "CCI_4D_Stoch" = list(
    fit = function(train_data, cfg) cci_fit_joint4_gauss(train_data),
    impute = function(fit_obj, data, cfg) cci_impute_joint4_gauss(fit_obj, data, m = 1, use_rW = "equal")
  ),
  "CCI_VAR_Det" = list(
    fit = function(train_data, cfg) cci_fit_var_gauss(train_data),
    impute = function(fit_obj, data, cfg) cci_impute_var_gauss(fit_obj, data, m = cfg$n_impute, use_rW = "var")
  ),
  "CCI_VAR_Stoch" = list(
    fit = function(train_data, cfg) cci_fit_var_gauss(train_data),
    impute = function(fit_obj, data, cfg) cci_impute_var_gauss(fit_obj, data, m = 1, use_rW = "equal")
  )
)

fit_imputation <- function(method, train_data, cfg) {
  handler <- imputation_registry[[method]]
  if (is.null(handler)) stop(sprintf("Unknown imputation method: %s", method))
  handler$fit(train_data, cfg)
}

apply_imputation <- function(method, fit_obj, data, cfg) {
  handler <- imputation_registry[[method]]
  if (is.null(handler)) stop(sprintf("Unknown imputation method: %s", method))
  handler$impute(fit_obj, data, cfg)
}

