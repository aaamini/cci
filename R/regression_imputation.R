# Student-t regression based imputation helpers.

fit_t_regression <- function(train_data, min_obs = 50L) {
  if (!requireNamespace("gamlss", quietly = TRUE) ||
      !requireNamespace("gamlss.dist", quietly = TRUE)) {
    stop("Packages 'gamlss' and 'gamlss.dist' are required for regression imputation.")
  }

  is_weekday <- !(weekdays(zoo::index(train_data)) %in% c("Saturday", "Sunday"))
  train_idx <- is_weekday & !is.na(train_data$r_eq) & !is.na(train_data$r_cr)

  df_train <- data.frame(
    r_eq = as.numeric(train_data$r_eq[train_idx]),
    r_cr = as.numeric(train_data$r_cr[train_idx])
  )

  if (nrow(df_train) < min_obs) {
    warning("Not enough weekday data to fit regression model; returning NULL fit.")
    return(NULL)
  }

  t_model <- gamlss::gamlss(
    r_eq ~ .,
    data = df_train,
    family = gamlss.dist::TF,
    trace = FALSE
  )
  list(model = t_model, train_df = df_train)
}

impute_with_t_model <- function(fit, impute_data,
                                stochastic = TRUE,
                                weekend_constraint = FALSE) {
  if (is.null(fit) || is.null(fit$model)) {
    warning("Invalid or NULL fit provided; returning unimputed data.")
    return(impute_data)
  }

  res <- impute_data
  idx <- zoo::index(res)

  is_weekend <- weekdays(idx) %in% c("Saturday", "Sunday")
  impute_idx <- which(is_weekend & is.na(res$r_eq) & !is.na(res$r_cr))

  if (length(impute_idx) == 0) {
    return(res)
  }

  df_predict <- data.frame(r_cr = as.numeric(res$r_cr[impute_idx]))

  if (!stochastic) {
    predicted_r_eq <- stats::predict(fit$model, what = "mu", newdata = df_predict, data = fit$train_df)
  } else {
    predicted_mu <- stats::predict(fit$model, what = "mu", newdata = df_predict, data = fit$train_df)
    predicted_sigma <- exp(stats::predict(fit$model, what = "sigma", newdata = df_predict, data = fit$train_df))
    predicted_nu <- exp(stats::predict(fit$model, what = "nu", newdata = df_predict, data = fit$train_df))
    predicted_r_eq <- gamlss.dist::rTF(
      n = length(predicted_mu),
      mu = predicted_mu,
      sigma = predicted_sigma,
      nu = predicted_nu
    )
  }

  res$r_eq[impute_idx] <- predicted_r_eq

  if (isTRUE(weekend_constraint)) {
    sat_positions <- which(weekdays(idx) == "Saturday")
    if (length(sat_positions) > 0) {
      for (sat_pos in sat_positions) {
        rW <- as.numeric(res$rW[sat_pos])
        if (is.na(rW)) next

        sun_date <- idx[sat_pos] + 1
        sun_pos <- match(sun_date, idx)
        if (is.na(sun_pos)) next

        res$r_eq[sun_pos] <- rW - as.numeric(res$r_eq[sat_pos])
      }
    }
  }

  res
}

impute_with_regression <- function(train_data, impute_data, stochastic = TRUE) {
  fit <- fit_t_regression(train_data)
  if (is.null(fit)) return(impute_data)
  impute_with_t_model(fit, impute_data, stochastic = stochastic, weekend_constraint = FALSE)
}

impute_with_regression_weekend <- function(train_data, impute_data, stochastic = TRUE) {
  fit <- fit_t_regression(train_data)
  if (is.null(fit)) return(impute_data)
  impute_with_t_model(fit, impute_data, stochastic = stochastic, weekend_constraint = TRUE)
}

