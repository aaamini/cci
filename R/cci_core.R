# Core Conditional Copula Imputation components.

pobs01 <- function(x) {
  r <- rank(x, ties.method = "average", na.last = "keep")
  n <- sum(!is.na(r))
  (r - 0.5) / n
}

cci_fit <- function(r_xts,
                    familyset = NA,
                    selcrit = "AIC") {
  stopifnot(requireNamespace("VineCopula", quietly = TRUE))

  wk <- weekdays(zoo::index(r_xts))
  idx <- which(
    wk %in% c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday") &
      !is.na(r_xts$r_eq) &
      !is.na(r_xts$r_cr)
  )
  if (length(idx) < 20) {
    stop("Not enough overlapping weekday observations to fit the copula.")
  }

  x1 <- as.numeric(r_xts$r_eq[idx])
  x2 <- as.numeric(r_xts$r_cr[idx])

  u1 <- pobs01(x1)
  u2 <- pobs01(x2)

  fit <- VineCopula::BiCopSelect(
    u1,
    u2,
    familyset = familyset,
    selectioncrit = selcrit,
    indeptest = TRUE
  )

  ecdf1 <- stats::ecdf(x1)
  ecdf2 <- stats::ecdf(x2)
  inv1 <- function(p) stats::quantile(x1, probs = p, type = 8, names = FALSE)

  list(copula = fit, F1 = ecdf1, F2 = ecdf2, F1_inv = inv1)
}

cci_impute <- function(fitobj, r_xts, m = 1000, seed = NULL, eps = 1e-6) {
  stopifnot(all(c("r_eq", "r_cr") %in% colnames(r_xts)))
  if (!is.null(seed)) set.seed(seed)

  res <- r_xts
  idx <- zoo::index(res)

  is_sat <- weekdays(idx) == "Saturday"
  sat_idx <- idx[is_sat & is.na(res$r_eq) & !is.na(res$r_cr)]
  if (length(sat_idx) == 0) {
    warning("No Saturday equity gaps found – nothing to impute.")
    return(res)
  }

  draws <- matrix(NA_real_, nrow = length(sat_idx), ncol = 2 * m)
  colnames(draws) <- c(paste0("sat", seq_len(m)), paste0("sun", seq_len(m)))
  draws <- xts::xts(draws, order.by = sat_idx)

  res$r_eq_sd <- NA_real_

  for (i in seq_along(sat_idx)) {
    sat_date <- sat_idx[i]

    r2 <- as.numeric(res$r_cr[sat_date])
    u2 <- fitobj$F2(r2)
    u2 <- min(max(u2, eps), 1 - eps)

    u1 <- VineCopula::BiCopHinv2(stats::runif(m), rep(u2, m), fitobj$copula)
    r1_draws_sat <- fitobj$F1_inv(u1)

    draws[sat_date, seq_len(m)] <- r1_draws_sat
    res$r_eq[sat_date] <- mean(r1_draws_sat)
    res$r_eq_sd[sat_date] <- stats::sd(r1_draws_sat)

    rW <- as.numeric(res$rW[sat_date])
    if (!is.na(rW)) {
      sun_date <- sat_date + 1
      if (sun_date %in% idx) {
        r1_draws_sun <- rW - r1_draws_sat
        draws[sat_date, m + seq_len(m)] <- r1_draws_sun
        res$r_eq[sun_date] <- mean(r1_draws_sun)
        res$r_eq_sd[sun_date] <- stats::sd(r1_draws_sun)
      }
    }
  }

  attr(res, "cci_draws") <- draws
  res
}

