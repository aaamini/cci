# Data acquisition and transformation helpers.

get_master_data <- function(eq_sym, cr_sym,
                            from = Sys.Date() - 365 * 3,
                            to = Sys.Date(),
                            tz = "UTC") {
  stopifnot(
    requireNamespace("quantmod", quietly = TRUE),
    requireNamespace("xts", quietly = TRUE),
    requireNamespace("lubridate", quietly = TRUE)
  )

  eq <- quantmod::getSymbols(
    eq_sym,
    from = from,
    to = to,
    src = "yahoo",
    auto.assign = FALSE,
    tz = tz
  )
  cr <- quantmod::getSymbols(
    cr_sym,
    from = from,
    to = to,
    src = "yahoo",
    auto.assign = FALSE,
    tz = tz
  )

  full_dates <- seq.Date(as.Date(from), as.Date(to), by = "day")
  master_xts <- xts::xts(order.by = full_dates)
  master_xts <- merge(master_xts, quantmod::Cl(eq), quantmod::Cl(cr), quantmod::Op(eq), all = TRUE)
  colnames(master_xts) <- c("p_eq", "p_cr", "p_eq_open")

  master_tbl <- xts_to_tibble_with_wday(master_xts) |>
    dplyr::mutate(r_eq = log(p_eq / dplyr::lag(p_eq))) |>
    dplyr::mutate(r_cr = log(p_cr / dplyr::lag(p_cr))) |>
    dplyr::mutate(
      next_open_bfill = dplyr::lead(
        rev(zoo::na.locf(rev(p_eq_open), na.rm = FALSE))
      ),
      rW = dplyr::if_else(
        wday == "Saturday",
        log(dplyr::lag(next_open_bfill) / dplyr::lag(p_eq)),
        NA_real_
      )
    ) |>
    dplyr::select(-"next_open_bfill")

  xts::xts(master_tbl |> dplyr::select(-"date", -"wday"), order.by = master_tbl$date)
}

get_returns <- function(eq_sym, cr_sym,
                        from = Sys.Date() - 365 * 3,
                        to = Sys.Date(),
                        tz = "UTC") {
  stopifnot(
    requireNamespace("quantmod", quietly = TRUE),
    requireNamespace("xts", quietly = TRUE),
    requireNamespace("lubridate", quietly = TRUE)
  )

  eq <- quantmod::getSymbols(eq_sym, from = from, to = to, src = "yahoo", auto.assign = FALSE, tz = tz)
  cr <- quantmod::getSymbols(cr_sym, from = from, to = to, src = "yahoo", auto.assign = FALSE, tz = tz)

  eq_close <- quantmod::Cl(eq)
  cr_close <- quantmod::Cl(cr)

  r_eq <- quantmod::dailyReturn(eq_close, type = "log")
  r_cr <- quantmod::dailyReturn(cr_close, type = "log")
  colnames(r_eq) <- "r_eq"
  colnames(r_cr) <- "r_cr"

  eq_open <- quantmod::Op(eq)
  fridays <- zoo::index(eq_close)[weekdays(zoo::index(eq_close)) == "Friday"]
  rW_vec <- rep(NA_real_, length(fridays))
  for (i in seq_along(fridays)) {
    fri <- fridays[i]
    mon_idx <- which(zoo::index(eq_open) > fri)[1]
    if (!is.na(mon_idx)) {
      p_fri <- as.numeric(eq_close[fri])
      p_mon <- as.numeric(eq_open[mon_idx])
      if (!anyNA(c(p_fri, p_mon))) {
        rW_vec[i] <- log(p_mon / p_fri)
      }
    }
  }
  rW_xts <- xts::xts(rW_vec, order.by = fridays + 1)
  colnames(rW_xts) <- "rW"

  Reduce(function(d1, d2) merge(d1, d2, join = "outer"), list(r_eq, r_cr, rW_xts))
}

get_prices_and_returns <- function(eq_sym, cr_sym,
                                   from = Sys.Date() - 365 * 3,
                                   to = Sys.Date(),
                                   tz = "UTC") {
  stopifnot(requireNamespace("quantmod", quietly = TRUE))

  eq <- quantmod::getSymbols(eq_sym, from = from, to = to, src = "yahoo", auto.assign = FALSE, tz = tz)
  cr <- quantmod::getSymbols(cr_sym, from = from, to = to, src = "yahoo", auto.assign = FALSE, tz = tz)

  eq_close <- quantmod::Cl(eq)
  cr_close <- quantmod::Cl(cr)
  colnames(eq_close) <- "p_eq"
  colnames(cr_close) <- "p_cr"

  r_eq <- quantmod::dailyReturn(eq_close, type = "log")
  r_cr <- quantmod::dailyReturn(cr_close, type = "log")
  colnames(r_eq) <- "r_eq"
  colnames(r_cr) <- "r_cr"

  eq_open <- quantmod::Op(eq)
  fridays <- zoo::index(eq_close)[weekdays(zoo::index(eq_close)) == "Friday"]
  rW_vec <- rep(NA_real_, length(fridays))
  for (i in seq_along(fridays)) {
    fri <- fridays[i]
    mon_idx <- which(zoo::index(eq_open) > fri)[1]
    if (!is.na(mon_idx)) {
      p_fri <- as.numeric(eq_close[fri])
      p_mon <- as.numeric(eq_open[mon_idx])
      if (!anyNA(c(p_fri, p_mon))) {
        rW_vec[i] <- log(p_mon / p_fri)
      }
    }
  }
  rW <- xts::xts(rW_vec, order.by = fridays + 1)
  colnames(rW) <- "rW"

  merge(r_eq, r_cr, rW, eq_close, cr_close, all = TRUE)
}

get_portfolio_returns <- function(all_data, type = "EW") {
  type <- toupper(type)
  if (type == "EW") {
    simple_returns <- exp(all_data[, c("r_eq", "r_cr")]) - 1
    pf_simple_returns <- 0.5 * simple_returns$r_eq + 0.5 * simple_returns$r_cr
    return(pf_simple_returns)
  }

  if (type == "EWBH") {
    stopifnot(requireNamespace("quantmod", quietly = TRUE))

    prices <- all_data[, c("p_eq", "p_cr")]
    first_day_prices <- na.omit(prices)[1, ]

    initial_investment <- 1000
    p_eq_start <- as.numeric(first_day_prices$p_eq)
    p_cr_start <- as.numeric(first_day_prices$p_cr)

    shares_eq <- (0.5 * initial_investment) / p_eq_start
    shares_cr <- (0.5 * initial_investment) / p_cr_start

    portfolio_value <- (shares_eq * prices$p_eq) + (shares_cr * prices$p_cr)

    return(quantmod::dailyReturn(portfolio_value, type = "log"))
  }

  stop("Only EW or EWBH portfolio types are supported.")
}

impute_weekend_prices_from_returns <- function(xts_data) {
  stopifnot("p_eq" %in% colnames(xts_data))
  result <- xts_data
  all_dates <- zoo::index(result)
  weekdays_vec <- weekdays(all_dates)
  n_total <- length(all_dates)

  sat_indices <- which(weekdays_vec == "Saturday")
  for (sat_idx in sat_indices) {
    if (!is.na(result$p_eq[sat_idx])) next
    if (sat_idx - 1 < 1) next
    if (weekdays_vec[sat_idx - 1] != "Friday") next

    p_fri <- as.numeric(result$p_eq[sat_idx - 1])
    r_sat <- as.numeric(result$r_eq[sat_idx])
    p_sat <- p_fri * exp(r_sat)
    result$p_eq[sat_idx] <- p_sat

    if (sat_idx == n_total) next
    if (weekdays_vec[sat_idx + 1] != "Sunday") next

    r_sun <- as.numeric(result$r_eq[sat_idx + 1])
    result$p_eq[sat_idx + 1] <- p_sat * exp(r_sun)
  }

  result
}
