# Utility helpers for working with time–series data and tibbles.

xts_to_tibble <- function(obj_xts) {
  tibble::as_tibble(obj_xts, rownames = "date") |>
    dplyr::mutate(date = as.Date(date))
}

append_weekday <- function(obj_tibble) {
  obj_tibble |>
    dplyr::mutate(wday = weekdays(date))
}

xts_to_tibble_with_wday <- function(obj_xts) {
  obj_xts |>
    xts_to_tibble() |>
    append_weekday()
}

is_weekend <- function(dates) {
  weekdays(dates) %in% c("Saturday", "Sunday")
}

tbl_to_xts <- function(tbl) {
  xts::xts(
    tbl |>
      dplyr::select(-dplyr::any_of(c("date", "wday"))),
    order.by = tbl$date
  )
}

impute_prices_from_returns <- function(tbl, max_iter = 10L) {
  for (i in seq_len(max_iter)) {
    previous <- tbl$p_eq
    tbl <- tbl |>
      dplyr::mutate(p_temp = dplyr::lag(p_eq) * exp(r_eq)) |>
      dplyr::mutate(p_eq = dplyr::coalesce(p_eq, p_temp)) |>
      dplyr::select(-"p_temp")
    if (identical(tbl$p_eq, previous)) break
  }
  tbl
}
