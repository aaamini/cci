# Gaussian copula extensions for joint weekend imputation.

to_u <- function(x, Fhat) {
  u <- Fhat(x)
  pmin(pmax(u, 1e-6), 1 - 1e-6)
}

to_z <- function(u) stats::qnorm(u)

from_z_to_x <- function(z, Finv) {
  u <- stats::pnorm(z)
  Finv(u)
}

rmvnorm2_chol <- function(n, mu, Sigma, jitter = 1e-12) {
  R <- try(stats::chol(Sigma), silent = TRUE)
  if (inherits(R, "try-error")) {
    Sigma <- Sigma + diag(jitter, 2L)
    R <- stats::chol(Sigma)
  }
  Z <- matrix(stats::rnorm(n * 2L), nrow = n, ncol = 2L)
  sweep(Z %*% R, 2, mu, `+`)
}

cci_fit_joint4_gauss <- function(r_xts) {
  stopifnot(all(c("r_eq", "r_cr") %in% colnames(r_xts)))

  idx <- zoo::index(r_xts)
  wk <- weekdays(idx)
  wday_ok <- wk %in% c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday")
  ok <- which(wday_ok & !is.na(r_xts$r_eq) & !is.na(r_xts$r_cr))
  ok <- ok[ok < length(idx)]
  ok <- ok[(ok + 1) %in% ok]

  if (length(ok) < 50) stop("Not enough weekday pairs to fit 4D model.")

  x_t <- as.numeric(r_xts$r_eq[ok])
  y_t <- as.numeric(r_xts$r_cr[ok])
  x_tp1 <- as.numeric(r_xts$r_eq[ok + 1])
  y_tp1 <- as.numeric(r_xts$r_cr[ok + 1])

  F1 <- stats::ecdf(as.numeric(r_xts$r_eq[wday_ok & !is.na(r_xts$r_eq)]))
  F2 <- stats::ecdf(as.numeric(r_xts$r_cr[!is.na(r_xts$r_cr)]))
  Finv1 <- function(p) stats::quantile(
    as.numeric(r_xts$r_eq[wday_ok & !is.na(r_xts$r_eq)]),
    probs = p,
    type = 8,
    names = FALSE
  )

  U <- cbind(to_u(x_t, F1), to_u(y_t, F2), to_u(x_tp1, F1), to_u(y_tp1, F2))
  Z <- stats::qnorm(U)
  mu <- colMeans(Z)
  S <- stats::cov(Z)

  list(mu = mu, Sigma = S, F1 = F1, F2 = F2, F1_inv = Finv1)
}

cci_impute_joint4_gauss <- function(fit4, r_xts, m = 1000,
                                    use_rW = c("none", "equal", "var")) {
  use_rW <- match.arg(use_rW)
  stopifnot(all(c("r_eq", "r_cr") %in% colnames(r_xts)))

  res <- r_xts
  idx <- zoo::index(res)
  is_sat <- weekdays(idx) == "Saturday"

  draws <- matrix(
    NA_real_,
    nrow = NROW(res),
    ncol = 2 * m,
    dimnames = list(
      as.character(idx),
      c(paste0("sat_", seq_len(m)), paste0("sun_", seq_len(m)))
    )
  )

  ix_x <- c(1L, 3L)
  ix_y <- c(2L, 4L)
  mu <- fit4$mu
  S <- fit4$Sigma
  Sxx <- S[ix_x, ix_x, drop = FALSE]
  Syy <- S[ix_y, ix_y, drop = FALSE]
  Sxy <- S[ix_x, ix_y, drop = FALSE]
  Syx <- t(Sxy)
  iSyy <- solve(Syy)

  res$r_eq_sd <- NA_real_

  for (k in which(is_sat)) {
    d_sat <- idx[k]
    d_sun <- d_sat + 1

    if (is.na(res$r_cr[k]) || !(d_sun %in% idx) || is.na(res$r_cr[d_sun])) next

    u_y1 <- to_u(as.numeric(res$r_cr[k]), fit4$F2)
    u_y2 <- to_u(as.numeric(res$r_cr[d_sun]), fit4$F2)
    z_y <- stats::qnorm(c(u_y1, u_y2))
    mu_y <- fit4$mu[ix_y]
    mu_x <- fit4$mu[ix_x]

    c_mean <- mu_x + Sxy %*% iSyy %*% (z_y - mu_y)
    c_S <- Sxx - Sxy %*% iSyy %*% Syx

    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Package 'MASS' is required for Gaussian draw generation.")
    }
    z_x_draws <- MASS::mvrnorm(n = m, mu = as.numeric(c_mean), Sigma = c_S)
    if (m == 1) {
      z_x_draws <- matrix(z_x_draws, nrow = 1)
    }

    r_sat_draws <- fit4$F1_inv(stats::pnorm(z_x_draws[, 1]))
    r_sun_draws <- fit4$F1_inv(stats::pnorm(z_x_draws[, 2]))

    if (use_rW != "none" && !is.na(res$rW[k])) {
      delta <- as.numeric(res$rW[k]) - (r_sat_draws + r_sun_draws)
      if (use_rW == "var") {
        if (m == 1) {
          z_tmp <- MASS::mvrnorm(n = 512, mu = as.numeric(c_mean), Sigma = c_S)
          r_sat_tmp <- fit4$F1_inv(stats::pnorm(z_tmp[, 1]))
          r_sun_tmp <- fit4$F1_inv(stats::pnorm(z_tmp[, 2]))
          v1 <- stats::var(r_sat_tmp)
          v2 <- stats::var(r_sun_tmp)
        } else {
          v1 <- stats::var(r_sat_draws)
          v2 <- stats::var(r_sun_draws)
        }
        w1 <- v1 / (v1 + v2 + 1e-12)
        w2 <- 1 - w1
      } else {
        w1 <- 0.5
        w2 <- 0.5
      }
      r_sat_draws <- r_sat_draws + w1 * delta
      r_sun_draws <- r_sun_draws + w2 * delta
    }

    res$r_eq[d_sat] <- mean(r_sat_draws)
    res$r_eq[d_sun] <- mean(r_sun_draws)
    res$r_eq_sd[d_sat] <- stats::sd(r_sat_draws)
    res$r_eq_sd[d_sun] <- stats::sd(r_sun_draws)

    draws[as.character(d_sat), seq_len(m)] <- r_sat_draws
    draws[as.character(d_sat), m + seq_len(m)] <- r_sun_draws
  }

  attr(res, "cci_draws") <- draws
  res
}

