# VAR(1) Gaussian copula based weekend imputation.

cci_fit_var_gauss <- function(r_xts) {
  stopifnot(all(c("r_eq", "r_cr") %in% colnames(r_xts)))

  idx <- zoo::index(r_xts)
  is_wkday <- weekdays(idx) %in% c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday")
  eq <- as.numeric(r_xts$r_eq)
  cr <- as.numeric(r_xts$r_cr)

  F1 <- stats::ecdf(eq[is_wkday & !is.na(eq)])
  F2 <- stats::ecdf(cr[!is.na(cr)])
  Finv1 <- function(p) stats::quantile(
    eq[is_wkday & !is.na(eq)],
    probs = p,
    type = 8,
    names = FALSE
  )

  U1 <- to_u(eq, F1)
  U2 <- to_u(cr, F2)
  Z1 <- stats::qnorm(U1)
  Z2 <- stats::qnorm(U2)

  ok <- which(is_wkday & !is.na(eq) & !is.na(cr))
  keep <- ok[(ok - 1L) %in% ok]
  if (length(keep) < 20) stop("Not enough weekday pairs for VAR(1).")

  Y <- cbind(Z1[keep], Z2[keep])
  X <- cbind(1, Z1[keep - 1L], Z2[keep - 1L])

  Bhat <- solve(t(X) %*% X, t(X) %*% Y)
  b <- as.numeric(Bhat[1, ])
  A <- t(Bhat[2:3, , drop = FALSE])

  E <- Y - X %*% Bhat
  S <- stats::cov(E)

  list(b = b, A = A, Sigma = S, F1 = F1, F2 = F2, F1_inv = Finv1)
}

cci_impute_var_gauss <- function(fitv, r_xts, m = 1000,
                                 use_rW = c("none", "equal", "var")) {
  use_rW <- match.arg(use_rW)
  res <- r_xts
  idx <- zoo::index(res)
  sat_idx <- which(weekdays(idx) == "Saturday")

  draws <- matrix(
    NA_real_,
    nrow = NROW(res),
    ncol = 2 * m,
    dimnames = list(
      as.character(idx),
      c(paste0("sat_", seq_len(m)), paste0("sun_", seq_len(m)))
    )
  )

  b <- fitv$b
  A <- fitv$A
  S <- fitv$Sigma

  eq <- r_xts$r_eq
  cr <- r_xts$r_cr

  res$r_eq_sd <- NA_real_

  for (k in sat_idx) {
    d_sat <- idx[k]
    d_fri <- d_sat - 1
    d_sun <- d_sat + 1

    if (!(d_fri %in% idx) || !(d_sun %in% idx)) next
    if (is.na(eq[d_fri]) || is.na(cr[d_fri])) next
    if (is.na(cr[d_sat]) || is.na(cr[d_sun])) next

    Z0 <- c(
      stats::qnorm(to_u(as.numeric(eq[d_fri]), fitv$F1)),
      stats::qnorm(to_u(as.numeric(cr[d_fri]), fitv$F2))
    )
    Zy1 <- stats::qnorm(to_u(as.numeric(cr[d_sat]), fitv$F2))
    Zy2 <- stats::qnorm(to_u(as.numeric(cr[d_sun]), fitv$F2))

    mu1 <- b + A %*% Z0
    mu2 <- b + A %*% mu1
    V11 <- S
    V22 <- A %*% S %*% t(A) + S
    V12 <- S %*% t(A)
    V21 <- t(V12)

    muW <- c(mu1[1], mu1[2], mu2[1], mu2[2])
    SW <- rbind(cbind(V11, V12), cbind(V21, V22))

    ix_x <- c(1L, 3L)
    ix_y <- c(2L, 4L)
    mu_x <- muW[ix_x]
    mu_y <- muW[ix_y]
    Sxx <- SW[ix_x, ix_x, drop = FALSE]
    Syy <- SW[ix_y, ix_y, drop = FALSE]
    Sxy <- SW[ix_x, ix_y, drop = FALSE]
    iSyy <- solve(Syy)

    z_y <- c(Zy1, Zy2)
    c_mean <- mu_x + Sxy %*% iSyy %*% (z_y - mu_y)
    c_S <- Sxx - Sxy %*% iSyy %*% t(Sxy)

    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Package 'MASS' is required for Gaussian draw generation.")
    }
    z_x_draws <- MASS::mvrnorm(n = m, mu = as.numeric(c_mean), Sigma = c_S)
    if (m == 1) {
      z_x_draws <- matrix(z_x_draws, nrow = 1)
    }

    r_sat_draws <- fitv$F1_inv(stats::pnorm(z_x_draws[, 1]))
    r_sun_draws <- fitv$F1_inv(stats::pnorm(z_x_draws[, 2]))

    if (use_rW != "none" && !is.na(res$rW[k])) {
      delta <- as.numeric(res$rW[k]) - (r_sat_draws + r_sun_draws)
      if (use_rW == "var") {
        if (m == 1) {
          z_x_tmp <- MASS::mvrnorm(n = 512, mu = as.numeric(c_mean), Sigma = c_S)
          r_sat_tmp <- fitv$F1_inv(stats::pnorm(z_x_tmp[, 1]))
          r_sun_tmp <- fitv$F1_inv(stats::pnorm(z_x_tmp[, 2]))
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

