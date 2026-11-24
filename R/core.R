impute_ic <- function(left, right, method = c("mid", "right")) {
  method <- match.arg(method)
  event <- is.finite(right)
  time <- switch(method,
                 mid   = ifelse(event, (left + right) / 2, left),
                 right = ifelse(event, right, left))
  data.frame(time = time, status = as.integer(event))
}
wmst_pseudo <- function(time, status, tau0, tau1) {
  validate_window(tau0, tau1)
  if (!is.numeric(time) || length(time) < 2L || any(!is.finite(time)) || any(time < 0) ||
      length(status) != length(time) || anyNA(status) || any(!status %in% c(0, 1)))
    stop("Supply finite nonnegative times and binary event indicators of equal length.", call. = FALSE)
  if (max(time) < tau1 && any(status == 0 & time == max(time)))
    stop("Follow-up does not support tau1: the final observation is censored before tau1.", call. = FALSE)
  n <- length(time)
  fit <- survfit(Surv(time, status) ~ 1)
  tk <- fit$time; S <- fit$surv; Y <- fit$n.risk; d <- fit$n.event
  left <- c(0, tk); right <- c(tk, Inf); Sv <- c(1, S)
  seg <- (pmin(pmax(right, tau0), tau1) - pmin(pmax(left, tau0), tau1)) * Sv
  theta <- sum(seg)
  tail_area <- rev(cumsum(rev(seg)))
  A <- tail_area[seq_along(tk) + 1]
  A[tk >= tau1] <- 0
  fac <- ifelse(Y > d, A / (Y - d), 0)
  term <- ifelse(Y > 0, fac * d / Y, 0)
  cumterm <- cumsum(term)
  idx <- findInterval(time, tk)
  part1 <- ifelse(idx > 0, cumterm[pmax(idx, 1)], 0)
  jump <- match(time, tk)
  part2 <- ifelse(status == 1 & !is.na(jump), fac[jump], 0)
  IF <- n * (part1 - part2)
  list(theta = theta, pseudo = theta + IF,
       se_greenwood = sqrt(sum(ifelse(Y > d, A^2 * d / (Y * (Y - d)), 0))))
}
pseudo_lm <- function(y, X) {
  X <- validate_design(X, length(y), "X")
  XtXi <- solve(crossprod(X))
  beta <- drop(XtXi %*% crossprod(X, y))
  e <- drop(y - X %*% beta)
  meat <- crossprod(X * e)
  V <- XtXi %*% meat %*% XtXi
  se <- sqrt(diag(V))
  data.frame(term = colnames(X), est = beta, se = se,
             lcl = beta - 1.96 * se, ucl = beta + 1.96 * se,
             p = 2 * pnorm(-abs(beta / se)), row.names = NULL)
}
wisp <- function(left, right, Z, tau0, tau1, impute = c("mid", "right"), strata = NULL) {
  validate_intervals(left, right)
  validate_window(tau0, tau1)
  Z <- as.data.frame(Z)
  if (nrow(Z) != length(left) || anyNA(Z))
    stop("Z must have one complete row per subject; missing values are not removed automatically.", call. = FALSE)
  if (!is.null(strata) && (length(strata) != length(left) || anyNA(strata)))
    stop("strata must have one nonmissing value per subject.", call. = FALSE)
  impute <- match.arg(impute)
  d <- impute_ic(left, right, impute)
  ps <- wmst_pseudo(d$time, d$status, tau0, tau1)
  if (!is.null(strata)) {
    pv <- numeric(length(left))
    for (s in unique(strata)) {
      k <- strata == s
      pv[k] <- wmst_pseudo(d$time[k], d$status[k], tau0, tau1)$pseudo
    }
    ps$pseudo <- pv
  }
  Z <- as.data.frame(Z)
  if (nrow(Z) == 0) stop("no rows in covariate data")
  X <- model.matrix(~ ., data = Z)
  out <- pseudo_lm(ps$pseudo, X)
  attr(out, "theta") <- ps$theta
  attr(out, "pseudo") <- ps$pseudo
  out
}
