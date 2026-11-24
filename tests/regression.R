library(WISP)
eq <- function(a, b, tol = 1e-8) stopifnot(isTRUE(all.equal(unname(a), unname(b), tolerance = tol, check.attributes = FALSE)))
fails <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
t <- c(0, 1, 2, 2, 3, 4, 5, 7)
z <- rep(1, length(t))
p <- wmst_pseudo(t, z, 1, 4)
eq(p$pseudo, pmin(pmax(t, 1), 4) - 1)
eq(p$theta, mean(p$pseudo))
eq(p$pseudo, wmst_pseudo(t, z, 0, 4)$pseudo - wmst_pseudo(t, z, 0, 1)$pseudo)
t <- c(1, 1, 2, 2, 3, 4, 5, 5)
s <- c(1, 0, 1, 1, 0, 1, 1, 0)
p <- wmst_pseudo(t, s, 1, 4)
weighted <- function(w) {
  fit <- survival::survfit(survival::Surv(t, s) ~ 1, weights = w)
  l <- c(0, fit$time)
  r <- c(fit$time, Inf)
  sum((pmin(pmax(r, 1), 4) - pmin(pmax(l, 1), 4)) * c(1, fit$surv))
}
eps <- 1e-5
for (i in seq_along(t)) {
  wp <- wm <- rep(1, length(t))
  wp[i] <- wp[i] + eps
  wm[i] <- wm[i] - eps
  eq(p$pseudo[i] - p$theta, length(t) * (weighted(wp) - weighted(wm)) / (2 * eps), 1e-6)
}
left <- c(0, 1, 1, 2, 2, 3, 4, 4)
right <- left + 1
Z <- data.frame(A = rep(c(0, 1), 4))
f <- wisp(left, right, Z, 1, 4)
y <- pmin(pmax((left + right) / 2, 1), 4) - 1
eq(f$est, coef(lm(y ~ A, Z)))
g <- wisp_ipcw(left, right, cbind(intercept = 1, A = Z$A), NULL, 1, 1, 4, rep(1, 8))
eq(f$est, g$est)
eq(f$se, g$se_hc0)
eq(sqrt(diag(crossprod(g$infl) / 64)), g$se_hc0)
fails(wisp(c(NA, 1), c(2, 3), data.frame(A = c(0, 1)), 0, 3))
fails(wisp(left, right, data.frame(A = rep(0, 8)), 1, 4))
fails(wmst_pseudo(c(1, 2, 3), c(1, 0, 0), 0, 4))
fails(wisp_ipcw(left, right + 0.5, cbind(1, Z$A), NULL, 1, 1, 4, rep(1, 8)))
fails(wisp_ipcw(left, right, cbind(1, Z$A), NULL, 1, 1, 4, rep(0, 8)))
ref <- readRDS("research-reference.rds")
d <- ref$data
f <- wisp(d$left, d$right, as.data.frame(d$W), 1, 4)
eq(f$est, ref$wisp$est)
eq(f$se, ref$wisp$se)
f <- wisp_ipcw(d$left, d$right, d$Z, d$W, 1, 1, 4)
eq(f$est, ref$ipcw$est, 1e-7)
eq(f$se_hc0, ref$ipcw$se_hc0, 1e-7)
eq(f$se_stack, ref$ipcw$se_stack, 1e-7)
eq(f$infl, ref$ipcw$infl, 1e-7)
eq(crossprod(f$infl) / nrow(d$Z)^2, f$vcov_stack)
f0 <- wisp_ipcw(d$left, d$right, d$Z, matrix(numeric(0), nrow(d$Z), 0), 1, 1, 4)
stopifnot(all(is.finite(f0$se_stack)))
fails(wisp_ipcw(d$left, d$right, d$Z, cbind(one = 1, d$W), 1, 1, 4))
