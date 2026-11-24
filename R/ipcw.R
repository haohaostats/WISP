dropout_records <- function(left, right, W, h, tau1) {
  n <- length(left); K <- round(tau1 / h)
  event <- is.finite(right)
  lastk <- ifelse(event, round(right / h) - 1, round(left / h))
  m <- pmin(lastk, K - 1) + 1
  m[m < 0] <- 0
  id <- rep(seq_len(n), m); k <- sequence(m) - 1
  D <- as.numeric(!event[id] & round(left[id] / h) == k)
  X <- cbind(outer(k, 0:(K - 1), "==") * 1, W[id, , drop = FALSE])
  colnames(X) <- c(paste0("k", 0:(K - 1)), colnames(W))
  list(id = id, k = k, D = D, X = X, K = K)
}
need_index <- function(left, right, h, tau1) {
  event <- is.finite(right)
  ifelse(event, pmin(round(right / h), round(tau1 / h)), round(tau1 / h))
}
complete_case <- function(left, right, tau1) is.finite(right) | left >= tau1 - 1e-9
wisp_ipcw <- function(left, right, Z, W, h, tau0, tau1, G_true = NULL, variance = TRUE) {
  validate_grid(left, right, h, tau0, tau1)
  Z <- validate_design(Z, length(left), "Z", intercept = TRUE)
  if (!is.logical(variance) || length(variance) != 1L || is.na(variance))
    stop("variance must be TRUE or FALSE.", call. = FALSE)
  if (!is.null(G_true) && (!is.numeric(G_true) || length(G_true) != length(left) ||
      any(!is.finite(G_true)) || any(G_true <= 0 | G_true > 1)))
    stop("G_true must contain one observation probability in (0, 1] per subject.", call. = FALSE)
  if (is.null(G_true)) {
    W <- as.matrix(W)
    if (!is.numeric(W) || nrow(W) != length(left) || any(!is.finite(W)))
      stop("W must be a finite numeric matrix with one row per subject and no intercept.", call. = FALSE)
    if (ncol(W) > 0L && (qr(W)$rank < ncol(W) || any(apply(W, 2, function(x) length(unique(x)) < 2L))))
      stop("W must have full rank and no constant columns.", call. = FALSE)
    if (ncol(W) > 0L && is.null(colnames(W))) colnames(W) <- paste0("w", seq_len(ncol(W)))
  }
  n <- length(left); K <- round(tau1 / h)
  Delta <- complete_case(left, right, tau1)
  Th <- ifelse(is.finite(right), (left + right) / 2, Inf)
  g <- pmin(pmax(Th, tau0), tau1) - tau0; g[!Delta] <- 0
  need <- need_index(left, right, h, tau1)
  if (!is.null(G_true)) {
    theta <- ifelse(Delta, g / G_true, 0)
    f <- pseudo_lm(theta, Z)
    Zm <- as.matrix(Z); nn <- nrow(Zm)
    psi <- Zm * drop(theta - Zm %*% f$est)
    return(list(est = f$est, se_hc0 = f$se, theta = theta, w = ifelse(Delta, 1 / G_true, 0),
                psi = psi,
                infl = psi %*% solve(crossprod(Zm) / nn)))
  }
  rec <- dropout_records(left, right, W, h, tau1)
  fit <- suppressWarnings(glm.fit(rec$X, rec$D, family = binomial(link = "cloglog")))
  if (!fit$converged) stop("drop-out model did not converge")
  alpha <- fit$coefficients
  if (any(!is.finite(alpha))) stop("drop-out model not identified")
  Wb <- drop(W %*% alpha[-(1:K)])
  eta <- outer(Wb, alpha[1:K], "+")
  lam <- -expm1(-exp(eta))
  use <- outer(need, 0:(K - 1), ">")
  hazards <- exp(eta)
  hazards[!use] <- 0
  logG <- -rowSums(hazards)
  G <- exp(logG)
  if (any(!is.finite(G[Delta]) | G[Delta] <= 0))
    stop("Estimated observation probabilities are numerically zero; inference is not supported.", call. = FALSE)
  theta <- ifelse(Delta, g / G, 0)
  f <- pseudo_lm(theta, Z)
  out <- list(est = f$est, se_hc0 = f$se, theta = theta, gamma = alpha[-(1:K)], maxw = max(1 / G[Delta]),
              w = ifelse(Delta, 1 / G, 0))
  out$G <- G
  out$complete <- Delta
  out$dropout_fit <- fit
  if (!variance) return(out)
  Z <- as.matrix(Z); p <- ncol(Z); q <- length(alpha)
  H <- crossprod(Z) / n
  U <- Z * drop(theta - Z %*% f$est)
  er <- drop(rec$X %*% alpha); mu <- 1 - exp(-exp(er)); mu <- pmin(pmax(mu, 1e-12), 1 - 1e-12)
  dmu <- exp(er - exp(er))
  S <- rowsum(rec$X * ((rec$D - mu) * dmu / (mu * (1 - mu))), rec$id, reorder = TRUE)
  Sfull <- matrix(0, n, q); Sfull[as.integer(rownames(S)), ] <- S
  Info <- crossprod(rec$X * sqrt(dmu^2 / (mu * (1 - mu)))) / n
  dm <- hazards
  a <- cbind(dm, rowSums(dm) * W)
  Dm <- crossprod(Z * theta, a) / n
  psi <- U + Sfull %*% t(Dm %*% solve(Info))
  Hi <- solve(H)
  V <- Hi %*% (crossprod(psi) / n^2) %*% Hi
  out$vcov_stack <- V
  out$se_stack <- sqrt(diag(V))
  out$psi <- psi
  out$infl <- psi %*% Hi
  out
}
