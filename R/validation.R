validate_window <- function(tau0, tau1) {
  if (!is.numeric(tau0) || length(tau0) != 1L || !is.finite(tau0) ||
      !is.numeric(tau1) || length(tau1) != 1L || !is.finite(tau1) ||
      tau0 < 0 || tau1 <= tau0) stop("Require finite 0 <= tau0 < tau1.", call. = FALSE)
}
validate_intervals <- function(left, right) {
  if (!is.numeric(left) || !is.numeric(right) || length(left) < 2L ||
      length(right) != length(left) || anyNA(left) || anyNA(right) ||
      any(!is.finite(left)) || any(left < 0) || any(right < left))
    stop("Supply equal-length interval vectors: finite nonnegative left and right >= left; use Inf for right censoring.", call. = FALSE)
}
validate_design <- function(X, n, name, intercept = FALSE) {
  X <- as.matrix(X)
  if (!is.numeric(X) || nrow(X) != n || ncol(X) < 1L || any(!is.finite(X)))
    stop(paste0(name, " must be a finite numeric design matrix with one row per subject."), call. = FALSE)
  if (nrow(X) <= ncol(X) || qr(X)$rank < ncol(X))
    stop(paste0(name, " must have full column rank and more rows than columns."), call. = FALSE)
  if (intercept && !any(vapply(seq_len(ncol(X)), function(j) all(X[, j] == 1), logical(1))))
    stop(paste0(name, " must include an intercept column of ones."), call. = FALSE)
  if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(ncol(X)))
  if (anyDuplicated(colnames(X))) stop(paste0(name, " must have unique column names."), call. = FALSE)
  X
}
validate_grid <- function(left, right, h, tau0, tau1) {
  validate_intervals(left, right)
  validate_window(tau0, tau1)
  if (!is.numeric(h) || length(h) != 1L || !is.finite(h) || h <= 0)
    stop("h must be a positive finite scalar.", call. = FALSE)
  v <- c(left, right[is.finite(right)], tau0, tau1) / h
  if (any(abs(v - round(v)) > 1e-7))
    stop("IPCW requires endpoints and window boundaries on the common grid.", call. = FALSE)
  ev <- is.finite(right)
  if (any(abs((right[ev] - left[ev]) / h - 1) > 1e-7))
    stop("IPCW requires event intervals of width h without missed assessments or exact events.", call. = FALSE)
}
