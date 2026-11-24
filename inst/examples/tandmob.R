if (!requireNamespace("icensBKL", quietly = TRUE)) stop("Install icensBKL before running this example.")
library(WISP)
e <- new.env()
utils::data("tandmob", package = "icensBKL", envir = e)
d <- e$tandmob
left <- d$CL16
right <- d$CR16
emergence <- ifelse(is.na(d$L16), 5, d$L16)
left <- ifelse(is.na(left) & !is.na(right), pmin(emergence, right - 0.01), left)
Z <- data.frame(brush_daily = d$BRUSHING, dmft = d$DMF_1,
                plaque_pits = as.integer(d$OH16o == 1), plaque_total = as.integer(d$OH16o == 2),
                sealant16 = d$SEAL16, girl = d$GENDER)
keep <- !(is.na(left) & is.na(right)) & complete.cases(Z)
left <- ifelse(is.na(left[keep]), 0, left[keep])
right <- ifelse(is.na(right[keep]), Inf, right[keep])
Z <- Z[keep, , drop = FALSE]
results <- do.call(rbind, lapply(c("mid", "right"), function(method) {
  rbind(cbind(model = "unadjusted", method = method, wisp(left, right, Z["brush_daily"], 8, 12, method)),
        cbind(model = "adjusted", method = method, wisp(left, right, Z, 8, 12, method)))
}))
windows <- do.call(rbind, lapply(c(7, 8, 9), function(t0) {
  fit <- wisp(left, right, Z, t0, 12)
  cbind(tau0 = t0, tau1 = 12, fit[fit$term == "brush_daily", ])
}))
print(c(n = nrow(Z), events = sum(is.finite(right))))
print(results, row.names = FALSE)
print(windows, row.names = FALSE)
