if (!requireNamespace("icensBKL", quietly = TRUE)) stop("Install icensBKL before running this example.")
library(WISP)
e <- new.env()
utils::data("breastCancer", package = "icensBKL", envir = e)
d <- e$breastCancer
left <- ifelse(is.na(d$low), 0, d$low)
right <- ifelse(is.na(d$upp), Inf, d$upp)
rad <- as.integer(!grepl("chemo", as.character(d$treat), ignore.case = TRUE))
results <- do.call(rbind, lapply(c(0, 12.5, 15, 17.5), function(t0) {
  do.call(rbind, lapply(c("mid", "right"), function(method) {
    fit <- wisp(left, right, data.frame(rad = rad), t0, 42, method, strata = rad)
    cbind(tau0 = t0, tau1 = 42, method = method, fit[fit$term == "rad", ])
  }))
}))
print(results, row.names = FALSE)
