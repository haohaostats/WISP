# WISP

**Window Mean Survival Time Regression for Interval-Censored Data Using Pseudo-Observations**

## Installation

```r
install.packages("remotes")
remotes::install_github("haohaostats/WISP")
library(WISP)
```

## Basic usage

Supply one row per subject. Use `right = Inf` for right censoring and `left = right` for exact events. Endpoints and window boundaries must use the same time unit. Supply complete covariates; an intercept is added automatically.

```r
left <- c(0, 1, 1, 2, 2, 3, 4, 4)
right <- left + 1
Z <- data.frame(treatment = rep(c(0, 1), 4))

fit <- wisp(left, right, Z, tau0 = 1, tau1 = 4)
fit
```

The output contains coefficient estimates (`est`), HC0 standard errors (`se`), 95% confidence limits (`lcl`, `ucl`) and p-values (`p`). Coefficients are differences in mean event-free time within the selected window.

```r
wisp(left, right, Z, 1, 4, impute = "right")
wisp(left, right, Z, 1, 4, strata = Z$treatment)
```

`impute = "right"` uses documented event dates. `strata` constructs pseudo-observations separately within each specified group.

## Covariate-dependent dropout

`wisp_ipcw()` supports a common assessment grid without missed visits. Supply an intercept in the outcome matrix `Z`, but not in the dropout matrix `W`. Use `se_stack` for inference with estimated weights.

```r
set.seed(42)
n <- 600
A <- rbinom(n, 1, 0.5)
X <- rnorm(n)
T <- rexp(n, rate = 0.15 * exp(-0.3 * A + 0.3 * X))
D <- rexp(n, rate = 0.08 * exp(0.2 * A + 0.3 * X))
V <- pmin(floor(D), 6)
R <- ceiling(T)
event <- R <= V

left <- ifelse(event, R - 1, V)
right <- ifelse(event, R, Inf)
Z <- cbind(intercept = 1, A = A, X = X)
W <- cbind(A = A, X = X)

fit <- wisp_ipcw(left, right, Z, W, h = 1, tau0 = 1, tau1 = 4)
cbind(estimate = fit$est, SE = fit$se_stack)
```

## Real-data examples

Run the breast-cosmesis and dental-caries analyses using data from `icensBKL`:

```r
install.packages("icensBKL")
source(system.file("examples", "bcos.R", package = "WISP"))
source(system.file("examples", "tandmob.R", package = "WISP"))
```

## Help

```r
?wisp
?wisp_ipcw
?wmst_pseudo
```

## Authors and license

Hao Chen and Clara Grazian · The University of Sydney  
Contact: hche2036@uni.sydney.edu.au · MIT License
