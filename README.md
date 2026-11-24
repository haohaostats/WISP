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

### Breast cosmetic deterioration

Compare radiotherapy alone with radiotherapy plus chemotherapy during months 15–42.

**1. Load the data.**

```r
install.packages("icensBKL")
library(WISP)
data("breastCancer", package = "icensBKL")
head(breastCancer)
```

**2. Prepare the event intervals and treatment indicator.**

`low` and `upp` are the event-interval endpoints in months. Missing upper endpoints indicate right censoring. `RT = 1` denotes radiotherapy alone; `RT = 0` denotes radiotherapy plus chemotherapy.

```r
left <- ifelse(is.na(breastCancer$low), 0, breastCancer$low)
right <- ifelse(is.na(breastCancer$upp), Inf, breastCancer$upp)
RT <- as.integer(breastCancer$treat == "radio only")
Z <- data.frame(RT = RT)
```

**3. Fit WISP and view the treatment contrast.**

```r
fit <- wisp(left, right, Z, tau0 = 15, tau1 = 42, strata = RT)
fit[fit$term == "RT", ]
```

Expected result (rounded):

| term | est | se | lcl | ucl | p |
|---|---:|---:|---:|---:|---:|
| RT | 7.533 | 2.281 | 3.062 | 12.005 | 0.000960 |

The radiotherapy-alone group had **7.53 more months of mean deterioration-free time within months 15–42** (95% CI: 3.06–12.00). This is an observed group difference, not a hazard ratio. `strata = RT` constructs pseudo-observations separately within each treatment group.

### Dental caries development

Estimate the association between daily brushing and caries-free time during ages 8–12, adjusting for baseline dental characteristics and sex.

**1. Load the data and define the covariates.**

```r
data("tandmob", package = "icensBKL")
d <- tandmob
Z <- data.frame(
  daily_brushing = d$BRUSHING,
  dmft = d$DMF_1,
  plaque_pits = as.integer(d$OH16o == 1),
  plaque_total = as.integer(d$OH16o == 2),
  sealant = d$SEAL16,
  girl = d$GENDER
)
```

**2. Prepare the intervals and retain complete cases.**

`CL16` and `CR16` record the caries-event interval for tooth 16 in years of age. For left-censored events, this example uses the tooth-emergence lower bound when available, otherwise a working age bound of 5 years. This is a choice for this dataset, not an automatic WISP default.

```r
left <- d$CL16
right <- d$CR16
lower_bound <- ifelse(is.na(d$L16), 5, d$L16)
left_censored <- is.na(left) & !is.na(right)
left[left_censored] <- pmin(
  lower_bound[left_censored], right[left_censored] - 0.01
)

keep <- !(is.na(left) & is.na(right)) & complete.cases(Z)
left <- left[keep]
right <- right[keep]
Z <- Z[keep, , drop = FALSE]
left[is.na(left)] <- 0
right[is.na(right)] <- Inf
nrow(Z)
```

This retains **369 children**. Both models below use these same children.

**3. Fit the unadjusted and adjusted models.**

```r
unadjusted <- wisp(left, right, Z["daily_brushing"], tau0 = 8, tau1 = 12)
adjusted <- wisp(left, right, Z, tau0 = 8, tau1 = 12)

rbind(
  unadjusted = unadjusted[unadjusted$term == "daily_brushing", ],
  adjusted = adjusted[adjusted$term == "daily_brushing", ]
)
adjusted
```

Expected daily-brushing results (rounded):

| Model | est | se | lcl | ucl | p |
|---|---:|---:|---:|---:|---:|
| Unadjusted | 0.665 | 0.252 | 0.170 | 1.159 | 0.0084 |
| Adjusted | 0.630 | 0.233 | 0.174 | 1.087 | 0.0068 |

After adjustment, daily brushing was associated with **0.630 additional caries-free years (about 7.6 months) within ages 8–12**. `adjusted` displays all covariate coefficients. Reference categories are less-than-daily brushing, no plaque, no sealant and boys; the dmft coefficient is per one-unit increase.

**4. Save your results.**

```r
write.csv(adjusted, "tandmob_wisp_results.csv", row.names = FALSE)
```

For your own data, replace the endpoint vectors and covariate data frame, and choose `tau0` and `tau1` in the same units as your event times. Keep endpoints and covariate rows aligned.

## Help

```r
?wisp
?wisp_ipcw
?wmst_pseudo
```

## Authors and license

Hao Chen and Clara Grazian · The University of Sydney  
Contact: hche2036@uni.sydney.edu.au · MIT License
