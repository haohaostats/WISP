\name{wisp}
\alias{wisp}
\title{Window Mean Survival Time Pseudo-Observation Regression}
\description{Fits an identity-link regression to infinitesimal jackknife pseudo-observations of a Kaplan--Meier window integral after event-time reconstruction.}
\usage{wisp(left, right, Z, tau0, tau1, impute = c("mid", "right"), strata = NULL)}
\arguments{
\item{left, right}{Numeric interval endpoints. Use finite nonnegative left endpoints and \code{Inf} right endpoints for right censoring. Equal finite endpoints denote exact events. Missing values are rejected. Specify any working lower bound for left-censored events explicitly.}
\item{Z}{Data frame or matrix of outcome-regression covariates, without an intercept. An intercept is added and factors are expanded using \code{model.matrix}. Rows must match endpoints. Missing values are rejected rather than silently omitted.}
\item{tau0, tau1}{Finite window boundaries satisfying \code{0 <= tau0 < tau1}, in the same units as the endpoints.}
\item{impute}{Midpoint reconstruction (default) or documented-date/right-endpoint reconstruction.}
\item{strata}{Optional nonmissing stratum labels. Pseudo-observations are constructed separately within each stratum; the same regression is fitted across all subjects. Each stratum must support the requested window.}
}
\details{The estimand is the linear projection of the reconstructed window response. For fixed assessment intervals this can differ from the latent event-time WMST projection. The coefficient units are units of time within the window, not hazard ratios. The default pooled estimator requires an appropriate independent-censoring representation; including censoring predictors in the outcome regression alone does not correct covariate-dependent dropout.

The product-limit derivative accounts for tied events. Inference uses HC0 and normal Wald intervals without a finite-sample multiplier. Follow-up ending in censoring strictly before the window upper boundary is rejected. No automatic missing-data imputation, left-censoring lower-bound choice, or censoring correction is performed.}
\value{A data frame with columns \code{term}, \code{est}, \code{se}, \code{lcl}, \code{ucl}, and \code{p}. Attributes \code{pseudo} and \code{theta} contain the pseudo-observations and pooled window integral. When \code{strata} is used, \code{theta} remains the pooled integral; it is not a covariate-adjusted mean.}
\examples{
left <- c(0, 1, 1, 2, 2, 3, 4, 4)
right <- left + 1
Z <- data.frame(treatment = rep(c(0, 1), 4))
wisp(left, right, Z, tau0 = 1, tau1 = 4)
}
\seealso{\code{\link{wmst_pseudo}}, \code{\link{wisp_ipcw}}}
