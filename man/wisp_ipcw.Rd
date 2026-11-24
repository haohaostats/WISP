\name{wisp_ipcw}
\alias{wisp_ipcw}
\title{Common-Grid IPCW Window Mean Survival Time Regression}
\description{Fits complete-response inverse-probability weighted regression targeting the midpoint-reconstructed WMST projection under baseline-dependent monotone dropout.}
\usage{wisp_ipcw(left, right, Z, W, h, tau0, tau1, G_true = NULL,
  variance = TRUE)}
\arguments{
\item{left, right}{Observed interval endpoints on the common grid. Finite event intervals must have width \code{h}. A right-censored subject has \code{right = Inf} and \code{left} equal to the last attended event-free grid assessment. No exact events or intermittent missed assessments are supported.}
\item{Z}{Finite full-rank numeric outcome design matrix INCLUDING an intercept column of ones. Unlike \code{wisp}, no intercept is added automatically.}
\item{W}{Finite numeric matrix of baseline dropout-model covariates, WITHOUT an intercept. Interval-specific intercepts are supplied internally. A zero-column matrix specifies no covariates. Ignored when \code{G_true} is supplied; pass \code{NULL} in that case.}
\item{h}{Positive common assessment interval.}
\item{tau0, tau1}{Window boundaries on the grid satisfying \code{0 <= tau0 < tau1}.}
\item{G_true}{Optional vector of KNOWN observation probabilities in (0,1], one per subject. For observed events evaluate at \code{min(right, tau1)}; for right-censored observations use the probability of attending through \code{tau1}. This is not a vector of inverse weights. Externally estimated probabilities passed here are treated as known; their estimation uncertainty is not included.}
\item{variance}{For estimated probabilities, whether to compute the stacked sandwich. HC0 is returned in either case. With known probabilities this argument does not disable HC0.}
}
\details{The method requires the potential last attended assessment to be conditionally independent of event time given baseline covariates, positivity, and an appropriate observation model. The fitted dropout model uses a complementary log-log link and one intercept per assessment period. This extension is for a common fixed grid without missed visits, not a general irregular-visit IPCW procedure.

Complete responses are observed events or subjects followed through the window endpoint. The response is the midpoint window contribution times its inverse observation probability, or zero if incomplete. All subjects enter the outcome regression; this is not weighted least squares restricted to complete cases. The target remains the reconstructed projection, not a bias-corrected latent-time projection.

For estimated weights, use \code{se_stack} for inference. \code{se_hc0} treats weights as fixed and is a diagnostic comparator. Influence rows in \code{infl} are coefficient influence functions, including the inverse regression bread, not raw estimating-equation contributions. Model failures and numerically zero complete-response probabilities cause errors. No truncation of weights is performed.}
\value{A list containing \code{est}, \code{se_hc0}, weighted responses \code{theta}, and inverse weights \code{w} (zero for incomplete responses). Known probabilities additionally return \code{psi} and \code{infl}. Estimated probabilities additionally return \code{gamma}, \code{maxw}, \code{G}, \code{complete}, and \code{dropout_fit}; when \code{variance = TRUE}, also \code{se_stack}, \code{vcov_stack}, \code{psi}, and \code{infl}. Matrix column order follows \code{Z}.}
\examples{
left <- c(0, 1, 1, 2, 2, 3, 4, 4)
right <- left + 1
Z <- cbind(intercept = 1, treatment = rep(c(0, 1), 4))
wisp_ipcw(left, right, Z, W = NULL, h = 1,
          tau0 = 1, tau1 = 4, G_true = rep(1, 8))$est
}
