\name{wmst_pseudo}
\alias{wmst_pseudo}
\title{Window Integral and Infinitesimal Jackknife Pseudo-Observations}
\description{Computes a Kaplan--Meier window mean and its infinitesimal jackknife pseudo-observations, including tied events.}
\usage{wmst_pseudo(time, status, tau0, tau1)}
\arguments{
\item{time}{Finite nonnegative observed times, after any reconstruction.}
\item{status}{Event indicators, 1 for events and 0 for right censoring.}
\item{tau0, tau1}{Window boundaries satisfying \code{0 <= tau0 < tau1}.}
}
\details{The functional integrates the Kaplan--Meier survival function from \code{tau0} to \code{tau1}. The supplied data must identify the window: a final censored observation before \code{tau1} is rejected. If the final event brings the curve to zero, integration beyond it is allowed.}
\value{A list with \code{theta} (window integral), \code{pseudo} (one pseudo-observation per subject), and \code{se_greenwood} (Greenwood-based standard error of the window integral).}
\examples{
wmst_pseudo(c(1, 2, 2, 4, 5), rep(1, 5), 1, 4)
}
