% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bomega_general.R
\name{bomega_general}
\alias{bomega_general}
\title{bomega_general: Bayesian Estimation of Coefficient Omega, General Form}
\usage{
bomega_general(lambda, psi, alpha, beta, CI)
}
\arguments{
\item{lambda}{vector of item loadings.}

\item{psi}{vector of item variances.}

\item{alpha}{Prior true score variance.}

\item{beta}{Prior error variance.}

\item{CI}{Credible interval quantile, as a decimal (ie, for 95 percent, 0.95).}
}
\value{
Returns estimated median and quantile based credible limits for omega.
}
\description{
This function estimates coefficient omega internal consistency reliability.
}
\examples{
lambda=c(0.7,0.5,0.6,0.7)
psi=c(0.2,0.4,0.3)
alpha=3.51
beta=1.75
bomega_general(lambda=lambda,psi=psi,alpha=alpha,beta=beta,CI=0.95)

}
