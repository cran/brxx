% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/brxx_Cor_general.R
\name{brxx_Cor_general}
\alias{brxx_Cor_general}
\title{brxx_Cor_general: Bayesian Estimation of Reliability from Correlation, General Form}
\usage{
brxx_Cor_general(cor, alpha, beta, CI, items)
}
\arguments{
\item{cor}{Correlation estimate.}

\item{alpha}{Prior true score variance.}

\item{beta}{Prior error variance.}

\item{CI}{Credible interval quantile, as a decimal (ie, for 95 percent, 0.95).}

\item{items}{Number of test items.}
}
\value{
Returns estimated median and quantile based credible limits for reliability.
}
\description{
This function estimates reliability from correlation given the correlation estimate.
}
\examples{
brxx_Cor_general(cor=0.85,alpha=3.51,beta=1.75,CI=0.95,items=10)

}
