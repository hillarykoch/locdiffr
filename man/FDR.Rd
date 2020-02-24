\name{False discovery control}
\alias{FDR}
\alias{FDX}
\title{Test for significantly differing interactions across scanning windows.}
\description{Based on MCMC output, determine differential genomic regions based on either the False Discovery Rate or False Discovery Exceedance criterion.}
\usage{
FDR(theta, alpha = .1, nthresh = 100)
FDX(theta, alpha = .1, beta = .9, nthresh = 100)
\arguments{
  \item{theta}{Boolean matrix describing how often each sampled value within each scan window was less than the mean process. This matrix is of dimension (number of predictions) by W.}
  \item{alpha}{Numeric between 0 and 1. Confidence level for making rejections.}
  \item{nthresh}{Integer. Number of tresholds to check for evaluating FDR/FDX.}
  \item{beta}{Exceedance level, for calculating FDX. Recommended to be set to 1-alpha.}
}
\value{
\item{level}{Confidence level used.}
\item{reject}{Boolean vector of length W detailing which windows to reject.}
\item{thresh}{Thresholds used in FDR/FDX calculations.}
\item{BFDR}{Bayes False Discovery Rate (may not be exactly alpha).}
\item{BFDX}{Bayes False Discovery Exceedance (may not be exactly alpha).}
}
\references{
Sun, W., Reich, B.J., Cai, T.T., Guindani, M. and Schwartzman, A., 2015. False discovery control in large‐scale spatial multiple testing. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 77(1), pp.59-83.
}
\author{\email{hbk5086@psu.edu}}

