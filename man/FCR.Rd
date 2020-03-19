\name{False cluster control}
\alias{FCR}
\title{Weighted test for significantly differing interactions across scanning windows, weighted by scanning window size.}
\description{Based on MCMC output, determine differential genomic regions based on False Cluster Rate criterion.}
\usage{
FCR(theta_list, alpha = .1, nthresh = 100)
}
\arguments{
  \item{theta_list}{List of Boolean matrices describing how often each sampled value within each scan window was less than the mean process. Each matrix is of dimension W by (number of predictions).}
  \item{alpha}{Numeric between 0 and 1. Confidence level for making rejections.}
  \item{nthresh}{Integer. Number of tresholds to check for evaluating FCR.}
}
\value{
\item{level}{Confidence level used.}
\item{reject}{Boolean vector of length W detailing which windows to reject.}
\item{thresh}{Thresholds used in FDR/FDX calculations.}
\item{BFCR}{Bayes False Cluster Rate (may not be exactly alpha).}
}
\references{
Sun, W., Reich, B.J., Cai, T.T., Guindani, M. and Schwartzman, A., 2015. False discovery control in large‐scale spatial multiple testing. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 77(1), pp.59-83.
}
\author{\email{hbk5086@psu.edu}}

