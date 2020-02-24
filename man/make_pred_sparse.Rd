\name{make_pred_sparse}
\alias{make_pred_sparse}
\title{Sample from a NNGP for given parameter settings}
\description{For a given output MCMC from \link{run_nnsgp}}
\usage{
make_pred_sparse(fit, y, s, X, stationary_iterations)
}
\arguments{
  \item{fit}{MCMC output from \link{run_nnsgp}.}
  \item{y}{W by n numeric matrix of n Gaussian processes coming from transformed SCC scans across W windows}
  \item{s}{W by 1 integer matrix of scanning window locations}
  \item{X}{W by p numeric design matrix for p regression coefficients (typically, p basis functions)}
  \item{stationary_iterations}{Iterations to retain for making predictions. Not using the first B iterations is akin to applying burning. Thinning can also be applied.}
}
\details{y, s, and X should all be the same as the y, s, and X used to generate the fit object.}
\value{A matrix with W rows and the number of columns equal to the length of stationary_iterations. Each column is a newly sampled nearest neighbor Gaussian process based on the current parameters at a given MCMC iteration.}
\author{\email{hbk5086@psu.edu}}
