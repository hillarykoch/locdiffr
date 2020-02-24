\name{extend_nn_mcmc}
\alias{extend_nn_mcmc}
\title{
Extend an existing MCMC chain for fitting and nearest neighbor Gaussian process.
}
\description{
This function adds iterations to the output of \link{run_nnsgp}.
}
\usage{
extend_nn_mcmc(fit, y, s, X, min_range = 1, max_max_range = 100, range_sd = sqrt(2), as = 2, bs = 2, tauinv = NULL, errvar = NULL, sd_beta = 3, nugget_sd = 1, iters = 500, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{Object output from \link{run_nnsgp}.}
  \item{y}{W by n numeric matrix of n Gaussian processes coming from transformed SCC scans across W windows}
  \item{s}{W by 1 integer matrix of scanning window locations}
  \item{X}{W by p numeric design matrix for p regression coefficients (typically, p basis functions)}
  \item{min_range}{Numeric. The minimum value of the spatial range parameter. Defaults to 1.}
  \item{max_max_range}{Numeric. Pre-determined maximum of the spatial range parameter, before subsequent shrinking. Defaults to diff(range(s))}
  \item{range_sd}{Numeric. Prior standard deviation on the spatial range parameter.}
  \item{as}{Numeric. Prior shape parameter for the Exponential variance. Defaults to 2.}
  \item{bs}{Numeric. Prior rate parameter for the Exponential variance. Defaults to 2.}
  \item{tauinv}{Numeric. Initial value of the Exponential variance. Defaults to estimation from an ordinal linear model.}
  \item{errvar}{Numeric. Initial value of the error variance. Defaults to estimation from an ordinal linear model.}
  \item{sd_beta}{Numeric. Prior standard deviation of the regression coefficients.}
  \item{nugget_sd}{Numeric. Prior standard deviation on the error variance. Defaults to 1.}
  \item{iters}{Integer. Number of MCMC iterations to tack on.}
}
\details{
ALL parameters (except fit and iters) should be the same as those used to run the original MCMC that produced the "fit" object.
}
\value{
    \item{beta}{iters by p matrix of MCMC samples of the regression coefficients.}
    \item{covar_params}{iters by 3 matrix of MCMC samples of the Exponential standard deviation, spatial range parameter, and error standard deviation}
    \item{neighbor_info}{Matrix containing information on nearest neighbors used in the analysis.}
}
\author{\email{hbk5086@psu.edu}}
\examples{
}
