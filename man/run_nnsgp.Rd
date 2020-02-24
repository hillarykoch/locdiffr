\name{run_nnsgp}
\alias{run_nnsgp}
\title{
Do inference on nearest neighbor Gaussian process from transformed Hi-C scan statistics
}
\description{
For a given window size, fit a nearest neighbor Gaussian process to a dataset coming from transformed SCC scans. Do inference with MCMC.
}
\usage{
run_nnsgp <- function(y, s, X, num_neighbors = 25, min_range = 1,
                    max_max_range = NULL, init_range = 5, range_sd = sqrt(2),
                    as = 2, bs = 2, tauinv = NULL, errvar = NULL,
                    sd_beta = 3, init_beta = NULL, nugget_sd = 1,
                    iters = 500)
}
\arguments{
  \item{y}{W by n numeric matrix of n Gaussian processes coming from transformed SCC scans across W windows}
  \item{s}{W by 1 integer matrix of scanning window locations}
  \item{X}{W by p numeric design matrix for p regression coefficients (typically, p basis functions)}
  \item{num_neighbors}{Integer. Number of neighbors used in the nearest neighbor Gaussian process}
  \item{min_range}{Numeric. The minimum value of the spatial range parameter. Defaults to 1.}
  \item{max_max_range}{Numeric. Pre-determined maximum of the spatial range parameter, before subsequent shrinking. Defaults to diff(range(s))}
  \item{init_range}{Numeric. Initial value of the spatial range parameter.}
  \item{range_sd}{Numeric. Prior standard deviation on the spatial range parameter.}
  \item{as}{Numeric. Prior shape parameter for the Exponential variance. Defaults to 2.}
  \item{bs}{Numeric. Prior rate parameter for the Exponential variance. Defaults to 2.}
  \item{tauinv}{Numeric. Initial value of the Exponential variance. Defaults to estimation from an ordinal linear model.}
  \item{errvar}{Numeric. Initial value of the error variance. Defaults to estimation from an ordinal linear model.}
  \item{sd_beta}{Numeric. Prior standard deviation of the regression coefficients.}
  \item{init_beta}{Numeric vector. Initial value of the regression coefficients. Defaults to estimation from an ordinal linear model.}
  \item{nugget_sd}{Numeric. Prior standard deviation on the error variance. Defaults to 1.}
  \item{iters}{Integer. Number of MCMC iterations to run.}
}
\details{The model provides a close approximation to an ordinary Gaussian process using nearest-neighbor information. The model has mean function following the design matrix X and regression coefficients beta, spatial random effects based on a Exponential covariance function, plus iid noise across all windows.}
\value{
    \item{beta}{iters by p matrix of MCMC samples of the regression coefficients.}
    \item{covar_params}{iters by 3 matrix of MCMC samples of the Exponential standard deviation, spatial range parameter, and error standard deviation}
    \item{neighbor_info}{Matrix containing information on nearest neighbors used in the analysis.}
}
\author{\email{hbk5086@psu.edu}}
\examples{
}
