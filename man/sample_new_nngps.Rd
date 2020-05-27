\name{sample_new_nngps}
\alias{sample_new_nngps}
\title{
Sample new nearest-neighbor Gaussian processes for doing inference.
}
\description{
Sample new nearest-neighbor Gaussian processes based on posterior MCMC samples.
}
\usage{
sample_new_nngps(scc_scan_file, mcmc_fit_file, outpath, stationary_iterations, parallel = FALSE, ncores = 10, return = FALSE)
}
\arguments{
  \item{scc_scan_file}{
Path to the rds file output by \link{run_scc_scan}
}
  \item{mcmc_fit_file}{
Path to the rds file output by \link{fit_nngp}
}
  \item{outpath}{
Path to where the output will be written. The output file should be given the .rds file extension.
}
  \item{stationary_iterations}{
Which MCMC iterations to use for posterior samples.
}
  \item{parallel}{
Logical. Should the MCMC fitting be parallelized
}
  \item{ncores}{
If parallel = TRUE, number of cores to parallelize over.
}
  \item{return}{
Logical. Should the output also be returned by the function?
}
}
\details{
This function samples new NNGPs based on the posterior MCMC samples, which are ultimately used for hypothesis testing. MCMC chains need some time to reach convergence, and posterior samples can produce quite large files. For that reason, it is suggested not to use all MCMC iterations as the stationary_iterations. For example, from an MCMC with 15,000 iterations, one could set stationary_iterations=5001:15000
}
\value{
If return = TRUE, this function outputs a list of nearest-neighbor Gaussian processes generated based on on posterior MCMC samples.
}
\examples{
# Use posterior samples to make new draws from the nearest-neighbor Gaussian process, for hypothesis testing
sample_new_nngps(scc_scan_file = "inst/extdata/output/z.rds", mcmc_fit_file = "inst/extdata/output/fit.rds", outpath = "inst/extdata/output/pred.rds", stationary_iterations = 101:150, parallel = FALSE)
}
