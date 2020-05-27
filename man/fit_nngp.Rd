\name{fit_nngp}
\alias{fit_nngp}
\title{
Fit the nearest-neighbor Gaussian process.
}
\description{
Run the MCMC to fit the nearest-neighbor Gaussian process for a collection of window sizes.
}
\usage{
fit_nngp(infile, outpath, num_neighbors = 15, iters = 15000, parallel = FALSE, ncores = 10, return = FALSE)
}
\arguments{
  \item{infile}{
Path to the rds file output by \link{run_scc_scan}
}
  \item{outpath}{
Path to where the output will be written. The output file should be given the .rds file extension.
}
  \item{num_neighbors}{
Integer number of nearest neighbors to use. Defaults to 15.
}
  \item{iters}{
Number of MCMC iterations to use. Defaults to 15,000
}
  \item{parallel}{
Logical. Should the MCMC fitting be parallelized?
}
  \item{ncores}{
If parallel = TRUE, number of cores to parallelize over.
}
  \item{return}{
Logical. Should the output also be returned by the function?
}
}
\details{
This function fits the nearest-neighbor Gaussian process for the window sizes used in the SCC scan associated with the infile.
}
\value{
If return = TRUE, this function outputs a list of posterior samples across the window sizes contained in the infile.
}

\examples{
# Fit the nearest-neighbor Gaussian process model via MCMC
fit_nngp(infile = "inst/extdata/output/z.rds", outpath = "inst/extdata/output/fit.rds", num_neighbors = 5, iters = 150, parallel = FALSE)
}
