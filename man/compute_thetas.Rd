\name{compute_thetas}
\alias{compute_thetas}
\title{
Compute theta indicators used for hypothesis testing.
}
\description{
Compute theta indicators used for hypothesis testing using any hesting quantity.
}
\usage{
compute_thetas(mcmc_fit_file, sampled_nngps_file, outpath, return = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcmc_fit_file}{
Path to the rds file output by \link{fit_nngp}
}
  \item{sampled_nngps_file}{
Path to the rds file output by \link{sample_new_nngps}
}
  \item{outpath}{
Path to where the output will be written. The output file should be given the .rds file extension.
}
  \item{return}{
Logical. Should the output also be returned by the function?
}
}
\value{
If return = TRUE, this function outputs a list of indicators indicating evidence for accepting or rejecting the null hypothesis.
}

\examples{
thetas <- compute_thetas(
    mcmc_fit_file =  "inst/extdata/output/fit.rds",
    sampled_nngps_file = "inst/extdata/output/pred.rds",
    outpath = "inst/extdata/output/thetas.rds",
    return = TRUE
    )

str(thetas)

}
