\name{test with weighted quantities}
\alias{test_by_wFDX}
\alias{test_by_wFDR}
\title{
Perform hypothesis testing using weighted false discovery controlling quantities wFDR and/or wFDX.
}
\description{
Perform hypothesis testing using weighted false discovery controlling quantities wFDR and/or wFDX.
}
\usage{
test_by_wFDR(scc_scan_file, sampled_nngps_file, outpath, alpha = .01, nthresh = 100, return = FALSE)
test_by_wFDX(scc_scan_file, sampled_nngps_file, outpath, alpha = .01, beta = 0.01, nthresh = 100, return = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scc_scan_file}{
Path to the rds file output by \link{run_scc_scan}
}
  \item{sampled_nngps_file}{
Path to the rds file output by \link{sample_new_nngps}
}
  \item{outpath}{
Path to where the output will be written. The output file should be given the .rds file extension.
}
  \item{alpha}{
Confidence level at which to control false discoveries for the hypothesis text.
}
  \item{beta}{
Additional exceedance threshold used for testing wFDX. Recommended to set to the same value as alpha.
}
  \item{nthresh}{
Number of thresholds used to approximate wFDR and wFDX. Defaults to 100.
}
  \item{return}{
Logical. Should the output also be returned by the function?
}
}
\value{
If return = TRUE, this function outputs a list of genomic loci and an indicator of whether or not the locus was rejected.
}

\examples{
# Compute rejections based on weighted false discovery rate
test_by_wFDR(
    scc_scan_file = "inst/extdata/output/z.rds",
    sampled_nngps_file = "inst/extdata/output/pred.rds",
    outpath = "inst/extdata/output/wfdr.rds",
    alpha = 0.1
    )

# Compute rejections based on weighted false discovery exceedence
test_by_wFDX(
    scc_scan_file = "inst/extdata/output/z.rds",
    sampled_nngps_file = "inst/extdata/output/pred.rds",
    outpath = "inst/extdata/output/wfdx.rds",
    alpha = 0.1,
    beta = 0.1
    )
}
