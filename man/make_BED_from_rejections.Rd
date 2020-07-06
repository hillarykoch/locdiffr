\name{make_BED_from_rejections}
\alias{make_BED_from_rejections}
\title{
Write BED files from testing output.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make_BED_from_rejections(infiles, chromosomes, null_outfile, rejected_outfile, resolution, region_size = 5000)
}
\arguments{
  \item{infiles}{
 Paths to output from one of test_by_wFDR or test_by_wFDX. For analyses of multiple chromosomes, this is a vector of paths.
}
  \item{chromosomes}{
    (Integer or character) vector of chromosome numbers, in the same order as the chromosomes examined in each infile.
}
  \item{null_outfile}{
path to file where BED file will be written for null regions.
}
  \item{rejected_outfile}{
path to file where BED file will be written for rejected regions.
}
  \item{resolution}{
Resolution of the Hi-C data.
}
  \item{region_size}{
Number of base pairs to fall into each region of the BED file. Must divide evenly into resolution. Defaults to 5,000.
}
}
\value{
Writes two outfiles in 3-column BED format, containing loci which are labeled as falling within either a null or rejected region based on the testing info supplied.
}

\examples{
make_BED_from_rejections(
    infiles = "inst/extdata/output/wfdr.rds",
    chromosomes = 19,
    null_outfile = "inst/extdata/output/null_wfdr.bed",
    rejected_outfile = "inst/extdata/output/rejected_wfdr.bed",
    resolution = 40000
    )
}
