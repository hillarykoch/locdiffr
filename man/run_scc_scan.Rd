\name{run_scc_scan}
\alias{run_scc_scan}
\title{
Compute the SCC scan statistic for a sequence of window sizes.
}
\description{
This function implements the SCC scan and transforms the scan to Z-scores. It can be paralllized across a grid of window sizes.
}
\usage{
run_scc_scan(infiles1, infiles2, outpath, resolution, winsizes = seq(10, 100, by = 10), min_count = 1, proportion_below_min_count = 0.25, parallel = FALSE, ncores = 10, offset = TRUE, return = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{infiles1}{
List or vector of paths to data files from the first condition.
}
  \item{infiles2}{
List or vector of paths to data files from the second condition.
}
  \item{outpath}{
Path to where the output will be written. The output file should be given the .rds file extension.
}
  \item{resolution}{
Resolution of the data being analysis, e.g., 50000 for 50 kilobase resolution.
}
  \item{winsizes}{
A vector of window sizes to use for the scan, in bins.
}
  \item{min_count}{
  A desired minimum number of counts within each bin. See details.
}
  \item{proportion_below_min_count}{
  If the proportion of bins within a window with counts less than "min_count" is greater than this value, the window is filtered. See details.
}
  \item{parallel}{
Logical. Should the function be parallelized?
}
  \item{ncores}{
If parallel=TRUE, how many cores to parallelize over.
}
  \item{offset}{
Should an offset of 1 be used? Use FALSE if the input files are one-indexed, but TRUE if the input files are zero-indexed.
}
\item{return}{
Logical. Should the output also be returned by the function?
}
}
\details{
This function writes the output to the provided outpath in an rds file.

Windows with insufficient reads can be filtered from the analysis using "min_count" and "proportion_below_min_count" options. The default sets min_count = 1 and proportion_below_min_count = 0.25. This means that if more than 25\% of the bins within a window have less than 1 measured interaction, the window is filtered from downstream analyses. To do no filtering, set proportion_below_min_count=1.
}
\value{
If return = TRUE, this function returns a list of length winsizes which contains the SCC scan for each window size.
}

\examples{
# Define paths to data from 2 conditions
all_data <- list.files(system.file("extdata", package = "sgp"))
ko_data <-
    system.file("extdata", all_data[grepl(pattern = "KO", x = all_data)], package = "sgp")
parent_data <-
    system.file("extdata", all_data[grepl(pattern = "parent", x = all_data)], package = "sgp")

# Compute the SCC scan statistics
run_scc_scan(
    infiles1 = ko_data,
    infiles2 = parent_data,
    outpath = "inst/extdata/output/z.rds",
    resolution = 40000,
    winsizes = c(10, 25),
    parallel = FALSE,
    offset = TRUE
)
}
