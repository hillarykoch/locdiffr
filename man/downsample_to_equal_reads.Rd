\name{downsample_to_equal_reads}
\alias{downsample_to_equal_reads}
\title{
Downsample replicate experiments to equal reads.
}
\description{
The SCC score will be biased if one experiment is much more deeply sequenced than another. This function accounts for that.
}
\usage{
downsample_to_equal_reads(data_list)
}
\arguments{
  \item{data_list}{List of lower-triangular count (sparse) matrices.}
}
\details{
The matrices MUST be lower-triangular.
}
\value{
A list of lower-triangular data matrices, downsampled to have equal numbers of read counts.
}
\author{\email{hbk5086@psu.edu}}
