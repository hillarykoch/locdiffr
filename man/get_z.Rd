\name{get_z}
\alias{get_z}
\title{
Get Gaussian process from SCC scan.
}
\description{
This function applies Fisher's z-transformation to an SCC scan.
}
\usage{
get_z(SCCs)
}
\arguments{
  \item{SCCs}{Numeric vector of SCC scores from a given scan.}
}
\value{
Data frame with 3 columns. crd is the scan locations, scc is the scc score, and z_s is the transformed z-score.
}
\author{\email{hbk5086@psu.edu}}
