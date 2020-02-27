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
\examples{
# make data lower-tri
brd2 <- map(brd2, ~replace(.x, upper.tri(.x), 0))

# downsample
downsampled_data <- downsample_to_equal_reads(brd2)
sccs <- map(1:2, ~get_loc_sim(downsampled_data[[.x]], downsampled_data[[.x * 2]], win_min = 10, win_max = 10))
zs <- map(sccs, get_z)

par(mfrow = c(1,2))
plot(sccs[[1]], type = "l", ylim = range(unlist(sccs)), ylab = "SCC scan statistic")
lines(sccs[[2]], col = "darkcyan")

plot(zs[[1]]$z_s, type = "l", ylim = range(unlist(map(zs, "z_s"))), ylab = "z-score")
lines(zs[[2]]$z_s, col = "darkcyan")
}
