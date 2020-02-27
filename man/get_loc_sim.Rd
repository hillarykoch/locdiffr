\name{get_loc_sim}
\alias{get_loc_sim}
\title{
Run an SCC scan (get the local similarity scores) between to Hi-C contact matrices.
}
\usage{
get_loc_sim(mat1, mat2, h = 0, win_min = 25, win_max = 25)
}
\arguments{
  \item{mat1}{lower-triangular Hi-C contact matrix 1.}
  \item{mat2}{lower-triangular Hi-C contact matrix 2.}
  \item{h}{Smoothing parameter used by HiCRep. Defaults to 0.}
  \item{win_min}{Integer. Minimum window size. Setting win_min = win_max will result in doing the SCC scan for one window size.}
  \item{win_max}{Integer. Maximum window size. Setting win_min = win_max will result in doing the SCC scan for one window size.}
}
\author{\email{hbk5086@psu.edu}}
\examples{
# make data lower-tri
brd2 <- map(brd2, ~replace(.x, upper.tri(.x), 0))

# downsample
downsampled_data <- downsample_to_equal_reads(brd2)
sccs <- map(1:2, ~get_loc_sim(downsampled_data[[.x]], downsampled_data[[.x * 2]], win_min = 10, win_max = 10))

plot(sccs[[1]], type = "l", ylim = range(unlist(sccs)), ylab = "SCC scan statistic")
lines(sccs[[2]], col = "darkcyan")
}
