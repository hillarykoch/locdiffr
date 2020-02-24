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
