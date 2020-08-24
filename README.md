# locdiffr (Analysis of LOCal DIFFerences in chromatin aRchitecture)

## Dependencies

1. **This package was developed in R version 3.6.2, back compatability is not guaranteed**

2. To use this package, you need a compiler that has support for C++11, such as
    *   GCC: [see here, for example](https://www.gnu.org/software/gcc/projects/cxx-status.html#cxx11)
    *   clang: [see here](http://clang.llvm.org/cxx_status.html)
    

## Installation

Installation is easy with the R package `devtools`. If you don't already have `devtools` installed, you will first need to type (from within R)
```r
install.packages("devtools")
```

With `devtools` installed, simply enter the following:
```r
devtools::install_github("hillarykoch/locdiffr", build_vignettes = TRUE)
```

Installation may take a couple of minutes due to building the vignette. The vignette contains all of the necessary instructions to run `locdiffr`. To view the vignette, enter the following within R:
```r
browseVignettes(package = "locdiffr")
```

## Input format

To begin the analysis, you need file paths to all of the data, split into two groups. The data should be in the following tab-delimited format:

| Loc1 | Loc2   | Counts |
|------|--------|--------|
| 0    | 0      | 100    |
| 0    | 50000  | 88     |
| 0    | 100000 | 40     |

Here, Loc1 is the start of the first bin, Loc2 the start of the second bin, and counts indicating the read counts for the interactions between those 2 bins. The value in Loc1 should always be less than or equal to the value in Loc2. _**An individual dataset is intended to be Hi-C data from an entire chromosome. Chromosomes should be analyzed separately.**_
