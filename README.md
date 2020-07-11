# Multi-Scale Differential Interactions

## Dependencies
**This package was developed in R version 3.6.2**

To use this package, you need a compiler that has support for C++11, such as
    *   GCC: [see here, for example](https://www.gnu.org/software/gcc/projects/cxx-status.html#cxx11)
    *   clang: [see here](http://clang.llvm.org/cxx_status.html)
    

Sometimes package won't load if the correct shared objects aren't found. So, if

```r
devtools::install_github("hillarykoch/sgp", build_vignettes = TRUE)
```

doesn't work, need to try

```console
git clone https://github.com/hillarykoch/sgp
R
```

```r
Sys.setenv("PKG_LIBS" = "/usr/lib64/R/modules/lapack.so")
devtools::build(pkg = "sgp/")
install.packages("sgp_0.0.99.tar.gz", repos = NULL, type = "source")
q()
n
```

```console
rm -rf sgp
rm sgp_0.0.99.tar.gz
```
