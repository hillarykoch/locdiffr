# Multi-Scale Differential Interactions

## Dependencies
1.  To use this package, you need to download the C++ library LEMON graph library; [download it here.](https://lemon.cs.elte.hu/trac/lemon/wiki/Downloads)

If you do not have them already, you do not need all of LEMON's dependencies to run the sgp package. To ease this during configuring LEMON, you can go into the INSTALL file after downloading LEMON and add options

```
-DLEMON_ENABLE_GLPK=NO
-DLEMON_ENABLE_COIN=NO
-DLEMON_ENABLE_ILOG=NO
```

<!---
LEMON citation:
Balázs Dezső, Alpár Jüttner, Péter Kovács. LEMON – an Open Source C++ Graph Template Library. Electronic Notes in Theoretical Computer Science, 264:23-45, 2011. Proceedings of the Second Workshop on Generative Technologies (WGT) 2010.
-->

2.  You need a compiler that has support for C++11, such as
    *   GCC: [see here, for example](https://www.gnu.org/software/gcc/projects/cxx-status.html#cxx11)
    *   clang: [see here](http://clang.llvm.org/cxx_status.html)
    

Sometimes package won't load if the correct shared objects aren't found. So, if

```r
devtools::install_github("hillarykoch/sgp")
```

doesn't work, need to try

```console
git clone https://github.com/hillarykoch/sgp
R CMD build sgp
R
```

```r
Sys.setenv("PKG_LIBS" = "/usr/lib64/R/modules/lapack.so")
install.packages("sgp_0.0.99.tar.gz", repos = NULL, type = "source")
q()
n
```

```console
rm -rf sgp
rm sgp_0.0.99.tar.gz
```
