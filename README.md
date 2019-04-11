# Multi-Scale Differential InteractionsSometimes package won't load if the correct shared objects aren't found. So, if

```{r}
devtools::install_github("hillarykoch/sgp")
```

doesn't work, need to try

```{bash}git clone https://github.com/hillarykoch/sgpR CMD build sgp
R
```

```{r}
Sys.setenv("PKG_LIBS" = "/usr/lib64/R/modules/lapack.so")install.packages("sgp_0.0.99.tar.gz", repos = NULL, type = "source")q()
n
```

```{bash}
rm -rf sgp
rm sgp_0.0.99.tar.gz
```