# Sparse functions, like in Finley et al. (2019)
library(Matrix)
library(tidyverse)
library(inline)
library(Rcpp)
library(RcppArmadillo)

get_neighbor_list <- function(y, num_neighbors) {
    neighbor_list <- list()
    neighbor_list[[1]] <- NA
    for(i in 2:nrow(y)) {
        if(i <= (num_neighbors+1)) {
            neighbor_list[[i]] <- 1:(i-1)
        } else {
            neighbor_list[[i]] <- (i-num_neighbors):(i-1)
        }
    }
    return(neighbor_list)
}

get_cor_sparse <- function(d, range = NULL, col_index) {
    # Compute the Exponential covariance function
    d[col_index] <- fields::Exponential(d[col_index], range = range)
    diag(d) <- 1
    d
}