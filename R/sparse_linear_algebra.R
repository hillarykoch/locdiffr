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
    # d[!(seq(d) %in% col_index)] <- 0
    d[col_index] <- fields::Exponential(d[col_index], range = range)
    diag(d) <- 1
    d
}

make_pred_sparse <- function(fit, y, s, X, stationary_iterations){
    stationary_beta <- fit$beta[stationary_iterations,]
    stationary_covar <- fit$covar_params[stationary_iterations,]
    Xb <- X %*% t(stationary_beta)
    
    np <- nrow(y)
    d <- fit$neighbor_info$d
    col_index <- fit$neighbor_info$neighbor_column_major_index
    cond_covs <-
        purrr::map2(
            .x = stationary_covar[, "sigma"],
            .y = stationary_covar[, "range_s"],
            ~ get_cor_sparse(
                d,
                range = .y,
                col_index = col_index
            ) * .x
        ) %>%
        purrr::map(as.matrix)
    
    preds <- matrix(NA, nrow = nrow(y), ncol = length(stationary_iterations))
    neighbor_list_mod <- purrr::map2(.x = fit$neighbor_info$neighbor_list,
                                     .y = seq_along(fit$neighbor_info$neighbor_list),
                                     ~ c(.x, .y))
    for(j in 1:ncol(preds)) {
        preds[,j] <- cmake_one_pred_sparse(neighbor_list_mod,
                                           y,
                                           s,
                                           X,
                                           cond_covs[[j]])
    }
    preds <- preds + Xb
    preds <- sweep(preds, 1, rowMeans(y), `+`)
    preds
}
