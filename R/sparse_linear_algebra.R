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

get_2d_neighbor_list <- function(y, num_neighbors_loc, num_neighbors_win) {
    neighbor_list <- list()
    neighbor_list[[1]] <- NA
    count <- 1
    windows <- unique(S$win)
    for(j in 2:(sum(dplyr::select(S, win) == windows[1]))) {
        count <- count + 1
        if(j <= (num_neighbors_loc+1)) {
            neighbor_list[[count]] <- 1:(count-1)
        } else {
            neighbor_list[[count]] <- (j - num_neighbors_loc):(count-1)
        }
    }

    for(i in (count+1):nrow(S)) {
        cur_loc <- S$crd[i]
        cur_win <- S$win[i]

        if(which(cur_win == windows) <= num_neighbors_win) { # low window
            if(cur_loc <= num_neighbors_loc) { # low window, low loc
                keepers <- which(((S$crd %in% seq(cur_loc)) & (S$win <= cur_win)))
            } else { # low window, high loc
                keepers <- which(((S$crd %in% ((cur_loc - num_neighbors_loc):cur_loc)) & (S$win <= cur_win)))
            }
            neighbor_list[[i]] <- keepers[-length(keepers)]
        } else { # high window
            if(cur_loc <= num_neighbors_loc) { # high window, low loc
                keepers <- which(((S$crd %in% seq(cur_loc)) & (S$win %in% tail(windows[windows <= cur_win], num_neighbors_win))))
            } else { # high window, high loc
                keepers <- which((S$crd %in% ((cur_loc - num_neighbors_loc):cur_loc)) & (S$win %in% tail(windows[windows <= cur_win], num_neighbors_win)))
            }
            neighbor_list[[i]] <- keepers[-length(keepers)]
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

make_pred_sparse <- function(fit, y, s, X, stationary_iterations, BOOT){
    # stationary_beta <- fit$beta[stationary_iterations,]
    # Xb <- X %*% t(stationary_beta)
    stationary_covar <- fit$covar_params[stationary_iterations,]

    np <- nrow(y)
    d <- fit$neighbor_info$d
    col_index <- fit$neighbor_info$neighbor_column_major_index

    # preds <- matrix(NA, nrow = nrow(y), ncol = length(stationary_iterations))
    preds <- array(NA, dim = c(nrow(y), length(stationary_iterations), BOOT))
    neighbor_list_mod <- purrr::map2(.x = fit$neighbor_info$neighbor_list,
                                     .y = seq_along(fit$neighbor_info$neighbor_list),
                                     ~ c(.x, .y))
    for(j in 1:ncol(preds)) {
        cond_cov <-
            as.matrix(get_cor_sparse(d,
                           range = stationary_covar[j, "range_s"],
                           col_index = col_index) *
            stationary_covar[j, "sigma"])
        preds[, j, ] <- cmake_one_pred_sparse(neighbor_list_mod,
                                              y,
                                              s,
                                              X,
                                              cond_cov,
                                              BOOT,
                                              round(runif(1, 1, .Machine$integer.max)))
    }
    preds <- sweep(preds, 1, rowMeans(y), `+`)
    preds
}