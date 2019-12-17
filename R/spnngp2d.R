run_2dnnsgp <- function(y,
                        s,
                        X,
                        num_neighbors_loc = 25, # Number of nearest neighbors within a given scan
                        num_neighbors_win = 2, # Number of nearest neighbors across different window sizes
                        min_range = 1,
                        max_max_range = NULL, # max value for the spatial range parameter
                        init_range = 5, # init value of the spatial range parameter
                        range_sd = sqrt(2), # prior sd for the spatial range
                        as = 2, # the prior for the exponential variance is InvG(as,bs)
                        bs = 2,
                        tauinv = NULL, # initial value of the variance term in exponential cov
                        errvar = NULL, # nugget variance
                        sd_beta = 3, # regression coefficients have N(0, sd_beta) priors
                        init_beta = NULL, # initial value of regression coefs
                        nugget_sd = 1, # prior sd for nugget factor
                        iters = 500) {
    # number of observation locations, number of process replicates (sample size)
    n <- nrow(S)
    if(typeof(y) == "list") {
        reps <- length(y)
    } else {
        reps <- ncol(y)
    }

    # number of predictors
    p <- length(X)
    precision_beta <- diag(p) / (sd_beta ^ 2)


    #----------------------------------------------------
    # Initial values
    #----------------------------------------------------
    if(is.null(max_max_range)) {
        max_max_range <- diff(range(S))
    }
    max_range <- 25 #get_max_range(y, max_max_range)

    # These all default to NULL
    beta <- init_beta
    rhos <- init_range

    # Get initial regression values from a simple linear model fit
    if (reps > 1) {
        X_unfolded <- sapply(X, function(Y)
            Y[keepidx])
        Y_unfolded <-
            as.matrix(as.vector(Reduce(`+`, y) / reps), ncol = 1)
        Y_bound <- bind_cols(y)
        lmfit <- lm(Y_unfolded ~ X_unfolded - 1)
    } else {
        X_unfolded <- sapply(X, as.vector)
        Y_unfolded <- Y_bound <- as.matrix(as.vector(y), ncol = 1)
        lmfit <- lm(Y_unfolded ~ X_unfolded - 1)
    }
    if (is.null(beta)) {
        beta <- lmfit$coef
    }
    if (is.null(tauinv)) {
        tauinv <- var(lmfit$res)
    }
    if (is.null(errvar)) {
        errvar <- 1
    }

    # init value of the log matern spatial range parameter
    if (is.null(rhos)) {
        rhos  <- (max_range - 1) / 2
    }

    # tau is precision of matern covariance function
    tau <- 1 / tauinv

    # errprec  is the precision of the nugget
    errprec <- 1/errvar
    rm(lmfit)

    tune_var_epsilon <- nugget_sd ^ 2 # proposal variance for nugget process
    tune_var <- range_sd ^ 2 # proposal variance for spatial range

    # cor cur is the sparse approximation to the
    #   exponential correlation plus the nugget (not quite a correlation actually)
    # dtemp <- cpwdist(s, s)
    dtemp <- fields::rdist(S)

    neighbor_list <- get_2d_neighbor_list(y, num_neighbors_loc = num_neighbors_loc, num_neighbors_win = num_neighbors_win)
    neighbor_idx <- cbind(unlist(purrr::map2(.x = seq_along(neighbor_list), .y = neighbor_list, ~ rep(.x, times = sum(!is.na(.y))))),
                          unlist(neighbor_list[-1]))

    # access non-zero matrix elements by the column major index
    neighbor_column_major_idx <- c(neighbor_idx[,1] + nrow(dtemp) * (neighbor_idx[,2] - 1),
                                   neighbor_idx[,2] + nrow(dtemp) * (neighbor_idx[,1] - 1))
    d <- Matrix::sparseMatrix(i = neighbor_idx[,1],
                              j = neighbor_idx[,2],
                              x = dtemp[neighbor_column_major_idx[1:nrow(neighbor_idx)]],
                              dims = dim(dtemp),
                              symmetric = TRUE)
    rm(dtemp)

    cor_cur <- get_cor_sparse(d, rhos, neighbor_column_major_idx) + diag(errvar, n)
    cov_cur <- cor_cur * tauinv

    # store regression coefficients and covariance params
    beta_chain <- matrix(0, iters, p)
    param_chain <- matrix(0, iters, 3)
    colnames(param_chain) <-
        c("sigma", "range_s", "epsilon")

    # Kick off solving for A and D
    A_and_D <- csolve_for_A_and_D_2d(cov_cur = cov_cur, neighbor_list = neighbor_list)

    pb <- txtProgressBar(min = 1,
                         max = iters,
                         style = 3)
    for (i in 1:iters) {
        #-------------------------------------------------------------
        # Update regression coefficients
        #-------------------------------------------------------------
        B_and_b <- csolve_for_B_and_b_2d(Y_unfolded, X_unfolded, A_and_D$A, A_and_D$D@x, neighbor_list, precision_beta)
        beta <- RandomFieldsUtils::solvex(B_and_b$B, B_and_b$b) +
            backsolve(t(chol(B_and_b$B)), rnorm(p), upper.tri = FALSE)
        Xb <- as.vector(X_unfolded %*% beta)

        #-------------------------------------------------------------
        # Update covariance parameters for the spatial signal process
        #-------------------------------------------------------------
        # tau ~ Gamma
        yminusXb <- sweep(Y_bound, 1, Xb)
        SS <-
            apply(yminusXb, 2, function(U)
                csparse_quadratic_form_symm_2d(
                    u = U,
                    A = A_and_D$A,
                    D = A_and_D$D@x,
                    neighbor_list = neighbor_list
                )) %>%
            sum
        SS_bare <- SS * tauinv

        tau <- rgamma(1, (n * reps) / 2 + as, SS_bare / 2 + bs)
        tauinv <- 1/tau

        cov_cur <- tauinv * cor_cur
        A_and_D <- csolve_for_A_and_D_2d(cov_cur, neighbor_list)

        #-----------------------------------------------------------------------
        # Sampling covariance parameters variable-at-a-time
        #-----------------------------------------------------------------------
        # Spatial range
        lb <- pnorm(-(rhos-min_range), mean = 0, sd = sqrt(tune_var))
        ub <- pnorm(max_range-rhos, mean = 0, sd = sqrt(tune_var))
        uuu <- runif(1, lb, ub)
        rhos_star <- rhos + qnorm(uuu, mean = 0, sd = sqrt(tune_var))

        # Update covariances
        cor_star <-
            get_cor_sparse(d, rhos_star, neighbor_column_major_idx) + diag(errvar, n)
        cov_star <- tauinv * cor_star
        A_and_D_star <- csolve_for_A_and_D_2d(cov_star, neighbor_list)

        # Already have SS from before
        SS_star <- apply(yminusXb, 2, function(U)
            csparse_quadratic_form_symm_2d(U,
                                        A = A_and_D_star$A,
                                        D = A_and_D_star$D@x,
                                        neighbor_list = neighbor_list)) %>%
            sum
        ldet_cur <- abs(sum(log(A_and_D$D@x)))
        ldet_star <- abs(sum(log(A_and_D_star$D@x)))

        R <- (reps/2) * (ldet_cur - ldet_star) + 0.5 * (SS - SS_star)
        if (runif(1) < exp(R)) {
            rhos <- rhos_star
            cor_cur <- cor_star
            cov_cur <- cov_star
            ldet_cur <- ldet_star
            A_and_D <- A_and_D_star
            SS <- SS_star
        }

        #-------------------------------------------------------------
        # Update error variance term
        #-------------------------------------------------------------
        # errprec ~ Gamma
        lb <- pnorm(-errprec, mean = 0, sd = sqrt(tune_var_epsilon))
        uuu <- runif(1, lb, 1)
        errprec_star <- errprec + qnorm(uuu, mean = 0, sd = sqrt(tune_var_epsilon))
        errvar_star <- 1/errprec_star

        cor_star <- cor_cur - diag(errvar, n) + diag(errvar_star, n)
        cov_star <- tauinv * cor_star

        # Already have SS from before
        A_and_D_star <- csolve_for_A_and_D_2d(cov_star, neighbor_list)
        SS_star <- apply(yminusXb, 2, function(U)
            csparse_quadratic_form_symm_2d(U,
                                        A = A_and_D_star$A,
                                        D = A_and_D_star$D@x,
                                        neighbor_list = neighbor_list)) %>%
            sum
        ldet_star <- abs(sum(log(A_and_D_star$D@x)))

        R <- dgamma(errprec_star, as, bs) - dgamma(errprec, as, bs) +
            0.5 * (SS - SS_star) +
            reps / 2 * (ldet_cur - ldet_star)

        if (runif(1) < exp(R)) {
            errvar <- errvar_star
            errprec <- errprec_star
            cor_cur <- cor_star
            A_and_D <- A_and_D_star
        }


        beta_chain[i,] <- beta
        param_chain[i,] <-
            c(1 / sqrt(tau), rhos, sqrt(errvar))

        setTxtProgressBar(pb, i)
    }
    close(pb)

    list("beta" = beta_chain,
         "covar_params" = param_chain,
         "neighbor_info" = list(
             "num_neighbors_loc" = num_neighbors_loc,
             "num_neighbors_win" = num_neighbors_win,
             "neighbor_list" = neighbor_list,
             "neighbor_column_major_index" = neighbor_column_major_idx,
             "d" = d)
    )
}



# extend_nn_mcmc <- function(fit,
#                            y,
#                            s,
#                            X,
#                            min_range = 1,
#                            max_max_range = 100,
#                            range_sd = sqrt(2),
#                            as = 2,
#                            bs = 2,
#                            tauinv = NULL,
#                            errvar = NULL,
#                            sd_beta = 3,
#                            nugget_sd = 1,
#                            iters = 500,
#                            ...) {
#     iters_so_far <- nrow(fit$covar_params)
#
#     errvar <- fit$covar_params[iters_so_far, "epsilon"] ^ 2
#     init_range <- fit$covar_params[iters_so_far, "range_s"]
#     tauinv <- fit$covar_params[iters_so_far, "sigma"] ^ 2
#     init_beta <- fit$beta[nrow(fit$beta),]
#
#     run_nnsgp(
#         y,
#         s,
#         X,
#         num_neighbors = fit$neighbor_info$num_neighbors,
#         min_range = min_range,
#         max_max_range = max_max_range,
#         init_range = init_range,
#         range_sd = range_sd,
#         as = as,
#         bs = bs,
#         tauinv = tauinv,
#         errvar = errvar,
#         sd_beta = sd_beta,
#         init_beta = init_beta,
#         nugget_sd = nugget_sd,
#         iters = iters
#     )
# }
#
#
