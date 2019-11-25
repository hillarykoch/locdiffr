run_sgp <- function(y,
                    s,
                    X,
                    min_range = 1,
                    max_max_range = NULL,
                    # prior sd of the log matern spatial range parameter
                    init_range = 5,
                    # init value of the log matern spatial range parameter
                    range_sd = sqrt(2),
                    # prior sd for the spatial range
                    as = 2,
                    # the prior for the variance is InvG(as,bs)
                    bs = 2,
                    tauinv = NULL,
                    # initial value of the variance term in matern cov
                    errvar = NULL,
                    # nugget variance
                    sd_beta = 3,
                    # regression coefficients have N(0, sd_beta) priors
                    init_beta = NULL,
                    # initial value of regression coefs
                    nugget_sd = 1,
                    # prior sd for white noise process factor
                    iters = 500) {
    # number of observation locations, number of process replicates (sample size)
    n <- nrow(y)
    reps <- ncol(y)

    # number of predictors
    p <- ncol(X)

    d <- cpwdist(s, s)
    npred <- nrow(s)

    tX <- t(X)
    precision_beta <- diag(p) / (sd_beta ^ 2)


    #----------------------------------------------------
    # Initial values
    #----------------------------------------------------
    if(is.null(max_max_range)) {
        max_max_range <- diff(range(s))
    }
    max_range <- get_max_range(y, max_max_range)

    # These all default to NULL
    beta <- init_beta
    rhos <- init_range

    # Get initial regression values from a simple linear model fit
    if (reps > 1) {
        lmfit <- lm(rowMeans(y) ~ X - 1)
    } else {
        lmfit <- lm(y ~ X - 1)
    }
    if (is.null(beta)) {
        beta <- lmfit$coef
    }
    if (is.null(tauinv)) {
        tauinv <- var(lmfit$res)
    }
    if (is.null(errvar)) {
        errvar <- var(lmfit$res)
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

    # cor cur is the exponential correlation plus the nugget (not quite a correlation actually)
    cor_cur <- get_cor(d, rhos) + diag(errvar, n)
    Xb   <- as.vector(X %*% beta)

    # chain of predicted values for y
    # If no prediction locations, just iters of zeroes that won't update
    ypred <- matrix(0, iters, npred)

    # store regression coefficients and covariance params
    beta_chain <- matrix(0, iters, p)
    param_chain <- matrix(0, iters, 3)
    colnames(param_chain) <-
        c("sigma", "range_s", "epsilon")

    pb <- txtProgressBar(min = 1,
                         max = iters,
                         style = 3)
    for (i in 1:iters) {
        #-------------------------------------------------------------
        # Update regression coefficients
        #-------------------------------------------------------------
        inv_cor_cur <- cinv(cor_cur)
        premultX <- reps * tau * tX %*% inv_cor_cur
        varterm <- cinv(precision_beta + premultX %*% X)
        muterm <- premultX %*% rowMeans(y)

        # Sample new betas from this full conditional
        beta <-
            varterm %*% muterm + t(RandomFieldsUtils::cholx(varterm)) %*% rnorm(p)
        Xb   <- as.vector(X %*% beta)

        #-------------------------------------------------------------
        # Update covariance parameters for the spatial signal process
        #-------------------------------------------------------------
        # tau ~ Gamma
        yminusXb <- sweep(y, 1, Xb)
        SS <- sum(apply(yminusXb, 2, function(X) emulator::quad.form(inv_cor_cur, X)))
        tau <- rgamma(1, (n * reps) / 2 + as, SS / 2 + bs)
        tauinv <- 1/tau

        #-----------------------------------------------------------------------
        # Sampling covariance parameters variable-at-a-time
        #-----------------------------------------------------------------------
        # Spatial range
        lb <- pnorm(-(rhos-min_range), mean = 0, sd = sqrt(tune_var))
        ub <- pnorm(max_range-rhos, mean = 0, sd = sqrt(tune_var))
        uuu <- runif(1, lb, ub)
        rhos_star <- rhos + qnorm(uuu, mean = 0, sd = sqrt(tune_var))

        cov_cur <- tauinv * cor_cur
        invcov_cur <- cinv(cov_cur)
        
        cor_cur_star <- get_cor(d, rhos_star) + diag(errvar, n)
        cov_star <- tauinv * cor_cur_star
        invcov_star <- cinv(cov_star)
        
        # ldet_cur <- get_log_det(cov_cur)
        # ldet_star <- get_log_det(cov_star)
        ldet_cur <- log(det(cov_cur))
        ldet_star <- log(det(cov_star))
        
        SS <- tau * SS
        SS_star <- sum(apply(yminusXb, 2, function(X) emulator::quad.form(invcov_star, X)))
        R <- (reps/2)*(ldet_cur - ldet_star) + 0.5 * (SS - SS_star)

        if (runif(1) < exp(R)) {
            rhos <- rhos_star
            cor_cur <- cor_cur_star
            cov_cur <- cov_star
            invcov_cur <- invcov_star
        }

        #-------------------------------------------------------------
        # Update error variance term
        #-------------------------------------------------------------
        # errprec ~ Gamma
        lb <- pnorm(-errprec, mean = 0, sd = sqrt(tune_var_epsilon))
        uuu <- runif(1, lb, 1)
        errprec_star <- errprec + qnorm(uuu, mean = 0, sd = sqrt(tune_var_epsilon))
        errvar_star <- 1/errprec_star

        cor_cur_star <- get_cor(d, rhos) + diag(errvar_star, n)
        cov_star <- tauinv * cor_cur_star
        invcov_star <- cinv(cov_star)
        
        R <- sum(dmvnorm(x = t(Y), mean = Xb, sigma = cov_star, log = TRUE)) -
            sum(dmvnorm(x = t(Y), mean = Xb, sigma = cov_cur, log = TRUE)) +
            dgamma(errprec_star, as, bs) -
            dgamma(errprec, as, bs)
        
        if (runif(1) < exp(R)) {
            errvar <- errvar_star
            errprec <- errprec_star
            cor_cur <- cor_cur_star
            # acpt_epsilon[(i+1) %% win_len] <- 1
        }
        
        #-----------------------------------------------------------------------
        # Update the tuning variance
        #-----------------------------------------------------------------------
        # if(i >= win_len) {
        #     gamma1 <- c0 / (i + tune_k) ^ (c1)
        #     acpt_rt_epsilon <- mean(acpt_epsilon, na.rm = TRUE)
        #     tune_var_epsilon <- update_var(tune_var_epsilon, acpt_rt_epsilon, .3, gamma1)
        #     acpt_chain[i-49] <- acpt_rt_epsilon
        #     tune_var_chain[i-49] <- tune_var_epsilon
        # }
        
        beta_chain[i,] <- beta
        param_chain[i,] <-
            c(1 / sqrt(tau), rhos, sqrt(errvar))

        #-----------------------------------------------------------------------
        # Draw a process from the current model
        #   might want to break this out, and do it on the MCMC after the fact
        #   That way, won't need to make predictions for every iteration and
        #   store all this info as the MCMC runs
        #-----------------------------------------------------------------------
        matern_cov <- fields::Exponential(d, range = rhos)
        ypred[i,] <-
            X %*% beta + make_pred(y - Xb, matern_cov, tau)

        setTxtProgressBar(pb, i)
    }
    close(pb)

    list("beta" = beta_chain,
         "covar_params" = param_chain,
         "pred" = ypred)#,
         # "acpt_chain" = acpt_chain,
         # "tune_chain" = tune_var_chain)
}
