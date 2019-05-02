run_sgp_nugget <- function(y,
                           s,
                           X,
                           cutoff = 0,
                           # null is signal values below this, alternative is above
                           init_nu = 0.5,
                           # fixed value of the log matern smoothness parameter
                           # mean_range = 0,
                           # # prior mean of the log matern spatial range parameter
                           # sd_range = 1,
                           # Try a uniform prior on the spatial range (1, max(crd))
                           min_range = 1,
                           max_range = NULL,
                           # prior sd of the log matern spatial range parameter
                           init_range = NULL,
                           # init value of the log matern spatial range parameter
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
                           iters = 500,
                           burnin = 100) {
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
    if(is.null(max_range)) {
        max_range <- diff(range(s))
    }

    # These all default to NULL
    beta <- init_beta
    rhos <- init_range
    nus <- init_nu

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
        rhos  <- unname(quantile(d, 0.1))
    }

    # tau is precision of matern covariance function
    tau <- 1 / tauinv

    # errprec  is the precision of the nugget
    errprec <- 1 / errvar
    rm(lmfit)

    # Precision matrix and log determinant for the spatial and error processes
    PLDs <- get_prec_and_det(d, 1, rhos, nus)
    Xb   <- as.vector(X %*% beta)

    # chain of predicted values for y
    # If no prediction locations, just iters of zeroes that won't update
    ypred <- matrix(0, iters, npred)

    # store regression coefficients and covariance params
    beta_chain <- matrix(0, iters, p)
    param_chain <- matrix(0, iters, 3)
    colnames(param_chain) <-
        c("sigma", "range_s", "err_sd")


    #-----------------------------------------------------------------------
    # Set up adaptive tuning
    #-----------------------------------------------------------------------
    c0 <- 10
    c1 <- 0.8
    tune_k <- 2
    win_len <- min(iters, 50)
    acpt_rhos <- c(1, rep(NA, win_len - 1))
    tune_var <- 1
    acpt_rt_rhos <- 1
    acpt_chain <- rep(NA,iters)
    tune_var_chain <- rep(NA,iters)

    pb <- txtProgressBar(min = 1,
                         max = iters,
                         style = 3)
    for (i in 1:iters) {
        #-------------------------------------------------------------
        # Update regression coefficients
        #-------------------------------------------------------------
        premultX <-
            reps * tau * tX %*% (PLDs$precision + diag(errprec, n))
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
        SS <- tau * sum(apply(yminusXb, 2, function(X)
            emulator::quad.form(PLDs$precision, X)))
        tau <- rgamma(1, (n * reps) / 2 + as, SS / 2 + bs)
        SS  <- tau * SS

        #-----------------------------------------------------------------------
        # Sampling covariance parameters variable-at-a-time
        #-----------------------------------------------------------------------

        # First do rho_s
        PLDs_star <- NA
        while(is.na(PLDs_star[1])) {
            # This takes the log of rhos, does a normal random walk proposal,
            #   and then exponentiates back
            # rhos_star <- exp(log(rhos) + rnorm(1, 0, sqrt(tune_var)))

            # Try instead to propose on the current scale,
            #   and just reject if rhos is negative (hopefully more stable)
            rhos_star <- rhos + rnorm(1, 0, sqrt(tune_var))
            while(rhos_star <= 0) {
                rhos_star <- rhos + rnorm(1, 0, sqrt(tune_var))
            }

            PLDs_star <- get_prec_and_det(d, 1, rhos_star, nus)
        }
        SS_star <- tau * sum(apply(yminusXb, 2, function(X)
            emulator::quad.form(PLDs_star$precision, X)))

        # R <- dnorm(log(rhos_star), mean_range, sd_range, log = T) -
        #     dnorm(log(rhos), mean_range, sd_range, log = T) +
        #     0.5 * (PLDs_star$ldeterminant - PLDs$ldeterminant) -
        #     0.5 * (SS_star - SS)
        R <- dunif(rhos_star, min_range, max_range, log = T) -
            dunif(rhos, min_range, max_range, log = T) +
            0.5 * (PLDs_star$ldeterminant - PLDs$ldeterminant) -
            0.5 * (SS_star - SS)
        if (!is.na(exp(R))) {
            if (runif(1) < exp(R)) {
                rhos <- rhos_star
                PLDs <- PLDs_star
                acpt_rhos[(i+1) %% win_len] <- 1
            } else {
                acpt_rhos[(i+1) %% win_len] <- 0
            }
        }

        #-----------------------------------------------------------------------
        # Update the tuning variance
        #-----------------------------------------------------------------------
        if(i >= 100) {
            gamma1 <- c0 / (i + tune_k) ^ (c1)
            acpt_rt_rhos <- mean(acpt_rhos, na.rm = TRUE)
            tune_var <- update_var(tune_var, acpt_rt_rhos, .3, gamma1)
            acpt_chain[i] <- acpt_rt_rhos
            tune_var_chain[i] <- tune_var
        }

        #-------------------------------------------------------------
        # Update error variance term
        #-------------------------------------------------------------

        # errprec ~ Gamma
        yminusXb <- sweep(y, 1, Xb)
        SSe <-
            errprec * sum(apply(yminusXb, 2, function(Y)
                t(Y) %*% Y))
        errprec <- rgamma(1, (n * reps) / 2 + as, SSe / 2 + bs)

        beta_chain[i,] <- beta
        param_chain[i,] <-
            c(1 / sqrt(tau), rhos, 1 / sqrt(errprec))

        if (i >= burnin) {
            matern_cov <- fields::Matern(d, range = rhos, smoothness = nus)
            ypred[i,] <-
                X %*% beta + make_pred(y - Xb, matern_cov, tau)
        }

        setTxtProgressBar(pb, i)
    }
    close(pb)

    list("beta" = beta_chain,
         "covar_params" = param_chain,
         "pred" = ypred[burnin:iters,],
         "acpt_chain" = acpt_chain,
         "tune_chain" = tune_var_chain)
}
