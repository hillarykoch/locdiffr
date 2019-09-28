run_sgp_correlated_errs <- function(y,
                        s,
                        X,
                        init_nus = 0.5,
                        init_nue = 0.5,
                        min_range = 1, # prior mean of the log matern spatial range parameter
                        max_max_range = NULL, # For computing the upper bound on the prior on the spatial range parameter
                        init_range_s = NULL, # init value of the log matern spatial range parameter
                        init_range_e = NULL, # init value of the log matern spatial range parameter
                        init_r = NULL, # how much of the error process is correlated vs iid noise
                        sd_r = 3, # standard deviation on prior for log(r/(1-r))
                        as = 2,  # the prior for the variance is InvG(as,bs)
                        bs = 2,
                        tauinv = NULL, # initial value of the variance term in matern cov
                        sd_beta = 3, # regression coefficients have N(0, sd_beta) priors
                        init_beta = NULL, # initial value of regression coefs
                        cpp = TRUE,
                        iters = 500,
                        burnin = 1) {
    # number of observation locations, number of process replicates (sample size)
    n <- nrow(y)
    reps <- ncol(y)

    # number of predictors
    p <- ncol(X)

    d <- cpwdist(s,s)
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
    r <- init_r
    rhos <- init_range_s
    rhoe <- init_range_e
    nus <- init_nus
    nue <- init_nue

    # Get initial regression values from a simple linear model fit
    if(reps > 1) {
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

    # init value of the log matern spatial range parameter
    if (is.null(rhos)) {
        rhos <- (max_range - 1)/2#unname(quantile(d, 0.1))
    }
    if (is.null(rhoe)) {
        rhoe <- (max_range - 1)/2#unname(quantile(d, 0.1))
    }

    # initialize with 80% of errors being correlated as opposed to the 20% iid noise
    if (is.null(r)) {
        r <- 0.8
    }

    # initial value of the log matern smoothness parameter
    if (is.null(nus)) {
        nus <- 0.5
    }

    if (is.null(nue)) {
        nue <- 0.5
    }

    # tau is precision of matern covariance function
    tau <- 1 / tauinv
    rm(lmfit)

    # Precision matrix and log determinant for the spatial and error processes
    PLDs <- get_prec_and_det(d, 1, rhos, nus)
    PLDe <- get_prec_and_det(d, r, rhoe, nue)
    Xb   <- as.vector(X %*% beta)

    # chain of predicted values for y
    # If no prediction locations, just iters of zeroes that won't update
    ypred <- matrix(0, iters, npred)

    # store regression coefficients and covariance params
    beta_chain <- matrix(0, iters, p)
    param_chain <- matrix(0, iters, 4)
    colnames(param_chain) <-
        c("sigma", "r", "range_s", "range_e")

    #-----------------------------------------------------------------------
    # Set up adaptive tuning
    #-----------------------------------------------------------------------
    c0 <- 10
    c1 <- 0.8
    tune_k <- 2
    win_len <- min(iters, 50)
    acpt_rhos <- acpt_rhoe <- c(1, rep(NA, win_len - 1))
    tune_vars <- tune_vare <- 1
    acpt_rt_rhos <- acpt_rt_rhoe <- 1
    acpt_chain <- tune_var_chain <- matrix(NA, nrow = iters - 99, ncol = 2)
    colnames(acpt_chain) <- colnames(tune_var_chain) <- c("rhos", "rhoe")

    if(cpp) {
        out <- crun_sgp_correlated(reps, n, p, tau, X, PLDs$precision, PLDs$ldeterminant,
                            PLDe$precision, PLDe$ldeterminant, errprec, precision_beta,
                            y, iters, as, bs, rhos, nus, rhoe, nue, r, sd_r, tune_vars,
                            tune_vare, min_range, max_range, win_len, d, acpt_rhos,
                            acpt_rhoe, c0, c1, tune_k, acpt_chain, tune_var_chain)
        colnames(out$covar_params) <- c("sigma", "r", "range_s", "range_e")
        colnames(out$acpt_chain) <- colnames(out$tune_chain) <- c("rhos", "rhoe")
        return(out)
    } else {
        pb <- txtProgressBar(min = 1, max = iters, style = 3)
        for (i in 1:iters) {
            #-------------------------------------------------------------
            # Update regression coefficients
            #-------------------------------------------------------------
            premultX <- reps * tau * tX %*% (PLDs$precision + PLDe$precision)
            varterm <- cinv(precision_beta + premultX %*% X)
            muterm <- premultX %*% rowMeans(y)

            # Sample new betas from this full conditional
            beta <- varterm %*% muterm + t(RandomFieldsUtils::cholx(varterm)) %*% rnorm(p)
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
            # Sampling range parameter
            #-----------------------------------------------------------------------
            PLDs_star <- NA
            while(is.na(PLDs_star[1])) {
                rhos_star <- rhos + rnorm(1, 0, sqrt(tune_vars))
                while(rhos_star < 1 | rhos_star > max_range) {
                    rhos_star <- rhos + rnorm(1, 0, sqrt(tune_vars))
                }
                PLDs_star <- get_prec_and_det(d, 1, rhos_star, nus)
            }

            # And the sum of squares star
            SS_star <- tau * sum(apply(yminusXb, 2, function(X)
                emulator::quad.form(PLDs_star$precision, X)))

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

            #-------------------------------------------------------------
            # Update covariance parameters for the error process
            #-------------------------------------------------------------

            # Propose the proportion of error which comes from the correlated
            #   error process, vs. iid noise

            # log(r/(1-r)) is in (-Inf, Inf), and we propose from here using
            #   a normal distribution whose result can be transformed back
            SS <- sum(apply(yminusXb, 2, function(X)
                emulator::quad.form(PLDe$precision, X)))
            lr <- log(r / (1 - r))

            PLDe_star <- NA
            while(is.na(PLDe_star[1])) {
                lr_star  <- rnorm(1, lr, 0.5)
                r_star  <- exp(lr_star) / (1 + exp(lr_star))
                PLDe_star <- get_prec_and_det(d, r_star, rhoe, nue)
            }

            SS_star <- sum(apply(yminusXb, 2, function(X)
                emulator::quad.form(PLDe_star$precision, X)))

            # The candidate lr and new zlr_star are evaluated at the prior on lr
            #   and then again the Multivariate normal likelihood
            R <- dnorm(lr_star, sd = sd_r, log = T) -
                dnorm(lr,  sd = sd_r, log = T) +
                0.5 * (PLDe_star$ldeterminant - PLDe$ldeterminant) -
                0.5 * (SS_star - SS)
            if (!is.na(exp(R))) {
                if (runif(1) < exp(R)) {
                    r <- r_star
                    PLDe <- PLDe_star
                    SS <- SS_star
                }
            }

            # Jointly propose new range and smoothness parameters for the error process
            # This procedure is exactly the same as the one used when updating
            #   nu and rho for the spatial signal Matern covariance

            # Sometimes the proposal is bad -- this prevents the program from dying
            PLDe_star <- NA
            while(is.na(PLDe_star[1])) {
                rhoe_star <- rhoe + rnorm(1, 0, sqrt(tune_vare))
                while(rhoe_star < 1 | rhoe_star > max_range) {
                    rhoe_star <- rhoe + rnorm(1, 0, sqrt(tune_vare))
                }

                PLDe_star <- get_prec_and_det(d, r, rhoe_star, nue)
            }

            SS_star <- sum(apply(yminusXb, 2, function(X)
                emulator::quad.form(PLDe_star$precision, X)))
            R <- dunif(rhoe_star, min_range, max_range, log = T) -
                dunif(rhoe, min_range, max_range, log = T) +
                0.5 * (PLDe_star$ldeterminant - PLDe$ldeterminant) -
                0.5 * (SS_star - SS)

            if (!is.na(exp(R))) {
                if (runif(1) < exp(R)) {
                    rhoe <- rhoe_star
                    PLDe <- PLDe_star
                    acpt_rhoe[(i+1) %% win_len] <- 1
                } else {
                    acpt_rhoe[(i+1) %% win_len] <- 0
                }
            }

            #-----------------------------------------------------------------------
            # Update the tuning variance
            #-----------------------------------------------------------------------
            if(i >= 100) {
                gamma1 <- c0 / (i + tune_k) ^ (c1)
                acpt_rt_rhos <- mean(acpt_rhos, na.rm = TRUE)
                acpt_rt_rhoe <- mean(acpt_rhoe, na.rm = TRUE)
                tune_vars <- update_var(tune_vars, acpt_rt_rhos, .3, gamma1)
                tune_vare <- update_var(tune_vare, acpt_rt_rhoe, .3, gamma1)
                acpt_chain[i-99,"rhos"] <- acpt_rt_rhos
                acpt_chain[i-99,"rhoe"] <- acpt_rt_rhoe
                tune_var_chain[i-99, "rhos"] <- tune_vars
                tune_var_chain[i-99, "rhoe"] <- tune_vare
            }


            beta_chain[i, ] <- beta
            param_chain[i, ] <- c(1 / sqrt(tau), r, rhos, rhoe)

            if (i >= burnin) {
                matern_cov <- fields::Matern(d, range = rhos, smoothness = nus)
                ypred[i, ] <- X %*% beta + make_pred(y - Xb, matern_cov, tau)
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
}
