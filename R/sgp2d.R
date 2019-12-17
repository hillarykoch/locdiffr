run_2dsgp <-
    function(y,
             # list of observations, where observations are matrices
             S,
             # Matrix of spatial locations
             X,
             # List of basis matrices
             min_range = 1,
             max_max_range = NULL,
             # prior sd of the log matern spatial range parameter
             init_range = 5,
             # init value of the log matern spatial range parameter
             range_sd = sqrt(2),
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
             iters = 500) {
        # number of observation locations, number of process replicates (sample size)
        n <- length(y[[1]])#prod(dim(y[[1]]))
        reps <- length(y)
        
        # number of predictors (number of basis functions)
        p <- length(X)
        
        d <- fields::rdist(S)
        # npred <- nrow(S)
        precision_beta <- diag(p) / (sd_beta ^ 2)
        
        
        #----------------------------------------------------
        # Initial values
        #----------------------------------------------------
        if (is.null(max_max_range)) {
            max_max_range <- diff(range(S))
        }
        max_range <-
            25# max(apply(Reduce(`+`, y), 2, function(xx) get_max_range(xx, max_max_range)))
        
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
            errvar <- var(lmfit$res)
        }
        
        # init value of the log matern spatial range parameter
        if (is.null(rhos)) {
            rhos  <- (max_range - 1) / 2
        }
        
        # tau is precision of matern covariance function
        tau <- 1 / tauinv
        
        # errprec  is the precision of the nugget
        errprec <- 1 / errvar
        rm(lmfit)
        
        tune_var_epsilon <- nugget_sd ^ 2 # proposal variance for nugget process
        tune_var <- range_sd ^ 2 # proposal variance for spatial range
        
        # cor cur is the exponential correlation plus the nugget (not quite a correlation actually)
        cor_cur <- get_cor(d, rhos) + diag(errvar, n)
        Xb   <- as.vector(X_unfolded %*% beta)
        tX_unfolded <- t(X_unfolded)
        
        # chain of predicted values for y
        # If no prediction locations, just iters of zeroes that won't update
        ypred <- matrix(0, nrow = iters, ncol = length(y[[1]]))
        
        # store regression coefficients and covariance params
        beta_chain <- matrix(0, iters, p)
        param_chain <- matrix(0, iters, 3)
        colnames(param_chain) <- c("sigma", "range_s", "epsilon")
        
        pb <- txtProgressBar(min = 1,
                             max = iters,
                             style = 3)
        for (i in 1:iters) {
            #-------------------------------------------------------------
            # Update regression coefficients
            #-------------------------------------------------------------
            inv_cor_cur <- cinv(cor_cur)
            premultX <- reps * tau * tX_unfolded %*% inv_cor_cur
            varterm <- cinv(precision_beta + premultX %*% X_unfolded)
            muterm <- premultX %*% Y_unfolded
            
            # Sample new betas from this full conditional
            beta <-
                varterm %*% muterm + t(RandomFieldsUtils::cholx(varterm)) %*% rnorm(p)
            Xb   <- as.vector(X_unfolded %*% beta)
            
            #-------------------------------------------------------------
            # Update covariance parameters for the spatial signal process
            #-------------------------------------------------------------
            # tau ~ Gamma
            yminusXb <- sweep(Y_bound, 1, Xb)
            SS <- sum(apply(yminusXb, 2, function(X) emulator::quad.form(inv_cor_cur, X)))
            tau <- rgamma(1, (n * reps) / 2 + as, SS / 2 + bs)
            tauinv <- 1/tau
            
            #-----------------------------------------------------------------------
            # Sampling covariance parameters variable-at-a-time
            #-----------------------------------------------------------------------
            # Spatial range
            lb <- pnorm(-(rhos - min_range), mean = 0, sd = sqrt(tune_var))
            ub <- pnorm(max_range - rhos, mean = 0, sd = sqrt(tune_var))
            uuu <- runif(1, lb, ub)
            rhos_star <- rhos + qnorm(uuu, mean = 0, sd = sqrt(tune_var))
            
            cov_cur <- tauinv * cor_cur
            invcov_cur <- cinv(cov_cur)
            
            cor_cur_star <- get_cor(d, rhos_star) + diag(errvar, n)
            cov_star <- tauinv * cor_cur_star
            invcov_star <- cinv(cov_star)
            ldet_cur <- abs(get_log_det(cov_cur))
            ldet_star <- abs(get_log_det(cov_star))
            
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
            
            R <- sum(dmvnorm(x = t(Y_bound), mean = Xb, sigma = cov_star, log = TRUE)) -
                sum(dmvnorm(x = t(Y_bound), mean = Xb, sigma = cov_cur, log = TRUE)) +
                dgamma(errprec_star, as, bs) -
                dgamma(errprec, as, bs)
            
            if (runif(1) < exp(R)) {
                errvar <- errvar_star
                errprec <- errprec_star
                cor_cur <- cor_cur_star
            }
            
            beta_chain[i, ] <- beta
            param_chain[i, ] <-
                c(1 / sqrt(tau), rhos, sqrt(errvar))
            
            matern_cov <- fields::Exponential(d, range = rhos)
            ypred[i, ] <-
                X_unfolded %*% beta + make_pred(Y_unfolded - Xb, matern_cov, tau)
            
            setTxtProgressBar(pb, i)
        }
        close(pb)
        
        list("beta" = beta_chain,
             "covar_params" = param_chain,
             "pred" = ypred)
    }
