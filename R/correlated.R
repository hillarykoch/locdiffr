run_sgp_correlated_errs <- function(y,
                        s,
                        X,
                        cutoff = 0, # null is signal values below this, alternative is above
                        mean_nu = log(1.5), # prior mean of the log matern smoothness parameter
                        sd_nu = 1, # prior sd of the log matern spatial range parameter
                        init_nu_s = NULL, # initial value of the log matern smoothness parameter
                        init_nu_e = NULL, # initial value of the log matern smoothness parameter
                        mean_range = -1, # prior mean of the log matern spatial range parameter
                        sd_range = 1, # prior sd of the log matern spatial range parameter
                        init_range_s = NULL, # init value of the log matern spatial range parameter
                        init_range_e = NULL, # init value of the log matern spatial range parameter
                        init_r = NULL, # how much of the error process is correlated vs iid noise
                        sd_r = 3, # standard deviation on prior for log(r/(1-r))
                        as = 5,  # the prior for the variance is InvG(as,bs)
                        bs = 2,
                        tauinv = NULL, # initial value of the variance term in matern cov
                        sd_beta = 10, # regression coefficients have N(0, sd_beta) priors
                        init_beta = NULL, # initial value of regression coefs
                        iters = 500,
                        burnin = 100) {
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

    # These all default to NULL
    beta <- init_beta
    r <- init_r
    rhos <- init_range_s
    rhoe <- init_range_e
    nus <- init_nu_s
    nue <- init_nu_e

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
        rhos <- unname(quantile(d, 0.1))
    }
    if (is.null(rhoe)) {
        rhoe <- unname(quantile(d, 0.1))
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
    param_chain <- matrix(0, iters, 6)
    colnames(param_chain) <-
        c("sigma", "r", "range_s", "nu_s", "range_e", "nu_e")

    # matern_proposal cov is 2 x 2 matrix with diagonal 1 and off-diagonal -0.5
    # This is for updating hyperparameters nu and rho in the MCMC
    matern_proposal_cov <- diag(2) * 1.5 - 0.5
    proposal_chol <- 0.1 * t(RandomFieldsUtils::cholx(matern_proposal_cov))

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

        # Jointly propose new range and smoothness parameters
        # Accept or reject based on R (log likelihood ratio)
        # The proposal takes the log of the range/smoothness (rho and nu both > 0),
        #   then proposes an update from a normal distribution, then exponentiates
        #   back. Actually they range and smoothness are proposed together, and are
        #   the proposal is from a MVNormal where range and smoothness are
        #   anti-correlated

        PLDs_star <- NA
        while(is.na(PLDs_star[1])) {
            epsilon <- proposal_chol %*% rnorm(2)
            rhos_star <- exp(log(rhos) + epsilon[1])
            nus_star <- exp(log(nus) + epsilon[2])
            PLDs_star <- get_prec_and_det(d, 1, rhos_star, nus_star)
        }

        # And the sum of squares star
        SS_star <- tau * sum(apply(yminusXb, 2, function(X)
                            emulator::quad.form(PLDs_star$precision, X)))

        # The acceptance ratio is the log proposal nu over the log current nu
        #   evaluated at the prior for nu (which is normal)
        # times the ratio of the log proposal rho over the log current rho
        #   evaluated at the prior for rho (which is also normal)
        # times the ratio of the data likelihoods
        # This is a symmetric proposal
        # NOTE: this can be adaptively tuned
        R <- dnorm(log(nus_star), mean_nu, sd_nu, log = T) -
            dnorm(log(nus), mean_nu, sd_nu, log = T) +
            dnorm(log(rhos_star), mean_range, sd_range, log = T) -
            dnorm(log(rhos), mean_range, sd_range, log = T) +
            0.5 * (PLDs_star$ldeterminant - PLDs$ldeterminant) -
            0.5 * (SS_star - SS)
        if (!is.na(exp(R))) {
            if (runif(1) < exp(R)) {
                nus <- nus_star
                rhos <- rhos_star
                PLDs <- PLDs_star
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
        lr     <- log(r / (1 - r))

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
            epsilon <- proposal_chol %*% rnorm(2)
            rhoe_star <- exp(log(rhoe) + epsilon[1])
            nue_star <- exp(log(nue) + epsilon[2])
            PLDe_star <- get_prec_and_det(d, r, rhoe_star, nue_star)
        }

        SS_star <- sum(apply(yminusXb, 2, function(X)
            emulator::quad.form(PLDe_star$precision, X)))
        R <- dnorm(log(nue_star), mean_nu, sd_nu, log = T) -
            dnorm(log(nue), mean_nu, sd_nu, log = T) +
            dnorm(log(rhoe_star), mean_range, sd_range, log = T) -
            dnorm(log(rhoe), mean_range, sd_range, log = T) +
            0.5 * (PLDe_star$ldeterminant - PLDe$ldeterminant) -
            0.5 * (SS_star - SS)
        if (!is.na(exp(R))) {
            if (runif(1) < exp(R)) {
                nue <- nue_star
                rhoe <- rhoe_star
                PLDe <- PLDe_star
            }
        }

        beta_chain[i, ] <- beta
        param_chain[i, ] <- c(1 / sqrt(tau), r, rhos, nus, rhoe, nue)

        if (i >= burnin) {
            matern_cov <- fields::Matern(d, range = rhos, smoothness = nus)
            ypred[i, ] <- X %*% beta + make_pred(y - Xb, matern_cov, tau)
        }

        setTxtProgressBar(pb, i)
    }
    close(pb)

    list("beta" = beta_chain, "covar_params" = param_chain, "pred" = ypred[burnin:iters, ])
}
