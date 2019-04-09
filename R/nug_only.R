run_sgp_nugget <- function(y,
                        s,
                        X,
                        sp = NULL, # prediction locations
                        Xp = NULL,
                        cutoff = 0, # null is signal values below this, alternative is above
                        mean_nu = log(1.5), # prior mean of the log matern smoothness parameter
                        sd_nu = 1, # prior sd of the log matern spatial range parameter
                        init_nu = NULL, # initial value of the log matern smoothness parameter
                        mean_range = -1, # prior mean of the log matern spatial range parameter
                        sd_range = 1, # prior sd of the log matern spatial range parameter
                        init_range = NULL, # init value of the log matern spatial range parameter
                        as = .01,  # the prior for the variance is InvG(as,bs)
                        bs = .01,
                        tauinv = NULL, # initial value of the variance term in matern cov
                        errvar = NULL, # nugget variance
                        sd_beta = 100, # regression coefficients have N(0, sd_beta) priors
                        init_beta = NULL, # initial value of regression coefs
                        iters = 500,
                        burnin = 100) {
    # number of observation locations, number of process replicates (sample size)
    n <- nrow(y)
    reps <- ncol(y)

    # number of predictors
    p <- ncol(X)

    d <- cpwdist(s,s)
    predictions <- !is.null(sp) & !is.null(Xp)
    npred <- 1
    if (predictions) {
        npred <- nrow(sp)
        d12 <- cpwdist(sp, s)
        d11 <- cpwdist(sp, sp)
    }
    tX <- t(X)
    precision_beta <- diag(p) / (sd_beta ^ 2)


    #----------------------------------------------------
    # Initial values
    #----------------------------------------------------

    # These all default to NULL
    beta <- init_beta
    rhos <- init_range
    nus <- init_nu

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
    if(is.null(errvar)) {
        errvar <- var(lmfit$res)
    }

    # init value of the log matern spatial range parameter
    if (is.null(rhos)) {
        rhos  <- unname(quantile(d, 0.1))
    }

    # initial value of the log matern smoothness parameter
    if (is.null(nus)) {
        nus <- 0.5
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
    param_chain <- matrix(0, iters, 4)
    colnames(param_chain) <-
        c("sigma", "range_s", "nu_s", "err_sd")

    # matern_proposal cov is 2 x 2 matrix with diagonal 1 and off-diagonal -0.5
    # This is for updating hyperparameters nu and rho in the MCMC
    matern_proposal_cov <- diag(2) * 1.5 - 0.5
    proposal_chol <- 0.1 * t(RandomFieldsUtils::cholx(matern_proposal_cov))

    pb <- txtProgressBar(min = 1, max = iters, style = 3)
    for (i in 1:iters) {
        #-------------------------------------------------------------
        # Update regression coefficients
        #-------------------------------------------------------------
        premultX <- reps * tau * tX %*% (PLDs$precision + diag(errprec, n))
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
        epsilon <- proposal_chol %*% rnorm(2)
        rhos_star <- exp(log(rhos) + epsilon[1])
        nus_star <- exp(log(nus) + epsilon[2])

        # Obtain the precision of the matern covariance based on the proposed parameters
        PLDs_star <- get_prec_and_det(d, 1, rhos_star, nus_star)

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
        # Update error variance term
        #-------------------------------------------------------------

        # errprec ~ Gamma
        yminusXb <- sweep(y, 1, Xb)
        SSe <- errprec * sum(apply(yminusXb, 2, function(Y) t(Y) %*% Y))
        errprec <- rgamma(1, (n * reps) / 2 + as, SSe / 2 + bs)

        beta_chain[i, ] <- beta
        param_chain[i, ] <- c(1 / sqrt(tau), rhos, nus, 1 / sqrt(errprec))

        if (i >= burnin & predictions) {
            S11 <- fields::Matern(d11, range = rhos, smoothness = nus) / tau
            S12 <- fields::Matern(d12, range = rhos, smoothness = nus)
            S22inv <- tau * PLDs$precision

            ypred[i, ] <- Xp %*% beta + proj(y - Xb, S12, S11, S22inv)
        }

        setTxtProgressBar(pb, i)
    }
    close(pb)

    list("beta" = beta_chain, "covar_params" = param_chain, "pred" = ypred[burnin:iters, ])
}
