get_prec_and_det <- function(d,
                             r,
                             range = NULL,
                             nu = NULL,
                             thresh = 1e-07) {
    # Compute the matern covariance function (r = 1 when signal effects)
    # cifelse returns a matrix with 1 where d = 0, and is 0 otherwise
    Q <-
        r * fields::Matern(d, range = range, smoothness = nu) + (1 - r) * cifelse(d)

    if (any(is.nan(Q))) {
        warning("nans in Q when calling get_prec_and_det.")
        return(NA)
    } else {
        # Get the eigenvalues/vectors
        eigvals <- cgeteigs(Q)

        # Adjustment for numerical stability (dont want, effectively, 1/0)
        D <- cifelse_eigvals(eigvals, thresh)

        # Inverse correlation and its log determinant
        return(list(
            "precision" = cinv(Q),
            "ldeterminant" = sum(log(D))
        ))
    }
}

# Treats as mean 0 here, and is added on to Xp %*% beta in the actual body of the code
make_pred <- function(y, matern_cov, tau) {
    np   <- nrow(matern_cov)
    cond_cov <- matern_cov / tau
    rowMeans(y) + t(RandomFieldsUtils::cholx(cond_cov)) %*% rnorm(np)
}

update_var <- function(cur_var, acpt_rt, opt_rt = .3, gamma1) {
    min(exp(log(cur_var) + gamma1 * (acpt_rt - opt_rt)), 25)
    #exp(log(cur_var) + gamma1 * (acpt_rt - opt_rt))
}

extend_mcmc <- function(fit,
                        y,
                        s,
                        X,
                        sp = NULL,
                        Xp = NULL,
                        corr_errs = FALSE,
                        cutoff = 0,
                        mean_nu = log(1.5),
                        sd_nu = 1,
                        mean_range = -1,
                        sd_range = 1,
                        as = 1.5,
                        bs = 1.5,
                        tauinv = NULL,
                        errvar = NULL,
                        sd_beta = 100,
                        iters = 500) {
    iters_so_far <- nrow(fit$covar_params)

    if (corr_errs) {
        init_nu_e <- fit$covar_params[iters_so_far, "nu_e"]
        init_range_e <- fit$covar_params[iters_so_far, "range_e"]
        init_r <- fit$covar_params[iters_so_far, "r"]
        errvar <- NULL
    } else {
        errvar <- fit$covar_params[iters_so_far, "err_sd"] ^ 2
        init_nu_e <- init_range_e <- init_r <- NULL
    }
    init_nu <- fit$covar_params[iters_so_far, "nu_s"]
    init_range <- fit$covar_params[iters_so_far, "range_s"]
    tauinv <- fit$covar_params[iters_so_far, "sigma"] ^ 2
    init_beta <- fit$beta[nrow(fit$beta),]

    run_sgp(
        y,
        s,
        X,
        sp = sp,
        Xp = Xp,
        corr_errs = corr_errs,
        cutoff = cutoff,
        mean_nu = mean_nu,
        sd_nu = sd_nu,
        init_nu_s = init_nu_s,
        init_nu_e = init_nu_e,
        mean_range = mean_range,
        sd_range = sd_range,
        init_range = init_range,
        init_range_e = init_range_e,
        init_r = init_r,
        as = as,
        bs = bs,
        tauinv = tauinv,
        errvar = errvar,
        sd_beta = sd_beta,
        init_beta = init_beta,
        iters = iters,
        burnin = 1
    )
}
