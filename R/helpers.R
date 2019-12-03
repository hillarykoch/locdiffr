get_cor <- function(d, range = NULL) {
    # Compute the Exponential covariance function
    fields::Exponential(d, range = range)
}

get_log_det <- function(covmat, thresh = 1e-07) {
    eigvals <- cgeteigs(covmat)

    # Adjustment for numerical stability (dont want, effectively, 1/0)
    D <- cifelse_eigvals(eigvals, thresh)
    return(sum(log(D)))
}

# Treats as mean 0 here, and is added on to Xp %*% beta in the actual body of the code
make_pred <- function(y, exp_cov, tau) {
    np   <- nrow(exp_cov)
    cond_cov <- exp_cov / tau
    rowMeans(y) + t(RandomFieldsUtils::cholx(cond_cov)) %*% rnorm(np)
}

update_var <- function(cur_var, acpt_rt, opt_rt = .3, gamma1) {
    #min(exp(log(cur_var) + gamma1 * (acpt_rt - opt_rt)), 25)
    exp(log(cur_var) + gamma1 * (acpt_rt - opt_rt))
}

# Compute upper bound on the prior for the spatial range parameter
get_max_range <- function(dat, max_max_range) {
    if(is.null(dim(dat))) {
        autocor <- as.vector(acf(dat, plot = FALSE, lag.max = max_max_range)$acf)
    } else {
        autocor <- as.vector(acf(rowMeans(dat), plot = FALSE, lag.max = max_max_range)$acf)
    }

    idx <- sapply(1:(max_max_range-1), function(X) autocor[X] < autocor[X+1]) & autocor[1:(max_max_range-1)] < .05
    if(any(idx)) {
        which(idx)[1]
    } else {
        max_max_range
    }
}

extend_mcmc <- function(fit,
                        y,
                        s,
                        X,
                        min_range = 1,
                        max_max_range = 100,
                        range_sd = sqrt(2),
                        as = 2,
                        bs = 2,
                        tauinv = NULL,
                        errvar = NULL,
                        sd_beta = 10,
                        nugget_sd = 1,
                        iters = 500,
                        ...) {
    iters_so_far <- nrow(fit$covar_params)

    errvar <- fit$covar_params[iters_so_far, "epsilon"] ^ 2
    init_range <- fit$covar_params[iters_so_far, "range_s"]
    tauinv <- fit$covar_params[iters_so_far, "sigma"] ^ 2
    init_beta <- fit$beta[nrow(fit$beta),]

    run_sgp(
        y,
        s,
        X,
        min_range = min_range,
        max_max_range = max_max_range,
        init_range = init_range,
        range_sd = range_sd,
        as = as,
        bs = bs,
        tauinv = tauinv,
        errvar = errvar,
        sd_beta = sd_beta,
        init_beta = init_beta,
        nugget_sd = nugget_sd,
        iters = iters,
        burnin = 0
    )
}
