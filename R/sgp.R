run_sgp <- function(y,
                    s,
                    X,
                    sp = NULL,
                    Xp = NULL,
                    corr_errs = FALSE,
                    cutoff = 0,
                    mean_nu = log(1.5),
                    sd_nu = 1,
                    init_nu = NULL,
                    mean_range = -1,
                    sd_range = 1,
                    init_range = NULL,
                    init_r = NULL,
                    as = .01,
                    bs = .01,
                    tauinv = NULL,
                    errvar = NULL,
                    sd_beta = 100,
                    init_beta = NULL,
                    iters = 500,
                    burnin = 100) {
    if (corr_errs) {
        if (ncol(y) < 2) {
            stop(
                "Replicates needed to estimate correlated error process. Try corr_errs = FALSE."
            )
        } else {
            run_sgp_correlated_errs(
                y,
                s,
                X,
                sp,
                Xp,
                cutoff,
                mean_nu,
                sd_nu,
                init_nu,
                mean_range,
                sd_range,
                init_range,
                init_r,
                as,
                bs,
                tauinv,
                sd_beta,
                init_beta,
                iters,
                burnin
            )
        }
    } else {
        run_sgp_nugget(
            y,
            s,
            X,
            sp,
            Xp,
            cutoff,
            mean_nu,
            sd_nu,
            init_nu,
            mean_range,
            sd_range,
            init_range,
            as,
            bs,
            tauinv,
            errvar,
            sd_beta,
            init_beta,
            iters,
            burnin
        )
    }
}
