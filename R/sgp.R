run_sgp <- function(y,
                    s,
                    X,
                    corr_errs = FALSE,
                    cutoff = 0,
                    mean_nu = 0,
                    sd_nu = 1,
                    init_nus = 0.5,
                    init_nue = NULL,
                    min_range = 1,
                    #max_range = NULL,
                    max_max_range = 100,
                    init_range_s = 5,
                    init_range_e = NULL,
                    init_r = NULL,
                    sd_r = 3,
                    as = 2,
                    bs = 2,
                    tauinv = NULL,
                    errvar = NULL,
                    sd_beta = 10,
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
                cutoff,
                init_nus,
                init_nue,
                min_range,
                #max_range,
                max_max_range,
                init_range_s,
                init_range_e,
                init_r,
                sd_r,
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
            cutoff,
            init_nus,
            min_range,
            #max_range,
            max_max_range,
            init_range_s,
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
