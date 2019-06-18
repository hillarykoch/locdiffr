get_rej_locs <- function(fit,
                         sccs,
                         h,
                         nb,
                         method = c("FDR", "FDX", "rank"),
                         burnin = 0,
                         alpha = 0.1,
                         resolution = 40000,
                         num_rejections = NULL) {
    #--------------------------------------------------
    # Filter unneeded windows
    if (length(sccs) > 1) {
        keepers <- Reduce(intersect, purrr::map(sccs, "crd"))
        for (i in seq_along(sccs)) {
            sccs[[i]] <- sccs[[i]][sccs[[i]]$crd %in% keepers,]
        }
    }
    
    #--------------------------------------------------
    # Reproduce basis functions
    s <- as.matrix(sccs$z1$crd, ncol = 1)
    X <-
        fda::getbasismatrix(s, fda::create.bspline.basis(range(s), nbasis = nb, norder = 4))
    
    #--------------------------------------------------
    # Compute indicator theta for testing
    theta <- matrix(NA, nrow(fit$pred) - burnin, ncol(fit$pred))
    for (i in 1:nrow(theta)) {
        theta[i,] <- fit$pred[i + burnin,] < X %*% fit$beta[i + burnin,]
    }
    
    if(method == "FDR") {
        rej_loc <- FDR(theta, alpha = alpha, nthresh = 100)$reject
    } else if(method == "FDX") {
        rej_loc <- FDX(theta,
            alpha = alpha,
            beta = 1 - alpha,
            nthresh = 100)$reject
    } else {
        rej_loc <- dplyr::between(rank(-colMeans(theta)), 1, num_rejections)
    }
    rej_loc
}
