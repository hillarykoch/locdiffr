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

#------------------------------------------------------------
# Compute BRD2 enrichment
#------------------------------------------------------------

library(sgp)
library(tidyverse)

get_background_enrichment <- function(num_locs, loc_range, peak_locs, win_size, multicount = TRUE) {
    locs <- sample(loc_range, size = num_locs, replace = TRUE)
    
    # are there peaks nearby? sometimes more than 1 peak in a given loc,
    #   this is not counted here if multicount = FALSE
    if(multicount) {
        bypeak <- sapply(locs, function(X) between(peak_locs, X, X+win_size-1)) %>%
            colSums
    } else {
        bypeak <- sapply(locs, function(X) between(peak_locs, X, X+win_size-1)) %>%
            apply(2, function(X) ifelse(any(X), 1, 0))
    }
    mean(bypeak)
}

get_rejection_enrichment <- function(rej_locs, peak_locs, win_size, multicount = TRUE) {
    # are there peaks nearby? sometimes more than 1 peak in a given loc,
    #   this is not counted here if multicount = FALSE
    if(multicount) {
        bypeak <- sapply(rej_locs, function(X) between(peak_locs, X, X+win_size-1)) %>%
            colSums
    } else {
        bypeak <- sapply(rej_locs, function(X) between(peak_locs, X, X+win_size-1)) %>%
            apply(2, function(X) ifelse(any(X), 1, 0))
    }
    mean(bypeak)
}


get_enrichment <- function(rej_locs, num_locs, loc_range, peak_locs, win_size, multicount = TRUE) {
    num <- get_rejection_enrichment(rej_locs, peak_locs, win_size, multicount)
    denom <- replicate(1000, get_background_enrichment(num_locs, loc_range, peak_locs, win_size, multicount))
    out <- c(num/mean(denom), 1 - pnorm(num, mean = mean(denom), sd = sd(denom)))
    names(out) <- c("ratio", "p-value")
    out
}
