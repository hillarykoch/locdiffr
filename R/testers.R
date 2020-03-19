FDR <- function(theta,
                alpha = .1,
                nthresh = 100) {
    # Theta is an indicator -- was the predicted value as location s
    #   less than X * beta at location s?
    # reject constructed such that E(mean(theta[reject])) < alpha
    
    rej_prob <- colMeans(theta)
    thresh <- seq(0, 1, length = nthresh)
    
    # The Bayes FDR at each threshold
    # (want to find one that is closest to and not greater than choice of alpha)
    BFDR <- rep(0, nthresh)
    for (j in 1:nthresh) {
        if (sum(rej_prob >= thresh[j]) > 0) {
            BFDR[j] <- 1 - mean(rej_prob[rej_prob >= thresh[j]])
        }
    }
    level <- min(thresh[BFDR < alpha])
    reject <- rej_prob > level
    list(
        level = level,
        reject = reject,
        thresh = thresh,
        BFDR = BFDR
    )
}


FDX <- function(theta,
                alpha = .1,
                beta = .9,
                nthresh = 100) {
    # reject constructed such that P(mean(theta[reject]) < alpha) < beta
    
    rej_prob <- colMeans(theta)
    thresh <- seq(0, max(rej_prob), length = nthresh)
    BFDX <- rep(0, nthresh)
    for (j in 1:nthresh) {
        xceeds <- rej_prob >= thresh[j]
        if (sum(xceeds) > 1) {
            BFDX[j] <- mean(1 - rowMeans(theta[, xceeds]) < alpha)
        }
    }
    
    # The level is the minimum value such that the BFDX exceeds bound beta
    level <- 1
    if (sum(BFDX > beta) > 0) {
        level <- min(thresh[BFDX > beta])
    }
    reject <- rej_prob > level
    list(
        level = level,
        reject = reject,
        thresh = thresh,
        BFDX = BFDX
    )
}

FCR <- function(theta_list,
                alpha = 0.1,
                nthresh = 100) {
    rej_prob_list <- map(theta_list, rowMeans)
    cluster_sizes <- # Area of the scanning windows
        (as.numeric(names(theta_list))) ^ 2 / 4
    cluster_size_vec <-
        # expanding the cluster sizes for ease of computing denominator of BFCR
        map2(cluster_sizes, rej_prob_list, ~ rep(.x, times = length(.y))) %>%
        unlist
    
    # the larger this value is, the higher the FCR
    weighted_rej <-
        map2(cluster_sizes, rej_prob_list, ~ .x * .y) %>% unlist
    rej_prob <- unlist(rej_prob_list)
    thresh <- seq(0, 1, length = nthresh)
    
    BFCR <- rep(0, nthresh)
    for (j in 1:nthresh) {
        idx <- rej_prob >= thresh[j]
        if (sum(idx) > 0) {
            BFCR[j] <- 1 - sum(weighted_rej[idx]) / sum(cluster_size_vec[idx])
        }
    }
    level <- min(thresh[BFCR < alpha])
    reject <- weighted_rej > level
    
    list(
        level = level,
        reject_list = split(reject, cluster_size_vec) %>% `names<-` (names(theta_list)),
        thresh = thresh,
        BFCR = BFCR
    )
}
