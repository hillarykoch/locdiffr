FDR <- function(theta,
                alpha = .1,
                nthresh = 100) {
    # Theta is an indicator -- was the predicted value as location s
    #   less than X * beta at location s?
    # reject constructed such that E(mean(theta[reject])) < alpha

    rej_prob <- rowMeans(theta)
    thresh <- seq(0, 1, length = nthresh)

    # The Bayes FDR at each threshold
    # (want to find one that is closest to and not greater than choice of alpha)
    BFDR <- rep(0, max(rej_prob))
    for (j in 1:nthresh) {
        if (any(rej_prob >= thresh[j])) {
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
                beta = .5,
                nthresh = 100) {
    # reject constructed such that P(mean(theta[reject]) > beta) <= alpha
    rej_prob <- rowMeans(theta)
    thresh <- seq(0, max(rej_prob), length = nthresh)
    BFDX <- rep(0, nthresh)
    for (j in 1:nthresh) {
        xceeds <- rej_prob >= thresh[j]
        if (sum(xceeds) > 1) {
            BFDX[j] <- mean( 1 - rowMeans(theta[xceeds,]) > beta)
        }
    }

    # The level is the minimum value such that the BFDX is less than alpha
    level <- 1
    if (any(BFDX <= alpha)) {
        level <- min(thresh[BFDX <= alpha])
    }

    reject <- rej_prob >= level
    list(
        level = level,
        reject = reject,
        thresh = thresh,
        BFDX = BFDX
    )
}


wFDR <- function(theta_list,
                 alpha = 0.1,
                 nthresh = 100) {
    rej_prob_list <- map(theta_list, rowMeans)

    cluster_sizes <- # Area of the scanning windows
        (as.numeric(names(theta_list))) ^ 2 / 4
    cluster_size_vec <-
        # expanding the cluster sizes for ease of computing denominator of BwFDR
        map2(cluster_sizes, rej_prob_list, ~ rep(.x, times = length(.y))) %>%
        unlist

    # the larger this value is, the higher the wFDR
    weighted_rej <-
        map2(cluster_sizes, rej_prob_list, ~ .x * .y) %>% unlist
    rej_prob <- unlist(rej_prob_list)
    thresh <- seq(0, max(rej_prob), length = nthresh)

    BwFDR <- ccompute_bwfdr(weighted_rej, rej_prob, thresh, cluster_size_vec)

    if(any(BwFDR <= alpha)) {
        reject <- rej_prob >= thresh[which(BwFDR <= alpha)[1]]
    } else {
        reject <- rep(FALSE, times = length(rej_prob))
    }
    reject_list <- split(unname(reject), cluster_size_vec) %>% `names<-` (names(theta_list))

    reject_list
}

wFDX <- function(theta_list,
                 alpha = .1,
                 beta = .1,
                 nthresh = 100,
                 bootstrap_replicates = 1000) {
    rej_prob_list <- map(theta_list, rowMeans)

    cluster_sizes <- # Area of the scanning windows
        (as.numeric(names(theta_list))) ^ 2 / 4
    cluster_size_vec <-
        # expanding the cluster sizes for ease of computing denominator of BwFDR
        map2(cluster_sizes, rej_prob_list, ~ rep(.x, times = length(.y))) %>%
        unlist

    # the larger this value is, the higher the wFDX
    weighted_rej <-
        map2(cluster_sizes, rej_prob_list, ~ .x * .y) %>% unlist
    rej_prob <- unlist(rej_prob_list)


    thresh <- seq(0, max(rej_prob), length = nthresh)
    BwFDX <- ccompute_bwfdx(weighted_rej, rej_prob, thresh, cluster_size_vec, bootstrap_replicates, beta)

    level <- 1
    if(any(BwFDX <= alpha)) {
        level <- min(thresh[BwFDX <= alpha], na.rm = TRUE)
    }


    reject <- rej_prob >= level
    reject_list <- split(unname(reject), cluster_size_vec) %>%
        `names<-` (names(theta_list))

    reject_list

}
