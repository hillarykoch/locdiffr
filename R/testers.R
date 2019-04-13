FDR <- function(theta,
                alpha = .1,
                nthresh = 100) {
    # Theta is an indicator -- was the predicted value as location s greater than mu0?
    #pick reject so that E(mean(theta[reject]))<confidence

    # Probability at each prediction location of an observation being greater than mu0
    inprob <- apply(theta, 2, mean)

    # A grid of thresholds to test
    thresh <- seq(0, 1, length = nthresh)

    # The Bayes FDR at each threshold
    # (want to find one that is closest to and not greater than choice of alpha)
    BFDR <- rep(0, nthresh)
    for (j in 1:nthresh) {
        if (sum(inprob >= thresh[j]) > 0) {
            BFDR[j] <- 1 - mean(inprob[inprob >= thresh[j]])
        }
    }
    level <- min(thresh[BFDR < alpha])
    reject <- inprob > level
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
    #pick reject so that P(mean(theta[reject])<alpha)<beta

    reject <- rep(0, nrow(theta))
    level <- 0

    inprob <- apply(theta, 2, mean) # proportion rejected at each spatial location
    thresh <- seq(0, max(inprob), length = nthresh)
    BFDX <- rep(0, nthresh)
    for (j in 1:nthresh) {
        # for fixed threshold in range 0 -- max(inprob), which sites have proportion rejected greater than that threshhold?
        these <- inprob >= thresh[j]
        if (sum(these) > 1) {
            # How many of these have propotion of rejections < 1-alpha
            BFDX[j] <- mean(1-apply(theta[, these], 1, mean) < alpha)
        }
    }

    # The level is the minimum value such that the BFDX exceeds bound beta
    level <- 1
    if (sum(BFDX > beta) > 0) {
        level <- min(thresh[BFDX > beta])
    }
    reject <- inprob > level
    list(
        level = level,
        reject = reject,
        thresh = thresh,
        BFDX = BFDX
    )
}
