
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
                 beta,
                 nthresh = 100,
                 bootstrap_replicates = 500) {
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
