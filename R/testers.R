
wFDR <- function(theta_list,
                 alpha = 0.01,
                 nthresh = 100) {
    rej_prob_list <- map(theta_list, rowMeans)

    cluster_sizes <- # Area of the scanning windows
        (as.numeric(names(theta_list))) ^ 2 / 4
    cluster_size_vec <-
        # expanding the cluster sizes for ease of computing denominator of BwFDR
        map2(cluster_sizes, rej_prob_list, ~ rep(.x, times = length(.y))) %>%
        unlist

    if(alpha == 0) {
        reject <- rep(0, length(unlist(rej_prob_list)))
        reject_list <- split(unname(reject), cluster_size_vec) %>%
            `names<-` (names(theta_list))

        return(reject_list)
    }

    # the larger this value is, the higher the wFDR
    weighted_rej <-
        map2(cluster_sizes, rej_prob_list, ~ .x * .y) %>% unlist
    rej_prob <- unlist(rej_prob_list)
    thresh <- seq(0, max(rej_prob), length = nthresh)

    BwFDR <- ccompute_bwfdr(weighted_rej, rej_prob, thresh, cluster_size_vec)

    if(sum(BwFDR <= alpha, na.rm = TRUE) > 0) {
        reject <- rej_prob > thresh[which(BwFDR <= alpha)[1]]
    } else {
        reject <- rep(FALSE, times = length(rej_prob))
    }
    reject_list <- split(unname(reject), cluster_size_vec) %>% `names<-` (names(theta_list))

    reject_list
}

wFDX <- function(theta_list,
                 alpha = .01,
                 beta = .01,
                 nthresh = 100,
                 bootstrap_replicates = 1000) {

    rej_prob_list <- map(theta_list, rowMeans)

    cluster_sizes <- # Area of the scanning windows
        (as.numeric(names(theta_list))) ^ 2 / 4
    cluster_size_vec <-
        # expanding the cluster sizes for ease of computing denominator of BwFDR
        map2(cluster_sizes, rej_prob_list, ~ rep(.x, times = length(.y))) %>%
        unlist

    if(alpha == 0 & beta == 0) {
        reject <- rep(0, length(unlist(rej_prob_list)))
        reject_list <- split(unname(reject), cluster_size_vec) %>%
            `names<-` (names(theta_list))

        return(reject_list)
    }

    # the larger this value is, the higher the wFDX
    weighted_rej <-
        map2(cluster_sizes, rej_prob_list, ~ .x * .y) %>% unlist
    rej_prob <- unlist(rej_prob_list)


    thresh <- seq(0, max(rej_prob), length = nthresh)
    BwFDX <- ccompute_bwfdx(weighted_rej, rej_prob, thresh, cluster_size_vec, bootstrap_replicates, beta)

    level <- 1
    if(sum(BwFDX <= alpha, na.rm = TRUE) > 0) {
        level <- min(thresh[BwFDX <= alpha], na.rm = TRUE)
    }


    reject <- rej_prob >= level
    reject_list <- split(unname(reject), cluster_size_vec) %>%
        `names<-` (names(theta_list))

    reject_list

}
