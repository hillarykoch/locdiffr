
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


get_enrichment <- function(rej_locs, num_locs, loc_range, peak_locs, win_size, nreps = 100, multicount = TRUE) {
    num <- get_rejection_enrichment(rej_locs, peak_locs, win_size, multicount)
    denom <- replicate(nreps, get_background_enrichment(num_locs, loc_range, peak_locs, win_size, multicount))
    out <- c(num/mean(denom), 1 - pnorm(num, mean = mean(denom), sd = sd(denom)))
    names(out) <- c("ratio", "p-value")
    out
}
