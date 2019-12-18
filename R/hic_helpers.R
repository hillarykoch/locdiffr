downsample_to_equal_reads <- function(data_list) {
    # total counts in each data matrix
    maxes <- purrr::map(data_list, sum) %>%
        purrr::map(max) %>%
        unlist
    
    # what level to downsample to
    tot_counts <- min(maxes)
    
    # only downsample data matrices who have more reads than tot_counts
    sample_idx <- which(!(maxes == tot_counts))
    out <- data_list
    for(i in sample_idx) {
        unlisted_data <- unlist(data_list[[i]])
        # data_list[[i]] <- 1
        expand <- rep.int(seq(length(unlisted_data)),
                          times = as.numeric(unlisted_data))
        samp <- sample(expand, size = tot_counts, replace = FALSE) %>% sort
        rles <- rle(samp)
        temp <- matrix(0, nrow = nrow(data_list[[i]]), ncol = ncol(data_list[[i]]))
        temp[rles$values] <- rles$lengths
        out[[i]] <- temp
    }
    out
}
