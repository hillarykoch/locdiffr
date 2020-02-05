downsample_to_equal_reads <- function(data_list, sparse = FALSE) {
    maxes <- purrr::map(data_list, sum) %>% purrr::map(max) %>% 
        unlist
    tot_counts <- min(maxes)
    sample_idx <- which(!(maxes == tot_counts))
    out <- data_list
    for (i in sample_idx) {
        if(sparse) {
            unlisted_data <- data_list[[i]]@x
        } else {
            unlisted_data <- unlist(data_list[[i]])    
        }
        
        expand <- rep.int(seq(length(unlisted_data)), times = as.numeric(unlisted_data))
        samp <- sample(expand, size = tot_counts, replace = FALSE) %>% 
            sort
        rles <- rle(samp)
        if(sparse) {
            temp <- matrix(0, nrow = data_list[[i]]@Dim[1], ncol = data_list[[i]]@Dim[1])
        } else {
            temp <- matrix(0, nrow = nrow(data_list[[i]]), ncol = ncol(data_list[[i]]))
        }
        temp[rles$values] <- rles$lengths
        out[[i]] <- temp
    }
    out <- purrr::map(out, as.matrix)
    out
}
