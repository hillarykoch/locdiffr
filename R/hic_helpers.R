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
    out <- purrr::map(out, as.matrix)
    out
}

get_permutation_data <- function(dat, modidx) {
    # dat is a Hi-C data matrix
    # modidx is the array of indices which will be replaced with permuted data
    molten <- reshape2::melt(unname(as.matrix(dat)))
    diffs <- molten[,2] - molten[,1]
    keepidx <- diffs > 0
    
    bkgrd <- data.frame("diffs" = diffs[keepidx], "value" = molten[keepidx,"value"])
    resampling <- tapply(bkgrd$value, INDEX = bkgrd$diffs, FUN = function(X) sample(X))
    for(i in seq_along(unique(bkgrd$diffs))[-length(unique(bkgrd$diffs))]) {
        curdiff <- unique(bkgrd$diffs)[i]
        bkgrd[bkgrd$diffs == curdiff,"value"] <- resampling[[curdiff]]
    }
    
    
    idx_list <- list()
    for (i in seq_along(unique(bkgrd$diffs))) {
        curdiff <- unique(bkgrd$diffs)[i]
        df <- data.frame("idx_start" = seq(sum(bkgrd$diffs == i)),
                         "idx_stop" = seq(sum(bkgrd$diffs == i)) + curdiff)
        idx_list[[curdiff]] <-
            dplyr::mutate(df, "val" = (
                dplyr::filter(bkgrd, diffs == curdiff) %>% select("value") %>% `[[` (1)
            ))
    }
    
    bind1 <- dplyr::bind_rows(idx_list)
    mirror <-
        data.frame(
            "idx_start" = bind1$idx_stop,
            "idx_stop" = bind1$idx_start,
            "val" = bind1$val
        )
    permuted_mat <- as.matrix(data.table::dcast(data = dplyr::bind_rows(list(bind1, mirror)), idx_start ~ idx_stop)[,-1])
    diag(permuted_mat) <- 0
    
    dat[modidx, modidx] <- permuted_mat[modidx, modidx]
    dat
}
