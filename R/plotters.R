plot_rejections_along_process <-
    function(scc_scan_file,
             mcmc_fit_file,
             sampled_nngps_file,
             rejection_files,
             rejection_names = NULL) {

    z <- readRDS(scc_scan_file)
    preds <- readRDS(sampled_nngps_file)$predictions
    fits <- readRDS(mcmc_fit_file)

    if(!is.null(rejection_names)) {
        rej <- purrr::map(rejection_files, readRDS) %>%
            setNames(rejection_names)
    } else {
        rej <- purrr::map(rejection_files, readRDS) %>%
            setNames(rejection_files)
    }

    winsizes <- names(z)
    p <- list()
    for(i in seq_along(z)) {
        winsize <- winsizes[i]

        #-----------------------------------------------------------------------
        # Prep data about observed scc scan statistics, estimated mean
        #   function, and average sampled nngps
        #-----------------------------------------------------------------------
        z_df <- purrr::imap(z[[i]], ~ mutate(.x, "process" = .y)) %>%
            bind_rows %>%
            mutate(scc = NULL) %>%
            setNames(c("crd", "val", "process")) %>%
            dplyr::filter(crd %in% rej[[1]][[i]]$crd)
        pred_df <-
            tidyr::tibble("crd" = z[[i]]$z1$crd,
                          "val" = rowMeans(preds[[i]]),
                          "process" = "z_star") %>%
            dplyr::filter(crd %in% rej[[1]][[i]]$crd)

        mean_process_df <- tidyr::tibble(
            "crd" = z[[i]]$z1$crd,
            "val" = as.vector(fits[[i]]$X %*% colMeans(fits[[i]]$chain$beta[preds$stationary_iterations,])),
            "process" = "mean_function"
        ) %>%
            dplyr::filter(crd %in% rej[[1]][[i]]$crd)

        line_df <- as_tibble(bind_rows(z_df, pred_df, mean_process_df)) %>%
            mutate(process = as.factor(process))

        #-----------------------------------------------------------------------
        # Prep data about rejection locations
        #-----------------------------------------------------------------------
        rej_df <- map(rej, i) %>%
            map2(.y = names(rej), ~ mutate(.x, "criterion" = .y)) %>%
            bind_rows() %>%
            dplyr::mutate(val = rep(pred_df$val, length(rej))) %>%
            dplyr::filter(reject) %>%
            mutate(reject = NULL)

        #-----------------------------------------------------------------------
        # Plot rejections along the process
        #-----------------------------------------------------------------------
        p[[i]] <- ggplot(line_df, aes(x = crd, y = val)) +
            geom_line(aes(color = process)) +
            geom_point(data = rej_df, aes(x = crd, y = val, shape = criterion), fill = "goldenrod3", size = 3) +
            scale_shape_manual(values = c(22, 2, 24, 8, 23, 21, 4, 3)) +
            ggtitle(paste0("Window size = ", winsize)) +
            labs(y = "z", x = "loc") +
            theme_minimal()
    }

    return(p)
}


plot_cond_vs_cond <-
    function(infiles1,
             infiles2,
             resolution,
             offset = TRUE,
             condition_names = NULL,
             sub_range = NULL) {

    # THIS FUNCTION DOWNSAMPLES THE DATA TO EQUAL READS FOR MORE ACCURATE COMPARISON
    # infiles1, infiles2: vector/list of paths to input files for both conditions
    # resolution: resolution of the data in base pairs, e.g. 50000
    # condition names: names for axes of the heatmap
    # offset: are the data 0-indexed?
    # sub_range: sub index to be plotted

    if(is.null(condition_names)) {
        condition_names <- c("condition 1", "condition 2")
    }

    d <- c(
        purrr::map(infiles1,
                   ~ readr::read_tsv(.x, col_names = FALSE)) %>%
            setNames(paste0("cond1_rep", seq_along(infiles1))),
        purrr::map(infiles2,
                   ~ readr::read_tsv(.x, col_names = FALSE)) %>%
            setNames(paste0("cond2_rep", seq_along(infiles2)))
    )

    cond1 <- grep("cond1", names(d))
    cond2 <- grep("cond2", names(d))

    # Ensure consistent dimensions across replicates
    offset <- as.numeric(offset)
    maxdim <-
        max(map_dbl(d, ~ max(.x$X1, .x$X2))) / resolution + offset

    if(is.null(sub_range)) {
        sub_range <- 1:maxdim
    } else {
        sub_range <- round(sub_range[1] / resolution):round(sub_range[2] / resolution)
    }

    # must be lower-tri
    dmat <-
        purrr::map(
            d,
            ~ Matrix::sparseMatrix(
                i = .x$X2 / resolution + offset,
                j = .x$X1 / resolution + offset,
                x = round(.x$X3),
                dims = rep(maxdim, 2)
            )
        )

    d1 <- Reduce(`+`, dmat[cond1])
    d2 <- Reduce(`+`, dmat[cond2])

    rm(d)
    rm(dmat)

    dd <- downsample_to_equal_reads(list(d1, d2))

    m1 <- reshape2::melt(as.matrix(dd[[1]])) %>%
        dplyr::filter(value != 0) %>%
        dplyr::filter(Var1 %in% sub_range & Var2 %in% sub_range)
    m2 <- reshape2::melt(as.matrix(dd[[2]])) %>%
        dplyr::filter(value != 0) %>%
        dplyr::filter(Var1 %in% sub_range & Var2 %in% sub_range)
    rm(d1)
    rm(d2)
    rm(dd)

    # Remove redundant diagonal information
    m1 <- dplyr::filter(m1, Var1 != Var2)


    # Make condition 1 lower-tri, condition 2 upper-tri
    colnames(m2) <- c("Var2", "Var1", "value")
    m2 <- m2[,c(2,1,3)]

    df <- bind_rows(m1, m2)

    p <- ggplot(data = df, aes(x = Var1, y = Var2, fill = log2(value + 1))) +
        geom_tile() +
        scale_fill_gradient(low = "white", high = "red", na.value = "white") +
        labs(x = condition_names[1], y = condition_names[2]) +
        theme_minimal() +
        theme(legend.position = "none")

    return(p)

}


plot_rej_vs_diffs <-
    function(infiles1,
             infiles2,
             rejections_file,
             resolution,
             condition_names = NULL,
             offset = TRUE,
             sub_range = NULL,
             absolute = FALSE) {
        # THIS FUNCTION DOWNSAMPLES THE DATA TO EQUAL READS FOR MORE ACCURATE COMPARISON
        # infiles1, infiles2: vector/list of paths to input files for both conditions
        # rejection_file: file output from test_with_FDR or test_with_FDX
        # resolution: resolution of the data in base pairs, e.g. 50000
        # condition names: names for axes of the heatmap
        # offset: are the data 0-indexed?
        # sub_range: sub index to be plotted

        if (is.null(condition_names)) {
            condition_names <- c("condition 1", "condition 2")
        }

        d <- c(
            purrr::map(infiles1,
                       ~ readr::read_tsv(.x, col_names = FALSE)) %>%
                setNames(paste0("cond1_rep", seq_along(infiles1))),
            purrr::map(infiles2,
                       ~ readr::read_tsv(.x, col_names = FALSE)) %>%
                setNames(paste0("cond2_rep", seq_along(infiles2)))
        )

        cond1 <- grep("cond1", names(d))
        cond2 <- grep("cond2", names(d))

        # Ensure consistent dimensions across replicates
        offset <- as.numeric(offset)
        maxdim <-
            max(map_dbl(d, ~ max(.x$X1, .x$X2))) / resolution + offset

        if (is.null(sub_range)) {
            sub_range <- 1:maxdim
        } else {
            sub_range <-
                round(sub_range[1] / resolution):round(sub_range[2] / resolution)
        }

        # must be lower-tri
        dmat <-
            purrr::map(
                d,
                ~ Matrix::sparseMatrix(
                    i = .x$X2 / resolution + offset,
                    j = .x$X1 / resolution + offset,
                    x = round(.x$X3),
                    dims = rep(maxdim, 2)
                )
            )

        d1 <- Reduce(`+`, dmat[cond1])
        d2 <- Reduce(`+`, dmat[cond2])

        rm(d)
        rm(dmat)

        dd <- downsample_to_equal_reads(list(d1, d2))
        
        if(absolute) {
            abs_diffs <- abs(dd[[1]] - dd[[2]])    
        } else {
            abs_diffs <- dd[[1]] - dd[[2]]
        }
        
        rejections <- readRDS(rejections_file)
        winsizes <- as.numeric(names(rejections))

        if(!any(map_lgl(map(rejections, "reject"), any))) {
            stop(paste0("No rejections are contained in file ", rejections_file, "."))
        } else {
            init_idx <- which(map_lgl(map(rejections, "reject"), any))[1]

            rej_mat <- Matrix::sparseMatrix(
                i = rejections[[init_idx]]$crd[which(rejections[[init_idx]]$reject)[1]],
                j = rejections[[init_idx]]$crd[which(rejections[[init_idx]]$reject)[1]] + 1,
                x = 1,
                dims = rep(maxdim, 2)
            )

            for(j in seq_along(rejections)) {
                if (any(rejections[[j]]$reject)) {
                    win_size <- winsizes[j]
                    rej_mat <- cpopulate_rejected_differences(
                        rej_mat,
                        rejections[[j]]$crd[rejections[[j]]$reject] - 1, win_size
                    )
                }
            }

            #-------------------------------------------------------------------
            # Make triangular (without diagonal) to make accurate comparisons
            # Lower tri is the absolute differences between conditions,
            #   upper tri are the rejections by given method
            rej_mat <- Matrix::triu(rej_mat, 1)

            diff_molten <- reshape2::melt(abs_diffs) %>%
                dplyr::filter(Var1 %in% sub_range & Var2 %in% sub_range)

            rej_molten <- reshape2::melt(as.matrix(rej_mat)) %>%
                dplyr::filter(Var1 %in% sub_range & Var2 %in% sub_range)

            #-------------------------------------------------------------------
            # plot
            if(absolute) {
                labels <- c(paste0("|", paste0(condition_names, collapse = "-"), "|"), "rejections")    
            } else {
                labels <- c(paste0(condition_names, collapse = "-"), "rejections")
            }
            
            
            if(absolute) {
                rej_p <- ggplot(data = dplyr::filter(rej_molten, Var1 < Var2), aes(
                    x = Var1,
                    y = Var2,
                    fill = value
                )) +
                    geom_tile() +
                    scale_fill_distiller(limits = c(0,1), palette = "Greys", direction = 1) +
                    ggnewscale::new_scale("fill") +
                    geom_tile(
                        mapping = aes(
                            x = Var1,
                            y = Var2,
                            fill = log2(value + 1)
                        ),
                        data = dplyr::filter(diff_molten, Var1 >= Var2)
                    ) +
                    coord_fixed() +
                    scale_fill_gradient(low = "white", high = "red", na.value = "white") +
                    theme_minimal() +
                    labs(x = labels[1], y = labels[2]) +
                    theme(legend.position = "none")
            } else {
                diff_molten[diff_molten$value < 0,"value"] <- -log2(abs(diff_molten[diff_molten$value < 0,"value"]))
                diff_molten[diff_molten$value > 0,"value"] <- log2(diff_molten[diff_molten$value > 0,"value"])
                rej_p <- ggplot(data = dplyr::filter(rej_molten, Var1 < Var2), aes(
                    x = Var1,
                    y = Var2,
                    fill = value
                )) +
                    geom_tile() +
                    scale_fill_distiller(limits = c(0,1), palette = "Greys", direction = 1) +
                    ggnewscale::new_scale("fill") +
                    geom_tile(
                        mapping = aes(
                            x = Var1,
                            y = Var2,
                            fill = value
                        ),
                        data = dplyr::filter(diff_molten, Var1 >= Var2)
                    ) +
                    coord_fixed() +
                    scale_fill_distiller(palette = "RdBu", na.value = "white") +
                    theme_minimal() +
                    labs(x = labels[1], y = labels[2]) +
                    theme(legend.position = "none")
            }
            return(rej_p)
        }
    }
