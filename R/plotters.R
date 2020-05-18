
plot_rejections_along_process <-
    function(scc_scan_file,
             mcmc_fit_file,
             sampled_nngps_file,
             rejection_files,
             rejection_names = NULL) {

    z <- readRDS(scc_scan_file)
    preds <- readRDS(sampled_nngps_file)
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
            setNames(c("crd", "val", "process"))
        pred_df <-
            tidyr::tibble("crd" = z[[i]]$z1$crd,
                          "val" = rowMeans(preds[[i]]),
                          "process" = "z_star")

        mean_process_df <- tidyr::tibble(
            "crd" = z[[i]]$z1$crd,
            "val" = as.vector(fits[[i]]$X %*% colMeans(fits[[i]]$chain$beta[preds$stationary_iterations,])),
            "process" = "mean_function"
        )

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
