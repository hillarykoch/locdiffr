#-------------------------------------------------------------------------------
# Wrapper functions
#-------------------------------------------------------------------------------

run_scc_scan <-
    function(infiles1,
             infiles2,
             outpath,
             resolution,
             winsizes = seq(10, 100, by = 10),
             parallel = FALSE,
             ncores = 10,
             offset = TRUE) {
        # cond1, cond2: list of paths to data from condition 1/condition 2 in BED-type format: LOC1, LOC2, COUNTS
        #   If different numbers of replicates, the function will automatically choose which to keep based on sequencing depth
        # outpath: where to write the output file
        # resolution: resolution of Hi-C matrix
        # winsizes: which window sizes (in bins) to use? A vector
        # parallel: do the scan in parallel?
        # ncores: how many cores to use, if parallel
        # offset: basically, if first bin is a location 0, set to TRUE to offset by 1 for matrix indexing

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



        # Don't make all datasets equal reads. Instead, match them by their rankings
        # E.g., take the lowest read control and lowest read treatment, and make them equal.
        # then the second lowest, etc.
        cond1_reord <-
            map(d[cond1], "X3") %>% map(~ sum(.x) * -1) %>% unlist %>% order
        cond2_reord <-
            map(d[cond2], "X3") %>% map(~ sum(.x) * -1) %>% unlist %>% order

        # This reorders the matrices by read count within treatment and within control
        split_dmat <-
            split(dmat, as.factor(grepl("cond2", names(d)))) %>%
            map2(.y = list(cond1_reord, cond2_reord), ~ .x[.y])

        # if different number of replicates, figure out how many to keep
        num_comparisons <- min(map_int(split_dmat, length))

        # Downsample in a read depth-aware manner.
        dd <-
            map2(.x = split_dmat[[1]][1:num_comparisons],
                 .y = split_dmat[[2]][1:num_comparisons],
                 ~ downsample_to_equal_reads(list(.x, .y))) %>%
            map( ~ setNames(.x, c("cond1", "cond2")))
        rm(dmat)



        if (parallel) {
            cluster <- makeCluster(ncores)
            registerDoParallel(cluster)
            clusterEvalQ(cluster, library(sgp))
            clusterEvalQ(cluster, library(tidyverse))

            zs <-
                foreach(i = seq_along(winsizes), .packages = "foreach") %dopar% {
                    winsize <- winsizes[i]
                    sccs <-
                        map2(
                            map(dd, "cond1"),
                            map(dd, "cond2"),
                            ~ sgp::get_loc_sim(
                                .x,
                                .y,
                                h = 0,
                                win_min = winsize,
                                win_max = winsize
                            )
                        )


                    # Get z-scores of these statistics
                    z <- map(sccs, get_z)

                    # Take out loci which are exactly equal (e.g., loci at the centromere)
                    oneidx <-
                        Reduce(intersect, map(sccs, ~ which(.x == 1)))

                    # save scc results
                    map(z, ~ dplyr::filter(.x, !(crd %in% oneidx))) %>%
                        setNames(paste0("z", seq_along(z)))
                }

            stopCluster(cluster)

        } else {
            zs <- list()
            for (i in seq_along(winsizes)) {
                winsize <- winsizes[i]
                sccs <-
                    map2(
                        map(dd, "cond1"),
                        map(dd, "cond2"),
                        ~ sgp::get_loc_sim(
                            .x,
                            .y,
                            h = 0,
                            win_min = winsize,
                            win_max = winsize
                        )
                    )


                # Get z-scores of these statistics
                z <- map(sccs, get_z)

                # Take out loci which are exactly equal (e.g., loci at the centromere)
                oneidx <-
                    Reduce(intersect, map(sccs, ~ which(.x == 1)))

                # save scc results
                zs[[i]] <-
                    map(z, ~ dplyr::filter(.x, !(crd %in% oneidx))) %>%
                    setNames(paste0("z", seq_along(z)))
            }
        }

        names(zs) <- winsizes

        saveRDS(zs, file = outpath)
    }

fit_nngp <-
    function(infile,
             outpath,
             num_neighbors = 15,
             iters = 15000,
             parallel = FALSE,
             ncores = 10) {
        # infile: path to .rds file output from run_scc_scan
        # outpath: where to write results
        # num_neighbors: how many neighbors to use for mcmc?
        # iters: number of mcmc iterations to do
        # parallel: do the scan in parallel?
        # ncores: how many cores to use, if parallel

        # Read in scanning data
        zs <- readRDS(infile)
        winsizes <- as.numeric(names(zs))


        if (parallel) {
            cluster <- makeCluster(ncores)
            registerDoParallel(cluster)
            clusterEvalQ(cluster, library(sgp))
            clusterEvalQ(cluster, library(tidyverse))
            clusterEvalQ(cluster, library(fda))

            fits <-
                foreach(i = seq_along(winsizes), .packages = "foreach") %dopar% {
                    s <- matrix(zs[[i]]$z1$crd, ncol = 1)
                    y <- map(zs[[i]], "z_s") %>%
                        dplyr::bind_cols() %>%
                        as.matrix
                    X <- fda::getbasismatrix(
                        evalarg = s,
                        fda::create.bspline.basis(
                            rangeval = range(s),
                            nbasis = max(round(nrow(s) / 300), 4),
                            norder = 4
                        )
                    )

                    list(
                        "X" = X,
                        "chain" = run_nnsgp(
                            y = y,
                            s = s,
                            X = X,
                            num_neighbors = num_neighbors,
                            iters = iters
                        )
                    )
                }

            stopCluster(cluster)

        } else {
            fits <- list()
            for (i in seq_along(winsizes)) {
                s <- matrix(zs[[i]]$z1$crd, ncol = 1)
                y <- map(zs[[i]], "z_s") %>%
                    dplyr::bind_cols() %>%
                    as.matrix
                X <- fda::getbasismatrix(
                    evalarg = s,
                    fda::create.bspline.basis(
                        rangeval = range(s),
                        nbasis = max(round(nrow(s) / 300), 4),
                        norder = 4
                    )
                )

                fits[[i]] <-
                    list(
                        "X" = X,
                        "chain" = run_nnsgp(
                            y = y,
                            s = s,
                            X = X,
                            num_neighbors = num_neighbors,
                            iters = iters
                        )
                    )
            }
        }

        names(fits) <- winsizes
        saveRDS(fits, file = outpath)
    }

sample_new_nngps <-
    function(scc_scan_file,
             mcmc_fit_file,
             outpath,
             stationary_iterations,
             parallel = FALSE,
             ncores = 10) {
        # scc_scan_file: path to rds file output from run_scc_scan
        # mcmc_fit_file: path to rds file output from fit_nngp
        # outpath: where to save newly sampled nngps
        # stationary_iterations: which iterations to draw new processes from
        # parallel: do the scan in parallel?
        # ncores: how many cores to use, if parallel

        #-------------------------------------------------------------------------------
        # Read in sgp output
        #-------------------------------------------------------------------------------
        zs <- readRDS(scc_scan_file)
        fits <- readRDS(mcmc_fit_file)

        #-------------------------------------------------------------------------------
        # prediction controls
        #-------------------------------------------------------------------------------
        iters <- nrow(fits[[1]]$chain$beta)
        winsizes <- as.numeric(names(zs))

        if (parallel) {
            cluster <- makeCluster(ncores)
            registerDoParallel(cluster)
            clusterEvalQ(cluster, library(sgp))
            clusterEvalQ(cluster, library(tidyverse))

            preds <-
                foreach(i = seq_along(winsizes), .packages = "foreach") %dopar% {
                    s <- matrix(zs[[i]]$z1$crd, ncol = 1)
                    y <- map(zs[[i]], "z_s") %>%
                        dplyr::bind_cols() %>%
                        as.matrix
                    X <- fits[[i]]$X

                    make_pred_sparse(fits[[i]]$chain, y, s, X, stationary_iterations)
                }

            stopCluster(cluster)
        } else {
            preds <- list()

            for (i in seq_along(winsizes)) {
                s <- matrix(zs[[i]]$z1$crd, ncol = 1)
                y <- map(zs[[i]], "z_s") %>%
                    dplyr::bind_cols() %>%
                    as.matrix
                X <- fits[[i]]$X

                preds[[i]] <-
                    make_pred_sparse(fits[[i]]$chain, y, s, X, stationary_iterations)
            }
        }

        names(preds) <- winsizes
        preds$stationary_iterations <- stationary_iterations

        saveRDS(preds, file = outpath)
    }
