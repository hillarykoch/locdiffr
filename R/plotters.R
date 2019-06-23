# setwd("~/Box Sync/School/research - Qunhua/taos_project/contrib/figs_etc/")
# 
# library(magrittr)
# library(readr)
# library(dplyr)
# library(fda)
# library(ggplot2)
# library(cowplot)
# library(purrr)
# 
# KO1 <-
#     as.matrix(
#         read.delim(
#             "../../received_progress/processedData/G1E-ER4-C8BRD2-KO-uninduced.Rep2.chr19.matrix"
#         )
#     )
# ui1 <-
#     as.matrix(
#         read.delim(
#             "../../received_progress/processedData/G1E-ER4-uninduced.Rep1.chr19.matrix"
#         )
#     )
# 
# scc_path <- "../fits/20/scc_C8_chr19_h2.Rdata"
# mcmc_path <- "../fits/20/C8_10bf_h2.Rdata"
# h <- 2
# nb <- 10
# load(mcmc_path)
# load(scc_path)
# 
# # BRD2 ChIP peaks (2-column matrix/df/tbl of peak locations)
# peaks <-
#     readr::read_tsv("GSM1532602_530.broadpeak", col_names = FALSE) %>%
#     dplyr::filter(X1 == "chr19") %>%
#     select(X2, X3)


plot_rejections <-
    function(fit,
             sccs,
             h,
             nb,
             peaks,
             burnin = 0,
             alpha = 0.1,
             peak_x = -2.5,
             resolution = 40000,
             main = "rejections",
             num_rejections = NULL) {
        #--------------------------------------------------
        # Filter unneeded windows
        if (length(sccs) > 1) {
            keepers <- Reduce(intersect, purrr::map(sccs, "crd"))
            for (i in seq_along(sccs)) {
                sccs[[i]] <- sccs[[i]][sccs[[i]]$crd %in% keepers,]
            }
        }

        #--------------------------------------------------
        # Reproduce basis functions
        s <- as.matrix(sccs$z1$crd, ncol = 1)
        X <-
            fda::getbasismatrix(s, fda::create.bspline.basis(range(s), nbasis = nb, norder = 4))

        #--------------------------------------------------
        # Compute indicator theta for testing
        theta <- matrix(NA, nrow(fit$pred) - burnin, ncol(fit$pred))
        for (i in 1:nrow(theta)) {
            theta[i,] <- fit$pred[i + burnin,] < X %*% fit$beta[i + burnin,]
        }
        
        if(is.null(num_rejections)) {
            #--------------------------------------------------
            # Compute rejections
            rej_FDR <- FDR(theta, alpha = alpha, nthresh = 100)$reject
            rej_FDX <-
                FDX(theta,
                    alpha = alpha,
                    beta = 1 - alpha,
                    nthresh = 100)$reject
            
            peakdf <- tibble("y" = rep(peak_x, times = nrow(peaks)),
                             "peak" = rowMeans(peaks) / resolution))
            
            ggdf <- tibble(
                crd = sccs$z1$crd,
                lb = apply(fit$pred[(burnin + 1):nrow(fit$pred),], 2, function(X)
                    quantile(X, .05)),
                ub = apply(fit$pred[(burnin + 1):nrow(fit$pred),], 2, function(X)
                    quantile(X, .95)),
                z_star = colMeans(fit$pred[(burnin + 1):nrow(fit$pred),]),
                FDR = rej_FDR,
                FDX = rej_FDX,
                mean_func = (X %*% colMeans(fit$beta[(burnin + 1):nrow(fit$beta),]))[, 1],
                z1 = sccs$z1$z_s
            )
            
            #--------------------------------------------------
            # Generalize function to any number of replicates
            if (length(sccs) > 1) {
                for (i in 2:length(sccs)) {
                    ggdf <- dplyr::mutate(ggdf, temp = sccs[[i]]$z_s)
                }
            }
            nm <-
                c("crd",
                  "lb",
                  "ub",
                  "z_star",
                  "FDR",
                  "FDX",
                  "mean_func",
                  paste0("z", seq_along(sccs)))
            names(ggdf) <- nm
            
            CIdf <- dplyr::select(ggdf, crd, lb, ub)
            line_molten <-
                dplyr::select(ggdf, crd, 7:ncol(ggdf), z_star) %>%
                data.table::melt(id.vars = "crd")
            meandf <- dplyr::select(ggdf, crd, mean_func) %>%
                mutate("mean func." = rep(1, nrow(ggdf)))
            fdrdf <-
                tibble(crd = ggdf$crd[ggdf$FDR],
                       y = ggdf$z_star[ggdf$FDR],
                       id = "FDR")
            fdxdf <-
                tibble(crd = ggdf$crd[ggdf$FDX],
                       y = ggdf$z_star[ggdf$FDX],
                       id = "FDX")
            rejdf <- rbind(fdrdf, fdxdf)
            
            ggplot(CIdf, aes(x = crd, ymin = lb, ymax = ub)) +
                geom_ribbon(color = "grey70", alpha = .4) +
                geom_line(
                    aes(x = crd,  y = value, color = variable),
                    data = line_molten,
                    inherit.aes = FALSE,
                    alpha = .75,
                    size = 1.5
                ) +
                geom_line(
                    aes(
                        x = crd,
                        y = mean_func,
                        linetype = "mean func."
                    ),
                    data = meandf,
                    color = "navyblue",
                    inherit.aes = FALSE,
                    size = 2
                ) +
                geom_point(
                    aes(
                        x = crd,
                        y = y,
                        shape = as.factor(id)
                    ),
                    size = c(rep(2, nrow(fdrdf)), rep(1.5, nrow(fdxdf))),
                    fill = c(rep("goldenrod3", nrow(fdrdf)), rep("hotpink", nrow(fdxdf))),
                    data = rejdf,
                    inherit.aes = FALSE
                ) +
                geom_point(
                    aes(x = peak, y = y, shape = "BRD2 peaks"),
                    data = peakdf,
                    color = "deeppink4",
                    size = 4,
                    inherit.aes = FALSE
                ) +
                scale_shape_manual(values = c(
                    "BRD2 peaks" = 1,
                    "FDR" = 23,
                    "FDX" = 24
                )) +
                theme_cowplot() +
                theme(legend.title = element_blank()) +
                ggtitle(main) +
                xlab("loc") + ylab("z")
        } else {
            rej <- dplyr::between(rank(-colMeans(theta)), 1, num_rejections)
            
            peakdf <- tibble("y" = rep(peak_x, times = nrow(peaks)),
                             "x" = rowMeans(peaks) / resolution,
                             "id" = "BRD2 peaks")
            
            ggdf <- tibble(
                crd = sccs$z1$crd,
                lb = apply(fit$pred[(burnin + 1):nrow(fit$pred),], 2, function(X)
                    quantile(X, .05)),
                ub = apply(fit$pred[(burnin + 1):nrow(fit$pred),], 2, function(X)
                    quantile(X, .95)),
                z_star = colMeans(fit$pred[(burnin + 1):nrow(fit$pred),]),
                rejections = rej,
                mean_func = (X %*% colMeans(fit$beta[(burnin + 1):nrow(fit$beta),]))[, 1],
                z1 = sccs$z1$z_s
            )
            
            #--------------------------------------------------
            # Generalize function to any number of replicates
            if (length(sccs) > 1) {
                for (i in 2:length(sccs)) {
                    ggdf <- dplyr::mutate(ggdf, temp = sccs[[i]]$z_s)
                }
            }
            nm <-
                c("crd",
                  "lb",
                  "ub",
                  "z_star",
                  "rejections",
                  "mean_func",
                  paste0("z", seq_along(sccs)))
            names(ggdf) <- nm
            
            CIdf <- dplyr::select(ggdf, crd, lb, ub)
            line_molten <-
                dplyr::select(ggdf, crd, 7:ncol(ggdf), z_star) %>%
                data.table::melt(id.vars = "crd")
            meandf <- dplyr::select(ggdf, crd, mean_func) %>%
                mutate("mean func." = rep(1, nrow(ggdf)))
            rejdf <-
                tibble(x = ggdf$crd[ggdf$rejections],
                       y = ggdf$z_star[ggdf$rejections],
                       id = "rejections")
            pointdf <- rbind(rejdf, peakdf)
            
            ggplot(CIdf, aes(x = crd, ymin = lb, ymax = ub)) +
                geom_ribbon(color = "grey70", alpha = .4) +
                geom_line(
                    aes(x = crd,  y = value, color = variable),
                    data = line_molten,
                    inherit.aes = FALSE,
                    alpha = .75,
                    size = 1.5
                ) +
                geom_line(
                    aes(
                        x = crd,
                        y = mean_func,
                        linetype = "mean func."
                    ),
                    data = meandf,
                    color = "navyblue",
                    inherit.aes = FALSE,
                    size = 2
                ) +
                geom_point(
                    aes(x = x, y = y, shape = id, color = id),
                    data = pointdf,
                    size = c(rep(2, nrow(rejdf)), rep(4, nrow(peakdf))),
                    color = c(rep("black", nrow(rejdf)), rep("deeppink4", nrow(peakdf))),
                    fill = "goldenrod3",
                    inherit.aes = FALSE
                ) +
                scale_shape_manual(values = c(
                    "BRD2 peaks" = 1,
                    "rejections" = 24
                )) +
                theme_cowplot() +
                theme(legend.title = element_blank()) +
                ggtitle(main) +
                xlab("loc") + ylab("z")
        }
    }

plot_heat <- function(hic, rang, cutoff = 1000) {
    subx <- hic[rang, rang]
    r <- quantile(as.vector(subx), prob = c(0.001, 0.999))
    colr <- c("white", "red")
    pal <- colorRampPalette(colr)(n = cutoff)
    defaultpal <- palette(pal)
    mycol <- round((as.vector(subx) - r[1]) * 100 / (r[2] - r[1]))
    mycol[mycol < 1] <- 1
    mycol[mycol > cutoff] <- cutoff

    l <- nrow(subx)

    mmycol <- matrix(mycol, nrow(subx), nrow(subx))
    molten <- reshape2::melt(mmycol)
    molten[, c("Var1", "Var2")] <-
        molten[, c("Var1", "Var2")] + min(rang)
    ggplot(molten, aes(x = Var1, y = Var2, fill = value)) +
        geom_tile() +
        scale_fill_continuous(low = "white", high = "red")
}

plot_heat_rejections <-
    function(hic1,
             hic2,
             fit,
             zoom,
             h,
             nb,
             win_size,
             sccs,
             method = c("FDR", "FDX", "rank"),
             burnin = 0,
             alpha = 0.1,
             num_rejections = NULL,
             main = "",
             rej_color = "#887711") {
        #--------------------------------------------------
        # Filter unneeded windows
        if (length(sccs) > 1) {
            keepers <- Reduce(intersect, purrr::map(sccs, "crd"))
            for (i in seq_along(sccs)) {
                sccs[[i]] <- sccs[[i]][sccs[[i]]$crd %in% keepers,]
            }
        }

        #--------------------------------------------------
        # Reproduce basis functions
        s <- as.matrix(sccs$z1$crd, ncol = 1)
        X <-
            fda::getbasismatrix(s, fda::create.bspline.basis(range(s), nbasis = nb, norder = 4))

        #--------------------------------------------------
        # Compute indicator theta for testing
        theta <- matrix(NA, nrow(fit$pred) - burnin, ncol(fit$pred))
        for (i in 1:nrow(theta)) {
            theta[i,] <- fit$pred[i + burnin,] < X %*% fit$beta[i + burnin,]
        }

        #--------------------------------------------------
        # Compute rejections
        if (method == "FDR") {
            rej_FDR <- FDR(theta, alpha = alpha, nthresh = 100)$reject
        } else if(method == "FDX") {
            rej_FDR <-
                FDX(
                    theta,
                    alpha = alpha,
                    beta = 1 - alpha,
                    nthresh = 100
                )$reject
        } else {
            rej_FDR <- dplyr::between(rank(-colMeans(theta)), 1, num_rejections)
        }

        smoo1 <- fastMeanFilter(hic1, h)
        smoo2 <- fastMeanFilter(hic2, h)
        fdr <- sccs$z1$crd[rej_FDR]
        sub <- fdr[fdr %in% zoom]
        segdf <-
            tibble(
                x1 = sub,
                y1 = sub,
                x2 = sapply(sub, function(X)
                    min(X + win_size, max(zoom))),
                y2 = sapply(sub, function(X)
                    min(X + win_size, max(zoom)))
            )


        ph1 <- plot_heat(log(smoo1 + 1), zoom, 1000) +
            xlab("") +
            ylab("") +
            ggtitle(main) +
            theme(legend.position = "none")
        ph2 <- plot_heat(log(smoo2 + 1), zoom, 1000) +
            xlab("") +
            ylab("") +
            ggtitle(main) +
            theme(legend.position = "none")

        diffplot <-
            plot_heat(log(abs(smoo1 - smoo2) + 1), zoom, 1000) +
            ggtitle("absolute difference + rejections") +
            xlab("") +
            ylab("") +
            theme(legend.position = "none") +
            geom_segment(
                aes(
                    x = x1,
                    xend = x2,
                    y = y1,
                    yend = y1
                ),
                data = segdf,
                inherit.aes = FALSE,
                color = rej_color,
                alpha = .8,
                size = 1.1
            ) +
            geom_segment(
                aes(
                    x = x2,
                    xend = x2,
                    y = y1,
                    yend = y2
                ),
                data = segdf,
                inherit.aes = FALSE,
                color = rej_color,
                alpha = .8,
                size = 1.1
            )

        list("hic1" = ph1,
             "hic2" = ph2,
             "abs_diff" = diffplot)
    }

add_rejections <-
    function(diffplot,
             fit,
             zoom,
             nb,
             win_size,
             sccs,
             method = c("FDR", "FDX", "rank"),
             burnin = 0,
             alpha = 0.1,
             num_rejections = NULL,
             rej_color = "#887711") {
        #--------------------------------------------------
        # Filter unneeded windows
        if (length(sccs) > 1) {
            keepers <- Reduce(intersect, purrr::map(sccs, "crd"))
            for (i in seq_along(sccs)) {
                sccs[[i]] <- sccs[[i]][sccs[[i]]$crd %in% keepers,]
            }
        }

        #--------------------------------------------------
        # Reproduce basis functions
        s <- as.matrix(sccs$z1$crd, ncol = 1)
        X <-
            fda::getbasismatrix(s, fda::create.bspline.basis(range(s), nbasis = nb, norder = 4))

        #--------------------------------------------------
        # Compute indicator theta for testing
        theta <- matrix(NA, nrow(fit$pred) - burnin, ncol(fit$pred))
        for (i in 1:nrow(theta)) {
            theta[i,] <- fit$pred[i + burnin,] < X %*% fit$beta[i + burnin,]
        }

        #--------------------------------------------------
        # Compute rejections
        if (method == "FDR") {
            rej_FDR <- FDR(theta, alpha = alpha, nthresh = 100)$reject
        } else if(method == "FDX") {
            rej_FDR <-
                FDX(
                    theta,
                    alpha = alpha,
                    beta = 1 - alpha,
                    nthresh = 100
                )$reject
        } else {
            rej_FDR <- dplyr::between(rank(-colMeans(theta)), 1, num_rejections)
        }

        fdr <- sccs$z1$crd[rej_FDR]
        sub <- fdr[fdr %in% zoom]
        segdf <-
            tibble(
                x1 = sub,
                y1 = sub,
                x2 = sapply(sub, function(X)
                    min(X + win_size, max(zoom))),
                y2 = sapply(sub, function(X)
                    min(X + win_size, max(zoom)))
            )
        diffplot + geom_segment(
            aes(
                x = x1,
                xend = x2,
                y = y1,
                yend = y1
            ),
            data = segdf,
            inherit.aes = FALSE,
            color = rej_color,
            alpha = .8,
            size = 1.1
        ) +
            geom_segment(
                aes(
                    x = x2,
                    xend = x2,
                    y = y1,
                    yend = y2
                ),
                data = segdf,
                inherit.aes = FALSE,
                color = rej_color,
                alpha = .8,
                size = 1.1
            )
    }



# #--------------------------------------------------
# # Test stuff
# plot_rejections(
#     fit,
#     sccs,
#     h = 2,
#     nb = nb,
#     peaks = peaks,
#     burnin = 0,
#     alpha = 0.2
# )
# ppp <-
#     plot_heat_rejections(
#         hic1 = KO1,
#         hic2 = ui1,
#         fit = fit1,
#         zoom = 200:400,
#         h = h,
#         nb = nb,
#         win_size = 20,
#         sccs = sccs,
#         method = "FDR",
#         alpha = 0.2
#     )
# p3 <- add_rejections(
#     diffplot = p2,#ppp[[3]],
#     fit = fit3,
#     zoom = 200:400,
#     nb = nb,
#     win_size = 5,
#     sccs = sccs,
#     method = "FDR",
#     burnin = 0,
#     alpha = 0.2,
#     rej_color = "black"
# )
