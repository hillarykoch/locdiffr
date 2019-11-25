# Get z-scores from SCCs
stripones <- function(sc) {
    # loc1 <- min(which(cumsum(sc - 1) != 0))
    r_sc <- rev(sc)
    loc2 <- length(r_sc) - min(which(cumsum(r_sc) != 0)) + 1

    # return(c(loc1, loc2))
    return(c(1,loc2))
}

get_z <- function(ys) {
    locs <- stripones(ys)
    data.frame(crd = locs[1]:locs[2], scc = ys[locs[1]:locs[2]]) %>%
        dplyr::mutate(z_s = fishT(scc))
}

fishT <- function(y) {
    y[y >= 1 - 1e-06] <- 1 - 1e-06
    y[y <= -(1 - 1e-06)] = -(1 - 1e-06)
    .5 * (log((1 + y) / (1 - y)))
}

# scan throught chromsome using different window sizes
get_loc_sim <-
    function(mat1,
             mat2,
             h = 5,
             resol = 40000,
             win_min = 6,
             win_max = 50,
             CI = 0) {
        smoo1 <- fastMeanFilter(mat1, h)
        smoo2 <- fastMeanFilter(mat2, h)
        nd <- nrow(mat1)
        scoremat <- stdmat <-  matrix(0, nd, win_max)

        for (i in win_min:win_max) {
            l <- nd - i + 1
            for (j in seq_len(l)) {
                win <- j:(j + i - 1)
                if (sum(smoo1[win, win]) == 0 |
                    sum(smoo2[win, win]) == 0) {
                    scoremat[j, i] <- 1
                    stdmat[j, i] <- 0
                } else {
                    LSIM <- hic_ld(smoo1, smoo2, win, resol)
                    scoremat[j, i] = LSIM$scc
                    stdmat[j, i] = LSIM$std
                }
            }
        }

        scoremat <- scoremat[, -(1:(win_min - 1))]
        stdmat <- stdmat[, -(1:(win_min - 1))]

        scoremat[which(is.na(scoremat), arr.ind = TRUE)] <- 1
        stdmat[which(is.na(stdmat), arr.ind = TRUE)] <- 0

        if (CI == 1) {
            lwb <- scoremat - 3 * stdmat
            list(scoremat, lwb)
        } else {
            scoremat
        }
    }

hic_ld <- function(smoo1, smoo2, rang, resol = 40000) {
    sub_smoo1 <- smoo1[rang, rang]
    sub_smoo2 <- smoo2[rang, rang]

    get_scc(
        sub_smoo1,
        sub_smoo2,
        resol,
        h = 0,
        lbr = resol,
        ubr = nrow(sub_smoo1) * resol
    )
}

# Add the procedure to train h, now assume 5
get_scc <- function(mat1, mat2, resol, h, lbr = 0, ubr = 5e+06) {
    if (h == 0) {
        smoo_R1 <- mat1
        smoo_R2 <- mat2
    }
    else {
        smoo_R1 <- fastMeanFilter(as.matrix(mat1), h)
        smoo_R2 <- fastMeanFilter(as.matrix(mat2), h)
    }
    rm(mat1)
    rm(mat2)
    lb <- floor(lbr/resol)
    ub <- floor(ubr/resol)

    nr <- nrow(smoo_R1)
    corr <- wei <- n <- array(ub - lb + 1)

    st <- sapply(seq(lb, (ub-1)), est_scc, ub, smoo_R1, smoo_R2, nr)
    corr0 <- unlist(st[1, ])
    wei0 <- unlist(st[2, ])
    corr <- corr0[!is.na(corr0)]
    wei <- wei0[!is.na(wei0)]
    scc <- corr %*% wei/sum(wei)
    std <- sqrt(sum(wei^2 * var(corr))/(sum(wei))^2)

    list(corr = corr, wei = wei, scc = scc, std = std)
}

est_scc <- function(dist, ub, smoo1, smoo2, nr){
    if (dist < ub - 1){
        ffd1 <- ffd2 <- NULL
        for (i in 1:(ncol(smoo1) - dist)) {
            ffd1 <- c(ffd1, smoo1[i + dist, i])
            ffd2 <- c(ffd2, smoo2[i + dist, i])
            filt <- which(ffd1 == 0 & ffd2 == 0)
            if (length(filt) == 0) {
                ffd <- cbind(ffd1, ffd2)
            } else{
                ffd <- cbind(ffd1[-filt], ffd2[-filt])
            }
        }
    } else if (dist == (ub - 1)){
        ffd1 <- c(smoo1[nr, 1], smoo1[nr-1, 1], smoo1[nr, 2])
        ffd2 <- c(smoo2[nr, 1], smoo2[nr-1, 1], smoo2[nr, 2])
        filt <- which(ffd1 == 0 & ffd2 == 0)
        if (length(filt) == 0) {
            ffd <- cbind(ffd1, ffd2)
        } else{
            ffd <- cbind(ffd1[-filt], ffd2[-filt])
        }
    }

    if (nrow(ffd) != 0) {
        n <- nrow(ffd)
        nd <- vstran(ffd)
        if (length(unique(ffd[, 1])) != 1 & length(unique(ffd[, 2])) != 1) {
            corr <- cor(ffd[, 1], ffd[, 2])
            wei <- sqrt(var(nd[, 1]) * var(nd[, 2])) * n
        } else {
            corr <- NA
            wei <- NA
        }
    } else {
        corr <- NA
        wei <- NA
    }
    list(corr = corr, wei = wei)
}

