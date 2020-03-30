#-------------------------------------------------------------------------------
# Compute error rate based off of overlap in area of rejection region with truth
#-------------------------------------------------------------------------------

compute_area_based_error_rate <- function(z_list, rej_list, true_differences) {
    L <- max(z_list[[1]]$z1$crd) + as.numeric(names(z_list)[1]) - 1
    truth_mat <-
        Matrix::sparseMatrix(
            i = true_differences$start[1],
            j = true_differences$start[1] + 1,
            x = 1,
            dims = c(L, L)
        )
    rej_mat <-
        Matrix::sparseMatrix(
            i = which(rej_list[[1]])[1],
            j = which(rej_list[[1]])[1] + 1,
            x = 1,
            dims = c(L, L)
        )

    #---------------------------------------------------------------------------
    # Loop over all true differences, making a boolean matrix where TRUE is a difference
    truth_mat <- cpopulate_true_differences(truth_rcpp, (as.matrix(true_differences)-1))

    # Look through rejections, making a boolean matrix where TRUE is a rejection
    for(i in seq_along(rej_list)) {
        if(any(rej_list[[i]])) {
            win_size <- as.numeric(names(z_list)[i])
            rej_mat <- cpopulate_rejected_differences(rej_mat, which(rej_list[[i]]) - 1, win_size)
        }
    }


    # Make triangular (without diagonal) to make accurate comparisons
    truth_mat <- Matrix::triu(truth_mat, 1)
    rej_mat <- Matrix::triu(rej_mat, 1)

    # compute error rate quantities
    magnitude_halftri <- (L^2 - L) / 2 + L
    true_discoveries <- rej_mat & truth_mat
    false_discoveries <- rej_mat & !truth_mat
    num_true_negatives <- sum(!rej_mat & !truth_mat) - magnitude_halftri
    total_negatives <- (sum(!rej_mat) - magnitude_halftri)
    total_discoveries <- sum(rej_mat)

    FPR <- sum(false_discoveries) / total_negatives
    FDR <- sum(false_discoveries) / (sum(false_discoveries) + num_true_negatives)
    list(
        "TPR" = sum(true_discoveries) / total_discoveries,
        "FPR" = FPR,
        "TNR" = 1 - FPR,
        "FDR" = FDR,
        "precision" = 1 - FDR
    )
}
