#-------------------------------------------------------------------------------
# Make adjacency matrices based on window sizes and bin locations
#-------------------------------------------------------------------------------
make_adj_matrix_from_z_list <- function(z_list) {
    locs <- map(z_list, "z1") %>% map("crd")
    winsizes <- as.numeric(names(locs))
    cumsum_nlocs <- c(0, map(locs, length) %>% unlist %>% cumsum)
    nlocs <- cumsum_nlocs[length(cumsum_nlocs)]
    
    adj_mat <- Matrix::sparseMatrix(i = 1, j = cumsum_nlocs[2] + 1, x = 1, dims = c(nlocs, nlocs))
    
    for(j in 1:(length(locs)-1)) {
        for(i in seq_along(locs[[j]])) {
            if(i <= diff(winsizes[j:(j+1)])) {
                # starting boundary effect
                adj_mat[cumsum_nlocs[j] + i, (cumsum_nlocs[j + 1] + 1):(cumsum_nlocs[j + 1] + i)] <- 1
            } else if(i >= ( length(locs[[j]]) - diff(winsizes[j:(j+1)]) + 1) ) {
                # ending boundary effect
                adj_mat[cumsum_nlocs[j] + i, (cumsum_nlocs[j + 1] + i - diff(winsizes[j:(j+1)])):cumsum_nlocs[j + 2] ] <- 1
            } else {
                # regular middle windows
                adj_mat[cumsum_nlocs[j] + i, (cumsum_nlocs[j + 1] + i - diff(winsizes[j:(j+1)])):(cumsum_nlocs[j + 1] + i)] <- 1
            }
        } 
    }
    list("adj" = adj_mat, "layers" = map2(seq_along(locs), locs, ~ rep(.x, times = length(.y))) %>% unlist)
}

#-------------------------------------------------------------------------------
# Write LGF file, for subsequent computation of effective node/leaf quantities
#-------------------------------------------------------------------------------
write_LGF <- function(z_list, path) {
    dag_info <- make_adj_matrix_from_z_list(z_list)
    dag <- graph_from_adjacency_matrix(dag_info$adj, mode = "directed")
    
    # For the @nodes portion of the LGF file
    node <-
        data.frame(
            "label" = seq_along(dag_info$layers),
            "layer" = dag_info$layers
        )
    
    # For the @arcs portion of the LGF file
    edges <- E(dag)
    arcs <- data.frame(
        "source" = as.numeric(tail_of(dag, edges)),
        "target" = as.numeric(head_of(dag, edges))
    )
    
    # For the @attributes portion of the LGF file
    attrib <- data.frame("type" = "max_depth",
                         "label" = max(node$layer))
    
    # Write the LGF file
    cat("Writing LGF file...")
    readr::write_tsv(data.frame("@nodes"),
                     path = path,
                     col_names = FALSE)
    readr::write_tsv(node,
                     path = path,
                     col_names = TRUE,
                     append = TRUE)
    cat("\n", file = path, append = TRUE)
    
    readr::write_tsv(
        data.frame("@arcs"),
        path = path,
        col_names = FALSE,
        append = TRUE
    )
    cat("\t\t -\n", file = path, append = TRUE)
    readr::write_tsv(
        arcs,
        path = path,
        col_names = FALSE,
        append = TRUE,
        na = ""
    )
    cat("\n", file = path, append = TRUE)
    readr::write_tsv(
        data.frame("@attributes"),
        path = path,
        col_names = FALSE,
        append = TRUE
    )
    readr::write_tsv(attrib,
                     path = path,
                     col_names = FALSE,
                     append = TRUE)
    cat("done!\n")
}