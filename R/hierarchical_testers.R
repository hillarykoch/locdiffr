#-------------------------------------------------------------------------------
# Functions for hierarchical testing with DAGGER (Ramdas et al. 2018)
#-------------------------------------------------------------------------------

make_adj_matrix_from_z_list <- function(z_list, irrep_list) {
    locs <- map(z_list, "z1") %>% map("crd")
    winsizes <- as.numeric(names(locs))
    cumsum_nlocs <- c(0, map(locs, length) %>% unlist %>% cumsum)
    nlocs <- cumsum_nlocs[length(cumsum_nlocs)]
    
    adj_mat <-
        Matrix::sparseMatrix(
            i = 1,
            j = cumsum_nlocs[2] + 1,
            x = 1,
            dims = c(nlocs, nlocs)
        )
    for (j in 1:(length(locs) - 1)) {
        adj_mat <-
            cmake_adj_mat(adj_mat, locs[[j]], winsizes, cumsum_nlocs, j)
    }
    
    # Remove irreproducible sites from the adjacency matrix
    irrep_idx <- map(irrep_list, which) %>%
        map2(.y = head(cumsum_nlocs,-1), ~ .x + .y) %>%
        unlist
    
    list(
        "adj" = adj_mat[-irrep_idx,-irrep_idx],
        "layers" = map2(seq_along(locs), irrep_list, ~ rep(.x, times = sum(!.y))) %>% unlist
    )
}

write_LGF <- function(z_list, irrep_list, path) {
    dag_info <- make_adj_matrix_from_z_list(z_list, irrep_list)
    dag <-
        igraph::graph_from_adjacency_matrix(dag_info$adj, mode = "directed")
    
    # For the @nodes portion of the LGF file
    node <-
        data.frame("label" = seq_along(dag_info$layers),
                   "layer" = dag_info$layers)
    
    # For the @arcs portion of the LGF file
    edges <- igraph::E(dag)
    arcs <- data.frame(
        "source" = as.numeric(igraph::tail_of(dag, edges)),
        "target" = as.numeric(igraph::head_of(dag, edges))
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