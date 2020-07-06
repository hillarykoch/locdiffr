make_BED_from_rejections <- function(infiles,
                                     chromosomes,
                                     null_outfile,
                                     rejected_outfile,
                                     resolution,
                                     region_size = 5000
                                     ) {
    # infiles: path to output from one of test_by_wFDR or test_by_wFDX
    # chromosomes: (integer or string) vector of chromosome numbers, matched in the same order as the infiles
    # null_outfile: path to file where BED file will be written for null regions
    # rejected_outfile: path to file where BED file will be written for rejected regions
    # resolution: resolution of the Hi-C data
    # region_size: number of base pairs to fall into each region of the BED file

    if(!(resolution %% region_size == 0)) {
        stop("resolution must be a multiple of region_size.")
    }

    rej_wfd <- map(infiles, readRDS) %>%
        setNames(chromosomes)

    rejected_crds <- null_crds <- list()
    for(i in seq_along(rej_wfd)) {
        # All loci falling wihtin a rejection region
        rejected_crds[[i]] <- purrr::map2(rej_wfd[[i]], as.numeric(names(rej_wfd[[i]])),
                                          ~ unique(as.vector(sapply(.x$crd[.x$reject], function(Z)
                                              Z:(Z + .y - 1))))) %>%
            unlist %>%
            unique %>%
            sort

        # All loci not falling within a rejection region
        null_crds[[i]] <- purrr::map2(rej_wfd[[i]], as.numeric(names(rej_wfd[[i]])),
                                      ~ unique(as.vector(sapply(.x$crd, function(Z)
                                          Z:(Z + .y - 1))))) %>%
            Reduce(`unique`, .) %>%
            setdiff(rejected_crds[[i]])
    }

    #-------------------------------------------------------------------------------
    # Sample null regions and rejected regions 1kb in size for plotting with deeptools
    #   This will just be points within rejected regions,
    #       as opposed to points near the boundaries
    #-------------------------------------------------------------------------------
    rej_kb_idx <-
        map(rejected_crds, ~ rep(round(
            seq(1, resolution - region_size + 1, length.out = resolution / region_size)
        ), times = length(.x)))

    null_kb_idx <- map(null_crds, ~ rep(round(
        seq(1, resolution - region_size + 1, length.out = resolution / region_size)
    ), times = length(.x)))

    rej_bed <- null_bed <- list()
    for(i in seq_along(rejected_crds)) {
        rej_bed[[i]] <- map2_dfr(
            rej_kb_idx[[i]],
            rep(rejected_crds[[i]], each = resolution / region_size),
            ~ tibble(
                "chr" = paste0("chr", chromosomes[i]),
                "start" = .y * resolution + .x,
                "stop" = .y * resolution + .x + region_size - 1
            )
        ) %>%
            dplyr::arrange(start)

        null_bed[[i]] <- map2_dfr(
            null_kb_idx[[i]],
            rep(null_crds[[i]], each = resolution / region_size),
            ~ tibble(
                "chr" = paste0("chr", chromosomes[i]),
                "start" = .y * resolution + .x,
                "stop" = .y * resolution + .x + region_size - 1
            )
        ) %>%
            dplyr::arrange(start)
    }

    readr::write_tsv(
        x = dplyr::bind_rows(null_bed),
        path = null_outfile,
        col_names = FALSE
    )
    readr::write_tsv(
        x = dplyr::bind_rows(rej_bed),
        path = rejected_outfile,
        col_names = FALSE
    )
}
