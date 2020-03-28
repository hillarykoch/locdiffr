// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>
#include <lemon/list_graph.h>
#include <lemon/lgf_reader.h>
#include <iostream>
#include <fstream>
#include <string>
#include "dagger.h"

using namespace Rcpp;
using namespace RcppArmadillo;
using namespace lemon;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// For making adjacency matrix, to then produce LGF
// [[Rcpp::export]]
arma::sp_mat cmake_adj_mat(arma::sp_mat adj_mat,
                            arma::colvec locsj,
                            arma::colvec winsizes,
                            arma::colvec cumsum_nlocs,
                            int j) {
    int diff_win = winsizes(j) - winsizes(j-1);
    for(int i = 1; i <= locsj.size(); i++) {
        if(i <= diff_win) {
            adj_mat.row(cumsum_nlocs(j-1) + i - 1).cols(cumsum_nlocs(j), cumsum_nlocs(j) + i - 1).ones();
        } else if(i >= (locsj.size() - diff_win + 1)) {
            adj_mat.row(cumsum_nlocs(j-1) + i - 1).cols(cumsum_nlocs(j) + i - 1 - diff_win, cumsum_nlocs(j+1) - 1).ones();
        } else {
            adj_mat.row(cumsum_nlocs(j-1) + i - 1).cols(cumsum_nlocs(j) + i - 1 - diff_win, cumsum_nlocs(j) + i - 1).ones();
        }
    }
    return adj_mat;
}

// Run DAGGER
// [[Rcpp::export]]
arma::mat ctest_hierarchically(std::string filepath, double alpha, arma::colvec prob_theta_equals_zero, arma::colvec rank_map) {
    ListDigraph gr;
    ListDigraph::NodeMap<int> layer(gr);
    ListDigraph::NodeMap<int> label(gr);
    ListDigraph::NodeMap<int> depth(gr);
    std::fstream infile;
    infile.open(filepath);
    int max_depth;

    // Read in Directed Graph from lgf.txt
    //  dim gives which "layer" the given node lies in
    //  max_depth is the number of layers in the graph
    digraphReader(gr, infile)
        .nodeMap("label", label)
        .nodeMap("layer", layer)
        .attribute("max_depth", max_depth)
        .run();

    // Create a dagger object from the graph and observed probabilites
    Dagger dagger(gr, layer, label, prob_theta_equals_zero, rank_map, alpha, max_depth);

    // Run the testing procedure
    int ncan;
    // For each layer, beginning at the root
    for(auto d = 0; d < max_depth; d++) {
        ncan = dagger.num_currently_candidates();
        std::cout << "\nncan = " << ncan << std::endl;
        if(ncan > 0) {
            dagger.test_hypothesis_at_current_layer();
        }
    }
    

    // Matrix whose first column is node label,
    //  and second column is whether or not that node is rejected
    arma::colvec out(dagger.total_nodes(), arma::fill::zeros);
    for(ListDigraph::NodeIt n(gr); n != INVALID; ++n){
        out(dagger.label(n)-1) = dagger.is_rejected(n);
    }
    return out;
}
