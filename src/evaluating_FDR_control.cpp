// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>


using namespace Rcpp;
using namespace RcppArmadillo;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// put 1's where there are differences (for truth)
// [[Rcpp::export]]
arma::sp_mat cpopulate_true_differences(arma::sp_mat& truth_mat, arma::umat& true_differences) {
    for(auto i = 0; i < true_differences.n_rows; i++) {
        truth_mat.submat(true_differences(i,0), true_differences(i,0), true_differences(i,1), true_differences(i,1)).ones();
    }
    return truth_mat;
}

// put 1's where there are differences (for rejections)
// [[Rcpp::export]]
arma::sp_mat cpopulate_rejected_differences(arma::sp_mat& rej_mat, arma::colvec& rej_list, uint win_size) {
    for(auto i = 0; i < rej_list.size(); i++) {
        rej_mat.submat(rej_list(i), rej_list(i), rej_list(i) + win_size - 1, rej_list(i) + win_size - 1).ones();
    }
    return rej_mat;
}