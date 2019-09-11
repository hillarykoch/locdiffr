// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace RcppArmadillo;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]


// [[Rcpp::export]]
arma::mat cpwdist(arma::colvec locs1, arma::colvec locs2) {
    int nrow = locs1.size();
    int ncol = locs2.size();

    arma::mat out(nrow, ncol, arma::fill::none);

    for(int i = 0; i < nrow; i++) {
        for(int j = 0; j < ncol; j++) {
            out(i,j) = std::abs(locs1(i) - locs2(j));
        }
    }

    return out;
}