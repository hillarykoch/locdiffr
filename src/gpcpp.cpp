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

// [[Rcpp::export]]
arma::mat cinv(arma::mat covmat) {
    try{
        return arma::inv_sympd(covmat);
    }
    catch(...){
	    return arma::inv(covmat);
    }
}

// [[Rcpp::export]]
arma::mat cifelse(const arma::mat& d) {
    arma::uvec diagidx = find(d == 0);
    const int n = d.n_rows;
    arma::mat out(n, n, arma::fill::zeros);
    out.elem(diagidx).ones();
    return out;
}

// [[Rcpp::export]]
arma::mat cifelse_eigvals(arma::colvec eigvals, const double thresh) {
    const int n = eigvals.size();
    for(auto i = 0; i < n; i++) {
        if(eigvals(i) < thresh) {
            eigvals(i) = 1 / thresh;
        } else {
            eigvals(i) = 1 / eigvals(i);
        }
    }

    return eigvals;
}

// [[Rcpp::export]]
arma::colvec cgeteigs(const arma::mat& Q) {
    arma::colvec eigval;
    arma::eig_sym(eigval, Q);
    return eigval;
}