// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace RcppArmadillo;

// [[Rcpp::export]]
arma::mat cinv(arma::mat covmat) {
    try{
        return arma::inv_sympd(covmat);
    }
    catch(...){
        return covmat.i();
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
arma::colvec cgeteigs(const arma::mat& Q) {
    arma::colvec eigval;
    arma::eig_sym(eigval, Q);
    return eigval;
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
double get_ldet(arma::mat Q) {
    arma::colvec eigvals = cgeteigs(Q);
    arma::mat D = cifelse_eigvals(eigvals, 1e-07);
    return(accu(log(D)));
}