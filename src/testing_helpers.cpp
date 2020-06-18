// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace RcppArmadillo;


// [[Rcpp::plugins(cpp11)]]


// For speeding up loop computing BwFDR at various thresholds
// [[Rcpp::export]]
arma::colvec ccompute_bwfdr(arma::colvec weighted_rej,
                            arma::colvec rej_prob,
                            arma::colvec thresh,
                            arma::colvec cluster_size_vec) {
    int nthresh = thresh.size();
    arma::uvec idx;
    arma::colvec bwfdr(nthresh, arma::fill::none);

    for(auto j = 0; j < nthresh; j++) {
        idx = find(rej_prob > thresh[j]);

        bwfdr[j] = 1 - accu(weighted_rej.elem(idx)) / accu(cluster_size_vec.elem(idx));
    }

    return(bwfdr);
}


// For speeding up loop computing BwFDX at various thresholds
// [[Rcpp::export]]
arma::colvec ccompute_bwfdx(arma::colvec rej_prob,
                            arma::colvec thresh,
                            arma::colvec cluster_size_vec,
                            double beta,
                            arma::mat bootstrapped_theta_mat) { // W x BOOT matrix of prob(theta=1)
    int nthresh = thresh.size();
    arma::uvec xceeds;
    arma::colvec bwfdx(nthresh, arma::fill::zeros);

    int bootstrap_replicates = bootstrapped_theta_mat.n_cols;
    arma::colvec boot(bootstrap_replicates, arma::fill::none);
    arma::mat subm;

    arma::uvec compare_to_beta;
    arma::colvec dbl_compare;

    for(auto j = 0; j < nthresh; j++) {
        xceeds = find(rej_prob > thresh[j]);

        if(xceeds.size() >= 1) {
            subm = bootstrapped_theta_mat.rows(xceeds);
            boot = (1 - (sum(subm, 0) / accu(cluster_size_vec.elem(xceeds)))).t();

            compare_to_beta = boot > beta;
            dbl_compare = arma::conv_to<arma::colvec>::from(compare_to_beta);
            bwfdx[j] = mean(dbl_compare);
        }
    }
    return(bwfdx);
}




