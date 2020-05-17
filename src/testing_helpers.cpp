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
        idx = find(rej_prob >= thresh[j]);

        bwfdr[j] = 1 - accu(weighted_rej.elem(idx)) / accu(cluster_size_vec.elem(idx));
    }

    return(bwfdr);
}


// For speeding up loop computing BwFDX at various thresholds
// [[Rcpp::export]]
arma::colvec ccompute_bwfdx(arma::colvec weighted_rej,
                            arma::colvec rej_prob,
                            arma::colvec thresh,
                            arma::colvec cluster_size_vec,
                            int bootstrap_replicates,
                            double beta) {
    int nthresh = thresh.size();
    int nsamp;
    arma::uvec xceeds;
    arma::colvec bwfdx(nthresh, arma::fill::ones);
    arma::colvec boot(bootstrap_replicates, arma::fill::none);

    arma::uvec samp;
    arma::uvec compare_to_beta;
    arma::colvec dbl_compare;


    for(auto j = 0; j < nthresh; j++) {
        xceeds = find(rej_prob >= thresh[j]);

        if(xceeds.size() > 1) {
            nsamp = round(xceeds.size() / 100);

            for(auto b = 0; b < bootstrap_replicates; b++) {
                samp = RcppArmadillo::sample(xceeds, nsamp, true);

                boot[b] = 1 - accu(weighted_rej.elem(samp)) / accu(cluster_size_vec.elem(samp));
            }

            compare_to_beta = boot > beta;
            dbl_compare = arma::conv_to<arma::colvec>::from(compare_to_beta);
            bwfdx[j] = mean(dbl_compare);

        }
    }

    return(bwfdx);

}




