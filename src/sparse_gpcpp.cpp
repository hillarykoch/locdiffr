// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>
#include <random>
#include <boost/math/distributions/normal.hpp>
#include "myrng.h"

using namespace Rcpp;
using namespace RcppArmadillo;
using namespace boost::math;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]

// accesses myrng.h to make a random uniform variable
// [[Rcpp::export]]
arma::colvec crunif(unsigned int n, unsigned int seed) {
    Uniform unif(seed);
    arma::colvec out(n, arma::fill::none);
    for(auto i = 0; i < n; i++) {
        out[i] = unif.draw();
    }
    return out;
}

// Pseudocode 2 of Finley et al. 2019
// [[Rcpp::export]]
Rcpp::List csolve_for_A_and_D(arma::sp_mat& cov_cur,
                              Rcpp::List& neighbor_list) {
    int s = cov_cur.n_rows;
    arma::sp_mat A = arma::sp_mat(s, s);
    arma::sp_mat D = arma::sp_mat(s, s);
    arma::uvec N;
    int sizen;
    D[0] = cov_cur[0];

    for(auto i = 0; i < (s-1); i++) {
        // minus 1 because 0 vs 1-based indexing
        N = as<arma::uvec>(neighbor_list[i+1]) - 1;
        sizen = N.size()-1;
        A.row(i+1).cols(N[0], N[sizen]) = arma::solve(arma::mat(cov_cur.cols(N[0], N[sizen]).rows(N[0], N[sizen])), arma::colvec(cov_cur.col(i+1).rows(N[0], N[sizen]))).t();
        D.row(i+1).col(i+1) = cov_cur.col(i+1).row(i+1) - dot(arma::rowvec(cov_cur.cols(N[0], N[sizen]).row(i+1)), A.cols(N[0], N[sizen]).row(i+1));
    }
    return Rcpp::List::create(
        Rcpp::Named("A") = A,
        Rcpp::Named("D") = D);
}

// [[Rcpp::export]]
double csparse_quadratic_form_symm(arma::colvec& u,
                                      arma::sp_mat& A,
                                      arma::colvec& D,
                                      Rcpp::List& neighbor_list) {
    double tmp = pow(u[0], 2) / D[0];
    int n = u.n_elem;
    arma::uvec N;
    int sizen;
    for(auto i = 1; i < n; i++) {
        N = as<arma::uvec>(neighbor_list[i]) - 1;
        sizen = N.size()-1;
        tmp += pow(u[i] - arma::dot(A.cols(N[0], N[sizen]).row(i), u.elem(N).t()), 2) / D[i];
    }
    return(tmp);
}

// [[Rcpp::export]]
double csparse_quadratic_form_asymm(arma::colvec u,
                                    arma::colvec v,
                                    arma::sp_mat A,
                                    arma::colvec D,
                                    Rcpp::List neighbor_list) {
    double tmp = u[0] * v[0] / D[0];
    int n = u.n_elem;
    arma::uvec N;
    int sizen;
    for(auto i = 1; i < n; i++) {
        N = as<arma::uvec>(neighbor_list[i]) - 1;
        sizen = N.size()-1;
        tmp += (u[i] - arma::dot(A.cols(N[0], N[sizen]).row(i), u.elem(N).t())) *
            (v[i] - arma::dot(A.cols(N[0], N[sizen]).row(i), v.elem(N).t())) / D[i];
    }
    return(tmp);
}

// [[Rcpp::export]]
Rcpp::List csolve_for_B_and_b(arma::mat& y,
                              arma::mat& X,
                              arma::sp_mat& A,
                              arma::colvec& D,
                              Rcpp::List& neighbor_list,
                              arma::mat& precision_beta) {
    int p = X.n_cols;
    arma::colvec b(p, arma::fill::ones);
    arma::mat B(p, p, arma::fill::ones);
    arma::colvec ybar = mean(y, 1);

    for(auto i = 0; i < p; i++) {
        b[i] = csparse_quadratic_form_asymm(X.col(i), ybar, A, D, neighbor_list);
        for(auto j = 0; j < p; j++) {
            B.col(j).row(i) = csparse_quadratic_form_asymm(X.col(i), X.col(j), A, D, neighbor_list);
        }
    }


    return(Rcpp::List::create(
            Rcpp::Named("B") = B + precision_beta,
            Rcpp::Named("b") = b)
               );
}

// [[Rcpp::export]]
arma::colvec cmake_one_pred_sparse(Rcpp::List& neighbor_list,
                                arma::mat& y,
                                arma::mat& s,
                                arma::mat& X,
                                arma::mat& cond_cov,
                                unsigned int SEED) {
    int n = y.n_rows;
    arma::colvec preds(n, arma::fill::none);
    arma::colvec z(n, arma::fill::none);
    arma::mat L;

    arma::colvec uniform_rv;
    uniform_rv = crunif(n, SEED);
    normal standard_normal;

    for(auto i = 0; i < n; i++) {
        z(i) = quantile(standard_normal, uniform_rv[i]);
    }

    arma::uvec Nidx;
    int sizen;
    preds[0] = sqrt(cond_cov[0]) * z[0];
    for(auto i = 1; i < n; i++) {
        Nidx = as<arma::uvec>(neighbor_list[i]) - 1;
        sizen = Nidx.size()-1;
        L = arma::chol(cond_cov.cols(Nidx[0], Nidx[sizen]).rows(Nidx[0], Nidx[sizen]), "upper");
        preds[i] = arma::dot(L.col(L.n_cols-1), z.subvec(Nidx[0], Nidx[sizen]));
    }

    return(preds);
}

// [[Rcpp::export]]
Rcpp::List csolve_for_A_and_D_2d(arma::sp_mat& cov_cur,
                                 Rcpp::List& neighbor_list) {
    int s = cov_cur.n_rows;
    arma::mat cov_cur_dense = arma::mat(cov_cur);

    arma::mat A = arma::mat(s, s);
    arma::sp_mat D = arma::sp_mat(s, s);
    arma::uvec N;

    D[0] = cov_cur[0];
    arma::uvec v = { 1 };
    for(auto i = 0; i < (s-1); i++) {
        // minus 1 because 0 vs 1-based indexing
        N = as<arma::uvec>(neighbor_list[i+1]) - 1;
        A.submat(v, N) = arma::solve(cov_cur_dense.submat(N, N), arma::colvec(cov_cur_dense.submat(N, v))).t();
        D.row(i+1).col(i+1) = cov_cur_dense.col(i+1).row(i+1) - dot(arma::rowvec(cov_cur_dense.submat(v, N)), A.submat(v, N));
        v = v+1;
    }

    return Rcpp::List::create(
        Rcpp::Named("A") = arma::sp_mat(A),
        Rcpp::Named("D") = D);
}

// [[Rcpp::export]]
double csparse_quadratic_form_symm_2d(arma::colvec& u,
                                   arma::sp_mat& A,
                                   arma::colvec& D,
                                   Rcpp::List& neighbor_list) {
    double tmp = pow(u[0], 2) / D[0];
    int n = u.n_elem;
    arma::uvec N;
    int sizen;

    arma::mat A_dense = arma::mat(A);
    arma::uvec v = { 1 };
    for(auto i = 1; i < n; i++) {
        N = as<arma::uvec>(neighbor_list[i]) - 1;
        sizen = N.size()-1;
        tmp += pow(u[i] - arma::dot(A_dense.submat(v, N), u.elem(N).t()), 2) / D[i];
        v = v+1;
    }
    return(tmp);
}

// [[Rcpp::export]]
double csparse_quadratic_form_asymm_2d(arma::colvec u,
                                    arma::colvec v,
                                    arma::sp_mat A,
                                    arma::colvec D,
                                    Rcpp::List neighbor_list) {
    double tmp = u[0] * v[0] / D[0];
    int n = u.n_elem;
    arma::uvec N;
    int sizen;

    arma::mat A_dense = arma::mat(A);
    arma::uvec vv = { 1 };
    for(auto i = 1; i < n; i++) {
        N = as<arma::uvec>(neighbor_list[i]) - 1;
        sizen = N.size()-1;
        tmp += (u[i] - arma::dot(A_dense.submat(vv, N), u.elem(N).t())) *
            (v[i] - arma::dot(A_dense.submat(vv, N), v.elem(N).t())) / D[i];
        vv = vv+1;
    }
    return(tmp);
}

// [[Rcpp::export]]
Rcpp::List csolve_for_B_and_b_2d(arma::mat& y,
                              arma::mat& X,
                              arma::sp_mat& A,
                              arma::colvec& D,
                              Rcpp::List& neighbor_list,
                              arma::mat& precision_beta) {
    int p = X.n_cols;
    arma::colvec b(p, arma::fill::ones);
    arma::mat B(p, p, arma::fill::ones);
    arma::colvec ybar = mean(y, 1);

    for(auto i = 0; i < p; i++) {
        b[i] = csparse_quadratic_form_asymm_2d(X.col(i), ybar, A, D, neighbor_list);
        for(auto j = 0; j < p; j++) {
            B.col(j).row(i) = csparse_quadratic_form_asymm_2d(X.col(i), X.col(j), A, D, neighbor_list);
        }
    }


    return(Rcpp::List::create(
            Rcpp::Named("B") = B + precision_beta,
            Rcpp::Named("b") = b)
    );
}
