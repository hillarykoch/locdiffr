// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// cinv
arma::mat cinv(arma::mat covmat);
RcppExport SEXP _sgp_cinv(SEXP covmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type covmat(covmatSEXP);
    rcpp_result_gen = Rcpp::wrap(cinv(covmat));
    return rcpp_result_gen;
END_RCPP
}
// cifelse
arma::mat cifelse(const arma::mat& d);
RcppExport SEXP _sgp_cifelse(SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(cifelse(d));
    return rcpp_result_gen;
END_RCPP
}
// cgeteigs
arma::colvec cgeteigs(const arma::mat& Q);
RcppExport SEXP _sgp_cgeteigs(SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(cgeteigs(Q));
    return rcpp_result_gen;
END_RCPP
}
// cifelse_eigvals
arma::mat cifelse_eigvals(arma::colvec eigvals, const double thresh);
RcppExport SEXP _sgp_cifelse_eigvals(SEXP eigvalsSEXP, SEXP threshSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type eigvals(eigvalsSEXP);
    Rcpp::traits::input_parameter< const double >::type thresh(threshSEXP);
    rcpp_result_gen = Rcpp::wrap(cifelse_eigvals(eigvals, thresh));
    return rcpp_result_gen;
END_RCPP
}
// get_ldet
double get_ldet(arma::mat Q);
RcppExport SEXP _sgp_get_ldet(SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(get_ldet(Q));
    return rcpp_result_gen;
END_RCPP
}
// fastMeanFilter
NumericMatrix fastMeanFilter(NumericMatrix mat, int h);
RcppExport SEXP _sgp_fastMeanFilter(SEXP matSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< int >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(fastMeanFilter(mat, h));
    return rcpp_result_gen;
END_RCPP
}
// cpwdist
arma::mat cpwdist(arma::colvec locs1, arma::colvec locs2);
RcppExport SEXP _sgp_cpwdist(SEXP locs1SEXP, SEXP locs2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type locs1(locs1SEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type locs2(locs2SEXP);
    rcpp_result_gen = Rcpp::wrap(cpwdist(locs1, locs2));
    return rcpp_result_gen;
END_RCPP
}
// crunif
arma::colvec crunif(unsigned int n, unsigned int seed);
RcppExport SEXP _sgp_crunif(SEXP nSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type n(nSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(crunif(n, seed));
    return rcpp_result_gen;
END_RCPP
}
// csolve_for_A_and_D
Rcpp::List csolve_for_A_and_D(arma::sp_mat& cov_cur, Rcpp::List& neighbor_list);
RcppExport SEXP _sgp_csolve_for_A_and_D(SEXP cov_curSEXP, SEXP neighbor_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type cov_cur(cov_curSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type neighbor_list(neighbor_listSEXP);
    rcpp_result_gen = Rcpp::wrap(csolve_for_A_and_D(cov_cur, neighbor_list));
    return rcpp_result_gen;
END_RCPP
}
// csparse_quadratic_form_symm
double csparse_quadratic_form_symm(arma::colvec& u, arma::sp_mat& A, arma::colvec& D, Rcpp::List& neighbor_list);
RcppExport SEXP _sgp_csparse_quadratic_form_symm(SEXP uSEXP, SEXP ASEXP, SEXP DSEXP, SEXP neighbor_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec& >::type u(uSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type D(DSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type neighbor_list(neighbor_listSEXP);
    rcpp_result_gen = Rcpp::wrap(csparse_quadratic_form_symm(u, A, D, neighbor_list));
    return rcpp_result_gen;
END_RCPP
}
// csparse_quadratic_form_asymm
double csparse_quadratic_form_asymm(arma::colvec u, arma::colvec v, arma::sp_mat A, arma::colvec D, Rcpp::List neighbor_list);
RcppExport SEXP _sgp_csparse_quadratic_form_asymm(SEXP uSEXP, SEXP vSEXP, SEXP ASEXP, SEXP DSEXP, SEXP neighbor_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type u(uSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type v(vSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type D(DSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type neighbor_list(neighbor_listSEXP);
    rcpp_result_gen = Rcpp::wrap(csparse_quadratic_form_asymm(u, v, A, D, neighbor_list));
    return rcpp_result_gen;
END_RCPP
}
// csolve_for_B_and_b
Rcpp::List csolve_for_B_and_b(arma::mat& y, arma::mat& X, arma::sp_mat& A, arma::colvec& D, Rcpp::List& neighbor_list, arma::mat& precision_beta);
RcppExport SEXP _sgp_csolve_for_B_and_b(SEXP ySEXP, SEXP XSEXP, SEXP ASEXP, SEXP DSEXP, SEXP neighbor_listSEXP, SEXP precision_betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type D(DSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type neighbor_list(neighbor_listSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type precision_beta(precision_betaSEXP);
    rcpp_result_gen = Rcpp::wrap(csolve_for_B_and_b(y, X, A, D, neighbor_list, precision_beta));
    return rcpp_result_gen;
END_RCPP
}
// cmake_one_pred_sparse
arma::colvec cmake_one_pred_sparse(Rcpp::List& neighbor_list, arma::mat& y, arma::mat& s, arma::mat& X, arma::mat& cond_cov, unsigned int SEED);
RcppExport SEXP _sgp_cmake_one_pred_sparse(SEXP neighbor_listSEXP, SEXP ySEXP, SEXP sSEXP, SEXP XSEXP, SEXP cond_covSEXP, SEXP SEEDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type neighbor_list(neighbor_listSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type s(sSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type cond_cov(cond_covSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type SEED(SEEDSEXP);
    rcpp_result_gen = Rcpp::wrap(cmake_one_pred_sparse(neighbor_list, y, s, X, cond_cov, SEED));
    return rcpp_result_gen;
END_RCPP
}
// csolve_for_A_and_D_2d
Rcpp::List csolve_for_A_and_D_2d(arma::sp_mat& cov_cur, Rcpp::List& neighbor_list);
RcppExport SEXP _sgp_csolve_for_A_and_D_2d(SEXP cov_curSEXP, SEXP neighbor_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type cov_cur(cov_curSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type neighbor_list(neighbor_listSEXP);
    rcpp_result_gen = Rcpp::wrap(csolve_for_A_and_D_2d(cov_cur, neighbor_list));
    return rcpp_result_gen;
END_RCPP
}
// csparse_quadratic_form_symm_2d
double csparse_quadratic_form_symm_2d(arma::colvec& u, arma::sp_mat& A, arma::colvec& D, Rcpp::List& neighbor_list);
RcppExport SEXP _sgp_csparse_quadratic_form_symm_2d(SEXP uSEXP, SEXP ASEXP, SEXP DSEXP, SEXP neighbor_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec& >::type u(uSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type D(DSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type neighbor_list(neighbor_listSEXP);
    rcpp_result_gen = Rcpp::wrap(csparse_quadratic_form_symm_2d(u, A, D, neighbor_list));
    return rcpp_result_gen;
END_RCPP
}
// csparse_quadratic_form_asymm_2d
double csparse_quadratic_form_asymm_2d(arma::colvec u, arma::colvec v, arma::sp_mat A, arma::colvec D, Rcpp::List neighbor_list);
RcppExport SEXP _sgp_csparse_quadratic_form_asymm_2d(SEXP uSEXP, SEXP vSEXP, SEXP ASEXP, SEXP DSEXP, SEXP neighbor_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type u(uSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type v(vSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type D(DSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type neighbor_list(neighbor_listSEXP);
    rcpp_result_gen = Rcpp::wrap(csparse_quadratic_form_asymm_2d(u, v, A, D, neighbor_list));
    return rcpp_result_gen;
END_RCPP
}
// csolve_for_B_and_b_2d
Rcpp::List csolve_for_B_and_b_2d(arma::mat& y, arma::mat& X, arma::sp_mat& A, arma::colvec& D, Rcpp::List& neighbor_list, arma::mat& precision_beta);
RcppExport SEXP _sgp_csolve_for_B_and_b_2d(SEXP ySEXP, SEXP XSEXP, SEXP ASEXP, SEXP DSEXP, SEXP neighbor_listSEXP, SEXP precision_betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type D(DSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type neighbor_list(neighbor_listSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type precision_beta(precision_betaSEXP);
    rcpp_result_gen = Rcpp::wrap(csolve_for_B_and_b_2d(y, X, A, D, neighbor_list, precision_beta));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_sgp_cinv", (DL_FUNC) &_sgp_cinv, 1},
    {"_sgp_cifelse", (DL_FUNC) &_sgp_cifelse, 1},
    {"_sgp_cgeteigs", (DL_FUNC) &_sgp_cgeteigs, 1},
    {"_sgp_cifelse_eigvals", (DL_FUNC) &_sgp_cifelse_eigvals, 2},
    {"_sgp_get_ldet", (DL_FUNC) &_sgp_get_ldet, 1},
    {"_sgp_fastMeanFilter", (DL_FUNC) &_sgp_fastMeanFilter, 2},
    {"_sgp_cpwdist", (DL_FUNC) &_sgp_cpwdist, 2},
    {"_sgp_crunif", (DL_FUNC) &_sgp_crunif, 2},
    {"_sgp_csolve_for_A_and_D", (DL_FUNC) &_sgp_csolve_for_A_and_D, 2},
    {"_sgp_csparse_quadratic_form_symm", (DL_FUNC) &_sgp_csparse_quadratic_form_symm, 4},
    {"_sgp_csparse_quadratic_form_asymm", (DL_FUNC) &_sgp_csparse_quadratic_form_asymm, 5},
    {"_sgp_csolve_for_B_and_b", (DL_FUNC) &_sgp_csolve_for_B_and_b, 6},
    {"_sgp_cmake_one_pred_sparse", (DL_FUNC) &_sgp_cmake_one_pred_sparse, 6},
    {"_sgp_csolve_for_A_and_D_2d", (DL_FUNC) &_sgp_csolve_for_A_and_D_2d, 2},
    {"_sgp_csparse_quadratic_form_symm_2d", (DL_FUNC) &_sgp_csparse_quadratic_form_symm_2d, 4},
    {"_sgp_csparse_quadratic_form_asymm_2d", (DL_FUNC) &_sgp_csparse_quadratic_form_asymm_2d, 5},
    {"_sgp_csolve_for_B_and_b_2d", (DL_FUNC) &_sgp_csolve_for_B_and_b_2d, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_sgp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
