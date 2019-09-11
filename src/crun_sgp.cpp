// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <random>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/gamma.hpp>

using namespace boost::math;
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
arma::mat cMatern(arma::mat d, double range, double smoothness) {
    // modeled after fields::Matern function
    d /= range;
    arma::uvec zeroidx = find(d == 0);
    d.elem(zeroidx) += 1e-10;
    double con = 1 / (pow(2, smoothness-1) * exp(lgamma(smoothness)));

    d.transform( [smoothness](double& dd) { return(pow(dd, smoothness) * cyl_bessel_k(smoothness, dd)); });
    return(d * con);
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
arma::mat getQ(arma::mat maternmat, arma::mat d, double r) {
    return(r * maternmat + (1-r)*cifelse(d));
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

// [[Rcpp::export]]
double mod(double a, double n) {
    return(a - trunc(a/n)*n);
}

double cupdate_var(double cur_var, double acpt_rt, double gamma1, double opt_rt = .3) {
    return(exp(log(cur_var) + gamma1 * (acpt_rt - opt_rt)));
}

// Treats as mean 0 here, and is added on to Xp %*% beta in the actual body of the code
arma::colvec cmake_pred(arma::mat y, arma::mat matern_cov, double tau) {
    arma::mat cond_cov;
    int np = matern_cov.n_rows;

    std::default_random_engine generator;
    std::uniform_real_distribution<double> runif(0.0,1.0);
    normal standard_normal;
    arma::colvec random_normals(np, arma::fill::none);

    for(auto pp = 0; pp < np; pp++) {
        random_normals(pp) = quantile(standard_normal, runif(generator));
    }

    cond_cov = matern_cov / tau;
    return(arma::mean(y,1) + arma::chol(cond_cov).t() * random_normals);
}

// [[Rcpp::export]]
Rcpp::List crun_sgp_nugget(
        double reps,
        int n,
        int p,
        double tau,
        arma::mat X,
        arma::mat PLDs_precision,
        double PLDs_ldet,
        double errprec,
        arma::mat precision_beta,
        arma::mat y,
        int iters,
        double aas,
        double bs,
        double rhos,
        double nus,
        double tune_var,
        double min_range,
        double max_range,
        double win_len,
        arma::mat d,
        arma::colvec acpt_rhos,
        double c0,
        double c1,
        double tune_k,
        arma::colvec acpt_chain,
        arma::colvec tune_var_chain) {
    arma::mat premultX;
    arma::mat varterm;
    arma::colvec muterm;
    arma::mat beta;
    arma::colvec Xb;
    arma::mat yminusXb;
    double SS;
    double SS_star;
    double SSe;
    double PLDs_ldet_star;
    double rhos_star;
    double R;
    // arma::mat Q;
    arma::mat maternmat;

    double gamma1;
    double acpt_rt_rhos;

    std::default_random_engine generator;
    std::uniform_real_distribution<double> runif(0.0,1.0);
    normal standard_normal;

    arma::colvec random_normals(p, arma::fill::none);
    arma::colvec temp1(reps, arma::fill::none);
    arma::colvec temp2(1, arma::fill::none);
    arma::mat PLDs_precision_star(n,n,arma::fill::none);

    arma::mat beta_chain(iters, p, arma::fill::none);
    arma::mat param_chain(iters, 3, arma::fill::none);

    arma::mat matern_cov;
    arma::mat ypred(iters, n, arma::fill::none);


    for(auto i = 1; i <= iters; i++) {
        //-------------------------------------------------------------
        // Update regression coefficients
        //-------------------------------------------------------------
        premultX = reps * tau * X.t() * (PLDs_precision + arma::diagmat(arma::ones(n) * errprec));
        varterm = cinv(precision_beta + premultX * X);
        muterm = premultX * mean(y, 1);

        for(auto pp = 0; pp < p; pp++) {
            random_normals(pp) = quantile(standard_normal, runif(generator));
        }

        beta = varterm * muterm + ((arma::chol(varterm)).t() * random_normals);
        Xb = vectorise(X * beta);

        //-------------------------------------------------------------
        // Update covariance parameters for the spatial signal process
        //-------------------------------------------------------------

        // tau ~ Gamma
        yminusXb = y.each_col() - Xb;
        for(auto rr = 0; rr < reps; rr++) {
            temp2 = (yminusXb.col(rr).t()) * PLDs_precision * yminusXb.col(rr);
            temp1(rr) = temp2(0);
        }
        SS = tau * arma::accu(temp1);

        gamma_distribution<> gamma_dist((n * reps) / 2 + aas, 1 / (SS / 2 + bs));
        tau = quantile(gamma_dist, runif(generator));
        SS  = tau * SS;

        //-----------------------------------------------------------------------
        // Sampling range parameter
        //-----------------------------------------------------------------------
        normal rhos_proposal(0.0, sqrt(tune_var));
        std::uniform_real_distribution<double> trunc_runif(cdf(rhos_proposal, -(rhos-min_range)), cdf(rhos_proposal, max_range - rhos));

        rhos_star = rhos + quantile(rhos_proposal, trunc_runif(generator));
        maternmat = cMatern(d, rhos_star, nus);
        // Q = getQ(maternmat, d, 1);

        PLDs_precision_star = cinv(maternmat);
        PLDs_ldet_star = get_ldet(maternmat);

        for(auto rr = 0; rr < reps; rr++) {
            temp2 = (yminusXb.col(rr).t()) * PLDs_precision_star * yminusXb.col(rr);
            temp1(rr) = temp2(0);
        }
        SS_star = tau * arma::accu(temp1);
        temp1.set_size(reps);

        R = 0.5 * (PLDs_ldet_star - PLDs_ldet) - 0.5 * (SS_star - SS);
        if(runif(generator) < exp(R)) {
            rhos = rhos_star;
            PLDs_precision = PLDs_precision_star;
            PLDs_ldet = PLDs_ldet_star;
            acpt_rhos(mod(i, win_len)) = 1;
        } else {
            acpt_rhos(mod(i, win_len)) = 0;
        }

        //-----------------------------------------------------------------------
        // Update the tuning variance
        //-----------------------------------------------------------------------
        if(i >= 100) {
            gamma1 = c0 / pow(i + tune_k, c1);
            acpt_rt_rhos = mean(acpt_rhos);
            tune_var = cupdate_var(tune_var, acpt_rt_rhos, gamma1, .3);
            acpt_chain(i-100) = acpt_rt_rhos;
            tune_var_chain(i-100) = tune_var;
        }

        //-------------------------------------------------------------
        // Update error variance term
        //-------------------------------------------------------------
        // errprec ~ Gamma
        yminusXb = y.each_col() - Xb;
        for(auto rr = 0; rr < reps; rr++) {
            temp2 = (yminusXb.col(rr).t()) * yminusXb.col(rr);
            temp1(rr) = temp2(0);
        }
        SSe = errprec * accu(temp1);

        gamma_distribution<> gamma_dist_e((n * reps) / 2 + aas, 1 / (SSe / 2 + bs));
        errprec = quantile(gamma_dist_e, runif(generator));

        //-------------------------------------------------------------
        // Update chain parameter log
        //-------------------------------------------------------------
        beta_chain.row(i-1) = beta.t();

        param_chain.row(i-1).col(0) = 1 / sqrt(tau);
        param_chain.row(i-1).col(1) = rhos;
        param_chain.row(i-1).col(2) = 1 / sqrt(errprec);

        //-------------------------------------------------------------
        // Make predictions
        //-------------------------------------------------------------
        matern_cov = cMatern(d, rhos, nus);
        ypred.row(i-1) = (X * beta + cmake_pred(yminusXb, matern_cov, tau)).t();
    }

    return (Rcpp::List::create(Rcpp::Named("beta") = beta_chain,
                              Rcpp::Named("covar_params") = param_chain,
                              Rcpp::Named("pred") = ypred,
                              Rcpp::Named("acpt_chain") = acpt_chain,
                              Rcpp::Named("tune_chain") = tune_var_chain));

}



// [[Rcpp::export]]
Rcpp::List crun_sgp_correlated(
        double reps,
        int n,
        int p,
        double tau,
        arma::mat X,
        arma::mat PLDs_precision,
        double PLDs_ldet,
        arma::mat PLDe_precision,
        double PLDe_ldet,
        double errprec,
        arma::mat precision_beta,
        arma::mat y,
        int iters,
        double aas,
        double bs,
        double rhos,
        double nus,
        double rhoe,
        double nue,
        double r,
        double sd_r,
        double tune_vars,
        double tune_vare,
        double min_range,
        double max_range,
        double win_len,
        arma::mat d,
        arma::colvec acpt_rhos,
        arma::colvec acpt_rhoe,
        double c0,
        double c1,
        double tune_k,
        arma::mat acpt_chain,
        arma::mat tune_var_chain) {

    arma::mat premultX;
    arma::mat varterm;
    arma::colvec muterm;
    arma::mat beta;
    arma::colvec Xb;
    arma::mat yminusXb;
    double SS;
    double SS_star;
    double PLDs_ldet_star;
    double rhos_star;
    double rhoe_star;
    double R;

    double lr;
    double lr_star;
    double r_star;
    arma::mat Q;
    arma::mat PLDe_precision_star(n,n,arma::fill::none);
    double PLDe_ldet_star;
    arma::mat maternmat;

    double gamma1;
    double acpt_rt_rhos;
    double acpt_rt_rhoe;

    std::default_random_engine generator;
    std::uniform_real_distribution<double> runif(0.0,1.0);
    normal lr_proposal(0.0, 0.5);
    normal standard_normal;

    arma::colvec random_normals(p, arma::fill::none);
    arma::colvec temp1(reps, arma::fill::none);
    arma::colvec temp2(1, arma::fill::none);

    arma::mat PLDs_precision_star(n,n,arma::fill::none);

    arma::mat beta_chain(iters, p, arma::fill::none);
    arma::mat param_chain(iters, 4, arma::fill::none);
    arma::mat matern_cov;
    arma::mat ypred(iters, n, arma::fill::none);

    for(auto i = 1; i <= iters; i++) {
        //-------------------------------------------------------------
        // Update regression coefficients
        //-------------------------------------------------------------
        premultX = reps * tau * X.t() * (PLDs_precision + PLDe_precision);
        varterm = cinv(precision_beta + premultX * X);
        muterm = premultX * mean(y, 1);

        for(auto pp = 0; pp < p; pp++) {
            random_normals(pp) = quantile(standard_normal, runif(generator));
        }

        beta = varterm * muterm + ((arma::chol(varterm)).t() * random_normals);
        Xb = vectorise(X * beta);

        //-------------------------------------------------------------
        // Update covariance parameters for the spatial signal process
        //-------------------------------------------------------------

        // tau ~ Gamma
        yminusXb = y.each_col() - Xb;
        for(auto rr = 0; rr < reps; rr++) {
            temp2 = (yminusXb.col(rr).t()) * PLDs_precision * yminusXb.col(rr);
            temp1(rr) = temp2(0);
        }
        SS = tau * arma::accu(temp1);

        gamma_distribution<> gamma_dist((n * reps) / 2 + aas, 1 / (SS / 2 + bs));
        tau = quantile(gamma_dist, runif(generator));
        SS  = tau * SS;

        //-----------------------------------------------------------------------
        // Sampling range parameter
        //-----------------------------------------------------------------------
        normal rhos_proposal(0.0, sqrt(tune_vars));
        std::uniform_real_distribution<double> trunc_runifs(cdf(rhos_proposal, -(rhos-min_range)), cdf(rhos_proposal, max_range - rhos));

        rhos_star = rhos + quantile(rhos_proposal, trunc_runifs(generator));
        maternmat = cMatern(d, rhos_star, nus);

        PLDs_precision_star = cinv(maternmat);
        PLDs_ldet_star = get_ldet(maternmat);

        for(auto rr = 0; rr < reps; rr++) {
            temp2 = (yminusXb.col(rr).t()) * PLDs_precision_star * yminusXb.col(rr);
            temp1(rr) = temp2(0);
        }
        SS_star = tau * arma::accu(temp1);
        temp1.set_size(reps);

        R = 0.5 * (PLDs_ldet_star - PLDs_ldet) - 0.5 * (SS_star - SS);
        if(runif(generator) < exp(R)) {
            rhos = rhos_star;
            PLDs_precision = PLDs_precision_star;
            PLDs_ldet = PLDs_ldet_star;
            acpt_rhos(mod(i, win_len)) = 1;
        } else {
            acpt_rhos(mod(i, win_len)) = 0;
        }

        //-------------------------------------------------------------
        // Update covariance parameters for the error process
        //-------------------------------------------------------------

        // Propose the proportion of error which comes from the correlated
        // error process, vs. iid noise
        for(auto rr = 0; rr < reps; rr++) {
            temp2 = (yminusXb.col(rr).t()) * PLDe_precision * yminusXb.col(rr);
            temp1(rr) = temp2(0);
        }
        SS = tau * arma::accu(temp1);
        lr = log(r / (1-r));
        lr_star = lr + quantile(lr_proposal, runif(generator));
        r_star = exp(lr_star) / (1 + exp(lr_star));

        Q = getQ(maternmat, d, r_star);
        PLDe_precision_star = cinv(Q);
        PLDe_ldet_star = get_ldet(Q);

        for(auto rr = 0; rr < reps; rr++) {
            temp2 = (yminusXb.col(rr).t()) * PLDe_precision_star * yminusXb.col(rr);
            temp1(rr) = temp2(0);
        }
        SS_star = arma::accu(temp1);
        temp1.set_size(reps);

        R = 0.5 * (PLDe_ldet_star - PLDe_ldet) - 0.5 * (SS_star - SS) +
            log(pdf(lr_proposal, lr_star)) - log(pdf(lr_proposal, lr));
        if(runif(generator) < exp(R)) {
            r = r_star;
            PLDe_precision = PLDe_precision_star;
            PLDe_ldet = PLDe_ldet_star;
            SS = SS_star;
        }

        normal rhoe_proposal(0.0, sqrt(tune_vare));
        std::uniform_real_distribution<double> trunc_runife(cdf(rhoe_proposal, -(rhoe-min_range)), cdf(rhoe_proposal, max_range - rhoe));
        rhoe_star = rhoe + quantile(rhoe_proposal, trunc_runife(generator));
        maternmat = cMatern(d, rhoe_star, nue);
        Q = getQ(maternmat, d, r);

        PLDe_precision_star = cinv(maternmat);
        PLDe_ldet_star = get_ldet(maternmat);

        for(auto rr = 0; rr < reps; rr++) {
            temp2 = (yminusXb.col(rr).t()) * PLDe_precision_star * yminusXb.col(rr);
            temp1(rr) = temp2(0);
        }
        SS_star = arma::accu(temp1);
        temp1.set_size(reps);

        R = 0.5 * (PLDe_ldet_star - PLDe_ldet) - 0.5 * (SS_star - SS);
        if(runif(generator) < exp(R)) {
            rhoe = rhoe_star;
            PLDe_precision = PLDe_precision_star;
            PLDe_ldet = PLDe_ldet_star;
            acpt_rhoe(mod(i, win_len)) = 1;
        } else {
            acpt_rhoe(mod(i, win_len)) = 0;
        }

        //-----------------------------------------------------------------------
        // Update the tuning variance
        //-----------------------------------------------------------------------
        if(i >= 100) {
            gamma1 = c0 / pow(i + tune_k, c1);
            acpt_rt_rhos = mean(acpt_rhos);
            acpt_rt_rhoe = mean(acpt_rhoe);
            tune_vars = cupdate_var(tune_vars, acpt_rt_rhos, gamma1, .3);
            tune_vare = cupdate_var(tune_vare, acpt_rt_rhoe, gamma1, .3);
            acpt_chain.row(i-100).col(0) = acpt_rt_rhos;
            acpt_chain.row(i-100).col(1) = acpt_rt_rhoe;
            tune_var_chain.row(i-100).col(0) = tune_vars;
            tune_var_chain.row(i-100).col(1) = tune_vare;
        }

        //-------------------------------------------------------------
        // Update chain parameter log
        //-------------------------------------------------------------
        beta_chain.row(i-1) = beta.t();

        param_chain.row(i-1).col(0) = 1 / sqrt(tau);
        param_chain.row(i-1).col(1) = r;
        param_chain.row(i-1).col(2) = rhos;
        param_chain.row(i-1).col(3) = rhoe;

        //-------------------------------------------------------------
        // Make predictions
        //-------------------------------------------------------------
        matern_cov = cMatern(d, rhos, nus);
        ypred.row(i-1) = (X * beta + cmake_pred(yminusXb, matern_cov, tau)).t();
    }
    return (Rcpp::List::create(Rcpp::Named("beta") = beta_chain,
                               Rcpp::Named("covar_params") = param_chain,
                               Rcpp::Named("pred") = ypred,
                               Rcpp::Named("acpt_chain") = acpt_chain,
                               Rcpp::Named("tune_chain") = tune_var_chain));
}
