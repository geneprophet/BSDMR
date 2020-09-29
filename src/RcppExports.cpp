// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// bpr_log_likelihood
double bpr_log_likelihood(const arma::vec& w, const arma::mat& X, const arma::mat& H, const double lambda, const bool is_nll);
RcppExport SEXP _BSDMR_bpr_log_likelihood(SEXP wSEXP, SEXP XSEXP, SEXP HSEXP, SEXP lambdaSEXP, SEXP is_nllSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type H(HSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const bool >::type is_nll(is_nllSEXP);
    rcpp_result_gen = Rcpp::wrap(bpr_log_likelihood(w, X, H, lambda, is_nll));
    return rcpp_result_gen;
END_RCPP
}
// bpr_gradient
Rcpp::NumericVector bpr_gradient(const arma::vec& w, const arma::mat& X, const arma::mat& H, const double lambda, const bool is_nll);
RcppExport SEXP _BSDMR_bpr_gradient(SEXP wSEXP, SEXP XSEXP, SEXP HSEXP, SEXP lambdaSEXP, SEXP is_nllSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type H(HSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const bool >::type is_nll(is_nllSEXP);
    rcpp_result_gen = Rcpp::wrap(bpr_gradient(w, X, H, lambda, is_nll));
    return rcpp_result_gen;
END_RCPP
}
// betareg_log_likelihood
double betareg_log_likelihood(const arma::vec& w, arma::mat& X, const arma::mat& H, const double lambda, const bool is_nll);
RcppExport SEXP _BSDMR_betareg_log_likelihood(SEXP wSEXP, SEXP XSEXP, SEXP HSEXP, SEXP lambdaSEXP, SEXP is_nllSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type H(HSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const bool >::type is_nll(is_nllSEXP);
    rcpp_result_gen = Rcpp::wrap(betareg_log_likelihood(w, X, H, lambda, is_nll));
    return rcpp_result_gen;
END_RCPP
}
// betareg_gradient
Rcpp::NumericVector betareg_gradient(const arma::vec& w, arma::mat& X, const arma::mat& H, const double lambda, const bool is_nll);
RcppExport SEXP _BSDMR_betareg_gradient(SEXP wSEXP, SEXP XSEXP, SEXP HSEXP, SEXP lambdaSEXP, SEXP is_nllSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type H(HSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const bool >::type is_nll(is_nllSEXP);
    rcpp_result_gen = Rcpp::wrap(betareg_gradient(w, X, H, lambda, is_nll));
    return rcpp_result_gen;
END_RCPP
}
// sum_weighted_bpr_lik
double sum_weighted_bpr_lik(const arma::vec& w, const Rcpp::List& X_list, const Rcpp::List& H_list, const arma::vec& r_nk, const double lambda, const bool is_nll);
RcppExport SEXP _BSDMR_sum_weighted_bpr_lik(SEXP wSEXP, SEXP X_listSEXP, SEXP H_listSEXP, SEXP r_nkSEXP, SEXP lambdaSEXP, SEXP is_nllSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type X_list(X_listSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type H_list(H_listSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type r_nk(r_nkSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const bool >::type is_nll(is_nllSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_weighted_bpr_lik(w, X_list, H_list, r_nk, lambda, is_nll));
    return rcpp_result_gen;
END_RCPP
}
// sum_weighted_bpr_grad
arma::rowvec sum_weighted_bpr_grad(const arma::vec& w, const Rcpp::List& X_list, const Rcpp::List& H_list, const arma::vec& r_nk, const double lambda, const bool is_nll);
RcppExport SEXP _BSDMR_sum_weighted_bpr_grad(SEXP wSEXP, SEXP X_listSEXP, SEXP H_listSEXP, SEXP r_nkSEXP, SEXP lambdaSEXP, SEXP is_nllSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type X_list(X_listSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type H_list(H_listSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type r_nk(r_nkSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const bool >::type is_nll(is_nllSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_weighted_bpr_grad(w, X_list, H_list, r_nk, lambda, is_nll));
    return rcpp_result_gen;
END_RCPP
}
// sum_weighted_betareg_lik
double sum_weighted_betareg_lik(const arma::vec& w, const Rcpp::List& X_list, const Rcpp::List& H_list, const arma::vec& r_nk, const double lambda, const bool is_nll);
RcppExport SEXP _BSDMR_sum_weighted_betareg_lik(SEXP wSEXP, SEXP X_listSEXP, SEXP H_listSEXP, SEXP r_nkSEXP, SEXP lambdaSEXP, SEXP is_nllSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type X_list(X_listSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type H_list(H_listSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type r_nk(r_nkSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const bool >::type is_nll(is_nllSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_weighted_betareg_lik(w, X_list, H_list, r_nk, lambda, is_nll));
    return rcpp_result_gen;
END_RCPP
}
// sum_weighted_betareg_grad
arma::rowvec sum_weighted_betareg_grad(const arma::vec& w, const Rcpp::List& X_list, const Rcpp::List& H_list, const arma::vec& r_nk, const double lambda, const bool is_nll);
RcppExport SEXP _BSDMR_sum_weighted_betareg_grad(SEXP wSEXP, SEXP X_listSEXP, SEXP H_listSEXP, SEXP r_nkSEXP, SEXP lambdaSEXP, SEXP is_nllSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type X_list(X_listSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type H_list(H_listSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type r_nk(r_nkSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const bool >::type is_nll(is_nllSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_weighted_betareg_grad(w, X_list, H_list, r_nk, lambda, is_nll));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _BSDMR_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// timesTwo
NumericVector timesTwo(NumericVector x);
RcppExport SEXP _BSDMR_timesTwo(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(timesTwo(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BSDMR_bpr_log_likelihood", (DL_FUNC) &_BSDMR_bpr_log_likelihood, 5},
    {"_BSDMR_bpr_gradient", (DL_FUNC) &_BSDMR_bpr_gradient, 5},
    {"_BSDMR_betareg_log_likelihood", (DL_FUNC) &_BSDMR_betareg_log_likelihood, 5},
    {"_BSDMR_betareg_gradient", (DL_FUNC) &_BSDMR_betareg_gradient, 5},
    {"_BSDMR_sum_weighted_bpr_lik", (DL_FUNC) &_BSDMR_sum_weighted_bpr_lik, 6},
    {"_BSDMR_sum_weighted_bpr_grad", (DL_FUNC) &_BSDMR_sum_weighted_bpr_grad, 6},
    {"_BSDMR_sum_weighted_betareg_lik", (DL_FUNC) &_BSDMR_sum_weighted_betareg_lik, 6},
    {"_BSDMR_sum_weighted_betareg_grad", (DL_FUNC) &_BSDMR_sum_weighted_betareg_grad, 6},
    {"_BSDMR_rcpp_hello_world", (DL_FUNC) &_BSDMR_rcpp_hello_world, 0},
    {"_BSDMR_timesTwo", (DL_FUNC) &_BSDMR_timesTwo, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_BSDMR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
