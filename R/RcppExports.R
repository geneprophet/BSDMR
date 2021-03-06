# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @name model_log_likelihood
#' 
#' @title model_log_likelihood
#' @param w a vector
#' @param X a matrix
#' @param H a matrix
#' @param lambda a number
#' @param is_nll boolean
#' @export
bpr_log_likelihood <- function(w, X, H, lambda, is_nll) {
    .Call(`_BSDMR_bpr_log_likelihood`, w, X, H, lambda, is_nll)
}

#' @rdname model_log_likelihood
#'
#' @export
bpr_gradient <- function(w, X, H, lambda, is_nll) {
    .Call(`_BSDMR_bpr_gradient`, w, X, H, lambda, is_nll)
}

#' @rdname model_log_likelihood
#'
#' @export
betareg_log_likelihood <- function(w, X, H, lambda, is_nll) {
    .Call(`_BSDMR_betareg_log_likelihood`, w, X, H, lambda, is_nll)
}

#' @rdname model_log_likelihood
#'
#' @export
betareg_gradient <- function(w, X, H, lambda, is_nll) {
    .Call(`_BSDMR_betareg_gradient`, w, X, H, lambda, is_nll)
}

#' @rdname model_log_likelihood
#' @param X_list  a list
#' @param H_list  a list
#' @param r_nk a vector
#' @export
sum_weighted_bpr_lik <- function(w, X_list, H_list, r_nk, lambda, is_nll) {
    .Call(`_BSDMR_sum_weighted_bpr_lik`, w, X_list, H_list, r_nk, lambda, is_nll)
}

#' @rdname model_log_likelihood
#'
#' @export
sum_weighted_bpr_grad <- function(w, X_list, H_list, r_nk, lambda, is_nll) {
    .Call(`_BSDMR_sum_weighted_bpr_grad`, w, X_list, H_list, r_nk, lambda, is_nll)
}

#' @rdname model_log_likelihood
#'
#' @export
sum_weighted_betareg_lik <- function(w, X_list, H_list, r_nk, lambda, is_nll) {
    .Call(`_BSDMR_sum_weighted_betareg_lik`, w, X_list, H_list, r_nk, lambda, is_nll)
}

#' @rdname model_log_likelihood
#'
#' @export
sum_weighted_betareg_grad <- function(w, X_list, H_list, r_nk, lambda, is_nll) {
    .Call(`_BSDMR_sum_weighted_betareg_grad`, w, X_list, H_list, r_nk, lambda, is_nll)
}

rcpp_hello_world <- function() {
    .Call(`_BSDMR_rcpp_hello_world`)
}

#' @name timesTwo
#' 
#' @title timesTwo
#' 
#' @param x A single integer.
#' 
#' @export
timesTwo <- function(x) {
    .Call(`_BSDMR_timesTwo`, x)
}

