// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// ratio_besselK
NumericVector cpp_ratio_besselK(arma::vec x, arma::vec nu);
RcppExport SEXP _cambridge_cpp_ratio_besselK(SEXP xSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_ratio_besselK(x, nu));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_cambridge_cpp_ratio_besselK", (DL_FUNC) &_cambridge_cpp_ratio_besselK, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_cambridge(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
