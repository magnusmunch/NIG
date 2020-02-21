// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// f_optim
double f_optim(arma::vec alpha, arma::vec lambda, double nu, double zeta, arma::mat Cmat, arma::mat Z, int n, int p, int D, int G, int H, arma::mat y, arma::mat x, arma::vec yty);
RcppExport SEXP _cambridge_f_optim(SEXP alphaSEXP, SEXP lambdaSEXP, SEXP nuSEXP, SEXP zetaSEXP, SEXP CmatSEXP, SEXP ZSEXP, SEXP nSEXP, SEXP pSEXP, SEXP DSEXP, SEXP GSEXP, SEXP HSEXP, SEXP ySEXP, SEXP xSEXP, SEXP ytySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< double >::type zeta(zetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Cmat(CmatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    Rcpp::traits::input_parameter< int >::type H(HSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type yty(ytySEXP);
    rcpp_result_gen = Rcpp::wrap(f_optim(alpha, lambda, nu, zeta, Cmat, Z, n, p, D, G, H, y, x, yty));
    return rcpp_result_gen;
END_RCPP
}
// f_optim_mult
double f_optim_mult(arma::vec alpha, arma::vec lambda, double nu, double zeta, arma::mat Cmat, arma::mat Z, int n, arma::vec p, int D, int G, int H, arma::mat y, List x, arma::vec yty);
RcppExport SEXP _cambridge_f_optim_mult(SEXP alphaSEXP, SEXP lambdaSEXP, SEXP nuSEXP, SEXP zetaSEXP, SEXP CmatSEXP, SEXP ZSEXP, SEXP nSEXP, SEXP pSEXP, SEXP DSEXP, SEXP GSEXP, SEXP HSEXP, SEXP ySEXP, SEXP xSEXP, SEXP ytySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< double >::type zeta(zetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Cmat(CmatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    Rcpp::traits::input_parameter< int >::type H(HSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< List >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type yty(ytySEXP);
    rcpp_result_gen = Rcpp::wrap(f_optim_mult(alpha, lambda, nu, zeta, Cmat, Z, n, p, D, G, H, y, x, yty));
    return rcpp_result_gen;
END_RCPP
}
// Sigma_unp
arma::mat Sigma_unp(double aold, arma::vec bold, arma::mat xu, arma::mat xr, int u, int r);
RcppExport SEXP _cambridge_Sigma_unp(SEXP aoldSEXP, SEXP boldSEXP, SEXP xuSEXP, SEXP xrSEXP, SEXP uSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type aold(aoldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type bold(boldSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xu(xuSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xr(xrSEXP);
    Rcpp::traits::input_parameter< int >::type u(uSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(Sigma_unp(aold, bold, xu, xr, u, r));
    return rcpp_result_gen;
END_RCPP
}
// Sigma
arma::mat Sigma(double aold, arma::vec bold, arma::mat x);
RcppExport SEXP _cambridge_Sigma(SEXP aoldSEXP, SEXP boldSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type aold(aoldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type bold(boldSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(Sigma(aold, bold, x));
    return rcpp_result_gen;
END_RCPP
}
// aux_var_unp
List aux_var_unp(double aold, arma::vec bold, arma::vec y, arma::mat xu, arma::mat xr, int u, int r);
RcppExport SEXP _cambridge_aux_var_unp(SEXP aoldSEXP, SEXP boldSEXP, SEXP ySEXP, SEXP xuSEXP, SEXP xrSEXP, SEXP uSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type aold(aoldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type bold(boldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xu(xuSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xr(xrSEXP);
    Rcpp::traits::input_parameter< int >::type u(uSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(aux_var_unp(aold, bold, y, xu, xr, u, r));
    return rcpp_result_gen;
END_RCPP
}
// aux_var
List aux_var(double aold, arma::vec bold, arma::vec y, arma::mat x, arma::rowvec ytx);
RcppExport SEXP _cambridge_aux_var(SEXP aoldSEXP, SEXP boldSEXP, SEXP ySEXP, SEXP xSEXP, SEXP ytxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type aold(aoldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type bold(boldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type ytx(ytxSEXP);
    rcpp_result_gen = Rcpp::wrap(aux_var(aold, bold, y, x, ytx));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_stan_fit4nig_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4nig_full_mod();

static const R_CallMethodDef CallEntries[] = {
    {"_cambridge_f_optim", (DL_FUNC) &_cambridge_f_optim, 14},
    {"_cambridge_f_optim_mult", (DL_FUNC) &_cambridge_f_optim_mult, 14},
    {"_cambridge_Sigma_unp", (DL_FUNC) &_cambridge_Sigma_unp, 6},
    {"_cambridge_Sigma", (DL_FUNC) &_cambridge_Sigma, 3},
    {"_cambridge_aux_var_unp", (DL_FUNC) &_cambridge_aux_var_unp, 7},
    {"_cambridge_aux_var", (DL_FUNC) &_cambridge_aux_var, 5},
    {"_rcpp_module_boot_stan_fit4nig_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4nig_mod, 0},
    {"_rcpp_module_boot_stan_fit4nig_full_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4nig_full_mod, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_cambridge(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
