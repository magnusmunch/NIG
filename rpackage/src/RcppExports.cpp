// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// ratio_besselK_cpp
NumericVector ratio_besselK_cpp(arma::vec x, int p);
RcppExport SEXP _cambridge_ratio_besselK_cpp(SEXP xSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(ratio_besselK_cpp(x, p));
    return rcpp_result_gen;
END_RCPP
}
// vb_update_cpp
List vb_update_cpp(arma::vec zetaold, arma::vec aold, double lambda, arma::vec theta, arma::mat x, arma::mat xxt, arma::mat xtx, arma::mat y);
RcppExport SEXP _cambridge_vb_update_cpp(SEXP zetaoldSEXP, SEXP aoldSEXP, SEXP lambdaSEXP, SEXP thetaSEXP, SEXP xSEXP, SEXP xxtSEXP, SEXP xtxSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type zetaold(zetaoldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type aold(aoldSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xxt(xxtSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xtx(xtxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(vb_update_cpp(zetaold, aold, lambda, theta, x, xxt, xtx, y));
    return rcpp_result_gen;
END_RCPP
}
// VarOneIter
Rcpp::List VarOneIter(arma::colvec myy, arma::mat myX, arma::colvec myP, arma::colvec aRand, arma::colvec bRand, arma::colvec bRandStarInit, double dSigmaStarInit, double cSigma, double dSigma, arma::mat lincomb);
RcppExport SEXP _cambridge_VarOneIter(SEXP myySEXP, SEXP myXSEXP, SEXP myPSEXP, SEXP aRandSEXP, SEXP bRandSEXP, SEXP bRandStarInitSEXP, SEXP dSigmaStarInitSEXP, SEXP cSigmaSEXP, SEXP dSigmaSEXP, SEXP lincombSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type myy(myySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type myX(myXSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type myP(myPSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type aRand(aRandSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type bRand(bRandSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type bRandStarInit(bRandStarInitSEXP);
    Rcpp::traits::input_parameter< double >::type dSigmaStarInit(dSigmaStarInitSEXP);
    Rcpp::traits::input_parameter< double >::type cSigma(cSigmaSEXP);
    Rcpp::traits::input_parameter< double >::type dSigma(dSigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type lincomb(lincombSEXP);
    rcpp_result_gen = Rcpp::wrap(VarOneIter(myy, myX, myP, aRand, bRand, bRandStarInit, dSigmaStarInit, cSigma, dSigma, lincomb));
    return rcpp_result_gen;
END_RCPP
}
// fixedPointIterEB
Rcpp::NumericVector fixedPointIterEB(arma::colvec initab, arma::colvec myallaRandStar, arma::colvec myallbRandStar, int mymaxiter, double myeps);
RcppExport SEXP _cambridge_fixedPointIterEB(SEXP initabSEXP, SEXP myallaRandStarSEXP, SEXP myallbRandStarSEXP, SEXP mymaxiterSEXP, SEXP myepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type initab(initabSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type myallaRandStar(myallaRandStarSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type myallbRandStar(myallbRandStarSEXP);
    Rcpp::traits::input_parameter< int >::type mymaxiter(mymaxiterSEXP);
    Rcpp::traits::input_parameter< double >::type myeps(myepsSEXP);
    rcpp_result_gen = Rcpp::wrap(fixedPointIterEB(initab, myallaRandStar, myallbRandStar, mymaxiter, myeps));
    return rcpp_result_gen;
END_RCPP
}
// BSEMVarOneIter
Rcpp::List BSEMVarOneIter(Rcpp::List ylist, Rcpp::List Xlist, Rcpp::List Plist, Rcpp::List alist, Rcpp::List blist, Rcpp::List bstarlist, double cSigma, double dSigma, arma::colvec dstarvec, Rcpp::List lincomblist);
RcppExport SEXP _cambridge_BSEMVarOneIter(SEXP ylistSEXP, SEXP XlistSEXP, SEXP PlistSEXP, SEXP alistSEXP, SEXP blistSEXP, SEXP bstarlistSEXP, SEXP cSigmaSEXP, SEXP dSigmaSEXP, SEXP dstarvecSEXP, SEXP lincomblistSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type ylist(ylistSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Xlist(XlistSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Plist(PlistSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type alist(alistSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type blist(blistSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type bstarlist(bstarlistSEXP);
    Rcpp::traits::input_parameter< double >::type cSigma(cSigmaSEXP);
    Rcpp::traits::input_parameter< double >::type dSigma(dSigmaSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type dstarvec(dstarvecSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type lincomblist(lincomblistSEXP);
    rcpp_result_gen = Rcpp::wrap(BSEMVarOneIter(ylist, Xlist, Plist, alist, blist, bstarlist, cSigma, dSigma, dstarvec, lincomblist));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_cambridge_ratio_besselK_cpp", (DL_FUNC) &_cambridge_ratio_besselK_cpp, 2},
    {"_cambridge_vb_update_cpp", (DL_FUNC) &_cambridge_vb_update_cpp, 8},
    {"_cambridge_VarOneIter", (DL_FUNC) &_cambridge_VarOneIter, 10},
    {"_cambridge_fixedPointIterEB", (DL_FUNC) &_cambridge_fixedPointIterEB, 5},
    {"_cambridge_BSEMVarOneIter", (DL_FUNC) &_cambridge_BSEMVarOneIter, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_cambridge(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}