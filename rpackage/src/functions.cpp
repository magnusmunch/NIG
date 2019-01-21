//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <R.h>
#include <Rcpp.h>

// To use functions in Rcpp library without the need of specifying "Rcpp::"
using namespace Rcpp;

// computes ratio besselK(x, nu - 1)/besselK(x, nu) (not tested)
// [[Rcpp::export]]
NumericVector ratio_besselK_cpp(arma::vec x, double nu) {

  // initialize variables
  int n = x.size();
  NumericVector value(n);
  double num;

  // check whether p/2 is an even and set the number of iteration accordingly
  int fnu = floor(nu);
  double rnu = nu - fnu;

  // iterate over the values in x and estimate the ratio iteratively
  for(int i=0; i<n; i++) {
    if(fnu==0) {
      value[i] = R::bessel_k(x[i], nu - 1, exp(1))/
        R::bessel_k(x[i], nu, exp(1));
    } else {
      num = R::bessel_k(x[i], rnu - 1, exp(1))/R::bessel_k(x[i], rnu, exp(1));
      for(int j=0; j<fnu; j++) {
        num = 1/(num + 2*(nu - fnu + j)/x[i]);
      }
      value[i] = num;
    }
  }

  // return a vector with the estimated values
  return value;
}

// one VB update with sigma^2 stochastic (not tested)
// [[Rcpp::export]]
List vb_update_cpp(arma::vec zetaold, arma::vec aold, double lambda, arma::vec theta, arma::mat x, 
                   arma::mat xxt, arma::mat xtx, arma::mat y) {
  
  // Initialize all constants and variables
  // int p = x.n_cols;
  // int D = zetaold.n_rows;
  // int n = x.n_rows;
  
  List Sigma;
  arma::mat mu;
  NumericVector delta;
  NumericVector zeta;
  NumericVector a;
  NumericVector v;
  
  
  
  return List::create(Named("Sigma") = Sigma, 
                      Named("mu")=mu, 
                      Named("delta")=delta, 
                      Named("zeta")=zeta, 
                      Named("a")=a, 
                      Named("v")=v);
}







