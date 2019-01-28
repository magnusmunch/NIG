//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <R.h>
#include <Rcpp.h>

// To use functions in Rcpp library without the need of specifying "Rcpp::"
using namespace Rcpp;

// computes ratio besselK(x, nu - 1)/besselK(x, nu) (tested)
// [[Rcpp::export]]
NumericVector ratio_besselK(arma::vec x, arma::vec nu) {

  // initialize variables
  int n = x.size();
  NumericVector value(n);
  double num;

  // check whether p/2 is an even and set the number of iteration accordingly
  

  // iterate over the values in x and estimate the ratio iteratively
  for(int i=0; i<n; i++) {
    int fnu = floor(nu[i]);
    double rnu = nu[i] - fnu;
    if(fnu==0) {
      value[i] = R::bessel_k(x[i], nu[i] - 1, exp(1))/
        R::bessel_k(x[i], nu[i], exp(1));
    } else {
      num = R::bessel_k(x[i], rnu - 1, exp(1))/R::bessel_k(x[i], rnu, exp(1));
      for(int j=0; j<fnu; j++) {
        num = 1/(num + 2*(nu[i] - fnu + j)/x[i]);
      }
      value[i] = num;
    }
  }

  // return a vector with the estimated values
  return value;
}






