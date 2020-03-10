// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <R.h>
#include <Rcpp.h>

// To use functions in Rcpp library without the need of specifying "Rcpp::"
using namespace Rcpp;
using namespace arma;

// function to optimise in EBridge with one x (not tested)
// [[Rcpp::export(".f.optim.mat")]]
double f_optim_mat(arma::vec alpha, arma::vec lambda, double nu, double zeta,
                   arma::mat Cmat, arma::mat Z, arma::vec n, int p, int D,
                   List idsel, int G, int H, List y, arma::mat x,
                   arma::vec yty) {

  arma::vec alphaf = alpha.head(G);
  arma::vec alphad = alpha.tail(H);

  arma::vec tau = exp(0.5*Z*alphad);
  arma::vec gamma = exp(0.5*Cmat*alphaf);
  double out=0.0;
  for(int d=0; d<D; d++) {
    arma::rowvec mat1(n(d));
    arma::mat mat2(n(d), n(d)); 
    arma::uvec cidsel = idsel[d];
    arma::mat cx = x.rows(cidsel - 1);
    arma::vec cy = y[d];

    mat2 = cx*trans(cx.each_row() % square(lambda(d)*tau(d)*
      gamma.subvec(d*p, (d + 1)*p - 1).t()));
    mat1 = cy.t()*mat2;
    mat2.diag() += 1;
    out += -0.5*real(log_det(mat2)) -
      (n(d)/2 + nu)*log(zeta + arma::conv_to<double>::from(
          0.5*yty(d) - 0.5*mat1*y[d] + 0.5*mat1*mat2.i()*mat1.t()));
  }
  return out;
}

// function to optimise in EBridge with multiple x (not tested)
// [[Rcpp::export(".f.optim.list")]]
double f_optim_list(arma::vec alpha, arma::vec lambda, double nu, double zeta,
                    arma::mat Cmat, arma::mat Z, arma::vec n, arma::vec p, 
                    int D, int G, int H, List y, List x, arma::vec yty) {

  arma::vec alphaf = alpha.head(G);
  arma::vec alphad = alpha.tail(H);

  arma::vec tau = exp(0.5*Z*alphad);
  arma::vec gamma = exp(0.5*Cmat*alphaf);
  arma::vec cp(D + 1);
  cp.subvec(1, D) = cumsum(p);
  cp(0) = 0;
  arma::vec out(D);
  for(int d=0; d<D; d++) {
    arma::rowvec mat1(n(d)); 
    arma::mat mat2(n(d), n(d));
    arma::mat cx = x[d];
    arma::vec cy = y[d];

    mat2 = cx*trans(cx.each_row() % square(lambda(d)*tau(d)*
      gamma.subvec(cp(d), cp(d + 1) - 1).t()));
    mat1 = cy.t()*mat2;
    mat2.diag() += 1;
    out(d) = -0.5*real(log_det(mat2)) -
      (n(d)/2 + nu)*log(zeta + arma::conv_to<double>::from(
          0.5*yty(d) - 0.5*mat1*cy + 0.5*mat1*mat2.i()*mat1.t()));

  }
  double val = sum(out);
  return val;
}

// calculate full covariance with unpenalized covariates (tested)
// [[Rcpp::export(".Sigma.unp")]]
arma::mat Sigma_unp(double aold, arma::vec bold, arma::mat xu, arma::mat xr,
                    int u, int r) {

  int p = u + r;
  arma::mat xrt = xr.t();
  arma::mat xut = xu.t();
  arma::vec hinv = 1/bold;
  arma::mat HinvXrt = xrt.each_col() % hinv;
  arma::mat XrHinvXrt, Om;
  XrHinvXrt = Om = xr*HinvXrt;
  Om.diag() += 1;
  Om = Om.i();
  arma::mat HinvXrtOm = HinvXrt*Om;
  arma::mat M11 = xut*Om*xu;
  M11 = M11.i();
  arma::mat M12 = - M11*xut*HinvXrt.t() + M11*xut*XrHinvXrt.t()*HinvXrtOm.t();
  arma::mat M22 = -HinvXrtOm*HinvXrt.t() - HinvXrt*xu*M12 +
    HinvXrtOm*XrHinvXrt*xu*M12;
  M22.diag() += hinv;
  arma::mat Sigma(p, p);
  Sigma.submat(0, 0, u - 1, u - 1) = M11;
  Sigma.submat(u, u, p - 1, p - 1) = M22;
  Sigma.submat(0, u, u - 1, p - 1) = M12;
  Sigma.submat(u, 0, p - 1, u - 1) = M12.t();
  Sigma = Sigma/aold;

  return Sigma;
}

// calculate full covariance (tested)
// [[Rcpp::export(".Sigma")]]
arma::mat Sigma(double aold, arma::vec bold, arma::mat x) {
  arma::vec hinv = 1/bold;
  arma::mat HinvXt = x.t();
  HinvXt = HinvXt.each_col() % hinv;
  arma::mat Sigma = x*HinvXt;
  Sigma.diag() += 1;
  Sigma = -HinvXt*Sigma.i()*HinvXt.t();
  Sigma.diag() += hinv;
  Sigma = Sigma/aold;
  return Sigma;
}

// calculates the auxiliary variables with unpenalized variables (tested)
// [[Rcpp::export(".aux.var.unp")]]
List aux_var_unp(double aold, arma::vec bold, arma::vec y, arma::mat xu, 
                 arma::mat xr, int u, int r) {

  int n = xu.n_rows;
  int p = u + r;
  arma::mat xrt = xr.t();
  arma::mat xut = xu.t();
  arma::vec hinv = 1/bold;
  arma::mat HinvXrt = xrt.each_col() % hinv;
  arma::mat XrHinvXrt, Om; 
  XrHinvXrt = Om = xr*HinvXrt;
  Om.diag() += 1;
  Om = Om.i();
  arma::mat HinvXrtOm = HinvXrt*Om;
  arma::mat M11 = xut*Om*xu;
  M11 = M11.i();
  arma::mat M12 = - M11*xut*HinvXrt.t() + M11*xut*XrHinvXrt.t()*HinvXrtOm.t();
  
  // auxiliaries
  arma::vec dSigma(p);
  dSigma.head(u) = M11.diag()/aold;
  dSigma.tail(r) = (hinv - sum(HinvXrt%HinvXrtOm, 1) -
    sum((HinvXrt*xu) % M12.t(), 1) +
    sum(HinvXrtOm % (M12.t()*xut*XrHinvXrt), 1))/aold;
  arma::mat mu(p, n);
  mu.head_rows(u) = M11*xut + M12*xrt;
  mu.tail_rows(r) = HinvXrt - HinvXrtOm*XrHinvXrt -
    (HinvXrt - HinvXrtOm*XrHinvXrt)*xu*M12*xrt;
  double trXtXSigma = (accu(xu % (xu*M11)) + 2*accu(xu % (xr*M12.t())) +
                       accu(xrt % mu.tail_rows(r)))/aold;
  mu.tail_rows(r) =  mu.tail_rows(r) + trans(xu*M12);
  mu = mu*y;
  arma::vec ytXmu(n);
  ytXmu = xu*mu.head_rows(u) + xr*mu.tail_rows(r);
  double mutXtXmu = arma::as_scalar(ytXmu.t()*ytXmu);
  ytXmu = ytXmu.t()*y;
  
  // log determinant of Sigma for ELBO calculations
  arma::mat amat = -xu*inv_sympd(xut*xu)*xut;
  amat.diag() += 1;
  amat = amat*XrHinvXrt;
  amat.diag() += 1;
  double val;
  double sign;
  log_det(val, sign, amat);
  double ldet = val;
  log_det(val, sign, xut*xu);
  double ldetSigma = -p*log(aold) - sum(log(bold)) - val - ldet;

  // changing type for easier R integration
  double dytXmu = arma::as_scalar(ytXmu);
  NumericVector vmu = wrap(mu);
  vmu.attr("dim") = R_NilValue;
  NumericVector vdSigma = wrap(dSigma);
  vdSigma.attr("dim") = R_NilValue;

  return List::create(Named("mu") = vmu,
                      Named("dSigma") = vdSigma,
                      Named("ytXmu") = dytXmu,
                      Named("trXtXSigma") = trXtXSigma,
                      Named("mutXtXmu") = mutXtXmu,
                      Named("ldetSigma") = ldetSigma);
}

// calculates the auxiliary variables (tested)
// [[Rcpp::export(".aux.var")]]
List aux_var(double aold, arma::vec bold, arma::vec y, arma::mat x, 
             arma::rowvec ytx) {
  
  arma::mat xt = x.t();
  arma::vec hinv = 1/bold;
  arma::mat HinvXt = xt.each_col() % hinv;
  arma::mat XHinvXt = x*HinvXt;
    
  // shared auxiliaries
  arma::mat dSigma = XHinvXt;
  dSigma.diag() += 1;
  double val;
  double sign;
  log_det(val, sign, dSigma);
  dSigma = HinvXt*dSigma.i();
  arma::mat mu, trXtXSigma, mutXt;
  mu = trXtXSigma = mutXt = HinvXt - dSigma*XHinvXt;
  dSigma = (hinv - sum(dSigma%HinvXt, 1))/aold;
  mu = mu*y;
  mutXt = ytx*mutXt;
  trXtXSigma = accu(trXtXSigma%xt)/aold;

  // separate auxiliaries
  double ytXmu = arma::as_scalar(ytx*mu);
  double mutXtXmu = arma::as_scalar(sum(square(mutXt), 1));
  
  // changing type for easier R integration
  double dtrXtXSigma = arma::as_scalar(trXtXSigma);
  NumericVector vmu = wrap(mu);
  vmu.attr("dim") = R_NilValue;
  NumericVector vdSigma = wrap(dSigma);
  vdSigma.attr("dim") = R_NilValue;
  
  // log determinant of Sigma for ELBO calculations
  int p = x.n_cols;
  double ldetSigma = -p*log(aold) - sum(log(bold)) - val;

  return List::create(Named("mu") = vmu,
                      Named("dSigma") = vdSigma,
                      Named("ytXmu") = ytXmu,
                      Named("trXtXSigma") = dtrXtXSigma,
                      Named("mutXtXmu") = mutXtXmu,
                      Named("ldetSigma") = ldetSigma);
}
