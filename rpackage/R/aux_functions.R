# ELBO (not tested)
.single.elbo <- function(p, n, zeta, yty, aux, g, b, delta, eta, lambdaf, 
                         lambdad, Zalphad, Calphaf) {
  
  chi <- switch(is.null(Zalphad) + 1, 1/Zalphad, NULL)
  phi <- switch(is.null(Calphaf) + 1, 1/Calphaf, NULL)
  elbo <- -0.5*(n + p + 1)*log(pi) + 0.5*(p + 1 - n)*log(2) + 0.5*(n + 1) + 
    p + lgamma(0.5*(n + p + 1)) - 0.5*(n + p + 1)*log(zeta) + 
    0.5*aux$ldetSigma - 0.25*(n + p + 1)/zeta*
    (yty - 2*aux$ytXmu + aux$trXtXSigma + aux$mutXtXmu + 
       g*sum(b*aux$dSigma) + g*sum(aux$mu^2*b)) + 
    switch(is.null(Zalphad) + 1, 0.25*(p + 3)/4*log(lambdad) + 
             0.5*(eta - lambdad)*g + lambdad/chi - 0.5*(p + 1)*log(chi) - 
             0.25*(p + 1)*log(eta) + 
             gsl::bessel_lnKnu(0.5*(p + 1), sqrt(lambdad*eta)/chi), 
           - 0.5*log(2) + 0.5*log(pi)) +
    switch(is.null(Calphaf) + 1, p*log(lambdaf) + sum(0.5*(delta - lambdaf)*b) + 
             sum(lambdaf/phi) - sum(log(phi)) - 0.5*sum(log(delta)) +
             sum(gsl::bessel_lnKnu(1, sqrt(lambdaf*delta)/phi)), 
           -p/2*log(2) + p/2*log(p))
  return(elbo)
  
}

# computes ratio besselK(x, nu - 1)/besselK(x, nu) (tested)
ratio_besselK <- function(x, nu) {
  res <- nu - floor(nu)
  val <- besselK(x, res - 1, expon.scaled=TRUE)/
    besselK(x, res, expon.scaled=TRUE)
  for(i in 0:(floor(nu) - 1)) {
    val <- 1/(val + 2*(res + i)/x)
  }
  return(val)
}
