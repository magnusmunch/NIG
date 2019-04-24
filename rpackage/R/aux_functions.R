# calculates the auxiliary variables in the VB step (tested)
.aux.var <- function(cold, aold, sv, uty, n, p) {
  
  # auxiliary variables involving mu and Sigma
  trSigma <- (sum(1/(sv^2 + cold)) + max(p - n, 0)/cold)/aold
  trXtXSigma <- sum(sv^2/(sv^2 + cold))/aold
  mutmu <- sum(sv^2*uty^2/(sv^2 + cold)^2)
  mutXtXmu <- sum(sv^4*uty^2/(sv^2 + cold)^2)
  logdetSigma <- -p*log(aold) - sum(log(sv^2 + cold)) - 
    max(p - n, 0)*log(cold)
  ytXmu <- sum(sv^2*uty^2/(sv^2 + cold))
  
  out <- list(trSigma=trSigma, trXtXSigma=trXtXSigma, mutmu=mutmu,
              mutXtXmu=mutXtXmu, logdetSigma=logdetSigma,
              ytXmu=ytXmu)
  return(out)
}

# calculates the auxiliary variables in the VB step with intercept (not tested)
.aux.var.int <- function(cold, aold, sv, uty, sumu, sumy, n, p) {
  
  s <- 1/(n - sum(sumu^2*sv^2/(sv^2 + cold)))
  
  aux1 <- sum(sv^2*sumu^2/(sv^2 + cold))
  aux2 <- sum(sv^2*uty^2/(sv^2 + cold))
  aux3 <- sum(sv^2*sumu*uty/(sv^2 + cold))
  aux4 <- sum(sv^4*sumu^2/(sv^2 + cold)^2)
  aux5 <- sum(sv^4*uty^2/(sv^2 + cold)^2)
  aux6 <- sum(sv^4*sumu*uty/(sv^2 + cold)^2)
  
  # auxiliary variables involving mu and Sigma
  trSigma <- (sum(1/(sv^2 + cold)) + s*sum(sv^2*sumu^2/(sv^2 + cold)^2) + 
                max(p - n, 0)/cold)/aold 
  trXtXSigma <- (sum(sv^2/(sv^2 + cold)) + s*(n - 2*aux1 + aux4))/aold
  mutmu <- s^2*(sumy^2 - 2*sumy*aux3 + aux3^2)*
    sum(sumu^2*sv^2/(sv^2 + cold)^2) + 
    2*s*(aux3 - sumy)*sum(sumu*uty*sv^2/(sv^2 + cold)^2) + 
    sum(uty^2*sv^2/(sv^2 + cold)^2)
  mutXtXmu <- n*sumy^2*s^2 + sumy*s*aux3*(2 - 2*n*s) - 2*sumy^2*s^2*aux1 + 
    s*aux3^2*(n*s - 2) + 4*sumy*s^2*aux1*aux3 - 2*s^2*aux1*aux3^2 +
    sumy^2*s^2*aux4 - 2*sumy*s*aux6 - 2*sumy*s^2*aux3*aux4 + 2*s*aux3*aux6 + 
    s^2*aux3^2*aux4 + aux5
  logdetSigma <- log(s) - (p + 1)*log(aold) - sum(log(sv^2 + cold)) - 
    max(p - n, 0)*log(cold)
  ytXmu <- sumy^2*s - 2*sumy*s*aux3 + aux2 + s*aux3^2
  
  out <- list(trSigma=trSigma, trXtXSigma=trXtXSigma, mutmu=mutmu,
              mutXtXmu=mutXtXmu, logdetSigma=logdetSigma,
              ytXmu=ytXmu, s=s, aux3=aux3)
  return(out)
}

# calculates the ELBO in the inverse Gaussian model (tested)
.elbo.inv.gauss <- function(aux, zeta, delta, a, b, c, e, lambda, theta, df, p, 
                            yty) {
  
  # calculate the elbo part that is constant after next eb update
  elbo.const <- 0.5*aux$logdetSigma -
    0.5*a*(yty - 2*aux$ytXmu + aux$mutXtXmu + aux$trXtXSigma) +
    0.5*b*(delta - (a + 1 - b/c)*(aux$mutmu + aux$trSigma)) - 
    0.5*(df + 1)*log(zeta) + 0.25*(p + 1)*log(lambda) - 
    0.5*(p + 1)*log(theta) - 0.25*(p + 1)*log(delta) + 
    gsl::bessel_lnKnu(0.5*(p + 1), sqrt(lambda*delta)/theta) + 
    0.5*lambda*e/theta^2
  
  # total elbo
  elbo <- elbo.const + lambda/theta + 0.5*log(lambda) - 0.5*lambda*e/theta^2 -
    0.5*lambda*b
  
  out <- list(elbo=elbo, elbo.const=elbo.const)
  return(out)
}

# calculates the ELBO in the inverse Gamma model (tested)
.elbo.inv.gamma <- function(aux, zeta, delta, a, b, c, e, lambda, eta, df, p, 
                            yty) {
  
  # calculate the elbo part that is constant after next eb update
  elbo.const <- 0.5*aux$logdetSigma -
    0.5*a*(yty - 2*aux$ytXmu + aux$mutXtXmu + aux$trXtXSigma) +
    0.5*b*(delta - (a + 1 - b/c)*(aux$mutmu + aux$trSigma)) - 
    0.5*(df + 1)*log(zeta) - 0.25*(p + 1)*log(delta) + 0.5*eta*e + 
    0.5*eta*log(2) + lgamma(0.5*(p + eta))
  
  # total elbo
  elbo <- elbo.const - lgamma(0.5*eta) - 0.5*eta*e - 0.5*eta*log(2) + 
    0.5*eta*log(lambda) - 0.5*lambda*b
  
  out <- list(elbo=elbo, elbo.const=elbo.const)
  return(out)
  
}

# ridge marginal log likelihood of ridge model (not tested)
.ridge.mll <- function(par, sv, uty, n) {
  sigma.sq <- par[1]
  gamma.sq <- par[2]
  mll <- -0.5*n*log(sigma.sq) - 0.5*n*log(gamma.sq) - 
    sum(log(sv^2 + 1/gamma.sq)) - 0.5*sum(uty^2/(sv^2 + 1/gamma.sq))/
    (sigma.sq*gamma.sq)
  return(mll)
}

# ridge marginal log lik for independent beta and sigma^2 prior (not tested)
.ridge.mll.ind <- function(par, sv, uty, n) {
  sigma.sq <- par[1]
  gamma.sq <- par[2]
  mll <- -0.5*n*log(sigma.sq) - 0.5*n*log(gamma.sq) - 
    sum(log(sv^2/sigma.sq + 1/gamma.sq)) - 
    0.5*sum(uty^2/(sv^2*gamma.sq + sigma.sq))
  return(mll)
}

# initial values (not tested)
.init.param.inv.gauss <- function(C, sv, uty, D, n, p) {
  # maximise ridge mll to find gamma_d^2 and sigma_d^2 estimates
  init.mll.est <- t(sapply(c(1:D), function(d) {
    constrOptim(c(1, 1), .ridge.mll, NULL, ui=diag(2), ci=c(0, 0),
                sv=sv, uty=uty[, d], n=n, control=list(fnscale=-1))$par}))
  
  # consider them samples from prior and estimate one inverse Gaussian prior
  theta <- mean(init.mll.est[, 2])
  lambda <- D/sum(1/init.mll.est[, 2] - 1/theta)
  
  # use empirical mode of gamma_d^2 to derive one delta
  d <- density(init.mll.est[, 2])
  mode <- d$x[which.max(d$y)]
  delta <- as.numeric(mode^2*lambda/theta^2 + (p + 3)*mode)
  a <- sqrt(lambda/(theta^2*delta))*
    ratio_besselK(sqrt(lambda*delta)/theta, p) + (p + 1)/delta
  
  # consider sigma_d^2 fixed and estimate inverse Gamma prior
  zeta <- 0.5*D*(n + p + 1)/sum(1/init.mll.est[, 1])
  
  out <- list(theta=rep(theta, D), alpha=c(theta, rep(0, ncol(C) - 1)), 
              lambda=lambda, zeta=rep(zeta, D), a=rep(a, D))
  return(out)
  
}

# initial values with independent beta and sigma^2 prior (not tested)
.init.param.inv.gauss.ind <- function(C, sv, uty, D, n, p) {
  # maximise ridge mll to find gamma_d^2 and sigma_d^2 estimates
  init.mll.est <- t(sapply(c(1:D), function(d) {
    constrOptim(c(1, 1), .ridge.mll.ind, NULL, ui=diag(2), ci=c(0, 0),
                sv=sv, uty=uty[, d], n=n, control=list(fnscale=-1))$par}))
  
  # consider them samples from prior and estimate one inverse Gaussian prior
  theta <- mean(init.mll.est[, 2])
  lambda <- D/sum(1/init.mll.est[, 2] - 1/theta)
  
  # use empirical mode of gamma_d^2 to derive one delta
  d <- density(init.mll.est[, 2])
  mode <- d$x[which.max(d$y)]
  delta <- as.numeric(mode^2*lambda/theta^2 + (p + 3)*mode)
  a <- sqrt(lambda/(theta^2*delta))*
    ratio_besselK(sqrt(lambda*delta)/theta, p) + (p + 1)/delta
  
  # consider sigma_d^2 fixed and estimate inverse Gamma prior
  zeta <- 0.5*D*(n + 1)/sum(1/init.mll.est[, 1])
  
  out <- list(theta=rep(theta, D), alpha=c(theta, rep(0, ncol(C) - 1)), 
              lambda=lambda, zeta=rep(zeta, D), a=rep(a, D))
  return(out)
}