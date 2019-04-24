# estimating equation for eta for root-finding (tested)
.froot.inv.gamma <- function(eta, lambda, esum, size) {
  log(lambda) - log(2) - digamma(eta/2) - esum/size
}

# inverse Gamma EB update (tested)
.eb.update.inv.gamma <- function(esum, bsum, b, delta, old.eta, elbo.const, 
                                 nclass, sclass, p, D, epsilon, maxit) {
  
  # obtaining initial values and bounds for root-finding
  old.eta <- eta.star <- 1/(log(bsum) + esum/sclass - log(sclass))
  old.lambda <- old.eta*sclass/bsum
  
  # iterate between root-finding alpha and computing lambda
  conv <- FALSE
  iter <- 0
  while(!(conv | iter >= maxit)) {
    
    # update iteration number and EB parameters
    iter <- iter + 1
    eta <- sapply(c(1:nclass), function(c) {
      uniroot(.froot.inv.gamma, c(eta.star[c], 2*eta.star[c]), 
              old.lambda[c], esum[c], sclass[c])$root})
    lambda <- eta*sclass/bsum
    conv <- all(abs(c(eta - old.eta, lambda - old.lambda)) < epsilon)
    old.lambda <- lambda
    old.eta <- eta
  }
  
  # elbo calculation involves constant part plus updated part
  elbo <- elbo.const - 0.5*(p + rep(eta, sclass))*log(delta) - 
    rep(lgamma(0.5*eta), sclass) + 
    0.5*rep(eta*digamma(0.5*(p + old.eta)), sclass) + 
    0.5*rep(eta*log(lambda), sclass) - 0.5*rep(lambda, sclass)*b
  
  out <- list(alpha=rep(0, nclass), eta=rep(eta, times=sclass), 
              theta=rep(Inf, D), lambda=rep(lambda, times=sclass), 
              elbo=elbo, conv=conv)
  return(out)
  
}

# one multiple lambda inverse Gaussian EB update (tested)
.eb.update.mult.lambda <- function(old.lambda, old.alpha, e, b, C, epsilon) {
  
  # create new IRLS updates
  alpha <- rowSums(solve(t(C*e*old.lambda) %*% C) %*% t(C*old.lambda))
  lambda <- 1/(b + e*as.numeric(C %*% alpha)^2 - 2*as.numeric(C %*% alpha))
  
  # check convergence of IRLS
  conv <- all(abs(c(alpha - old.alpha, lambda - old.lambda)) < epsilon)
  
  out <- list(alpha=unname(alpha), lambda=unname(lambda), conv=unname(conv))
  return(out)
}

# inverse Gaussian EB update (tested)
.eb.update.inv.gauss <- function(e, b, elbo.const, C, D, mult.lambda=FALSE, 
                                 epsilon, maxit) {
  
  # eb parameters
  alpha <- rowSums(solve(t(C*e) %*% C) %*% t(C))
  lambda <- D/(sum(b) - sum(t(C)*alpha))
  
  # IRLS
  if(mult.lambda) {
    new.eb <-  list(alpha=alpha, lambda=rep(lambda, D), conv=FALSE)
    iter <- 0
    while(!(new.eb$conv | iter >= maxit)) {
      
      # update iteration number and EB parameters
      iter <- iter + 1
      new.eb <- .eb.update.mult.lambda(new.eb$lambda, new.eb$alpha, e, b, C, 
                                       epsilon)
    }
    alpha <- new.eb$alpha
    lambda <- new.eb$lambda
    conv <- new.eb$conv
  } else {
    conv <- NA
  }
  
  # compute theta
  theta <- 1/as.numeric(C %*% alpha)
  
  # elbo calculation involves constant part plus updated part
  elbo <- elbo.const + lambda/theta + 0.5*log(lambda) - 0.5*lambda*e/theta^2 -
    0.5*lambda*b
  
  out <- list(alpha=alpha, eta=rep(1, D), theta=theta, lambda=lambda, elbo=elbo, 
              conv=conv)
  return(out)
}