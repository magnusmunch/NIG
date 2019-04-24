# single VB update (tested)
.single.vb.update <- function(aold, bold, eta, theta, lambda, sv, n, p, uty, 
                              yty, conjugate, hyperprior) {
  
  # variables that change with conjugacy
  cold <- bold/(aold^(1 - conjugate))
  df <- n + conjugate*p
  
  # auxiliary variables involving mu and Sigma
  aux <- .aux.var(cold, aold, sv, uty, n, p)
  
  # vb parameters
  delta <- unname((bold/cold)*(aux$mutmu + aux$trSigma) + lambda)
  if(hyperprior=="inv. Gaussian") {
    b <- unname(sqrt(lambda/(theta^2*delta))*
                  ratio_besselK(sqrt(lambda*delta/theta^2), 
                                0.5*(p + eta)) + (p + eta)/delta)
    
    # needed in EB step
    e <- (b - (p + eta)/delta)*delta*theta^2/lambda
  } else {
    b <- unname((p + eta)/delta)
    e <- log(delta) - digamma(0.5*(p + eta)) - log(2)
  }
  zeta <- unname(0.5*(yty - 2*aux$ytXmu + aux$trXtXSigma + aux$mutXtXmu + 
                        ifelse(conjugate, b, 0)*(aux$trSigma + aux$mutmu)))
  a <- unname((df + 1)/(2*zeta))
  c <- b/(a^(1 - conjugate))
  
  # calculate the elbo separately for the inverse Gaussian and Gamma models
  if(hyperprior=="inv. Gaussian") {
    elbo <- .elbo.inv.gauss(aux, zeta, delta, a, b, c, e, lambda, theta, df, p, 
                            yty)
  } else {
    elbo <- .elbo.inv.gamma(aux, zeta, delta, a, b, c, e, lambda, eta, df, p, 
                            yty)
  }
  
  
  out <- c(delta=delta, zeta=zeta, a=a, b=b, e=e, elbo=unname(elbo$elbo), 
           elbo.const=unname(elbo$elbo.const))
  return(out)
}

# single VB update with intercept (not tested)
.single.vb.update.int <- function(aold, bold, eta, theta, lambda, sv, n, p, uty, 
                                  yty, sumu, sumy, conjugate, hyperprior) {
  
  # variables that change with conjugacy
  cold <- bold/(aold^(1 - conjugate))
  df <- n + conjugate*p
  
  # auxiliary variables involving mu and Sigma
  aux <- .aux.var.int(cold, aold, sv, uty, sumu, sumy, n, p)
  
  # vb parameters
  delta <- unname((bold/cold)*(aux$mutmu + aux$trSigma) + lambda)
  if(hyperprior=="inv. Gaussian") {
    b <- unname(sqrt(lambda/(theta^2*delta))*
                  ratio_besselK(sqrt(lambda*delta/theta^2), 
                                0.5*(p + eta)) + (p + eta)/delta)
    
    # needed in EB step
    e <- (b - (p + eta)/delta)*delta*theta^2/lambda
  } else {
    b <- unname((p + eta)/delta)
    e <- log(delta) - digamma(0.5*(p + eta)) - log(2)
  }
  zeta <- unname(0.5*(yty - 2*aux$ytXmu + aux$trXtXSigma + aux$mutXtXmu + 
                        ifelse(conjugate, b, 0)*(aux$trSigma + aux$mutmu)))
  a <- unname((df + 1)/(2*zeta))
  c <- b/(a^(1 - conjugate))
  
  # calculate the elbo separately for the inverse Gaussian and Gamma models
  if(hyperprior=="inv. Gaussian") {
    elbo <- .elbo.inv.gauss(aux, zeta, delta, a, b, c, e, lambda, theta, df, p, 
                            yty)
  } else {
    elbo <- .elbo.inv.gamma(aux, zeta, delta, a, b, c, e, lambda, eta, df, p, 
                            yty)
  }
  
  
  out <- c(delta=delta, zeta=zeta, a=a, b=b, e=e, elbo=unname(elbo$elbo), 
           elbo.const=unname(elbo$elbo.const))
  return(out)
}

# estimate variational Bayes posterior (not tested)
est.vb <- function(x, y, theta, lambda, intercept=TRUE, init=NULL, 
                   posterior=FALSE,
                   control=list(epsilon.vb=1e-3, maxit.vb=100, trace=TRUE)) {
  
  n <- nrow(x)
  p <- ncol(x)
  D <- ncol(y)
  
  svd.x <- svd(x)
  sv <- svd.x$d
  v <- svd.x$v
  uty <- t(svd.x$u) %*% y
  yty <- colSums(y^2)
  sumu <- colSums(svd.x$u)
  sumy <- colSums(y)
  
  # set initial values if none given
  if(is.null(init)) {
    init$a <- rep(1, D)
    init$b <- rep(1, D)
  }
  
  if(intercept) {
    old.vb <- sapply(c(1:D), function(d) {
      .single.vb.update.int(init$a[d], init$b[d], 1, theta,
                            lambda, sv, n, p, uty[, d], yty[d], sumu, sumy[d],
                            TRUE, "inv. Gaussian")})
  } else {
    old.vb <- sapply(c(1:D), function(d) {
      .single.vb.update(init$a[d], init$b[d], 1, theta,
                        lambda, sv, n, p, uty[, d], yty[d], TRUE,
                        "inv. Gaussian")})
  }
  
  
  # set all convergence check parameters
  conv.vb <- FALSE
  iter.vb <- 0
  check.vb <- FALSE
  while(!check.vb) {
    
    # increase iteration number by one
    iter.vb <- iter.vb + 1
    
    # update the VB parameters and elbo
    new.vb <- sapply(c(1:D), function(d) {
      .single.vb.update.int(old.vb[3, d], old.vb[4, d],  1, theta, lambda, sv, 
                            n, p, uty[, d], yty[d], sumu, sumy[d], TRUE, 
                            "inv. Gaussian")})
    
    # check convergence of the VB iterations
    conv.vb <- all(abs(new.vb - old.vb) <= control$epsilon.vb)
    check.vb <- conv.vb | (iter.vb >= control$maxit.vb)
    old.vb <- new.vb
  }
  
  # bold is same as cold here
  aux <- lapply(c(1:D), function(d) {
    .aux.var.int(new.vb[4, d], new.vb[3, d], sv, uty[, d], sumu, 
                 sumy[d], n, p)})
  if(D==1) {
    aux <- as.matrix(unlist(aux))
  } else {
    aux <- Reduce(function(s1, s2) {cbind(unlist(s1), unlist(s2))}, aux)
  }
  
  # creating Sigma and mu 
  mu <- sapply(c(1:D), function(d) {
    mat <- svd.x$v %*% (t(svd.x$u)*sv/(sv^2 + new.vb[4, d]))
    unname(c(aux[7, d]*(sumy[d] - aux[8, d]),
             aux[7, d]*(aux[8, d] - sumy[d])*rowSums(mat) + mat %*% y[, d]))})
  
  if(posterior) {
    Sigma <- sapply(c(1:D), function(d) {
      mat <- sumu %*% (t(v)*sv/(sv^2 + new.vb[4, d]))
      if(p > n) {
        v <- svd(x, nu=nrow(x), nv=ncol(x))$v
      } 
      out <- matrix(NA, nrow=p + intercept, ncol=p + intercept)
      out[1, 1] <- aux[7, d]/new.vb[3, d]
      out[1, 2:(p + 1)] <- out[2:(p + 1), 1] <- -aux[7, d]*as.numeric(mat)/
        new.vb[3, d]
      out[2:(p + 1), 2:(p + 1)] <- (v %*% (t(v)/(c(
        sv^2, rep(0, max(p - n, 0))) + new.vb[4, d])) + aux[7, d]*t(mat) %*% 
          mat)/new.vb[3, d]
      return(out)}, simplify=FALSE)
  } else {
    Sigma <- sapply(c(1:D), function(d) {
      unname(c(aux[7, d], colSums((t(v2)/sqrt(sv2^2 + new.vb[4, d]))^2) +
                 aux[7, d]*rowSums(svd.x$v %*% 
                                     (t(svd.x$u)*sv/(sv^2 + new.vb[4, d])))^2)/
               new.vb[3, d])})
  }
  
  # preparing the output
  out <- list(vb.post=list(mu=mu, Sigma=Sigma, 
                           delta=new.vb[1, ], zeta=new.vb[2, ]), 
              conv=conv.vb, iter=iter.vb)
  return(out)
}