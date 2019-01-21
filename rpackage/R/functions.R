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

# single VB update, drugs covariate model only (tested)
.single.vb.update <- function(aold, bold, eta, theta, lambda, sv, n, p, uty, 
                              yty, conjugate, inv.gauss) {
  
  # variables that change with conjugacy
  cold <- bold/(aold^(1 - conjugate))
  df <- n + conjugate*p
  
  # auxiliary variables involving mu and Sigma
  aux <- .aux.var(cold, aold, sv, uty, n, p)
  
  # vb parameters
  delta <- (bold/cold)*aux$mutmu + lambda
  if(inv.gauss) {
    b <- sqrt(lambda/(theta^2*delta))*
      ratio_besselK_cpp(sqrt(lambda*delta/theta^2), 0.5*(p + eta)) + 
      (p + eta)/delta
    
    # needed in EB step
    e <- (b - (p + eta)/delta)*delta*theta^2/lambda
  } else {
    b <- (p + eta)/delta
    e <- log(delta) - digamma(0.5*(p + eta)) - log(2)
  }
  zeta <- 0.5*(yty - 2*aux$ytXmu + aux$trXtXSigma + aux$mutXtXmu + 
                 ifelse(conjugate, b, 0)*(aux$trSigma + aux$mutmu))
  a <- (df + 1)/(2*zeta)
  c <- b/(a^(1 - conjugate))
  
  # calculate the elbo separately for the inverse Gaussian and Gamma models
  if(inv.gauss) {
    elbo <- .elbo.inv.gauss(aux, zeta, delta, a, b, c, e, lambda, theta, df, p, 
                            yty)
  } else {
    elbo <- .elbo.inv.gamma(aux, zeta, delta, a, b, c, e, lambda, eta, df, p, 
                            yty)
  }
  
  
  out <- c(delta=delta, zeta=zeta, a=a, b=b, e=e, elbo=elbo$elbo, 
           elbo.const=elbo$elbo.const)
  return(out)
}

# estimating equation for eta for root-finding (tested)
.froot.inv.gamma <- function(eta, lambda, esum, size) {
  log(lambda) - log(2) - digamma(eta/2) - esum/size
}

# inverse Gamma EB update (tested)
.eb.update.inv.gamma <- function(esum, bsum, b, delta, old.eta, elbo.const, 
                                 nclass, sclass, p, epsilon, maxit) {
  
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
                              
  out <- list(eta=eta, lambda=lambda, elbo=elbo, conv=conv)
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
  
  out <- list(theta=theta, alpha=alpha, lambda=lambda, elbo=elbo, conv=conv)
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
.init.param <- function(C, sv, uty, D, n, p) {
  # maximise ridge mll to find gamma_d^2 and sigma_d^2 estimates
  init.mll.est <- t(sapply(c(1:D), function(d) {
      constrOptim(c(1, 1), ridge.mll, NULL, ui=diag(2), ci=c(0, 0),
                  sv=sv, uty=uty[, d], n=n, control=list(fnscale=-1))$par}))
    
  # consider them samples from prior and estimate one inverse Gaussian prior
  theta <- mean(init.mll.est[, 2])
  lambda <- D/sum(1/init.mll.est[, 2] - 1/theta)
  
  # use empirical mode of gamma_d^2 to derive one delta
  d <- density(init.mll.est[, 2])
  mode <- d$x[which.max(d$y)]
  delta <- as.numeric(mode^2*lambda/theta^2 + (p + 3)*mode)
  a <- sqrt(lambda/(theta^2*delta))*
    ratio_besselK_cpp(sqrt(lambda*delta)/theta, p) + (p + 1)/delta
  
  # consider sigma_d^2 fixed and estimate inverse Gamma prior
  zeta <- 0.5*D*(n + p + 1)/sum(1/init.mll.est[, 1])
  
  out <- list(theta=rep(theta, D), alpha=c(theta, rep(0, ncol(C) - 1)), 
              lambda=lambda, zeta=rep(zeta, D), a=rep(a, D))
  return(out)
  
}

# initial values with independent beta and sigma^2 prior (not tested)
.init.param.ind <- function(C, sv, uty, D, n, p) {
  # maximise ridge mll to find gamma_d^2 and sigma_d^2 estimates
  init.mll.est <- t(sapply(c(1:D), function(d) {
        constrOptim(c(1, 1), ridge.mll.ind, NULL, ui=diag(2), ci=c(0, 0),
                    sv=sv, uty=uty[, d], n=n, control=list(fnscale=-1))$par}))
  
  # consider them samples from prior and estimate one inverse Gaussian prior
  theta <- mean(init.mll.est[, 2])
  lambda <- D/sum(1/init.mll.est[, 2] - 1/theta)
  
  # use empirical mode of gamma_d^2 to derive one delta
  d <- density(init.mll.est[, 2])
  mode <- d$x[which.max(d$y)]
  delta <- as.numeric(mode^2*lambda/theta^2 + (p + 3)*mode)
  a <- sqrt(lambda/(theta^2*delta))*
    ratio_besselK_cpp(sqrt(lambda*delta)/theta, p) + (p + 1)/delta
  
  # consider sigma_d^2 fixed and estimate inverse Gamma prior
  zeta <- 0.5*D*(n + 1)/sum(1/init.mll.est[, 1])
  
  out <- list(theta=rep(theta, D), alpha=c(theta, rep(0, ncol(C) - 1)), 
              lambda=lambda, zeta=rep(zeta, D), a=rep(a, D))
  return(out)
}

# fit model without tissue effect and molecular feature groups (not tested)
est.model <- function(x, y, C, inv.gauss, conjugate,
                      control=list(epsilon.eb=1e-3, epsilon.vb=1e-3, 
                                    maxit.eb=20, maxit.vb=2, trace=TRUE), 
                      init=NULL) {
  # init is either NULL or list(alpha, lambda, a, zeta)
  # control$epsilon.vb=Inf is equivalent to control$maxit.vb=1
  # inv.gauss indicates inv. Gaussian model (TRUE) or inv. Gamma (FALSE)
  # conjugate indicates independent beta and sigma^2 (FALSE) or dependent (TRUE)
  
  # fixed parameters
  n <- nrow(x)
  p <- ncol(x)
  D <- ncol(y)
  dcov <- ncol(C)
  svd.x <- svd(x)
  sv <- svd.x$d
  uty <- t(svd.x$u) %*% y
  yty <- colSums(y^2)
  
  # create initial values if none given
  if(is.null(init)) {
    if(test$ind) {
      init <- init.param.ind(C, sv, uty, D, n, p)
    } else {
      init <- init.param(C, sv, uty, D, n, p)
    }
  } else {
    init$theta <- 1/as.numeric(C %*% init$alpha)
  }

  # set initial values for VB, EB, and ELBO estimates
  old.eb <- list(theta=init$theta, alpha=init$alpha, lambda=init$lambda, 
                 elbo=NA)
  if(test$ind) {
    old.vb <- sapply(c(1:D), function(d) {
      single.vb.update.ind(init$zeta[d], init$a[d], init$lambda, init$theta[d], 
                           sv, n, p, uty[, d], yty[d])})
  } else {
    old.vb <- sapply(c(1:D), function(d) {
      single.vb.update(init$zeta[d], init$a[d], init$lambda, init$theta[d], 
                       sv, n, p, uty[, d], yty[d])})
  }
  old.elbo <- old.vb[5, ]
  
  # prepare objects to store EB and ELBO iterations
  seq.eb <- list(theta=numeric(0), alpha=numeric(0), lambda=numeric(0))
  seq.elbo <- numeric(0)
  
  # set all convergence check parameters
  conv.eb <- FALSE
  conv.vb <- NULL
  iter.eb <- 0
  iter.vb <- NULL
  check.eb <- FALSE
  check.vb <- FALSE
  
  # outer EB loop
  while(!check.eb) {
    
    # increase iteration number by one
    iter.eb <- iter.eb + 1
    
    # possibly print iteration number
    if(control$trace) {cat("\r", "iteration", iter.eb)}
    
    # update the EB parameters and store result in iteration track object
    new.eb <- eb.update(old.vb[4, ], old.vb[3, ], old.vb[6, ], C, D)
    new.elbo <- new.eb$elbo
    seq.eb <- sapply(c(1:length(seq.eb)), function(s) {
      rbind(seq.eb[[s]], new.eb[[s]])})
    seq.elbo <- rbind(seq.elbo, new.elbo)
    
    # check convergence of EB estimates
    conv.eb <- all((abs(c(new.eb$lambda - old.eb$lambda, 
                          new.eb$theta - old.eb$theta)) <= control$epsilon.eb))
    check.eb <- conv.eb | (iter.eb >= control$maxit.eb)
    old.eb <- new.eb
    
    
    # set convergence check parameters for VB
    conv.vb <- c(conv.vb, FALSE)
    check.vb <- FALSE
    iter.vb <- c(iter.vb, 0)
    
    # inner VB loop
    while(!check.vb) {
      
      # increase iteration number by one
      iter.vb[iter.eb] <- iter.vb[iter.eb] + 1
      
      # update the VB parameters and elbo
      if(test$ind) {
        new.vb <- sapply(c(1:D), function(d) {
          single.vb.update.ind(old.vb[2, d], old.vb[3, d], old.eb$lambda, 
                               old.eb$theta[d], sv, n, p, uty[, d], yty[d])})
      } else {
        new.vb <- sapply(c(1:D), function(d) {
          single.vb.update(old.vb[2, d], old.vb[3, d], old.eb$lambda, 
                           old.eb$theta[d], sv, n, p, uty[, d], yty[d])})
      }
      new.elbo <- new.vb[5, ]
      seq.elbo <- rbind(seq.elbo, new.elbo)

      # check convergence of the VB iterations
      conv.vb[iter.eb] <- all(abs(new.elbo - old.elbo) <= control$epsilon.vb)
      check.vb <- conv.vb[iter.eb] | (iter.vb[iter.eb] >= control$maxit.vb)
      old.elbo <- new.elbo
      old.vb <- new.vb
    }
  }
  names(seq.eb) <- names(old.eb)[-4]
  
  # creating Sigma and mu 
  mu <- sapply(c(1:D), function(d) {
    svd.x$v %*% (t(svd.x$u)*sv/(sv^2 + new.vb[3, d])) %*% y[, d]})
  if(p > n) {
    svd.x <- svd(x, nu=nrow(x), nv=ncol(x))
  }
  dSigma <- sapply(c(1:D), function(d) {
    2*new.vb[2, d]*colSums(t(svd.x$v)*(t(svd.x$v)/(c(sv^2, rep(
      0, max(p - n, 0))) + new.vb[3, d])))/(n + p + 1)})
  
  # preparing the output
  out <- list(vb.post=list(mu=mu, dSigma=dSigma, 
                           delta=new.vb[1, ], zeta=new.vb[2, ]), 
              seq.eb=seq.eb,
              seq.elbo=seq.elbo,
              conv=list(eb=conv.eb, vb=conv.vb),
              iter=list(eb=iter.eb, vb=iter.vb))
  return(out)
}

# fits Gwen's independent inverse Gamma priors model (tested)
est.gwen <- function(x, y, eqid,
                     control=list(epsilon.eb=1e-3, epsilon.vb=1e-3, 
                                  epsilon.opt=1e-6, maxit.eb=20, maxit.vb=2, 
                                  maxit.opt=20, conv.vb=c("elbo", "param"), 
                                  trace=TRUE), 
                     init=list(aprior=0.001, bprior=0.001)) {
  
  p <- ncol(x)
  D <- ncol(y)
  nclass <- length(unique(eqid))
  listX <- rep(list(x), D)
  
  # Generate matrix P (3 global shrinkage priors)
  listP <- lapply(c(1:D), function(d) {rep(eqid[d], p)})
  
  # Generate outcomes
  listy <- unname(split(y, rep(1:ncol(y), each=nrow(y))))
  
  # Input arguments
  EBid <- 1:nclass # indexes of gamma priors for which EB is desired. If no EB then non-informative. 
  
  # Intialization
  priors <- sort(unique(unlist(listP)))
  nbpriors <- length(priors)
  aprior <- matrix(init$aprior, nrow=1, ncol=nbpriors)
  bprior <- matrix(init$bprior, nrow=1, ncol=nbpriors)
  colnames(aprior) <- paste("a", 1:nbpriors, sep="")
  colnames(bprior) <- paste("b", 1:nbpriors, sep="")
  idxPriorList <- lapply(listP, function(x){sort(unique(x))})
  allMLs <- numeric(0)
  
  # set initial values for convergence check reasons
  old.elbo <- -Inf
  old.vb <- list(tau=lapply(c(1:D), matrix, data=Inf, nrow=1, ncol=2),
                 sigma=lapply(c(1:D), matrix, data=Inf, nrow=2, ncol=1),
                 beta=lapply(c(1:D), matrix, data=Inf, nrow=p, ncol=2))
  
  # outer EB loop
  conv.eb <- FALSE
  conv.vb <- NULL
  iter.eb <- 0
  iter.vb <- NULL
  check.eb <- FALSE
  check.vb <- FALSE
  while(!check.eb) {
    iter.eb <- iter.eb + 1
    if(control$trace) {cat("\r", "iteration", iter.eb)}
    
    # Prior as lists
    inputaList <- lapply(idxPriorList, function(x){aprior[iter.eb, x]})
    inputbList <- lapply(idxPriorList, function(x){bprior[iter.eb, x]})
    
    if(iter.eb==1) {
      mydstarvec <- rep(0, length(inputaList))
      mybstarlist <- lapply(idxPriorList, function(x){rep(0, length(x))})
    } else {
      mydstarvec <- unlist(lapply(new.res$postSigList, function(x) {x[2]}))
      mybstarlist <- lapply(new.res$postRandList, function(x) {x[, 2]})
    }
    
    # Fit BSEM
    conv.vb <- c(conv.vb, FALSE)
    check.vb <- FALSE
    iter.vb <- c(iter.vb, 0)
    while(!check.vb) {
      iter.vb[iter.eb] <- iter.vb[iter.eb] + 1
      
      # one iteration of the BSEM VB step
      new.res <- BSEMVarOneIter(ylist=listy, Xlist=listX, Plist=listP, 
                                alist=inputaList, blist=inputbList, 
                                bstarlist=mybstarlist, cSigma=0.001, 
                                dSigma=0.001, dstarvec=mydstarvec, 
                                lincomblist=list())
      
      # storing the new ELBO
      new.elbo <- as.numeric(new.res$allmargs)
      allMLs <- rbind(allMLs, new.elbo)
      
      # storing the new variational parameters
      new.vb <- list(tau=new.res$postRandList, sigma=new.res$postSigList, 
                     beta=lapply(new.res$postBetaList, function(s) {
                       cbind(s[, 1], s[, 2]^2)}))
      
      # check convergence of the VB, either by the ELBO or by the parameters
      if(control$conv.vb[1]=="elbo") {
        conv.vb[iter.eb] <- sum(abs(new.elbo - old.elbo)) < 
          control$epsilon.vb
      } else {
        conv.vb[iter.eb] <- all(abs(unlist(new.vb) - unlist(old.vb)) < 
                                  control$epsilon.vb)
      }
      
      check.vb <- conv.vb[iter.eb] | (iter.vb[iter.eb] >= control$maxit.vb)
      
      old.elbo <- new.elbo
      old.vb <- new.vb
    }
    
    # Empirical Bayes
    aprior <- rbind(aprior, NA)
    bprior <- rbind(bprior, NA)
    for(ii in 1:nbpriors){
      if(ii%in%EBid){
        # Get posterior shape and rate parameters
        allaRandStar <- sapply(1:length(idxPriorList), function(x) {
          new.res$postRandList[[x]][idxPriorList[[x]]==ii,1]}, simplify=TRUE)
        if(is.list(allaRandStar)){
          allaRandStar <- unlist(allaRandStar)
        }
        allbRandStar <- sapply(1:length(idxPriorList), function(x) {
          new.res$postRandList[[x]][idxPriorList[[x]]==ii,2]}, simplify=TRUE)
        if(is.list(allbRandStar)){
          allbRandStar <- unlist(allbRandStar)
        }
        
        # Variational Empirical Bayes using fixed-point iteration as in Valpola and Honkela (2006)
        ab <- c(aprior[iter.eb, ii], bprior[iter.eb, ii])
        ab <- fixedPointIterEB(initab=ab, myallaRandStar=allaRandStar, 
                               myallbRandStar=allbRandStar, 
                               mymaxiter=control$maxit.opt, 
                               myeps=control$maxit.opt)
        aprior[iter.eb + 1, ii] <- ab[1]
        bprior[iter.eb + 1, ii] <- ab[2]
        
      }
    }
    
    # to make sure the iteration numbers match with other models
    allMLs <- rbind(allMLs, new.elbo)
    
    # check convergence of EB iteration
    conv.eb <- all(c(abs(aprior[iter.eb + 1, ] - aprior[iter.eb, ]),
                     abs(bprior[iter.eb + 1, ] - bprior[iter.eb, ])) <= 
                     control$epsilon.eb)
    
    check.eb <- conv.eb | (iter.eb >= control$maxit.eb)
  }
  
  # create an output list
  out <- list(vb.post=list(mu=sapply(new.vb$beta, function(s) {s[, 1]}),
                           dSigma=sapply(new.vb$beta, function(s) {s[, 2]}),
                           apost=sapply(new.vb$tau, function(s) {s[, 1]}),
                           bpost=sapply(new.vb$tau, function(s) {s[, 2]}),
                           cpost=sapply(new.vb$sigma, function(s) {s[1, ]}),
                           dpost=sapply(new.vb$sigma, function(s) {s[2, ]})),
              seq.eb=list(aprior=aprior[-1, ], bprior=bprior[-1, ]),
              seq.elbo=allMLs,
              conv=list(eb=conv.eb, vb=conv.vb),
              iter=list(eb=iter.eb, vb=iter.vb))
  return(out)
  
}




