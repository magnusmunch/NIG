# one fast VB update with sigma^2 stochastic (tested)
single.vb.update <- function(zetaold, aold, lambda, theta, sv, n, p, uty, yty) {
  
  # auxiliary variables involving mu and Sigma
  trSigma <- 2*zetaold*(sum(1/(sv^2 + aold)) + max(p - n, 0)/aold)/(n + p + 1)
  mutmu <- sum(sv^2*uty^2/(sv^2 + aold)^2)
  trXtXSigma <- 2*zetaold*sum(sv^2/(sv^2 + aold))/(n + p + 1)
  mutXtXmu <- sum(sv^4*uty^2/(sv^2 + aold)^2)
  ytXmu <- sum(sv^2*uty^2/(sv^2 + aold))
  logdetSigma <- p*log(zetaold) - sum(log(sv^2 + aold)) - 
    max(p - n, 0)*log(aold)
  
  # vb parameters and auxiliaries
  delta <- 0.5*(n + p + 1)*(trSigma + mutmu)/zetaold + lambda
  ratio <- ratio_besselK_cpp(sqrt(lambda*delta/theta^2), p)
  a <- sqrt(lambda/(theta^2*delta))*ratio + (p + 1)/delta
  zeta <- 0.5*(yty + a*trSigma + a*mutmu + trXtXSigma + mutXtXmu - 2*ytXmu)
  v <- sqrt(delta*theta^2/lambda)*ratio
  
  # calculate the elbo part that is constant after next eb update
  elbo.const <- 0.5*logdetSigma - 0.5*(n + p + 1)*log(zeta) - 
    0.25*(n + p + 1)*(yty + mutXtXmu + trXtXSigma - 2*ytXmu)/zeta -
    0.25*(p + 1)*log(delta) + 0.25*(p + 1)*log(lambda) - 
    0.5*(p + 1)*log(theta) + 
    bessel_lnKnu(0.5*(p + 1), sqrt(lambda*delta)/theta) + 
    0.5*(delta - 0.5*(n + p + 1)*(mutmu + trSigma)/zeta)*a +
    0.5*lambda*v/theta^2
  
  # total elbo calculation
  elbo <- elbo.const  + 0.5*log(lambda) - 0.5*lambda*v/theta^2 + lambda/theta -
    0.5*lambda*a
  
  out <- c(delta=unname(delta), zeta=unname(zeta), a=unname(a), v=unname(v), 
           elbo=unname(elbo), elbo.const=unname(elbo.const))
  return(out)
}

# one fast VB update with independent beta and sigma^2 prior (tested)
single.vb.update.ind <- function(zetaold, aold, lambda, theta, sv, n, p, uty, 
                                 yty) {
  
  # auxiliary variables involving mu and Sigma
  val <- 2*zetaold/(n + 1)
  trSigma <- val*sum(1/(sv^2 + val*aold)) + max(p - n, 0)/aold
  mutmu <- sum(sv^2*uty^2/(sv^2 + val*aold)^2)
  trXtXSigma <- val*sum(sv^2/(sv^2 + val*aold))
  mutXtXmu <- sum(sv^4*uty^2/(sv^2 + val*aold)^2)
  ytXmu <- sum(sv^2*uty^2/(sv^2 + val*aold))
  logdetSigma <- p*log(zetaold) - sum(log(sv^2 + val*aold)) - 
    max(p - n, 0)*(log(aold) + log(val))
  
  # vb parameters and auxiliaries
  delta <- trSigma + mutmu + lambda
  ratio <- ratio_besselK_cpp(sqrt(lambda*delta/theta^2), p)
  a <- sqrt(lambda/(theta^2*delta))*ratio + (p + 1)/delta
  zeta <- 0.5*(yty + trXtXSigma + mutXtXmu - 2*ytXmu)
  v <- sqrt(delta*theta^2/lambda)*ratio
  
  # calculate the elbo part that is constant after next eb update
  elbo.const <- 0.5*logdetSigma - 0.5*(n + 1)*log(zeta) - 
    0.25*(n + 1)*(yty + mutXtXmu + trXtXSigma - 2*ytXmu)/zeta -
    0.25*(p + 1)*log(delta) + 0.25*(p + 1)*log(lambda) - 
    0.5*(p + 1)*log(theta) + 
    bessel_lnKnu(0.5*(p + 1), sqrt(lambda*delta)/theta) + 
    0.5*(delta - mutmu - trSigma)*a + 0.5*lambda*v/theta^2
  
  # total elbo calculation
  elbo <- elbo.const  + 0.5*log(lambda) - 0.5*lambda*v/theta^2 + lambda/theta -
    0.5*lambda*a
  
  out <- c(delta=unname(delta), zeta=unname(zeta), a=unname(a), v=unname(v), 
           elbo=unname(elbo), elbo.const=unname(elbo.const))
  return(out)
  
}

# one EB update (tested)
eb.update <- function(v, a, elbo.const, C, D) {
  
  # eb parameters
  alpha <- rowSums(solve(t(C*v) %*% C) %*% t(C))
  lambda <- D/(sum(a) - sum(t(C)*alpha))
  theta <- 1/as.numeric(C %*% alpha)
  
  # elbo calculation involves constant part plus updated part
  elbo <- elbo.const + 0.5*log(lambda) - 0.5*lambda*v/theta^2 + lambda/theta -
    0.5*lambda*a
  
  out <- list(theta=theta, alpha=alpha, lambda=lambda, elbo=elbo)
  return(out)
}

# ridge marginal log likelihood of ridge model (not tested)
ridge.mll <- function(par, sv, uty, n) {
  sigma.sq <- par[1]
  gamma.sq <- par[2]
  mll <- -0.5*n*log(sigma.sq) - 0.5*n*log(gamma.sq) - 
    sum(log(sv^2 + 1/gamma.sq)) - 0.5*sum(uty^2/(sv^2 + 1/gamma.sq))/
    (sigma.sq*gamma.sq)
  return(mll)
}

# ridge marginal log lik for independent beta and sigma^2 prior (not tested)
ridge.mll.ind <- function(par, sv, uty, n) {
  sigma.sq <- par[1]
  gamma.sq <- par[2]
  mll <- -0.5*n*log(sigma.sq) - 0.5*n*log(gamma.sq) - 
    sum(log(sv^2/sigma.sq + 1/gamma.sq)) - 
    0.5*sum(uty^2/(sv^2*gamma.sq + sigma.sq))
  return(mll)
}

# initial values (not tested)
init.param <- function(C, sv, uty, D, n, p) {
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
init.param.ind <- function(C, sv, uty, D, n, p) {
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
est.igauss <- function(x, y, C,
                       control=list(epsilon.eb=1e-3, epsilon.vb=1e-3, 
                                    maxit.eb=20, maxit.vb=2, trace=TRUE), 
                       init=NULL,
                       test=list(ind=FALSE)) {
  # init is either NULL or list(alpha, lambda, a, zeta)
  # control$epsilon.vb=Inf is equivalent to control$maxit.vb=1
  # test$ind=TRUE estimates the model with independent beta and sigma^2
  
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