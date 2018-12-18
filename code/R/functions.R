# computes ratio besselK(x, p/2 - 1/2)/besselK(x, p/2 + 1/2) (tested)
ratio.besselK <- function(x, p) {
  out <- sapply(x, function(s) {
    if((p %% 2)==0) {
      value <- 1
      for(i in 1:(p/2)) {
        value <- 1/(value + 2*(i - 0.5)/s)
      }
    } else {
      value <- besselK(s, 1, expon.scaled=TRUE)/besselK(s, 0, expon.scaled=TRUE)
      for(i in 1:(p/2 + 0.5)) {
        value <- 1/(value + 2*(i - 1)/s)
      }
    }
    return(value)})
  return(out)
}

# one VB update with sigma^2 stochastic (not tested)
vb.update <- function(zetaold, aold, lambda, theta, x, xxt, xtx, y) {
  p <- ncol(x)
  D <- length(length(zetaold))
  n <- nrow(y)
  
  # sort of efficient matrix inversion
  mat <- lapply(c(1:D), function(d) {
    mat <- xxt
    diag(mat) <- diag(mat) + aold[d]
    mat <- -(t(x) %*% solve(mat) %*% x)/aold[d]
    diag(mat) <- diag(mat) + 1/aold[d]
    return(mat)})
  
  # vb parameters
  Sigma <- lapply(c(1:D), function(d) {2*zetaold[d]*mat[[d]]/(n + p + 1)})
  mu <- sapply(c(1:D), function(d) {mat[[d]] %*% t(x) %*% y[, d]})
  delta <- sapply(c(1:D), function(d) {
    0.5*(n + p + 1)*(sum(diag(Sigma[[d]])) + sum(mu[, d]^2))/
      zetaold[d] + lambda})
  zeta <- sapply(c(1:D), function(d) {
    0.5*aold[d]*(sum(y[, d]^2) - 2*t(y[, d]) %*% x %*% mu[, d] + 
                   sum((xtx + diag(p))*Sigma[[d]]) + 
                   t(mu[, d]) %*% (xtx + diag(p)) %*% mu[, d])})
  
  # calculate the ratio of two BesselK functions
  ratio <- ratio_besselK_cpp(sqrt(lambda*delta/theta^2), p)
  
  # additional VB and EB parameters used in next iteration
  a <- sqrt(lambda/(theta^2*delta))*ratio + (p + 1)/delta
  v <- sqrt(delta*theta^2/lambda)*ratio
  
  out <- list(Sigma=Sigma, mu=mu, delta=delta, zeta=zeta, a=a, v=v)
  return(out)
}

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
  zeta <- 0.5*a*(yty + trSigma + mutmu + trXtXSigma + mutXtXmu - 2*ytXmu)
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

# naive updates to check the results
single.naive.update <- function(zetaold, aold, lambda, theta, n, p, x, y) {
  
  # these updates are computationally expensive
  Sigma <- 2*zetaold*solve(t(x) %*% x + aold*diag(p))/(n + p + 1)
  mu <- solve(t(x) %*% x + aold*diag(p)) %*% t(x) %*% y
  delta <- 0.5*(n + p + 1)*(sum(diag(Sigma)) + sum(mu^2))/zetaold + lambda
  a <- sqrt(lambda/(theta^2*delta))*
    ratio_besselK_cpp(sqrt(lambda*delta)/theta, p) + (p + 1)/delta
  zeta <- 0.5*a*(sum(y^2) - 2*as.numeric(t(y) %*% x %*% mu) + 
                   sum(diag((t(x) %*% x + diag(p)) %*% Sigma)) + 
                   as.numeric(t(mu) %*% (t(x) %*% x + diag(p)) %*% mu))
  v <- sqrt(delta*theta^2/lambda)*
    ratio_besselK_cpp(sqrt(lambda*delta)/theta, p)
  
  out <- c(delta=unname(delta), zeta=unname(zeta), a=unname(a), v=unname(v), 
           elbo=NA, elbo.const=NA)
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

# one VB update with independent beta and sigma^2 prior (not tested)
vb.update2 <- function(zetaold, aold, lambda, theta, x, xxt, xtx, y) {
  p <- ncol(x)
  D <- length(zetaold)
  n <- nrow(y)
  val <- (2*zetaold)/(n + 1)
  
  # vb parameters
  Sigma <- lapply(c(1:D), function(d) {
    mat <- xxt
    diag(mat) <- diag(mat) + val[d]*aold[d]
    mat <- -(t(x) %*% solve(mat) %*% x)/aold[d]
    diag(mat) <- diag(mat) + 1/aold[d]
    return(mat)})
  mu <- sapply(c(1:D), function(d) {Sigma[[d]] %*% t(x) %*% y[, d]/val[d]})
  delta <- sapply(c(1:D), function(d) {
    sum(diag(Sigma[[d]])) + sum(mu[, d]^2) + lambda})
  zeta <- sapply(c(1:D), function(d) {
    0.5*(sum(y[, d]^2) - 2*t(y[, d]) %*% x %*% mu[, d] + 
           sum(xtx*Sigma[[d]]) + t(mu[, d]) %*% xtx %*% mu[, d])})
  
  # calculate the ratio of two BesselK functions
  ratio <- ratio_besselK_cpp(sqrt(lambda*delta/theta^2), p)
  
  # additional VB and EB parameters used in next iteration
  a <- sqrt(lambda/(theta^2*delta))*ratio + (p + 1)/delta
  v <- sqrt(delta*theta^2/lambda)*ratio
  
  out <- list(Sigma=Sigma, mu=mu, delta=delta, zeta=zeta, a=a, v=v)
  return(out)
}

# fits model (not tested)
est.igauss <- function(x, y, C,
                       control=list(epsilon.eb=1e-3, epsilon.vb=1e-3, 
                                    maxit.eb=100, maxit.vb=100, trace=TRUE), 
                       init=NULL,
                       test=list(ind.var=FALSE)) {
  # init is either NULL or a list(alpha, lambda, a, zeta)
  # epsilon.vb=Inf means one vb iteration per eb iteration
  
  # fixed parameters
  n <- nrow(x)
  p <- ncol(x)
  D <- length(sigma)
  dcov <- ncol(C)
  xxt <- x %*% t(x)
  xtx <- t(x) %*% x
  
  # create initial values if none given
  if(is.null(init)) {
    svd.x <- svd(x)
    
    # maximise ridge mll to find gamma_d^2 and sigma_d^2 estimates
    if(test$ind.var) {
      init.mll.est <- t(sapply(c(1:D), function(d) {
        constrOptim(c(1, 1), init.mll2, NULL, ui=diag(2), ci=c(0, 0),
                    sv=svd.x$d, uty=t(svd.x$u) %*% y[, d], n=n,
                    control=list(fnscale=-1))$par}))
    } else {
      init.mll.est <- t(sapply(c(1:D), function(d) {
        constrOptim(c(1, 1), init.mll, NULL, ui=diag(2), ci=c(0, 0),
                    sv=svd.x$d, uty=t(svd.x$u) %*% y[, d], n=n,
                    control=list(fnscale=-1))$par}))
    }
    
    # consider them samples from prior and estimate one inverse Gaussian prior
    init.theta <- mean(init.mll.est[, 2])
    init.lambda <- D/sum(1/init.mll.est[, 2] - 1/init.theta)
    
    # use empirical mode of gamma_d^2 to derive one delta
    d <- density(init.mll.est[, 2])
    mode <- d$x[which.max(d$y)]
    init.delta <- as.numeric(mode^2*init.lambda/init.theta^2 + (p + 3)*mode)
    init.a <- sqrt(init.lambda/(init.theta^2*init.delta))*
      ratio_besselK_cpp(sqrt(init.lambda*init.delta)/init.theta, p) +
      (p + 1)/init.delta
    
    # consider sigma_d^2 fixed and estimate inverse Gamma prior
    if(test$ind.var) {
      init.zeta <- 0.5*D*(n + 1)/sum(1/init.mll.est[, 1])
    } else {
      init.zeta <- 0.5*D*(n + p + 1)/sum(1/init.mll.est[, 1])
    }
    
    
    init <- list(alpha=c(init.theta, rep(0, ncol(C) - 1)), lambda=init.lambda,
                 a=rep(init.a, D), zeta=rep(init.zeta, D))
    
  }
  
  # set initial values for VB, EB, and ELBO estimates
  old.eb <- list(alpha=init$alpha, lambda=init$lambda, 
                 theta=as.numeric(C %*% init$alpha))
  if(test$ind.var) {
    old.vb <- vb.update2(init$zeta, init$a, old.eb$lambda, old.eb$theta, x, xxt, 
                         xtx, y)
  } else {
    old.vb <- vb.update(init$zeta, init$a, old.eb$lambda, old.eb$theta, x, xxt, 
                        xtx, y)
  }
  old.elbo <- -Inf
  
  # prepare objects to store EB and ELBO iterations
  eb.seq <- list(alpha=NULL, lambda=NULL, theta=NULL)
  elbo.seq <- numeric(0)
  
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
    new.eb <- eb.update(old.vb$v, old.vb$a, -Inf, C, D)
    eb.seq <- sapply(c(1:length(eb.seq)), function(s) {
      rbind(eb.seq[[s]], new.eb[[s]])})
    
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
      
      # update the VB parameters
      if(test$ind.var) {
        new.vb <- vb.update2(old.vb$zeta, old.vb$a, old.eb$lambda, old.eb$theta, 
                             x, xxt, xtx, y)
      } else {
        new.vb <- vb.update(old.vb$zeta, old.vb$a, old.eb$lambda, old.eb$theta, 
                            x, xxt, xtx, y)
      }
      
      # ELBO calculation is still problematic
      # new.elbo <- elbo.eval(new.vb$mu, new.vb$Sigma, new.vb$delta,
      #                       new.eb$lambda, new.eb$theta, x, xtx, y)
      new.elbo <- rep(-Inf, D)
      elbo.seq <- rbind(elbo.seq, new.elbo)
      # 
      # conv.vb[iter.eb] <- ((new.elbo - old.elbo) <= control$epsilon.vb) |
      #   (iter.vb[iter.eb] > maxit.vb)
      # old.elbo <- new.elbo
      
      # check convergence of the VB iterations
      conv.vb[iter.eb] <- all(abs(unlist(new.vb) - unlist(old.vb)) < 
                                control$epsilon.vb)
      check.vb <- conv.vb[iter.eb] | (iter.vb[iter.eb] >= control$maxit.vb)
      old.vb <- new.vb
    }
  }
  names(eb.seq) <- names(old.eb)[-4]
  
  # preparing the output
  out <- list(vb.post=list(mu=old.vb$mu, Sigma=old.vb$Sigma, 
                           delta=old.vb$delta, zeta=old.vb$zeta), 
              eb.seq=eb.seq,
              elbo.seq=elbo.seq,
              conv=list(eb=conv.eb, vb=conv.vb),
              iter=list(eb=iter.eb, vb=iter.vb))
  return(out)
}

# fits independent inverse Gamma prior (Gwen's) model (not tested)
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
    
    # check convergence of EB iteration
    conv.eb <- all(c(abs(aprior[iter.eb + 1, ] - aprior[iter.eb, ]),
                     abs(bprior[iter.eb + 1, ] - bprior[iter.eb, ])) <= 
                     control$epsilon.eb)
    
    check.eb <- conv.eb | (iter.eb >= control$maxit.eb)
  }
  
  # create an output list
  out <- list(vb.post=list(mu=t(sapply(new.vb$beta, function(s) {s[, 1]})),
                           dSigma=t(sapply(new.vb$beta, function(s) {s[, 2]})),
                           apost=sapply(new.vb$tau, function(s) {s[, 1]}),
                           bpost=sapply(new.vb$tau, function(s) {s[, 2]}),
                           cpost=sapply(new.vb$sigma, function(s) {s[1, ]}),
                           dpost=sapply(new.vb$sigma, function(s) {s[2, ]})),
              eb.seq=list(aprior=aprior, bprior=bprior),
              elbo.seq=allMLs,
              conv=list(eb=conv.eb, vb=conv.vb),
              iter=list(eb=iter.eb, vb=iter.vb))
  return(out)
  
}

# MLL in regular ridge model (not tested)
init.mll <- function(par, sv, uty, n) {
  sigma.sq <- par[1]
  gamma.sq <- par[2]
  mll <- -0.5*n*log(sigma.sq) - 0.5*n*log(gamma.sq) - 
    sum(log(sv^2 + 1/gamma.sq)) - 0.5*sum(uty^2/(sv^2 + 1/gamma.sq))/
    (sigma.sq*gamma.sq)
  return(mll)
}

# MLL in regular ridge model for independent beta and sigma^2 prior (not tested)
init.mll2 <- function(par, sv, uty, n) {
  sigma.sq <- par[1]
  gamma.sq <- par[2]
  mll <- -0.5*n*log(sigma.sq) - 0.5*n*log(gamma.sq) - 
    sum(log(sv^2/sigma.sq + 1/gamma.sq)) - 
    0.5*sum(uty^2/(sv^2*gamma.sq + sigma.sq))
  return(mll)
}


