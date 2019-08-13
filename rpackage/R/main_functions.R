# x <- rep(list(x), D)
# mult.lambda=TRUE
# intercept.eb=TRUE
# fixed.eb=c("none")
# full.post=FALSE
# init=NULL
# control=list(conv.post=TRUE, trace=TRUE,
#              epsilon.eb=1e-3, epsilon.vb=1e-3,
#              epsilon.opt=sqrt(.Machine$double.eps),
#              maxit.eb=10, maxit.vb=2, maxit.opt=100,
#              maxit.post=100)

# estimate ENIG model (not tested)
enig <- function(x, y, C, mult.lambda=FALSE, intercept.eb=TRUE,
                 fixed.eb=c("none", "lambda", "both"), full.post=FALSE, 
                 init=NULL,
                 control=list(conv.post=TRUE, trace=TRUE,
                              epsilon.eb=1e-3, epsilon.vb=1e-3, 
                              epsilon.opt=sqrt(.Machine$double.eps),
                              maxit.eb=10, maxit.vb=2, maxit.opt=100,
                              maxit.post=100)) {
  
  # save the arguments
  cl <- match.call()
  
  # fixed parameters
  n <- nrow(y)
  p <- sapply(x, ncol)
  D <- ncol(y)
  ncd <- ifelse(is.null(C), 0, ncol(C[[1]]))
  
  # add intercept and transform co-data to matrix
  if(intercept.eb & is.null(C)) {
    C <- lapply(1:D, function(d) {
      m <- matrix(1, nrow=p[d]); colnames(m) <- "intercept"; return(m)})
  } else if(intercept.eb) {
    C <- lapply(C, function(cc) {cbind(intercept=1, cc)})
  }
  Cmat <- Reduce("rbind", C)
  
  # fixed data functions
  yty <- colSums(y^2)
  ytx <- lapply(1:D, function(d) {as.numeric(y[, d] %*% x[[d]])})
  
  # set initial values if not given
  if(is.null(init)) {
    init$aold <- rep(1, D)
    init$bold <- lapply(1:D, function(d) {rep(1, p[d])})
    init$lambda <- rep(1, ifelse(mult.lambda, D, 1))
    if(intercept.eb) {
      init$alpha <- c(1, rep(0, ncd))
      init$Calpha <- lapply(1:D, function(d) {
        as.numeric(C[[d]] %*% init$alpha)})
    } else {
      init$alpha <- rep(0, ncd)
      init$Calpha <- lapply(1:D, function(d) {rep(1, p[d])})
    }
  }
  init$mprior <- lapply(1:D, function(d) {1/(init$Calpha[[d]])})
  init$vprior <- lapply(1:D, function(d) {
    1/(init$Calpha[[d]]*init$lambda[ifelse(mult.lambda, d, 1)])})
  
  # set initial values for VB, EB, and ELBO estimates
  old.eb <- list(alpha=init$alpha, lambda=init$lambda, mprior=init$mprior,
                 vprior=init$vprior, Calpha=init$Calpha, conv=FALSE, iter=0)
  old.vb <- lapply(c(1:D), function(d) {
    .single.vb.update.enig(init$aold[d], init$bold[[d]], init$Calpha[[d]], 
                           init$lambda[ifelse(mult.lambda, d, 1)], y[, d], 
                           x[[d]], ytx[[d]], yty[d], n, p[d], D)})
  old.vb <- setNames(lapply(names(old.vb[[1]]), function(var) {
    lapply(old.vb, "[[", var)}), names(old.vb[[1]]))
  
  # prepare objects to store EB and ELBO iterations
  seq.eb <- list(alpha=old.eb$alpha, lambda=old.eb$lambda, 
                 mprior=lapply(1:D, function(d) {1/old.eb$Calpha[[d]]}), 
                 vprior=lapply(1:D, function(d) {
                   1/(old.eb$Calpha[[d]]^3*
                        old.eb$lambda[ifelse(mult.lambda, d, 1)])}))
  
  # set all convergence check parameters
  conv.eb <- FALSE
  conv.vb <- logical(0)
  conv.opt <- logical(0)
  iter.eb <- 0
  iter.vb <- numeric(0)
  iter.opt <- numeric(0)
  check.eb <- FALSE
  check.vb <- FALSE
  
  # outer EB loop
  while(!check.eb) {
    
    # increase iteration number by one
    iter.eb <- iter.eb + 1
    
    # if EB is required
    if(fixed.eb!="both") {
      
      # possibly print iteration number
      if(control$trace & fixed.eb!="both") {cat("\r", "iteration", iter.eb)}
      
      # update the EB parameters
      new.eb <- .eb.update.enig(old.vb$e, old.vb$b, old.eb$lambda, Cmat, p, D, 
                                mult.lambda, fixed.eb, control$epsilon.opt, 
                                control$maxit.opt) 
      
      # paste new hyperparameters to previous
      seq.eb <- lapply(1:length(seq.eb), function(s) {
        if(!is.list(seq.eb[[s]])) {
          rbind(seq.eb[[s]], new.eb[[s]])
        } else {
          lapply(1:D, function(d) {rbind(seq.eb[[s]][[d]], new.eb[[s]][[d]])})
        }})
      iter.opt <- c(iter.opt, new.eb$iter)

      # check convergence of EB estimates
      conv.eb <- all(abs(unlist(new.eb[c("mprior", "vprior")]) -
                           unlist(old.eb[c("mprior", "vprior")]))/
                       unlist(old.eb[c("mprior", "vprior")]) <= 
                       control$epsilon.eb)
      conv.opt <- c(conv.opt, new.eb$conv)
      check.eb <- conv.eb | (iter.eb >= control$maxit.eb)
      old.eb <- new.eb
    } else {
      check.eb <- TRUE
    }
    
    
    # set convergence check parameters for VB
    conv.vb <- c(conv.vb, FALSE)
    check.vb <- FALSE
    iter.vb <- c(iter.vb, 0)
    
    # inner VB loop
    while(!check.vb) {
      
      # increase iteration number by one
      iter.vb[iter.eb] <- iter.vb[iter.eb] + 1
      
      # update the VB parameters and elbo
      new.vb <- lapply(c(1:D), function(d) {
        .single.vb.update.enig(old.vb$a[[d]], old.vb$b[[d]], old.eb$Calpha[[d]], 
                               old.eb$lambda[ifelse(mult.lambda, d, 1)], y[, d], 
                               x[[d]], ytx[[d]], yty[d], n, p[d], D)})
      new.vb <- setNames(lapply(names(new.vb[[1]]), function(var) {
        lapply(new.vb, "[[", var)}), names(new.vb[[1]]))
      
      # check convergence of the VB iterations
      conv.vb[iter.eb] <- all(abs(unlist(new.vb[c("delta", "zeta")]) - 
                                    unlist(old.vb[c("delta", "zeta")]))/
                                unlist(old.vb[c("delta", "zeta")]) <= 
                                control$epsilon.vb)
      
      # last vb iterations until convergence or not
      if(control$conv.post & check.eb) {
        check.vb <- conv.vb[iter.eb] | (iter.vb[iter.eb] >= control$maxit.post)
      } else {
        check.vb <- conv.vb[iter.eb] | (iter.vb[iter.eb] >= control$maxit.vb)
      }
      old.vb <- new.vb
    }
  }
  
  # computing mu and the diagonal of Sigma
  aux <- lapply(1:D, function(d) {.aux.var.enig(old.vb$a[[d]], old.vb$b[[d]], 
                                                y[, d], x[[d]], ytx[[d]])})
  mu <- lapply(aux, "[[", "mu")
  if(full.post) {
    Sigma <- lapply(1:D, function(d) {.Sigma.enig(old.vb$a[[d]], old.vb$b[[d]], 
                                                  x[[d]])})
  } else {
    Sigma <- lapply(aux, "[[", "dSigma")  
  }
  
  
  # computing posterior means and variances
  mpost <- list(beta=mu, gamma.sq=sapply(1:D, function(d) {
    old.eb$mprior[[d]]*sqrt(old.vb$delta[[d]]/old.eb$lambda[
      ifelse(mult.lambda, d, 1)])*ratio_besselK(sqrt(
        old.vb$delta[[d]]*old.eb$lambda[ifelse(mult.lambda, d, 1)])/
          old.eb$mprior[[d]], 1)}),
    sigma.sq=2*unlist(old.vb$zeta)/(n + p - 1))
  vpost <- list(beta=Sigma, gamma.sq=sapply(1:D, function(d) {
    old.eb$mprior[[d]]^2*
      old.vb$delta[[d]]/old.eb$lambda[ifelse(mult.lambda, d, 1)]*
      (1 - ratio_besselK(sqrt(old.vb$delta[[d]]*old.eb$lambda[ifelse(
        mult.lambda, d, 1)])/old.eb$mprior[[d]], 1)^2)}),
    sigma.sq=8*unlist(old.vb$zeta)^2/((n + p - 1)^2*(n + p - 3)))
  names(seq.eb) <- c("alpha", "lambda", "mprior", "vprior")
  
  # creating output object
  out <- list(call=cl,
              vb=list(mu=mu, Sigma=Sigma, delta=old.vb$delta, 
                      zeta=unlist(old.vb$zeta), mpost=mpost, vpost=vpost), 
              eb=old.eb[!(names(old.eb) %in% c("Calpha", "conv", "iter"))],
              seq.eb=seq.eb,
              conv=list(eb=conv.eb, vb=conv.vb, opt=conv.opt),
              iter=list(eb=iter.eb, vb=iter.vb, opt=iter.opt))
  return(out)
}

# fit model without tissue effect and molecular feature groups (not tested)
est.model <- function(x, y, C, hyperprior=c("inv. Gaussian", "inv. Gamma"), 
                      conjugate, init=NULL,
                      control=list(epsilon.eb=1e-3, epsilon.vb=1e-3, 
                                   epsilon.opt=sqrt(.Machine$double.eps),
                                   maxit.eb=20, maxit.vb=2, maxit.opt=100,
                                   trace=TRUE)) {
  # hyperprior indicates inv. Gaussian model (TRUE) or inv. Gamma (FALSE)
  # conjugate indicates independent beta and sigma^2 (FALSE) or dependent (TRUE)
  # init is NULL, list(a, b, theta, lambda) if hyperprior="inv. Gaussian",
  # or list(a, b, eta, lambda) if hyperprior="inv. Gamma"
  # control$epsilon.vb=Inf is equivalent to control$maxit.vb=1
  
  # save the arguments
  cl <- match.call()
  
  # fixed parameters
  n <- nrow(x)
  p <- ncol(x)
  D <- ncol(y)
  if(hyperprior=="inv. Gaussian") {
    dcov <- ncol(C)
    init$eta <- rep(1, D)
    init$alpha <- c(1, rep(0, dcov - 1))
  } else if(hyperprior=="inv. Gamma") {
    dcov <- length(unique(C))
    sclass <- rle(sort(C))$lengths
    init$theta <- rep(Inf, D)
    init$alpha <- rep(0, dcov)
  } else {
    stop('hyperprior is either "inv. Gaussian" or "inv. Gamma"')
  }
  svd.x <- svd(x)
  sv <- svd.x$d
  uty <- t(svd.x$u) %*% y
  yty <- colSums(y^2)
  
  # create initial values if none given
  if(is.null(init)) {
    if(conjugate & hyperprior=="inv. Gaussian") {
      init <- .init.param.inv.gauss(C, sv, uty, D, n, p)
    } else if(!conjugate & hyperprior=="inv. Gaussian") {
      init <- .init.param.inv.gauss.ind(C, sv, uty, D, n, p)
    } else {
      stop("currently only initial values for inv. Gaussian supported")
    }
  }
  mult.lambda <- length(init$lambda) > 1
  if(!mult.lambda) {init$lambda <- rep(init$lambda, D)}

  # set initial values for VB, EB, and ELBO estimates
  old.eb <- list(alpha=init$alpha, eta=init$eta, theta=init$theta, 
                 lambda=init$lambda, elbo=NA, conv=FALSE)
  old.vb <- sapply(c(1:D), function(d) {
    .single.vb.update(init$a[d], init$b[d], old.eb$eta[d], old.eb$theta[d],
                      old.eb$lambda[d], sv, n, p, uty[, d], yty[d], conjugate,
                      hyperprior)})
  old.elbo <- old.vb[6, ]
  
  # prepare objects to store EB and ELBO iterations
  seq.eb <- list(alpha=numeric(0), eta=numeric(0), theta=numeric(0), 
                 lambda=numeric(0))
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
    if(hyperprior=="inv. Gaussian") {
      new.eb <- .eb.update.inv.gauss(old.vb[5, ], old.vb[4, ], old.vb[7, ], C, 
                                     D, mult.lambda, control$epsilon.opt, 
                                     control$maxit.opt)
    } else {
      esum <- sapply(unique(C), function(g) {sum(old.vb[5, C==g])})
      bsum <- sapply(unique(C), function(g) {sum(old.vb[4, C==g])})
      new.eb <- .eb.update.inv.gamma(esum, bsum, old.vb[4, ], old.vb[1, ], 
                                     old.eb$eta, old.vb[7, ], dcov, sclass, p, 
                                     D, control$epsilon.opt, control$maxit.opt)  
    }
    if(!mult.lambda) {new.eb$lambda <- rep(new.eb$lambda, D)}
    new.elbo <- new.eb$elbo
    seq.eb <- sapply(c(1:length(seq.eb)), function(s) {
      rbind(seq.eb[[s]], new.eb[[s]])})
    seq.elbo <- rbind(seq.elbo, new.elbo)
    
    # check convergence of EB estimates
    if(hyperprior=="inv. Gaussian") {
      conv.eb <- all(unlist(sapply(c(2:4), function(s) {
        abs(new.eb[[s]] - old.eb[[s]])})) <= control$epsilon.eb)
    } else {
      conv.eb <- all(unlist(sapply(c(1, 4), function(s) {
        abs(new.eb[[s]] - old.eb[[s]])})) <= control$epsilon.eb)
    }
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
      new.vb <- sapply(c(1:D), function(d) {
        .single.vb.update(old.vb[3, d], old.vb[4, d],  old.eb$eta[d], 
                          old.eb$theta[d], old.eb$lambda[d], sv, n, p, 
                          uty[, d], yty[d], conjugate, hyperprior)})
        
      new.elbo <- new.vb[6, ]
      seq.elbo <- rbind(seq.elbo, new.elbo)

      # check convergence of the VB iterations
      conv.vb[iter.eb] <- all(abs(new.elbo - old.elbo) <= control$epsilon.vb)
      check.vb <- conv.vb[iter.eb] | (iter.vb[iter.eb] >= control$maxit.vb)
      old.elbo <- new.elbo
      old.vb <- new.vb
    }
  }
  names(seq.eb) <- names(old.eb)[-c(5:6)]
  
  # creating Sigma and mu 
  cold <- new.vb[4, ]/new.vb[3, ]^(conjugate - 1)
  mu <- sapply(c(1:D), function(d) {
    svd.x$v %*% (t(svd.x$u)*sv/(sv^2 + cold[d])) %*% y[, d]})
  if(p > n) {
    svd.x <- svd(x, nu=nrow(x), nv=ncol(x))
  }
  dSigma <- sapply(c(1:D), function(d) {
    colSums(t(svd.x$v)*(t(svd.x$v)/(c(sv^2, rep(0, max(p - n, 0))) + cold[d])))/
      new.vb[3, d]})
  
  # preparing the output
  out <- list(call=cl,
              vb.post=list(mu=mu, dSigma=dSigma, 
                           delta=new.vb[1, ], zeta=new.vb[2, ]), 
              seq.eb=seq.eb,
              seq.elbo=seq.elbo,
              conv=list(eb=conv.eb, vb=conv.vb),
              iter=list(eb=iter.eb, vb=iter.vb))
  return(out)
}
