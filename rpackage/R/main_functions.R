# x <- list(x, x); y <- cbind(y, y); C=NULL; unpenalized=NULL;
# standardize=FALSE; intercept=FALSE; intercept.eb=TRUE;
# mult.lambda=FALSE; fixed.eb=c("both"); full.post=TRUE;
# init=list(aold=rep(1, 2), bold=list(p, p),
#           lambda=lambda, alpha=1,
#           Calpha=rep(1/ctalphainv, 2));
# control=list(conv.post=TRUE, trace=TRUE,
#              epsilon.eb=1e-3, epsilon.vb=1e-3,
#              epsilon.opt=sqrt(.Machine$double.eps),
#              maxit.eb=10, maxit.vb=2, maxit.opt=100,
#              maxit.post=100)
# detach("package:cambridge", unload=TRUE)
# library(cambridge)

# estimate SEMNIG model
semnig <- function(x, y, C, unpenalized=NULL, standardize=TRUE, intercept=TRUE, 
                   intercept.eb=TRUE, mult.lambda=FALSE, 
                   fixed.eb=c("none", "lambda", "both"), 
                   full.post=FALSE, init=NULL,
                   control=list(conv.post=TRUE, trace=TRUE,
                                epsilon.eb=1e-3, epsilon.vb=1e-3, 
                                epsilon.opt=sqrt(.Machine$double.eps),
                                maxit.eb=10, maxit.vb=2, maxit.opt=100,
                                maxit.post=100)) {
  
  # save the arguments
  cl <- match.call()
  
  # check for unpenalized covariates
  unp <- intercept | !is.null(unpenalized)
  
  # fixed parameters
  D <- ncol(y)
  n <- nrow(y)
  if(is.null(unpenalized)) {
    u <- rep(0, D)
  } else {
    u <- sapply(unpenalized, function(s) {ncol(s) + intercept})
  }
  r <- sapply(x, ncol)
  p <- sapply(1:D, function(d) {u[d] + r[d]})
  ncd <- ifelse(is.null(C), 0, ncol(C[[1]]))
  
  # standardize data and add intercept if need be
  xr <- x
  xu <- unpenalized
  if(standardize) {
    xr <- lapply(xr, scale)
  }
  if(intercept) {
    if(is.null(xu)) {
      xu <- rep(list(matrix(1, nrow=1, ncol=1)), D)
    } else {
      xu <- lapply(xu, function(s) {cbind(1, s)})  
    }
  }
  
  # add intercept and transform co-data to matrix
  if(intercept.eb & is.null(C)) {
    C <- lapply(1:D, function(d) {
      m <- matrix(1, nrow=r[d]); colnames(m) <- "intercept"; return(m)})
  } else if(intercept.eb) {
    C <- lapply(C, function(cc) {cbind(intercept=1, cc)})
  }
  Cmat <- Reduce("rbind", C)
  
  # fixed data functions
  yty <- colSums(y^2)
  if(!unp) {
    ytx <- lapply(1:D, function(d) {as.numeric(y[, d] %*% x[[d]])})  
  }
  
  # set initial values if not given
  if(is.null(init)) {
    init$aold <- rep(1, D)
    init$bold <- lapply(1:D, function(d) {rep(1, r[d])})
    init$lambda <- rep(1, ifelse(mult.lambda, D, 1))
    if(intercept.eb) {
      init$alpha <- c(1, rep(0, ncd))
      init$Calpha <- lapply(1:D, function(d) {
        as.numeric(C[[d]] %*% init$alpha)})
    } else {
      init$alpha <- rep(0, ncd)
      init$Calpha <- lapply(1:D, function(d) {rep(1, r[d])})
    }
  }
  init$mprior <- lapply(1:D, function(d) {1/(init$Calpha[[d]])})
  init$vprior <- lapply(1:D, function(d) {
    1/(init$Calpha[[d]]*init$lambda[ifelse(mult.lambda, d, 1)])})
  
  # set initial values for VB, EB, and ELBO estimates
  old.eb <- list(alpha=init$alpha, lambda=init$lambda, mprior=init$mprior,
                 vprior=init$vprior, Calpha=init$Calpha, conv=FALSE, iter=0)
  if(unp) {
    old.vb <- lapply(c(1:D), function(d) {
      .single.vb.update.unp(init$aold[d], init$bold[[d]], init$Calpha[[d]], 
                            init$lambda[ifelse(mult.lambda, d, 1)], y[, d], 
                            xu[[d]], xr[[d]], yty[d], n, u[d], r[d])})  
  } else {
    old.vb <- lapply(c(1:D), function(d) {
      .single.vb.update(init$aold[d], init$bold[[d]], init$Calpha[[d]], 
                        init$lambda[ifelse(mult.lambda, d, 1)], y[, d], 
                        x[[d]], ytx[[d]], yty[d], n, r[d])})  
  }
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
      new.eb <- .eb.update(old.vb$e, old.vb$b, old.eb$lambda, Cmat, r, D, 
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
      if(unp) {
        new.vb <- lapply(c(1:D), function(d) {
          .single.vb.update.unp(old.vb$a[[d]], old.vb$b[[d]], 
                                old.eb$Calpha[[d]], 
                                old.eb$lambda[ifelse(mult.lambda, d, 1)], 
                                y[, d], xu[[d]], xr[[d]], yty[d], n, u[d], 
                                r[d])})
      } else {
        new.vb <- lapply(c(1:D), function(d) {
          .single.vb.update(old.vb$a[[d]], old.vb$b[[d]], old.eb$Calpha[[d]], 
                            old.eb$lambda[ifelse(mult.lambda, d, 1)], y[, d], 
                            x[[d]], ytx[[d]], yty[d], n, r[d])})  
      }
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
  if(unp) {
    aux <- lapply(1:D, function(d) {
      .aux.var.unp(old.vb$a[[d]], old.vb$b[[d]], y[, d], xu[[d]], xr[[d]], u[d], 
                   r[d])})
  } else {
    aux <- lapply(1:D, function(d) {
      .aux.var(old.vb$a[[d]], old.vb$b[[d]], y[, d], x[[d]], ytx[[d]])})
  }
  
  mu <- lapply(aux, "[[", "mu")
  if(full.post) {
    if(unp) {
      Sigma <- lapply(1:D, function(d) {
        .Sigma.unp(old.vb$a[[d]], old.vb$b[[d]], xu[[d]], xr[[d]], u[d], r[d])}) 
    } else {
      Sigma <- lapply(1:D, function(d) {
        .Sigma(old.vb$a[[d]], old.vb$b[[d]], x[[d]])})  
    }
    
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
  class(out) <- c("semnig")
  return(out)
}
