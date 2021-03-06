# estimate NIG model
nig <- function(x, y, C, Z, unpenalized=NULL, standardize=FALSE, 
                intercept=FALSE, fixed.eb="none", var.scale=1,
                full.post=FALSE, init=NULL,
                control=list(conv.post=TRUE, trace=TRUE,
                             epsilon.eb=1e-3, epsilon.vb=1e-3,
                             maxit.eb=100, maxit.vb=2, maxit.post=100,
                             maxit.block=0)) {
  
  # save the arguments
  cl <- match.call()
  
  # check for unpenalized covariates
  unp <- intercept | !is.null(unpenalized)
  
  # fixed parameters
  D <- ifelse(is.matrix(y), ncol(y), length(y))
  n <- sapply(y, length)
  if(is.null(unpenalized)) {
    u <- rep(intercept, D)
  } else {
    u <- sapply(unpenalized, function(s) {ncol(s) + intercept})
  }
  r <- sapply(x, ncol)
  p <- sapply(1:D, function(d) {u[d] + r[d]})
  ncd <- ifelse(is.null(C), 0, ncol(C[[1]]))
  nzd <- ifelse(is.null(Z), 0, ncol(Z))
  if(length(var.scale)==1) {
    var.scale <- rep(var.scale, D)
  }
  
  # standardize data and add intercept if need be
  xr <- x
  xu <- unpenalized
  if(standardize) {
    xr <- lapply(xr, scale)
  }
  if(intercept) {
    if(is.null(xu)) {
      xu <- lapply(1:D, function(d) {matrix(1, nrow=n[d], ncol=1)})
    } else {
      xu <- lapply(xu, function(s) {cbind(1, s)})  
    }
  }
  
  # transform co-data to matrix
  Cmat <- Reduce("rbind", C)
  
  # fixed data functions
  yty <- sapply(y, function(s) {sum(s^2)})
  if(!unp) {
    ytx <- lapply(1:D, function(d) {colSums(y[[d]]*x[[d]])})  
  }
  
  # set initial values if not given
  if(is.null(init$aold)) {
    init$aold <- rep(1, D)
  } 
  if(is.null(init$bold)) {
    init$bold <- lapply(1:D, function(d) {rep(1, r[d])})
  }
  if(is.null(init$gold)) {
      init$gold <- rep(1, D)
  }  
  if(is.null(init$lambdaf)) {
    init$lambdaf <- 1
  }
  if(is.null(init$lambdad)) {
    init$lambdad <- 1
  }
  if(is.null(init$fold)) {
    init$fold <- rep(1, D)
  }
  if(is.null(init$alphaf)) {
    init$alphaf <- switch(is.null(C[[1]]) + 1, c(1, rep(0, ncd - 1)), NULL)
  }
  init$Calphaf <- lapply(1:D, function(d) {
    switch(is.null(C[[d]]) + 1, as.numeric(C[[d]] %*% init$alphaf), NULL)})
  if(is.null(init$alphad)) {
    init$alphad <- switch(is.null(Z) + 1, c(1, rep(0, nzd - 1)), NULL)  
  }
  init$Zalphad <- switch(is.null(Z) + 1, as.numeric(Z %*% init$alphad), NULL)
  init$mpriorf <- lapply(1:D, function(d) {1/(init$Calphaf[[d]])})
  init$vpriorf <- lapply(1:D, function(d) {
    1/(init$Calphaf[[d]]^3*init$lambdaf)})  
  init$mpriord <- 1/init$Zalphad
  init$vpriord <- 1/(init$Zalphad^3*init$lambdad)
  
  # set initial values for VB, EB, and ELBO estimates
  old.eb <- list(alphaf=init$alphaf, lambdaf=init$lambdaf, mpriorf=init$mpriorf,
                 vpriorf=init$vpriorf, Calphaf=init$Calphaf, alphad=init$alphad, 
                 lambdad=init$lambdad, mpriord=init$mpriord, 
                 vpriord=init$vpriord, Zalphad=init$Zalphad)
  if(unp) {
    old.vb <- lapply(c(1:D), function(d) {
      .single.vb.update.unp(init$aold[d], init$bold[[d]], init$gold[[d]], 
                            init$Calphaf[[d]], init$Zalphad[[d]], init$lambdaf, 
                            init$lambdad, 
                            switch(is.list(y) + 1, y[, d], y[[d]]), xu[[d]], 
                            xr[[d]], yty[d], n[d], u[d], r[d], var.scale[d])})
  } else {
    old.vb <- lapply(c(1:D), function(d) {
      .single.vb.update(init$aold[d], init$bold[[d]], init$gold[[d]], 
                        init$Calphaf[[d]], init$Zalphad[[d]], init$lambdaf, 
                        init$lambdad, switch(is.list(y) + 1, y[, d], y[[d]]), 
                        x[[d]], ytx[[d]], yty[d], n[d], r[d], var.scale[d])})  
  }
  aux <- lapply(old.vb, function(s) {s$aux})
  old.vb <- setNames(lapply(names(old.vb[[1]])[
    names(old.vb[[1]])!="aux"], function(var) {lapply(old.vb, "[[", var)}), 
    names(old.vb[[1]])[names(old.vb[[1]])!="aux"])
  
  # calculate ELBO
  old.elbo <- sapply(1:D, function(d) {
    .single.elbo(p[d], n[d], old.vb$zeta[[d]], yty[[d]], aux[[d]], 
                 old.vb$g[[d]], old.vb$b[[d]], old.vb$delta[[d]], 
                 old.vb$eta[[d]], old.eb$lambdaf, old.eb$lambdad, 
                 old.eb$Zalphad[[d]], old.eb$Calphaf[[d]])})
  
  # prepare objects to store EB and ELBO iterations
  seq.eb <- old.eb
  seq.elbo <- old.elbo
  
  # set all convergence check parameters
  conv.eb <- FALSE
  conv.vb <- logical(0)
  conv.block <- setNames(c(FALSE, FALSE), c("feat", "drug"))
  iter.eb <- 0
  iter.vb <- numeric(0)
  iter.block <- 0
  check.eb <- FALSE
  check.vb <- FALSE
  if(control$maxit.block!=0) {
    block <- "feat"
    old.eb.block <- old.eb
  } else {
    block <- c("feat", "drug")
  }
  
  srt <- proc.time()[3]
  # outer EB loop
  while(!check.eb) {
    
    # increase iteration number by one
    iter.eb <- iter.eb + 1
    if(control$maxit.block!=0) {
      iter.block <- iter.block + 1  
    }
    
    # if EB is required
    if(fixed.eb!="all") {
      
      # possibly print iteration number
      if(control$trace) {
        cat("\r", "Estimating EB parameters, iteration ", iter.eb, ", alphaf: ", 
            paste(switch(
              is.null(old.eb$alphaf) + 1,
              round(old.eb$alphaf[1:min(5, length(old.eb$alphaf))], 2), NULL),
              collapse=" "), ", alphad: ", 
            paste(switch(
              is.null(old.eb$alphad) + 1, 
              round(old.eb$alphad[1:min(5, length(old.eb$alphad))], 2), NULL),
              collapse=" "), sep="")
      }
      
      # update the EB parameters
      if(!is.null(C[[1]]) & ("feat" %in% block)) {
        new.ebf <- .eb.updatef(old.vb$e, old.vb$b, Cmat, r, D,
                               switch((fixed.eb=="lambda") + 1, NULL, 
                                      init$lambdaf))   
      } else {
        new.ebf <- old.eb[c("alphaf", "lambdaf", "mpriorf", "vpriorf", 
                            "Calphaf")]
      }
      if(!is.null(Z) & ("drug" %in% block)) {
        new.ebd <- .eb.updated(unlist(old.vb$f), unlist(old.vb$g), Z, D,
                               switch((fixed.eb=="lambda") + 1, NULL, 
                                      init$lambdad))  
      } else {
        new.ebd <- old.eb[c("alphad", "lambdad", "mpriord", "vpriord", 
                            "Zalphad")]
      }
      new.eb <- c(new.ebf, new.ebd)
      if(fixed.eb=="lambda") {
        new.eb["lambdad"] <- old.eb["lambdad"]
        new.eb["lambdaf"] <- old.eb["lambdaf"]
      }
      
      # paste new hyperparameters to previous
      seq.eb <- sapply(names(seq.eb), function(s) {
        if(!is.list(seq.eb[[s]])) {
          rbind(seq.eb[[s]], new.eb[[s]])
        } else {
          lapply(1:D, function(d) {rbind(seq.eb[[s]][[d]], new.eb[[s]][[d]])})
        }}, USE.NAMES=TRUE, simplify=FALSE)
      
      # check convergence of EB estimates
      if(control$maxit.block!=0 & (iter.block < control$maxit.block)) {
        conv.eb <- FALSE
      } else if(control$maxit.block!=0 & iter.block==control$maxit.block) {
        cur.par <- paste0(c("mprior", "vprior"), substr(block, 1, 1))
        conv.block[block] <- all(abs(unlist(new.eb[cur.par]) - 
                                       unlist(old.eb.block[cur.par])) <=
          control$epsilon.eb)
        conv.eb <- all(conv.block)
        old.eb.block <- new.eb
        iter.block <- 0
        block <- c("feat", "drug")[c("feat", "drug")!=block]
      } else {
        conv.eb <- all(abs(unlist(new.eb[
          c("mpriorf", "vpriorf", "mpriord", "vpriord")]) - 
            unlist(old.eb[c("mpriorf", "vpriorf", "mpriord", "vpriord")])) <= 
            control$epsilon.eb)    
      }
      
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
          .single.vb.update.unp(
            old.vb$a[[d]], ifelse(is.null(C), rep(1, r[d]), old.vb$b[[d]]), 
            old.vb$g[[d]], old.eb$Calphaf[[d]], old.eb$Zalphad[[d]],
            old.eb$lambdaf, old.eb$lambdad, 
            switch(is.list(y) + 1, y[, d], y[[d]]), xu[[d]], xr[[d]], yty[d], 
            n[d], u[d], r[d], var.scale[d])})  
      } else {
        new.vb <- lapply(c(1:D), function(d) {
          .single.vb.update(old.vb$a[[d]], old.vb$b[[d]], old.vb$g[[d]],
                            old.eb$Calphaf[[d]], old.eb$Zalphad[[d]],
                            old.eb$lambdaf, old.eb$lambdad, 
                            switch(is.list(y) + 1, y[, d], y[[d]]), x[[d]], 
                            ytx[[d]], yty[d], n[d], r[d], var.scale[d])})  
      }
      aux <- lapply(new.vb, function(s) {s$aux})
      new.vb <- setNames(lapply(names(new.vb[[1]])[
        names(new.vb[[1]])!="aux"], function(var) {
          lapply(new.vb, "[[", var)}), 
        names(new.vb[[1]])[names(new.vb[[1]])!="aux"])
      
      # check convergence of the VB iterations
      conv.vb[iter.eb] <- all(abs(unlist(new.vb[c("delta", "zeta", "eta")]) - 
                                    unlist(old.vb[c("delta", "zeta", "eta")]))/
                                unlist(old.vb[c("delta", "zeta", "eta")]) <= 
                                control$epsilon.vb)
      
      # save times
      if(check.eb & iter.vb[iter.eb]==1) {
        eb.time <- proc.time()[3] - srt
        srt <- proc.time()[3]
        if(control$trace) {
          cat("\n", "Estimated EB parameters in ", round(eb.time), " seconds", 
              sep="")  
          if(control$conv.post) {
            cat("\n", "Estimating VB parameters", sep="")  
          }
        }
      }
      
      # last vb iterations until convergence or not
      if(control$conv.post & check.eb) {
        check.vb <- conv.vb[iter.eb] | (iter.vb[iter.eb] >= control$maxit.post)
      } else {
        check.vb <- conv.vb[iter.eb] | (iter.vb[iter.eb] >= control$maxit.vb)
      }
      old.vb <- new.vb
    }
   
    # calculate ELBO
    new.elbo <- sapply(1:D, function(d) {
      .single.elbo(p[d], n[d], new.vb$zeta[[d]], yty[[d]], aux[[d]],
                   new.vb$g[[d]], new.vb$b[[d]], new.vb$delta[[d]],
                   new.vb$eta[[d]], new.eb$lambdaf, new.eb$lambdad,
                   new.eb$Zalphad[[d]], new.eb$Calphaf[[d]])})
    seq.elbo <- rbind(seq.elbo, new.elbo)
    old.elbo <- new.elbo
    
  }
  
  # calculate posterior
  mu <- lapply(aux, "[[", "mu")
  if(full.post) {
    if(unp) {
      Sigma <- lapply(1:D, function(d) {
          .Sigma.unp(old.vb$a[[d]], old.vb$b[[d]]*old.vb$g[[d]]/var.scale[d], 
                     xu[[d]], xr[[d]], u[d], r[d])}) 
    } else {
      Sigma <- lapply(1:D, function(d) {
        .Sigma(old.vb$a[[d]], old.vb$b[[d]]*old.vb$g[[d]]/var.scale[d], 
               x[[d]])})  
    }
  } else {
    Sigma <- lapply(aux, "[[", "dSigma")  
  }
  
  # calculating VB time
  if(control$conv.post) {
    post.time <- proc.time()[3] - srt
    if(control$trace) {
      cat("\n", "Estimated VB parameters in ", round(post.time), " seconds", 
          sep="")
    }
  } else {
    post.time <- NULL
  }
  
  # create posterior mean and variance objects
  mpost <- list(beta=mu, gamma.sq=lapply(1:D, function(d) {
    old.eb$mpriorf[[d]]*sqrt(old.vb$delta[[d]]/old.eb$lambdaf)*
      ratio_besselK(sqrt(old.vb$delta[[d]]*old.eb$lambdaf)/
                      old.eb$mpriorf[[d]], 1)}),
    sigma.sq=2*unlist(old.vb$zeta)/(n + p - 1),
    tau.sq=sapply(1:D, function(d) {
      old.eb$mpriord[d]*sqrt(old.vb$eta[[d]]/old.eb$lambdad)*
        ratio_besselK(sqrt(old.vb$eta[[d]]*old.eb$lambdad)/
                        old.eb$mpriord[d], (r[d] + 1)/2)}))
  vpost <- list(beta=Sigma, gamma.sq=lapply(1:D, function(d) {
    old.eb$mpriorf[[d]]^2*old.vb$delta[[d]]/old.eb$lambdaf*
      (1 - ratio_besselK(sqrt(old.vb$delta[[d]]*old.eb$lambdaf)/
                           old.eb$mpriorf[[d]], 1)^2)}),
    sigma.sq=8*unlist(old.vb$zeta)^2/((n + p - 1)^2*(n + p - 3)),
    tau.sq=sapply(1:D, function(d) {old.eb$mpriord[d]^2*
        old.vb$eta[[d]]/old.eb$lambdad*
        (1 - ratio_besselK(sqrt(old.vb$eta[[d]]*old.eb$lambdad)/
                             old.eb$mpriord[d], (r[d] + 1)/2)^2)}))
  
  # creating and returning output object
  out <- list(call=cl,
              vb=list(mu=mu, Sigma=Sigma, delta=old.vb$delta, 
                      eta=unlist(old.vb$eta),
                      zeta=unlist(old.vb$zeta), mpost=mpost, vpost=vpost), 
              eb=old.eb[names(old.eb)!="Calpha"], 
              seq.eb=seq.eb,
              seq.elbo=seq.elbo,
              conv=list(eb=conv.eb, vb=conv.vb),
              iter=list(eb=iter.eb, vb=iter.vb),
              time=c(eb.time=unname(eb.time), post.time=unname(post.time)))
  class(out) <- c("nig")
  return(out)
  
}

# estimate ELBO on new data
new.elbo <- function(object, newx, newy) {
  D <- ncol(newy)
  p <- sapply(newx, ncol)
  n <- nrow(newy)
  out <- sapply(1:D, function(d) {
    aux <- list(ldetSigma=determinant(object$vb$Sigma[[d]])$modulus, 
                ytXmu=as.numeric(t(newy[, d]) %*% newx[[d]] %*% 
                                  object$vb$mu[[d]]), 
                trXtXSigma=sum(diag(t(newx[[d]]) %*% newx[[d]] %*% 
                                      object$vb$Sigma[[d]])), 
                mutXtXmu=sum((t(object$vb$mu[[d]]) %*% t(newx[[d]]))^2), 
                dSigma=sum(diag(object$vb$Sigma[[d]])), 
                mu=object$vb$mu[[d]])
    g <- object$vb$mpost$tau.sq[[d]]*object$eb$lambdad*
      object$eb$Zalphad[[d]]^2/object$vb$eta[[d]] + (p[d] + 1)/
      object$vb$eta[[d]]
    g <- switch(as.numeric(length(g)==0) + 1, g, 1)
    b <- object$vb$mpost$gamma.sq[[d]]*object$eb$lambdaf*
      object$eb$Calphaf[[d]]^2/object$vb$delta[[d]] + 2/object$vb$delta[[d]]
    b <- switch(as.numeric(length(b)==0) + 1, b, 1)
    .single.elbo(p[d], n, object$vb$zeta[[d]], sum(newy[, d]^2),
                 aux, g, b, object$vb$delta[[d]], 
                 object$vb$eta[[d]], object$eb$lambdaf,
                 object$eb$lambdad, object$eb$Zalphad[[d]], 
                 object$eb$Calphaf[[d]])}) 
  return(out)
}

# calculate log cpo on test data
logcpo <- function(xtest, ytest, ntrain, fit) {
  D <- length(ytest)
  ntest <- sapply(ytest, length)
  out <- lapply(ntest, numeric)
  for(d in 1:D) {
    xtSigmax <- rowSums((xtest[[d]] %*% fit$vb$Sigma[[d]])*xtest[[d]])
    xtmu <- as.numeric(xtest[[d]] %*% fit$vb$mu[[d]])
    p <- length(fit$vb$mu[[d]])
    zeta <- fit$vb$zeta[d]
    for(i in 1:ntest[d]) {
      int <- integrate(.f.int.cpo, -Inf, Inf, xtSigmax=xtSigmax[i], n=ntrain[d], 
                       p=p, zeta=zeta, y=ytest[[d]][i], xtmu=xtmu[i],
                       stop.on.error=FALSE)
      int.val <- ifelse(int$message=="the integral is probably divergent", NA, 
                        int$value) 
      out[[d]][i] <- -log(2) - log(pi) - 0.5*log(xtSigmax[i]) - 0.5*log(zeta) +
        lgamma((ntrain[d] + p + 2)/2) - lgamma((ntrain[d] + p + 1)/2) + 
        log(int.val)
    }
  }
  return(out)
}
