# x=rep(list(xtrain), D); y=ytrain
# unpenalized=NULL; standardize=FALSE; 
# intercept=FALSE; fixed.eb="none"; var.scale=1;
# full.post=FALSE; init=NULL;
# library(Rcpp)
# sourceCpp("rpackage/src/aux_functions.cpp")
# source("rpackage/R/vb_functions.R")
# source("rpackage/R/aux_functions.R")
# source("rpackage/R/eb_functions.R")
# detach("package:cambridge", unload=TRUE)
# library(cambridge)

# estimate NIG model (not tested)
semnig <- function(x, y, C, Z, unpenalized=NULL, standardize=FALSE, 
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
  class(out) <- c("semnig")
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

# calculate log cpo on test data (not tested)
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

# cross-validate the hyperparameters of the MAP estimator
cv.semnig <- function(x, y, nfolds=10, foldid=NULL, seed=NULL, phi=phi,
                      chi=chi, lambdaf=lambdaf, lambdad=lambdad,
                      type.measure="mse", control=list(trace=FALSE)) {
  
  D <- ifelse(is.matrix(y), ncol(y), 1)
  fit <- lapply(1:D, function(d) {
    if(control$trace) {cat("\r", "drug ", d)}
    .cv.single.semnig(x[[d]], y[, d], nfolds=nfolds, foldid=foldid, seed=seed,
                      phi=phi, chi=chi, lambdaf=lambdaf, lambdad=lambdad, 
                      type.measure=type.measure, control=control)})
  
  return(fit)
}


# x=xtrain; y=ytrain; C=lapply(C, function(s) {matrix(s[, -1])}); Z=NULL;
# mult.lambda=TRUE; fix.lambda=TRUE;
# nfolds=10; foldid=foldid;
# hyper=list(lambda=1/sqrt(n*sapply(fit.ridge1, "[[", "lambda.min")), 
#            zeta=0, nu=0);
# control=control.ebridge
# library(Rcpp)
# sourceCpp("rpackage/src/aux_functions.cpp")
# EBridge estimation
ebridge <- function(x, y, C, Z, mult.lambda=TRUE, fix.lambda=TRUE, 
                    nfolds=10, foldid=NULL,
                    hyper=list(lambda=NULL, zeta=0, nu=0),
                    control=list(epsilon=sqrt(.Machine$double.eps), 
                                 maxit=500, trace=FALSE, glmnet.fit2=FALSE,
                                 beta2=FALSE)) {
  
  yty <- sapply(y, function(s) {sum(s^2)})
  Cmat <- Reduce("rbind", C)
  n <- sapply(y, length)
  D <- length(y)
  if(is.matrix(x)) {
    p <- ncol(x)
    idsel <- lapply(y, function(s) {
      match(names(s), rownames(x))})
  } else {
    p <- sapply(x, ncol)
  }
  H <- ifelse(is.null(Z), 0, ncol(Z))
  G <- ncol(Cmat)
  
  if(is.null(foldid)) {
    foldid <- lapply(1:D, function(d) {
      sample(c(rep(1:nfolds, each=n[[d]] %/% nfolds),
               rep(1:(n[[d]] %% nfolds), (n[[d]] %% nfolds)!=0)))})
  }
  
  if(is.null(hyper$lambda) & fix.lambda) {
    srt <- proc.time()[3]
    if(control$trace) {
      cat("\r", "Cross-validating penalty parameters")
    }
    if(is.matrix(x)) {
      cv.fit1 <- lapply(1:D, function(d) {
        cv.glmnet(x[idsel[[d]], ], y[[d]], alpha=0, intercept=FALSE, 
                  foldid=foldid[[d]])})  
    } else {
      cv.fit1 <- lapply(1:D, function(d) {
        cv.glmnet(x[[d]], y[[d]], alpha=0, intercept=FALSE, 
                  foldid=foldid[[d]])})  
    }
    if(mult.lambda) {
      hyper$lambda <- 1/sqrt(n*sapply(cv.fit1, "[[", "lambda.min"))
    } else {
      hyper$lambda <- 1/sqrt(
        exp(mean(log(n) + log(sapply(cv.fit1, "[[", "lambda.min")))))
    }
    cv.time <- proc.time()[3] - srt
  } else if(!is.null(hyper$lambda)) {
    if(length(hyper$lambda)==1) {
      hyper$lambda <- rep(hyper$lambda, D)  
    }
    cv.time <- 0
  } else {
    if(mult.lambda) {
      Cmat <- cbind(model.matrix(
        ~ as.character(rep(c(1:D), times=rep(p, ifelse(length(p)==1, D, 1))))), 
        Cmat)
    } else {
      Cmat <- cbind(1, Cmat)    
    }
    G <- ncol(Cmat)
    cv.time <- 0
    hyper$lambda <- rep(1, D)
  }
  if(control$glmnet.fit2 & fix.lambda) {
    if(is.matrix(x)) {
      cv.fit2 <- lapply(1:D, function(d) {
        glmnet(x[idsel[[d]], ], y[[d]], alpha=0, 
               lambda=1/(exp(2*mean(log(hyper$lambda)))*n[d]), 
               intercept=FALSE)})
    } else {
      cv.fit2 <- lapply(1:D, function(d) {
        glmnet(x[[d]], y[[d]], alpha=0, 
               lambda=1/(exp(2*mean(log(hyper$lambda)))*n[d]), 
               intercept=FALSE)})
    }
  }
  
  method <- ifelse((H + G)==1, "Brent", "Nelder-Mead")
  lower <- ifelse(method=="Brent", -20, -Inf)
  upper <- ifelse(method=="Brent", 20, Inf)
  control$fnscale=-1
  names(control)[1] <- "reltol"
  srt <- proc.time()[3]
  if(control$trace) {
    cat("\r", "Estimating EB parameters")
  }
  
  if(is.matrix(x)) {
    opt <- optim(par=rep(0, H + G), fn=.f.optim.mat, lambda=hyper$lambda, 
                 nu=hyper$nu, zeta=hyper$zeta, Cmat=Cmat, 
                 Z=switch(as.numeric(is.null(Z)) + 1, Z, matrix(0)), 
                 n=n, p=p, D=D, idsel=idsel, G=G, H=H, y=y, x=x, yty=yty, 
                 Zpres=!is.null(Z), method=method, lower=lower, 
                 upper=upper,
                 control=control[!(names(control) %in% 
                                     c("glmnet.fit2", "beta2"))])
  } else {
    opt <- optim(par=rep(0, H + G), fn=.f.optim.list, lambda=hyper$lambda, 
                 nu=hyper$nu, zeta=hyper$zeta, Cmat=Cmat, 
                 Z=switch(as.numeric(is.null(Z)) + 1, Z, matrix(0)), 
                 n=n, p=p, D=D, G=G, H=H, y=y, x=x, yty=yty, Zpres=!is.null(Z),
                 method=method, lower=lower, 
                 upper=upper,
                 control=control[!(names(control) %in% 
                                     c("glmnet.fit2", "beta2"))])
  }
  eb.time <- proc.time()[3] - srt
  conv <- opt$convergence
  niter <- unname(opt$counts[1])
  alphaf <- opt$par[1:G]
  if(fix.lambda) {
    gamma <- lapply(C, function(s) {exp(colSums(alphaf*t(s))/2)})  
  } else {
    if(mult.lambda) {
      hyper$lambda <- exp((alphaf[1] + c(0, alphaf[2:D]))/2)
      gamma <- lapply(C, function(s) {exp(colSums(alphaf[-c(1:D)]*t(s))/2)})
    } else {
      hyper$lambda <- exp(rep(alphaf[1], D)/2)
      gamma <- lapply(C, function(s) {exp(colSums(alphaf[-1]*t(s))/2)})
    }
  }
  if(is.null(Z)) {
    alphad <- NULL
    tau <- rep(1, D)
  } else {
    alphad <- opt$par[-c(1:G)]
    tau <- exp(colSums(alphad*t(Z))/2)
  }
  
  # create estimates
  beta1 <- lapply(1:D, function(d) {
    h <- tau[d]*gamma[[d]]
    if(is.matrix(x)) {
      xh <- t(t(x[idsel[[d]], ])*h)
    } else {
      xh <- t(t(x[[d]])*h)
    }
    fit <- glmnet(xh, y[[d]], alpha=0, lambda=1/(hyper$lambda[d]^2*n[d]),
                  standardize=FALSE, intercept=FALSE)
    unname(coef(fit)[-1, 1])*h})
  
  if(control$beta2 & fix.lambda) {
    beta2 <- lapply(1:D, function(d) {
      h <- tau[d]*gamma[[d]]
      if(is.matrix(x)) {
        xh <- x[idsel[[d]], ]
      } else {
        xh <- x[[d]]
      }
      fit <- glmnet(xh, y[[d]], alpha=0, lambda=1/(hyper$lambda[d]^2*n[d]),
                    intercept=FALSE, penalty.factor=1/h^2)
      unname(coef(fit)[-1, 1])})
  } else if(control$beta2) {
    beta2 <- lapply(1:D, function(d) {
      h <- 1/(n[d]*hyper$lambda[d]^2*tau[d]*gamma[[d]])^2
      s <- mean(h)
      if(is.matrix(x)) {
        xh <- x[idsel[[d]], ]
      } else {
        xh <- x[[d]]
      }
      fit <- glmnet(xh, y[[d]], alpha=0, intercept=FALSE)
      unname(coef(fit, x=xh, y=y[[d]], exact=TRUE, s=s, penalty.factor=h, 
                  intercept=FALSE)[-1, 1])})
  } else {
    beta2 <- NULL
  }
  
  out <- list(beta1=beta1, beta2=beta2, alphaf=alphaf, alphad=alphad, 
              lambda=hyper$lambda, tau=tau, 
              gamma=gamma, time=list(cv.time=cv.time, eb.time=eb.time))
  
  if(exists("cv.fit1")) {
    out$glmnet.fit1 <- cv.fit1
  } else {
    out$glmnet.fit1 <- NULL
  }
  if(exists("cv.fit2")) {
    out$glmnet.fit2 <- cv.fit2  
  } else {
    out$glmnet.fit2 <- NULL
  }
  if(exists("beta2")) {
    out$beta2 <- beta2
  } else {
    out$beta2 <- NULL
  }
  
  return(out)
  
}

################################### ebridge2 ###################################
# theta update
.single.theta.update <- function(mat, x, y, n, lambda, nu, zeta) {
  theta <- lambda^2*colSums(x*(mat %*% x))/2
  phi <- (n/2 + nu)/(lambda^2*(2*zeta + sum(y*(mat %*% y))))
  return(c(phi, theta))
}

theta.update <- function(mat, x, y, n, p, D, lambda, nu, zeta, idsel,
                         fix.sigma) {
  theta <- lapply(1:D, function(d) {
    .single.theta.update(mat[[d]], x[idsel[[d]], ], y[[d]], n[d], lambda[d], nu,
                         zeta)})
  if(fix.sigma[1]!=FALSE) {
    out <- list(theta=unlist(lapply(theta, "[", -1)), phi=1/fix.sigma^2)
  } else {
    out <- list(theta=unlist(lapply(theta, "[", -1)), phi=sapply(theta, "[", 1))
  }
  return(out)
}

# alpha update
alpha.update <- function(delta, theta, phi, alpha.old, tCmat, p, D, lambda,
                         method, maxit) {
  phi <- rep(phi, times=p)
  deltatildesq <- delta^2*phi
  opt <- optim(alpha.old, .f.opt.alpha, tCmat=tCmat, theta=theta,
               deltatildesq=deltatildesq, method=method,
               control=list(maxit=maxit))
  alpha <- opt$par
  h <- exp(colSums(tCmat*alpha))*rep(lambda^2, times=p)
  return(list(alpha=alpha, h=h, conv=opt$convergence, iter=opt$counts))
}

.f.opt.alpha <- function(alpha, tCmat, theta, deltatildesq) {
  vec <- exp(colSums(tCmat*alpha))
  return(sum(theta*vec) + sum(deltatildesq/vec))
}

# update delta
.single.delta.update <- function(h, x, y) {
  fit <- glmnet(t(t(x)*sqrt(h)), y, alpha=0, lambda=1, intercept=FALSE, 
                standardize=FALSE)
  return(sqrt(h)*coef(fit, s=1, exact=TRUE)[-1, ])
  # temp <- t(x)*h
  # mat <- x %*% temp
  # diag(mat) <- diag(mat) + 1
  # mat <- solve(mat)
  # delta <- colSums(t(temp %*% mat)*y)
  # return(list(delta=delta, mat=mat))
}

delta.update <- function(h, x, y, idsel, D, p) {
  h <- split(h, rep(1:D, p))
  delta <- lapply(1:D, function(d) {
    .single.delta.update(h[[d]], x[idsel[[d]], ], y[[d]])})
  return(list(delta=unlist(lapply(delta, "[[", "delta")),
              mat=lapply(delta, "[[", "mat")))
}

# x=xtrain; y=ytrain; C=C; Z=Z; foldid=foldid;
# varsel=NULL
# hyper=list(lambda=fit1.ebridge$lambda, zeta=0, nu=0);
# control=list(epsilon=sqrt(.Machine$double.eps),
#              maxit=list(eb.outer=500, eb.inner=100,
#                         opt=500),
#              fix.sigma=FALSE, trace=TRUE,
#              opt.method="Nelder-Mead",
#              standardize.C=FALSE,
#              standardize.Z=FALSE)
# ebridge2 function
ebridge2 <- function(x, y, C, Z, varsel=NULL, nfolds=10, foldid=NULL,
                     hyper=list(lambda=NULL, zeta=0, nu=0),
                     control=list(epsilon=sqrt(.Machine$double.eps),
                                  maxit=list(eb.outer=500, eb.inner=100,
                                             opt=500),
                                  fix.sigma=FALSE, trace=TRUE,
                                  opt.method="Nelder-Mead",
                                  standardize.C=TRUE,
                                  standardize.Z=TRUE)) {
  
  # auxiliary variables
  D <- length(y)
  if(is.null(varsel)) {
    p <- rep(ncol(x), D)
  } else {
    p <- sapply(varsel, length)
  }
  n <- sapply(y, length)
  H <- ncol(Z)
  G <- ncol(C[[1]])
  idsel <- lapply(y, function(s) {
    match(names(s), rownames(x))})
  if(control$standardize.C) {
    C <- lapply(C, function(s) {t(t(s) - colMeans(s))})
  }
  if(control$standardize.Z) {
    Z <- t(t(Z) - colMeans(Z))
  }
  Ctilde <- lapply(1:D, function(d) {
    cbind(C[[d]], matrix(rep(Z[d, ], each=p[d]), ncol=H))})
  tCmat <- t(Reduce("rbind", Ctilde))

  # cross-validate overall penalty
  if(is.null(foldid)) {
    foldid <- lapply(1:D, function(d) {
      sample(c(rep(1:nfolds, each=n[[d]] %/% nfolds),
               rep(1:(n[[d]] %% nfolds), (n[[d]] %% nfolds)!=0)))})
  }
  if(is.null(hyper$lambda)) {
    check.cv <- TRUE
    srt <- proc.time()[3]
    if(control$trace) {
      cat("Estimating penalty parameters by CV", "\n", sep="")
    }
    cv.fit <- lapply(1:D, function(d) {
      cv.glmnet(x[idsel[[d]], ], y[[d]], alpha=0, intercept=FALSE,
                foldid=foldid[[d]])})
    hyper$lambda <- 1/sqrt(n*sapply(cv.fit, "[[", "lambda.min"))
    cv.time <- proc.time()[3] - srt
    if(control$trace) {
      printtext <- round(hyper$lambda, 2)
      printtext <- paste(paste(printtext[-length(printtext)], collapse=", "),
                         printtext[length(printtext)], sep=", and ")
      if(nchar(printtext) < 50) {
        cat("\r", "Penalty parameters estimated at ", printtext,
            " in ", round(cv.time, 2), " seconds ", sep="")
      } else {
        cat("\r", "Penalty parameters estimated at ", substr(printtext, 1, 47),
            "..., in ", round(cv.time, 2), " seconds ", sep="")
      }
    }
  } else {
    check.cv <- FALSE
    if(length(hyper$lambda)==1) {
      hyper$lambda <- rep(hyper$lambda, D)
    }
    cv.time <- 0
  }

  # starting values
  alpha.old <- alpha.seq <- rep(0, G + H)
  delta.old <- rep(0, sum(p))
  theta.old <- rep(hyper$lambda^2/2, p)
  if(control$fix.sigma[1]!=FALSE) {
    control$fix.sigma <- as.numeric(rep(
      control$fix.sigma, ifelse(length(control$fix.sigma)==1, D, 1)))
    phi.old <- 1/control$fix.sigma^2
  } else {
    phi.old <- (n/2 + hyper$nu)/(hyper$lambda^2*(2*hyper$zeta + 1))
  }
  par2.new <- list(mat=lapply(1:D, function(d) {0.5*diag(n[d])}))

  # outer loop
  iter.inner <- iter.opt <- numeric(0)
  conv.inner <- conv.opt <- logical(0)
  iter.outer <- 0
  check.outer <- FALSE
  if(control$trace) {
    if(check.cv) {
      cat("\n", "Estimating hyperparameters by EB", "\n", sep="")
    } else {
      cat("Estimating hyperparameters by EB", "\n", sep="")
    }
  }
  srt <- proc.time()[3]
  while(!check.outer) {

    # track iterations
    iter.outer <- iter.outer + 1
    if(control$trace) {
      printtext <- round(alpha.old, 2)
      printtext <- paste(paste(printtext[-length(printtext)], collapse=", "),
                         printtext[length(printtext)], sep=", and ")
      if(nchar(printtext) < 50) {
        cat("\r", "Iteration ", iter.outer, ", estimated EB parameters: ",
            printtext, "      ", sep="")
      } else {
        cat("\r", "Iteration ", iter.outer, ", estimated EB parameters: ",
            substr(printtext, 1, 47), "...", "      ", sep="")
      }
    }

    # update theta and phi
    par1.new <- theta.update(par2.new$mat, x, y, n, p, D, hyper$lambda,
                             hyper$nu, hyper$zeta, idsel, control$fix.sigma)
    theta.new <- par1.new$theta
    phi.new <- par1.new$phi

    # check convergence of outer loop
    conv.outer <- all(abs((theta.new - theta.old)/theta.old) < control$epsilon) &
      all(abs((phi.new - phi.old)/phi.old) < control$epsilon)
    check.outer <- conv.outer | (iter.outer >= control$maxit$eb.outer)

    # inner loop
    iter.inner <- c(iter.inner, 0)
    conv.inner <- c(conv.inner, FALSE)
    check.inner <- FALSE
    while(!check.inner) {
      iter.inner[iter.outer] <- iter.inner[iter.outer] + 1
      print(iter.inner[iter.outer])

      # update parameters
      eb.new <- alpha.update(delta.old, theta.new, phi.new, alpha.old, tCmat,
                             p, D, hyper$lambda, control$opt.method,
                             control$maxit$opt)
      par2.new <- delta.update(eb.new$h, x, y, idsel, D, p)
      alpha.new <- eb.new$alpha
      delta.new <- par2.new$delta
      conv.opt <- c(conv.opt, eb.new$conv)
      iter.opt <- rbind(iter.opt, eb.new$iter)

      # check convergence
      conv.inner[iter.outer] <-
        all(abs((alpha.old - alpha.new)/alpha.old) < control$epsilon)
      check.inner <- conv.inner[iter.outer] | (iter.inner[iter.outer] >=
                                                 control$maxit$eb.inner)

      # assign old parameters
      alpha.seq <- rbind(alpha.seq, alpha.new)
      alpha.old <- alpha.new
      delta.old <- delta.new
    }

    # assign old theta
    theta.old <- theta.new
    phi.old <- phi.new

  }
  eb.time <- proc.time()[3] - srt
  if(control$trace) {
    printtext <- round(alpha.new, 2)
    printtext <- paste(paste(printtext[-length(printtext)], collapse=", "),
                       printtext[length(printtext)], sep=", and ")
    if(nchar(printtext) < 50) {
      cat("\r", "Hyperparameters estimated at ", printtext,
          " in ", round(eb.time, 2), " seconds ", sep="")
    } else {
      cat("\r", "Hyperparameters estimated at ", substr(printtext, 1, 47), "...,",
          " in ", round(eb.time, 2), " seconds ", sep="")
    }
  }

  # fitting final betas
  beta <- lapply(1:D, function(d) {
    h <- exp(colSums(t(Ctilde[[d]])*alpha.new))
    xh <- t(t(x[idsel[[d]], ])*h)
    fit <- glmnet(xh, y[[d]], alpha=0, lambda=1/(hyper$lambda[d]^2*n[d]),
                  standardize=FALSE, intercept=FALSE)
    unname(coef(fit)[-1, 1])*h})

  # creating and returning output object
  if(!exists("cv.fit")) {
    cv.fit <- NULL
  }
  out <- list(beta=beta, alphaf=alpha.new[1:G], alphad=alpha.new[-c(1:G)],
              lambda=hyper$lambda, glmnet=cv.fit,
              time=list(cv=cv.time, eb=eb.time),
              conv=list(eb.outer=conv.outer, eb.inner=conv.inner,
                        opt=conv.opt==0),
              iter=list(eb.outer=iter.outer, eb.inner=iter.inner,
                        opt=iter.opt))
  return(out)
}

################################### ebridge3 ###################################
# alpha <- alpha.old; lambda <- hyper$lambda;
# nu <- hyper$nu; zeta <- hyper$zeta; fix.sigma <- control$fix.sigma
# theta update
par.update <- function(alpha, x, y, idsel, D, n, Ctilde, lambda, nu, zeta,
                       fix.sigma) {
  par <- lapply(1:D, function(d) {
    mat <- lambda[d]^2*x[idsel[[d]], ] %*% 
      (t(x[idsel[[d]], ])*exp(colSums(t(Ctilde[[d]])*alpha)))
    diag(mat) <- diag(mat) + 1
    mat <- solve(mat)
    theta <- lambda[d]^2*colSums(x[idsel[[d]], ]*(mat %*% x[idsel[[d]], ]))
    if(fix.sigma[d]==FALSE) {
      phi <- (n[d] + 2*nu)/(2*zeta + sum(colSums(mat*y[[d]])*y[[d]]))
    } else {
      phi <- 1/fix.sigma[d]^2
    }
    return(c(phi, theta))})
  return(list(theta=unlist(lapply(par, "[", -1)), phi=sapply(par, "[", 1)))
}

# theta <- theta.new; phi <- phi.new; lambda <- hyper$lambda; 
# method <- control$opt.method; maxit <- control$maxit$opt
# alpha update
eb.update <- function(theta, phi, alpha.old, tCmat, x, y, idsel, cp, D, 
                      lambda, method, maxit) {
  opt <- optim(alpha.old, .f.opt.alpha, theta=theta, phi=phi, tCmat=tCmat, x=x,
               y=y, idsel=idsel, cp=cp, D=D, lambda=lambda, 
               method=method, control=list(maxit=maxit))
  alpha <- opt$par
  return(list(alpha=alpha, conv=opt$convergence, iter=opt$counts))
}
.f.opt.alpha <- function(alpha, theta, phi, tCmat, x, y, idsel, cp, D, lambda) {
  vec <- exp(colSums(alpha*tCmat))
  val <- sapply(1:D, function(d) {
    mat <- lambda[d]^2*x[idsel[[d]], ] %*% 
      (t(x[idsel[[d]], ])*vec[(cp[d] + 1):cp[d + 1]])
    diag(mat) <- diag(mat) + 1
    sum(colSums(solve(mat)*y[[d]])*y[[d]])})
  return(sum(theta*vec) + sum(phi*val))
}


# x=xtrain; y=ytrain; C=C; Z=Z; foldid=foldid;
# varsel=NULL
# hyper=list(lambda=fit1.ebridge$lambda, zeta=0, nu=0);
# control=list(epsilon=sqrt(.Machine$double.eps),
#              maxit=list(eb=500, opt=500),
#              fix.sigma=FALSE, trace=TRUE,
#              opt.method="Nelder-Mead",
#              standardize.C=FALSE,
#              standardize.Z=FALSE)
# ebridge2 function
ebridge3 <- function(x, y, C, Z, varsel=NULL, nfolds=10, foldid=NULL,
                     hyper=list(lambda=NULL, zeta=0, nu=0),
                     control=list(epsilon=sqrt(.Machine$double.eps),
                                  maxit=list(eb=500, opt=500),
                                  fix.sigma=FALSE, trace=TRUE,
                                  opt.method="Nelder-Mead",
                                  standardize.C=TRUE,
                                  standardize.Z=TRUE)) {
 
  # auxiliary variables
  D <- length(y)
  if(is.null(varsel)) {
    p <- rep(ncol(x), D)
  } else {
    p <- sapply(varsel, length)
  }
  cp <- c(0, cumsum(p))
  n <- sapply(y, length)
  H <- ncol(Z)
  G <- ncol(C[[1]])
  idsel <- lapply(y, function(s) {
    match(names(s), rownames(x))})
  if(control$standardize.C) {
    C <- lapply(C, function(s) {t(t(s) - colMeans(s))})
  }
  if(control$standardize.Z) {
    Z <- t(t(Z) - colMeans(Z))
  }
  Ctilde <- lapply(1:D, function(d) {
    cbind(C[[d]], matrix(rep(Z[d, ], each=p[d]), ncol=H))})
  tCmat <- t(Reduce("rbind", Ctilde))
  
  # cross-validate overall penalty
  if(is.null(foldid)) {
    foldid <- lapply(1:D, function(d) {
      sample(c(rep(1:nfolds, each=n[[d]] %/% nfolds),
               rep(1:(n[[d]] %% nfolds), (n[[d]] %% nfolds)!=0)))})
  }
  if(is.null(hyper$lambda)) {
    check.cv <- TRUE
    srt <- proc.time()[3]
    if(control$trace) {
      cat("Estimating penalty parameters by CV", "\n", sep="")
    }
    cv.fit <- lapply(1:D, function(d) {
      cv.glmnet(x[idsel[[d]], ], y[[d]], alpha=0, intercept=FALSE,
                foldid=foldid[[d]])})
    hyper$lambda <- 1/sqrt(n*sapply(cv.fit, "[[", "lambda.min"))
    cv.time <- proc.time()[3] - srt
    if(control$trace) {
      printtext <- round(hyper$lambda, 2)
      printtext <- paste(paste(printtext[-length(printtext)], collapse=", "),
                         printtext[length(printtext)], sep=", and ")
      if(nchar(printtext) < 50) {
        cat("\r", "Penalty parameters estimated at ", printtext,
            " in ", round(cv.time, 2), " seconds ", sep="")
      } else {
        cat("\r", "Penalty parameters estimated at ", substr(printtext, 1, 47),
            "..., in ", round(cv.time, 2), " seconds ", sep="")
      }
    }
  } else {
    check.cv <- FALSE
    if(length(hyper$lambda)==1) {
      hyper$lambda <- rep(hyper$lambda, D)
    }
    cv.time <- 0
  }
  control$fix.sigma <- rep(control$fix.sigma, 
                           ifelse(length(control$fix.sigma)==1, D, 1))
  
  # starting values
  alpha.old <- alpha.seq <- rep(0, G + H)
  theta.old <- rep(1, sum(p))
  phi.old <- ifelse(control$fix.sigma==FALSE, 
                    (n/2 + hyper$nu)/(hyper$lambda^2*(2*hyper$zeta + 1)),
                    1/control$fix.sigma^2)
  
  # outer loop
  iter.eb <- 0
  iter.opt <- numeric(0)
  conv.eb <- FALSE
  conv.opt <- numeric(0)
  check.eb <- FALSE
  if(control$trace) {
    if(check.cv) {
      cat("\n", "Estimating hyperparameters by EB", "\n", sep="")
    } else {
      cat("Estimating hyperparameters by EB", "\n", sep="")
    }
  }
  srt <- proc.time()[3]
  while(!check.eb) {
    
    # track iterations
    iter.eb <- iter.eb + 1
    if(control$trace) {
      printtext <- round(alpha.old, 2)
      printtext <- paste(paste(printtext[-length(printtext)], collapse=", "),
                         printtext[length(printtext)], sep=", and ")
      if(nchar(printtext) < 50) {
        cat("\r", "Iteration ", iter.eb, ", estimated EB parameters: ",
            printtext, "      ", sep="")
      } else {
        cat("\r", "Iteration ", iter.eb, ", estimated EB parameters: ",
            substr(printtext, 1, 47), "...", "      ", sep="")
      }
    }
    
    # update theta and phi
    par.new <- par.update(alpha.old, x, y, idsel, D, n, Ctilde, hyper$lambda,
                          hyper$nu, hyper$zeta, control$fix.sigma)
    theta.new <- par.new$theta
    phi.new <- par.new$phi
    
    # update alpha
    eb.new <- eb.update(theta.new, phi.new, alpha.old, tCmat, x, y, idsel, cp,
                        D, hyper$lambda, control$opt.method, 
                        control$maxit$opt)
    alpha.new <- eb.new$alpha
    
    # check convergence of EB loop
    conv.eb <- all(abs(c((theta.new - theta.old)/theta.old, 
                         (phi.new - phi.old)/phi.old,
                         (alpha.new - alpha.old)/alpha.old)) < control$epsilon)
    check.eb <- conv.eb | (iter.eb >= control$maxit$eb)
    
    # track convergence and iterations of optimisation
    conv.opt <- c(conv.opt, eb.new$conv)
    iter.opt <- c(iter.opt, eb.new$iter)
    
    # assign old parameters
    theta.old <- theta.new
    phi.old <- phi.new
    alpha.old <- alpha.new
    alpha.seq <- rbind(alpha.seq, alpha.new)
    
  }
  eb.time <- proc.time()[3] - srt
  if(control$trace) {
    printtext <- round(alpha.new, 2)
    printtext <- paste(paste(printtext[-length(printtext)], collapse=", "),
                       printtext[length(printtext)], sep=", and ")
    if(nchar(printtext) < 50) {
      cat("\r", "Hyperparameters estimated at ", printtext,
          " in ", round(eb.time, 2), " seconds ", sep="")
    } else {
      cat("\r", "Hyperparameters estimated at ", substr(printtext, 1, 47), "...,",
          " in ", round(eb.time, 2), " seconds ", sep="")
    }
  }
  
  # fitting final betas
  beta <- lapply(1:D, function(d) {
    h <- exp(colSums(t(Ctilde[[d]])*alpha.new))
    xh <- t(t(x[idsel[[d]], ])*h)
    fit <- glmnet(xh, y[[d]], alpha=0, lambda=1/(hyper$lambda[d]^2*n[d]),
                  standardize=FALSE, intercept=FALSE)
    unname(coef(fit)[-1, 1])*h})
  
  # creating and returning output object
  if(!exists("cv.fit")) {
    cv.fit <- NULL
  }
  out <- list(beta=beta, alphaf=alpha.new[1:G], alphad=alpha.new[-c(1:G)],
              lambda=hyper$lambda, glmnet=cv.fit,
              time=list(cv=cv.time, eb=eb.time),
              conv=list(eb=conv.eb, opt=conv.opt),
              iter=list(eb=iter.eb, opt=iter.opt))
  return(out)
}