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
