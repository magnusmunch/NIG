################################# VB functions #################################
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
                              yty, conjugate, hyperprior) {
  
  # variables that change with conjugacy
  cold <- bold/(aold^(1 - conjugate))
  df <- n + conjugate*p
  
  # auxiliary variables involving mu and Sigma
  aux <- .aux.var(cold, aold, sv, uty, n, p)
  
  # vb parameters
  delta <- unname((bold/cold)*aux$mutmu + lambda)
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
.init.param.inv.gauss <- function(C, sv, uty, D, n, p) {
  # maximise ridge mll to find gamma_d^2 and sigma_d^2 estimates
  init.mll.est <- t(sapply(c(1:D), function(d) {
      constrOptim(c(1, 1), .ridge.mll, NULL, ui=diag(2), ci=c(0, 0),
                  sv=sv, uty=uty[, d], n=n, control=list(fnscale=-1))$par}))
    
  # consider them samples from prior and estimate one inverse Gaussian prior
  theta <- mean(init.mll.est[, 2])
  lambda <- D/sum(1/init.mll.est[, 2] - 1/theta)
  
  # use empirical mode of gamma_d^2 to derive one delta
  d <- density(init.mll.est[, 2])
  mode <- d$x[which.max(d$y)]
  delta <- as.numeric(mode^2*lambda/theta^2 + (p + 3)*mode)
  a <- sqrt(lambda/(theta^2*delta))*
    ratio_besselK(sqrt(lambda*delta)/theta, p) + (p + 1)/delta
  
  # consider sigma_d^2 fixed and estimate inverse Gamma prior
  zeta <- 0.5*D*(n + p + 1)/sum(1/init.mll.est[, 1])
  
  out <- list(theta=rep(theta, D), alpha=c(theta, rep(0, ncol(C) - 1)), 
              lambda=lambda, zeta=rep(zeta, D), a=rep(a, D))
  return(out)
  
}

# initial values with independent beta and sigma^2 prior (not tested)
.init.param.inv.gauss.ind <- function(C, sv, uty, D, n, p) {
  # maximise ridge mll to find gamma_d^2 and sigma_d^2 estimates
  init.mll.est <- t(sapply(c(1:D), function(d) {
        constrOptim(c(1, 1), .ridge.mll.ind, NULL, ui=diag(2), ci=c(0, 0),
                    sv=sv, uty=uty[, d], n=n, control=list(fnscale=-1))$par}))
  
  # consider them samples from prior and estimate one inverse Gaussian prior
  theta <- mean(init.mll.est[, 2])
  lambda <- D/sum(1/init.mll.est[, 2] - 1/theta)
  
  # use empirical mode of gamma_d^2 to derive one delta
  d <- density(init.mll.est[, 2])
  mode <- d$x[which.max(d$y)]
  delta <- as.numeric(mode^2*lambda/theta^2 + (p + 3)*mode)
  a <- sqrt(lambda/(theta^2*delta))*
    ratio_besselK(sqrt(lambda*delta)/theta, p) + (p + 1)/delta
  
  # consider sigma_d^2 fixed and estimate inverse Gamma prior
  zeta <- 0.5*D*(n + 1)/sum(1/init.mll.est[, 1])
  
  out <- list(theta=rep(theta, D), alpha=c(theta, rep(0, ncol(C) - 1)), 
              lambda=lambda, zeta=rep(zeta, D), a=rep(a, D))
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

################################# Gibbs sampler ################################
# sample one beta in the conjugate case
.sample.beta.conj <- function(gamma, sigma, v, sv, svvt, y) {
  uvec <- rnorm(p, 0, gamma*sigma)
  dvec <- rnorm(n, 0, 1)
  beta <- uvec + as.numeric(t(t(v)*(sv/(sv^2 + 1/gamma^2))) %*% 
                              (sigma*y - sigma*svvt %*% uvec + dvec))
  return(beta)
}

# sample one gamma in the conjugate case
.sample.gamma.conj <- function(beta, sigma, eta, lambda, theta, p) {
  gamma <- sqrt(rgig(1, sum(beta^2)/sigma^2 + lambda, 
                     lambda/theta^2, -0.5*(p + eta)))
  return(gamma)
}

# sample one sigma in the conjugate case
.sample.sigma.conj <- function(beta, gamma, yty, ytx, x, p, n) {
  alpha <- 0.5*(p + n + 1)
  beta <- 0.5*(yty + sum(beta^2)/gamma^2 + sum((x %*% beta)^2) - 2*ytx %*% beta)
  gamma <- 1/sqrt(rgamma(1, alpha, beta))
  return(gamma)
}

nsamp <- control$nsamp
init <- list(beta=beta, sigma=sigma, gamma=gamma)
# full sample from the posterior
.sample.post.conj <- function(eta, lambda, theta, nsamp, init, v, sv, svvt, y, 
                              yty, ytx, x, p, n) {
  
  # create the output and initial values
  out <- list(beta=matrix(NA, nrow=p, ncol=nsamp), sigma=rep(NA, nsamp),
              gamma=rep(NA, nsamp))
  beta <- init$beta
  sigma <- init$sigma
  
  # run the Gibbs sampler
  for(m in 1:nsamp) {
    gamma <- .sample.gamma.conj(beta, sigma, eta, lambda, theta, p)
    beta <- .sample.beta.conj(gamma, sigma, v, sv, svvt, y)
    sigma <- .sample.sigma.conj(beta, gamma, yty, ytx, x, p, n)
    out$beta[, m] <- beta
    out$sigma[m] <- sigma
    out$gamma[m] <- gamma
  }
  
  return(out)
}

D <- 20
n <- 20
p <- 30
x <- matrix(rnorm(n*p), nrow=n, ncol=p)
y <- matrix(rnorm(n*D), ncol=D, nrow=n)
eta <-rchisq(D, 1)
theta <- rchisq(D, 1)
lambda <- rchisq(D, 1)
init <- list(beta=matrix(rnorm(p*D), nrow=p, ncol=D), 
             gamma=rchisq(D, 1), sigma=rchisq(D, 1))
control <- list(nsamp=1000, parallel=TRUE, ncores=2)
gibbs.model <- function(x, y, hyperprior=c("inv. Gaussian", "inv. Gamma"), 
                        conjugate, init=NULL,
                        control=list(nsamp, parallel=FALSE, ncores=1)) {
  
  # set the number of cores
  ncores <- min(control$ncores, detectCores())
  
  # set fixed parameters
  D <- ncol(y)
  svd.x <- svd(x)
  v <- svd.x$v
  sv <- svd.x$d
  svvt <- t(v)*sv
  yty <- colSums(y^2)
  ytx <- t(y) %*% x
  
  # divide the equations over the cores
  # if(parallel) {
  #   coreid <- sort(rep(1:ncores, times=round(c(rep(
  #     D %/% ncores + as.numeric((D %% ncores)!=0), times=D %% ncores),
  #     rep(D %/% ncores, times=ncores - D %% ncores)))))
  # }
  
  out <- mclapply(c(1:D), function(d) {
    cinit <- list(beta=init$beta[, d], gamma=init$gamma[d], sigma=init$sigma[d])
    .sample.post.conj(eta[d], lambda[d], theta[d], control$nsamp, cinit, v, sv, 
                      svvt, y[, d], yty[d], ytx[d, ], x, p, n)}, 
    mc.cores=ncores)
  
  
}

mclapply(starts, fx, mc.cores = numCores)



library(parallel)




