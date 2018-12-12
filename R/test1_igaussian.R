### paths
path.code <- "/Users/magnusmunch/Documents/OneDrive/PhD/project_cambridge/code/"

### libraries
library(Rcpp)
library(RcppArmadillo)
library(statmod)
library(gsl) 

### Compile cpp functions
sourceCpp(paste(path.code, "myfunctions.cpp", sep=""))

### functions
# computes the ratio besselK(x, p/2 - 1/2)/besselK(x, p/2 + 1/2)
# tested and works
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

# c++ version of ratio computation
# tested and works
cppFunction("
NumericVector ratio_besselK_cpp(NumericVector x, int p) {
  int n = x.size();
  NumericVector value(n);
  bool even = (p % 2)==0;
  int niter = floor(p/2) + 1 - even;
  double num;
  if(even) {
    num = 0.5;
  } else {
    num = 1;
  }
  for(int i=0; i<n; i++) {
    if(even) {
      value[i] = 1;
    } else {
      value[i] = R::bessel_k(x[i], 1, exp(1))/R::bessel_k(x[i], 0, exp(1));
    }
    for(int j=0; j<niter; j++) {
      value[i] = 1/(value[i] + 2*(j + 1 - num)/x[i]);
    }
  }
  return value;
}")

# updates the variational parameters
# not tested
vb.update.fixed.sigma <- function(aold, lambda, theta, sigma, x, xxt, y) {
  p <- ncol(x)
  D <- length(sigma)
  
  # Sigma estimation checked and correct
  Sigma <- lapply(c(1:D), function(d) {
    mat <- xxt
    diag(mat) <- diag(mat) + aold[d]
    mat <- -(t(x) %*% solve(mat) %*% x)/aold[d]
    diag(mat) <- diag(mat) + 1/aold[d]
    return(sigma[d]^2*mat)})
  
  
  mu <- sapply(c(1:D), function(d) {Sigma[[d]] %*% t(x) %*% y[, d]/sigma[d]^2})
  delta <- sapply(c(1:D), function(d) {
    (sum(diag(Sigma[[d]])) + sum(mu[, d]^2))/sigma[d]^2 + lambda})
  a <- sqrt(lambda/theta^2*delta)*
    ratio_besselK_cpp(sqrt(lambda*delta/theta^2), p) + (p + 1)/delta
  return(list(Sigma=Sigma, mu=mu, delta=delta, a=a))
}

# VB update with sigma stochastic
vb.update <- function(zetaold, aold, lambda, theta, x, xxt, xtx, y) {
  p <- ncol(x)
  D <- length(sigma)
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

# VB update under independent prior beta and sigmas
vb.update2 <- function(zetaold, aold, lambda, theta, x, xxt, xtx, y) {
  p <- ncol(x)
  D <- length(sigma)
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

# updates the hyper-parameters
eb.update <- function(v, a, C) {
  D <- length(v)
  alpha <- rowSums(solve(t(C) %*% diag(v) %*% C) %*% t(C))
  lambda <- D/(sum(a) - sum(t(alpha) %*% t(C)))
  theta <- 1/as.numeric(C %*% alpha)
  return(list(alpha=alpha, lambda=lambda, theta=theta))
}


# mu=new.vb$mu
# Sigma=new.vb$Sigma
# delta=new.vb$delta
# lambda=old.eb$lambda
# theta=old.eb$theta

# evaluate the evidence lower bound
# doesn't work, numerical overflow for log of BesselK 
elbo.eval <- function(mu, Sigma, delta, lambda, theta, sigma, x, xtx, y) {
  D <- length(sigma)
  p <- nrow(mu)
  
  # alternative for log BesselK function:
  logK <- bessel_lnKnu(0.5*(p + 1), sqrt(lambda*delta)/theta)
  # logK <- log(besselK(sqrt(lambda*delta)/theta, 0.5*(p + 1)))
  
  elbo <- as.numeric(sapply(c(1:D), function(d) {
    t(y[, d]) %*% x %*% mu[, d]/sigma[d]^2 
    - 0.5*sum(diag(xtx %*% Sigma[[d]]))/sigma[d]^2
    - 0.5*t(mu[, d]) %*% xtx %*% mu[, d]/sigma[d]^2
    + 0.5*determinant(Sigma[[d]], log=TRUE)$modulus
    + lambda/theta[d]
    - 0.5*(p + 1)*log(theta[d])
    - 0.25*(p + 1)*log(delta[d])
    + logK[d]
    +0.5*(delta[d] - lambda - t(mu[, d]) %*% mu[, d]/sigma[d]^2 
          - sum(diag(Sigma[[d]]))/sigma[d]^2)*
      (sqrt(lambda/delta[d])*ratio.besselK(sqrt(lambda*delta[d])/theta[d], p)/
         theta[d] + (p + 1)/delta[d])})) + 0.25*(p + 3)*log(lambda)
  return(elbo)
}

# epsilon.vb=Inf means one vb iteration per eb iteration
est.igauss <- function(x, y, C,
                       control=list(epsilon.eb=1e-3, epsilon.vb=1e-3, 
                                    maxit.eb=100, maxit.vb=100, trace=TRUE), 
                       init=list(alpha, lambda, a, zeta),
                       test=list(ind.var=FALSE)) {
  
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
    init.mll.est <- t(sapply(c(1:D), function(d) {
      constrOptim(c(1, 1), marg.lik1, NULL, ui=diag(2), ci=c(0, 0),
                  sv=svd.x$d, uty=t(svd.x$u) %*% y[, d], n=n,
                  control=list(fnscale=-1))$par}))
    
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
    init.zeta <- 0.5*D*(n + p + 1)/sum(1/init.mll.est[, 1])
    
    init <- list(alpha=c(init.theta, rep(0, ncol(C) - 1)), lambda=init.lambda,
                 a=init.a, zeta=init.zeta)
    
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
    new.eb <- eb.update(old.vb$v, old.vb$a, C)
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
  names(eb.seq) <- names(old.eb)
  
  # preparing the output
  out <- list(vb.post=list(mu=old.vb$mu, Sigma=old.vb$Sigma, 
                           delta=old.vb$delta, zeta=old.vb$zeta), 
              eb.seq=eb.seq,
              elbo.seq=elbo.seq,
              conv=list(eb=conv.eb, vb=conv.vb),
              iter=list(eb=iter.eb, vb=iter.vb))
  return(out)
}

# Gwen's model
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


set.seed(123)
n <- 100
p <- 200
D <- 100
nclass <- 5

alpha <- c(1:nclass)
C <- model.matrix(~ 1 + factor(rep(c(1:nclass), each=D/nclass)))
theta <- as.numeric(1/(C %*% alpha))
gamma <- sqrt(theta)
sigma <- rep(1, D)
SNR <- 10

beta <- sapply(1:D, function(d) {rnorm(p, 0, sigma[d]*gamma[d])})
x <- matrix(rnorm(n*p, 0, sqrt(SNR/(mean(gamma^2)*p))), nrow=n, ncol=p)
y <- sapply(1:D, function(d) {rnorm(n, x %*% beta[, d], sigma[d])})

fit1.igauss1 <- est.igauss(x, y, C, 
                           control=list(epsilon.eb=1e-3, epsilon.vb=1e-3, 
                                        maxit.eb=20, maxit.vb=2, trace=TRUE), 
                           init <- list(lambda=1, alpha=rep(1, ncol(C)), 
                                        a=rep(0.1, D), zeta=rep(0.1, D)),
                           test=list(ind.var=FALSE))

fit1.igauss2 <- est.igauss(x, y, C, 
                           control=list(epsilon.eb=1e-3, epsilon.vb=1e-3, 
                                        maxit.eb=20, maxit.vb=2, trace=TRUE), 
                           init <- list(lambda=1, alpha=rep(1, ncol(C)), 
                                        a=rep(0.1, D), zeta=rep(0.1, D)),
                           test=list(ind.var=TRUE))

fit1.gwen <- est.gwen(x, y, eqid=rep(c(1:nclass), each=D/nclass),
                      control=list(epsilon.eb=1e-3, epsilon.vb=1e-3, 
                                   epsilon.opt=1e-6, maxit.eb=20, maxit.vb=2, 
                                   maxit.opt=20, conv.vb="elbo", 
                                   trace=TRUE), 
                      init=list(aprior=0.001, bprior=0.001))


# ### convergence
# # ELBOs
# plot(rowSums(fit1.igauss1$elbo.seq), type="l")
# plot(rowSums(fit1.gwen$elbo.seq), type="l")
# 
# # igauss alpha convergences
# plot(fit1.igauss1$eb.seq$alpha[, 1], type="l", ylim=range(fit1.igauss1$eb.seq$alpha))
# lines(fit1.igauss1$eb.seq$alpha[, 2], col=2)
# lines(fit1.igauss1$eb.seq$alpha[, 3], col=3)
# lines(fit1.igauss1$eb.seq$alpha[, 4], col=4)
# lines(fit1.igauss1$eb.seq$alpha[, 5], col=5)
# plot(fit1.igauss2$eb.seq$alpha[, 1], type="l", ylim=range(fit1.igauss2$eb.seq$alpha))
# lines(fit1.igauss2$eb.seq$alpha[, 2], col=2)
# lines(fit1.igauss2$eb.seq$alpha[, 3], col=3)
# lines(fit1.igauss2$eb.seq$alpha[, 4], col=4)
# lines(fit1.igauss2$eb.seq$alpha[, 5], col=5)
# 
# # theta convergence
# plot(fit1.igauss1$eb.seq$theta[, 1], type="l", ylim=range(fit1.igauss1$eb.seq$theta))
# lines(fit1.igauss1$eb.seq$theta[, D/nclass + 1], col=2)
# lines(fit1.igauss1$eb.seq$theta[, 2*D/nclass + 1], col=3)
# lines(fit1.igauss1$eb.seq$theta[, 3*D/nclass + 1], col=4)
# lines(fit1.igauss1$eb.seq$theta[, 4*D/nclass + 1], col=5)
# plot(fit1.igauss2$eb.seq$theta[, 1], type="l", ylim=range(fit1.igauss2$eb.seq$theta))
# lines(fit1.igauss2$eb.seq$theta[, D/nclass + 1], col=2)
# lines(fit1.igauss2$eb.seq$theta[, 2*D/nclass + 1], col=3)
# lines(fit1.igauss2$eb.seq$theta[, 3*D/nclass + 1], col=4)
# lines(fit1.igauss2$eb.seq$theta[, 4*D/nclass + 1], col=5)
# 
# # igauss lambda convergence
# plot(fit1.igauss1$eb.seq$lambda, type="l")
# plot(fit1.igauss2$eb.seq$lambda, type="l")

### posterior mean estimates against true values
# sigma^2
omar <- par()$mar
par(mfrow=c(1, 3), mar=omar + c(0, 1, 0, 0))
plot(sigma^2, 2*fit1.igauss1$vb.post$zeta/(n + p - 1), xlab=expression(sigma^2),
     ylab=expression(E(sigma^2 ~ "|" ~ y)), main="a)")
plot(sigma^2, 2*fit1.igauss2$vb.post$zeta/(n - 2.5), xlab=expression(sigma^2),
     ylab=expression(E(sigma^2 ~ "|" ~ y)), main="b)")
plot(sigma^2, fit1.gwen$vb.post$dpost/(fit1.gwen$vb.post$cpost - 1), 
     xlab=expression(sigma^2), ylab=expression(E(sigma^2 ~ "|" ~ y)), main="c)")
par(mfrow=c(1, 1), mar=omar)

# gamma^2
est1 <- sqrt(fit1.igauss1$vb.post$delta*
               fit1.igauss1$eb.seq$theta[fit1.igauss1$iter$eb, ]^2/
               fit1.igauss1$eb.seq$lambda[fit1.igauss1$iter$eb, ])*
  ratio_besselK_cpp(sqrt(fit1.igauss1$vb.post$delta*fit1.igauss1$eb.seq$lambda[
    fit1.igauss1$iter$eb, ]/fit1.igauss1$vb.post$delta*
      fit1.igauss1$eb.seq$theta[fit1.igauss1$iter$eb, ]^2), p)
est2 <- sqrt(fit1.igauss2$vb.post$delta*
               fit1.igauss2$eb.seq$theta[fit1.igauss2$iter$eb, ]^2/
               fit1.igauss2$eb.seq$lambda[fit1.igauss2$iter$eb, ])*
  ratio_besselK_cpp(sqrt(fit1.igauss2$vb.post$delta*fit1.igauss2$eb.seq$lambda[
    fit1.igauss2$iter$eb, ]/fit1.igauss2$vb.post$delta*
      fit1.igauss2$eb.seq$theta[fit1.igauss2$iter$eb, ]^2), p)
omar <- par()$mar
par(mfrow=c(1, 3), mar=omar + c(0, 1, 0, 0))
plot(gamma^2, est1, xlab=expression(gamma^2),
     ylab=expression(E(gamma^2 ~ "|" ~ y)), main="a)")
abline(a=0, b=1, lty=2)
plot(gamma^2, est2, xlab=expression(gamma^2),
     ylab=expression(E(gamma^2 ~ "|" ~ y)), main="b)")
abline(a=0, b=1, lty=2)
plot(gamma^2, fit1.gwen$vb.post$bpost/(fit1.gwen$vb.post$apost - 1),
     ylab=expression(E(gamma^2 ~ "|" ~ y)), main="c)")
abline(a=0, b=1, lty=2)
par(mfrow=c(1, 1), mar=omar)

# beta
omar <- par()$mar
par(mfrow=c(1, 3), mar=omar + c(0, 1, 0, 0))
plot(beta, fit1.igauss1$vb.post$mu, xlab=expression(beta),
     ylab=expression(E(beta ~ "|" ~ y)), main="a)")
abline(a=0, b=1, lty=2)
plot(beta, fit1.igauss2$vb.post$mu, xlab=expression(beta),
     ylab=expression(E(beta ~ "|" ~ y)), main="b)")
abline(a=0, b=1, lty=2)
plot(beta, t(fit1.gwen$vb.post$mu), xlab=expression(beta),
     ylab=expression(E(beta ~ "|" ~ y)), main="c)")
abline(a=0, b=1, lty=2)
par(mfrow=c(1, 1), mar=omar)

### EB estimation
# prior mean estimates of gamma^2 against true values
omar <- par()$mar
par(mfrow=c(1, 3), mar=omar + c(0, 1, 0, 0))
plot(theta, 1/(C %*% fit1.igauss1$eb.seq$alpha[fit1.igauss1$iter$eb, ]),
     xlab=expression(E(gamma^2)), ylab=expression(hat(E(gamma^2))), main="a)")
abline(a=0, b=1, lty=2)
plot(theta, 1/(C %*% fit1.igauss2$eb.seq$alpha[fit1.igauss2$iter$eb, ]),
     xlab=expression(E(gamma^2)), ylab=expression(hat(E(gamma^2))), main="b)")
abline(a=0, b=1, lty=2)
plot(theta, rep(fit1.gwen$eb.seq$bprior[fit1.gwen$iter$eb + 1, ]/
                  (fit1.gwen$eb.seq$aprior[fit1.gwen$iter$eb + 1, ] - 1), 
                each=D/nclass),
     xlab=expression(E(gamma^2)), ylab=expression(hat(E(gamma^2))), main="c)")
abline(a=0, b=1, lty=2)
par(mfrow=c(1, 1), mar=omar)





# f.root <- function(x, vary, q) {
#   gamma(3/2)*(1 - q - pgamma(0.5*x*3/vary, 3/2))
# }
# sigma.squared.hat <- (3/5)*apply(y, 2, function(c) {
#   uniroot(f.root, interval=c(0.0001, var(c)), vary=var(c), q=0.95)$root})
#
# f.root <- function(x, vary, n, p, prob) {
#   pgamma(1/vary, 0.5*(n + p + 1), x) + prob - 1
# }
# 
# curve(f.root(x, var(y[, 1]), n, p, 0.9), 0.0001, 100)
# apply(y, 2, function(c) {uniroot(f.root, c(0.0001, var(c)), vary=var(c), n=n, 
#                                  p=p, prob=0.9)$root})
# 
# zeta.start <- 
# library(invgamma)
# curve(dinvgamma(x, 0.5*(p + n + 1), 1), 0.001, 1)
# 
# 
# theta.start <- mean(sigma.squared.hat)
# lambda.start <- length(sigma.squared.hat)/sum(1/sigma.squared.hat - 1/theta.start)
# curve(dinvgauss(x, theta.start, lambda.start), 0.01, 10)

# check log besselK function
x <- seq(0.1, 1000, 0.1)
nu <- seq(1, 1000, 0.5)
grid <- expand.grid(x, nu)
library(gsl)
bessel_lnKnu


# starting values problem
init.mll <- function(par, sv, uty, n) {
  sigma.sq <- par[1]
  gamma.sq <- par[2]
  mll <- -0.5*n*log(sigma.sq) - 0.5*n*log(gamma.sq) - 
    sum(log(sv^2 + 1/gamma.sq)) - 0.5*sum(uty^2/(sv^2 + 1/gamma.sq))/
    (sigma.sq*gamma.sq)
  return(mll)
}

# gives correct gamma^2*sigma^2, but unstable in terms of separatio of the two
# marg.lik2 <- function(par, sv, uty, n) {
#   lsigma.sq <- par[1]
#   lgamma.sq <- par[2]
#   mll <- -0.5*n*lsigma.sq - 0.5*n*lgamma.sq - 
#     sum(log(sv^2 + exp(-lgamma.sq))) - 
#     0.5*exp(-lsigma.sq - lgamma.sq)*
#     sum(uty^2*exp(lgamma.sq)/(sv^2*exp(lgamma.sq) + 1))
#   return(mll)
# }









