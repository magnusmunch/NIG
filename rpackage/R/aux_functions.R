# function to optimise in EBridge (not tested)
.f.optim <- function(alpha, lambda, nu, zeta, Cmat, Z, n, p, D, G, y, x, yty) {
  
  alphaf <- alpha[1:G]
  alphad <- alpha[-c(1:G)]
  tau <- exp(colSums(alphad*t(Z))/2)
  gamma <- exp(colSums(alphaf*t(Cmat))/2)
  cp <- c(0, cumsum(p))
  out <- sum(sapply(1:D, function(d) {
    mat1 <- mat2 <- x[[d]] %*% (t(x[[d]])*lambda[d]^2*tau[d]^2*
                                  gamma[(cp[d] + 1):cp[d + 1]]^2);
    mat1 <- t(y[, d]) %*% mat1;
    diag(mat2) <- diag(mat2) + 1;
    -determinant(mat2, log=TRUE)$modulus/2 -
      (n/2 + nu)*log(zeta + yty[d]/2 - as.numeric(mat1 %*% y[, d])/2 +
                       as.numeric(mat1 %*% solve(mat2) %*% t(mat1))/2)}))
  return(out)
}

# function to integrate for cpo calculation (not tested)
.f.int.cpo <- function(x, xtSigmax, n, p, zeta, y, xtmu) {
  ((y - x)^2/(2*zeta) + 1)^(-(n + p + 1)/2)*exp(-(x - xtmu)^2/(2*xtSigmax))
}

# ELBO (not tested)
.single.elbo <- function(p, n, zeta, yty, aux, g, b, delta, eta, lambdaf, 
                         lambdad, Zalphad, Calphaf) {
  
  chi <- switch(is.null(Zalphad) + 1, 1/Zalphad, NULL)
  phi <- switch(is.null(Calphaf) + 1, 1/Calphaf, NULL)
  elbo <- -0.5*(n + p + 1)*log(pi) + 0.5*(p + 1 - n)*log(2) + 0.5*(n + 1) + 
    p + lgamma(0.5*(n + p + 1)) - 0.5*(n + p + 1)*log(zeta) + 
    0.5*aux$ldetSigma - 0.25*(n + p + 1)/zeta*
    (yty - 2*aux$ytXmu + aux$trXtXSigma + aux$mutXtXmu + 
       g*sum(b*aux$dSigma) + g*sum(aux$mu^2*b)) + 
    switch(is.null(Zalphad) + 1, 0.25*(p + 3)/4*log(lambdad) + 
             0.5*(eta - lambdad)*g + lambdad/chi - 0.5*(p + 1)*log(chi) - 
             0.25*(p + 1)*log(eta) + 
             gsl::bessel_lnKnu(0.5*(p + 1), sqrt(lambdad*eta)/chi), 
           - 0.5*log(2) + 0.5*log(pi)) +
    switch(is.null(Calphaf) + 1, p*log(lambdaf) + sum(0.5*(delta - lambdaf)*b) + 
             sum(lambdaf/phi) - sum(log(phi)) - 0.5*sum(log(delta)) +
             sum(gsl::bessel_lnKnu(1, sqrt(lambdaf*delta)/phi)), 
           -p/2*log(2) + p/2*log(p))
  return(elbo)
  
}

# computes ratio besselK(x, nu - 1)/besselK(x, nu) (tested)
ratio_besselK <- function(x, nu) {
  res <- nu - floor(nu)
  val <- besselK(x, res - 1, expon.scaled=TRUE)/
    besselK(x, res, expon.scaled=TRUE)
  for(i in 0:(floor(nu) - 1)) {
    val <- 1/(val + 2*(res + i)/x)
  }
  return(val)
}

# cross-validate the hyperparameters of a single equation's MAP
.cv.single.semnig <- function(x, y, nfolds=10, foldid=NULL, seed=NULL,
                              phi, chi, lambdaf, lambdad, type.measure="mse",
                              control=list(trace=FALSE)) {
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  n <- nrow(x)
  p <- ncol(x)
  if(is.null(foldid)) {
    foldid <- sample(c(rep(1:nfolds, each=n %/% nfolds), 
                       rep(1:(n %% nfolds), (n %% nfolds)!=0)))
  }
  par.grid <- expand.grid(phi=phi, chi=chi, lambdaf=lambdaf, lambdad=lambdad)
  
  cv.grid <- matrix(NA, ncol=nrow(par.grid), nrow=nfolds)
  for(l in 1:nrow(par.grid)) {
    if(control$trace) {
      cat("\r", "iteration ", l, " of ", nrow(par.grid))
    }
    for(r in 1:nfolds) {
      
      # split data
      ntrain <- sum(foldid!=r)
      xtrain <- scale(as.matrix(x[foldid!=r, ], nrow=ntrain))
      ytrain <- scale(y[foldid!=r], scale=FALSE)[, 1]
      ntest <- sum(foldid==r)
      xtest <- scale(as.matrix(x[foldid==r, ], nrow=ntest))
      ytest <- scale(y[foldid==r], scale=FALSE)[, 1]
      
      # fit model to training data
      fit <- rstan::optimizing(stanmodels$nig, 
                               data=list(p=p, n=ntrain, x=xtrain, y=ytrain, 
                                         phi=rep(par.grid[l, "phi"], p), 
                                         lambdaf=par.grid[l, "lambdaf"], 
                                         chi=par.grid[l, "chi"], 
                                         lambdad=par.grid[l, "lambdad"]))
      best <- fit$par[names(fit$par) %in% paste0("beta[", 1:p, "]")]
      pred <- as.numeric(xtest %*% best)
      if(type.measure=="mse") {
        cv.grid[r, l] <- mean((ytest - pred)^2)
      }
    }  
  }
  cv.mean <- colMeans(cv.grid)
  id.min <- which.min(cv.mean)
  par.min <- setNames(as.numeric(par.grid[id.min, ]), 
                      c("phi", "chi", "lambdaf", "lambdad"))
  fit <- rstan::optimizing(stanmodels$nig, 
                           data=list(p=p, n=n, x=x, y=y, 
                                     phi=rep(par.min["phi"], p), 
                                     lambdaf=par.min["lambdaf"], 
                                     chi=par.min["chi"], 
                                     lambdad=par.min["lambdad"]))
  out <- list(cvm=cv.mean, cvmat=cv.grid, id.min=id.min, par.min=par.min,
              par.grid=par.grid, fit=fit)
  return(out)
}