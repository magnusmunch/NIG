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

# calculate full covariance with unpenalized covariates
.Sigma.unp <- function(aold, bold, xu, xr, u, r) {
  
  hinv <- 1/bold
  HinvXrt <- t(xr)*hinv
  XrHinvXrt <- Om <- xr %*% HinvXrt
  diag(Om) <- diag(Om) + 1
  Om <- solve(Om)
  HinvXrtOm <- HinvXrt %*% Om
  M11 <- solve(t(xu) %*% Om %*% xu)
  M12 <- - M11 %*% t(xu) %*% t(HinvXrt) + M11 %*% t(xu) %*% t(XrHinvXrt) %*% 
    t(HinvXrtOm)
  M22 <- -HinvXrtOm %*% t(HinvXrt) - HinvXrt %*% xu %*% M12 + HinvXrtOm %*%
    XrHinvXrt %*% xu %*% M12
  diag(M22) <- diag(M22) + hinv
  
  Sigma <- matrix(NA, nrow=u + r, ncol=u + r)
  Sigma[c(1:u), c(1:u)] <- M11/aold
  Sigma[c(1:u), -c(1:u)] <- M12/aold
  Sigma[-c(1:u), c(1:u)] <- t(M12)/aold
  Sigma[-c(1:u), -c(1:u)] <- M22/aold
  
  return(Sigma)
}

# calculate full covariance
.Sigma <- function(aold, bold, x) {
  hinv <- 1/bold
  HinvXt <- t(x)*hinv
  Sigma <- x %*% HinvXt
  diag(Sigma) <- diag(Sigma) + 1
  Sigma <- -HinvXt %*% solve(Sigma) %*% t(HinvXt)
  diag(Sigma) <- diag(Sigma) + hinv
  Sigma <- Sigma/aold
  return(Sigma)
}

# calculates the auxiliary variables with unpenalized variables (not tested)
.aux.var.unp <- function(aold, bold, y, xu, xr, u, r) {
  
  hinv <- 1/bold
  HinvXrt <- t(xr)*hinv
  XrHinvXrt <- Om <- xr %*% HinvXrt
  diag(Om) <- diag(Om) + 1
  Om <- solve(Om)
  HinvXrtOm <- HinvXrt %*% Om
  M11 <- solve(t(xu) %*% Om %*% xu)
  M12 <- - M11 %*% t(xu) %*% t(HinvXrt) + M11 %*% t(xu) %*% t(XrHinvXrt) %*% 
    t(HinvXrtOm)
  
  # auxiliaries
  dSigma <- c(diag(M11), hinv - rowSums(HinvXrt*HinvXrtOm) - 
                rowSums((HinvXrt %*% xu)*t(M12)) + 
                rowSums(HinvXrtOm*(t(M12) %*% t(xu) %*% XrHinvXrt)))
  mu <- rbind(M11 %*% t(xu) + M12 %*% t(xr),
              HinvXrt - HinvXrtOm %*% XrHinvXrt -
                HinvXrt %*% xu %*% M12 %*% t(xr) + HinvXrtOm %*% XrHinvXrt %*%
                xu %*% M12 %*%t(xr))
  trXtXSigma <- (sum(xu*(xu %*% M11)) + 2*sum(xu*(xr %*% t(M12))) + 
                   sum(t(xr)*mu[-c(1:u), ]))/aold
  mu[-c(1:u), ] <- mu[-c(1:u), ] + t(xu %*% M12)
  mu <- colSums(t(mu)*y)
  ytXmu <- cbind(xu, xr) %*% mu
  mutXtXmu <- as.numeric(t(ytXmu) %*% ytXmu)
  ytXmu <- sum(as.numeric(ytXmu)*y)
  trHSigma <- sum(dSigma[-c(1:u)]*bold)
  mutHmu <- sum(mu[-c(1:u)]^2*bold)
  
  out <- list(mu=mu, dSigma=dSigma, ytXmu=ytXmu, trXtXSigma=trXtXSigma,
              mutXtXmu=mutXtXmu, trHSigma=trHSigma, mutHmu=mutHmu)
  return(out)
}

# calculates the auxiliary variables (not tested)
.aux.var <- function(aold, bold, y, x, ytx) {
  
  hinv <- 1/bold
  HinvXt <- t(x)*hinv
  XHinvXt <- x %*% HinvXt
  
  # shared auxiliaries
  dSigma <- XHinvXt
  diag(dSigma) <- diag(dSigma) + 1
  dSigma <- HinvXt %*% solve(dSigma)
  mu <- trXtXSigma <- mutXt <-
    HinvXt - dSigma %*% XHinvXt
  dSigma <- (hinv - rowSums(dSigma*HinvXt))/aold
  mu <- as.numeric(mu %*% y)
  mutXt <- as.numeric(ytx %*% mutXt)
  trXtXSigma <- sum(trXtXSigma*t(x))/aold
  
  # separate auxiliaries
  ytXmu <- as.numeric(ytx %*% mu)
  mutXtXmu <- sum(mutXt^2)
  trHSigma <- sum(dSigma*bold)
  mutHmu <- sum(mu^2*bold)

  out <- list(mu=mu, dSigma=dSigma, ytXmu=ytXmu, trXtXSigma=trXtXSigma,
              mutXtXmu=mutXtXmu, trHSigma=trHSigma, mutHmu=mutHmu)
  return(out)
}

# to integrate in CPO calculation
.f.int.cpo <- function(val, y, zeta, n, p, xtmu, xtSigmax) {
  (val^2/(2*zeta) + 1)^(-(n - length(xtmu) + p + 2)/2)*
    exp(-(val - xtmu + y)^2/(2*xtSigmax))
}

# calculate the conditional predictive ordinates
.cpo <- function(model, x, y, n, p, D) {
  out <- lapply(1:D, function(d) {
    zeta <- model$vb$zeta[d]
    xtmu <- as.numeric(x[[d]] %*% model$vb$mu[[d]])
    xtSigmax <- rowSums((x[[d]] %*% model$vb$Sigma[[d]])*x[[d]])
    int <- sapply(1:nrow(matrix(y, ncol=D)), function(i) {
      unlist(integrate(.f.int.cpo, -Inf, Inf, y=matrix(y, ncol=D)[i, d], 
                       zeta=zeta, n=n, p=p[d], xtmu=xtmu[i], 
                       xtSigmax=xtSigmax[i])[c(1:2)])})
    
    
    list((2*pi)^(-1)*exp(lgamma((n- nrow(matrix(y, ncol=D)) + p[d] + 2)/2) - 
                           lgamma((n - nrow(matrix(y, ncol=D)) + p[d] + 1)/2))*
           (zeta*xtSigmax)^(-1/2)*int[1, ], int[2, ])})
  
  return(list(value=unname(sapply(out, "[[", 1)), 
              abs.error=unname(sapply(out, "[[", 2))))
}

# calculate log pseudo marginal likelihood
lpml <- function(x, y, C, unpenalized=NULL, standardize=TRUE, intercept=TRUE, 
                 intercept.eb=TRUE, mult.lambda=FALSE, 
                 fixed.eb=c("none", "lambda", "both"), full.post=FALSE, 
                 init=NULL,
                 control=list(conv.post=TRUE, trace=TRUE,
                              epsilon.eb=1e-3, epsilon.vb=1e-3, 
                              epsilon.opt=sqrt(.Machine$double.eps),
                              maxit.eb=10, maxit.vb=2, maxit.opt=100,
                              maxit.post=100),
                 nfolds, foldid=NULL) {
  
  # create the folds
  p <- sapply(x, ncol)
  D <- ncol(y)
  n <- nrow(y)
  if(is.null(foldid)) {
    foldid <- sample(rep(1:nfolds, times=round(c(rep(
      n %/% nfolds + as.numeric((n %% nfolds)!=0), times=n %% nfolds),
      rep(n %/% nfolds, times=nfolds - n %% nfolds)))))  
  } else {
    nfolds <- length(unique(foldid))
  }
  
  # prepare object to store CPOi
  val <- err <- matrix(NA, ncol=D, nrow=n)
  
  # loop over the folds
  for(k in 1:nfolds) {
    print(paste0("fold ", k))
    xtrain <- lapply(x, function(x) {x[foldid!=k, ]})
    ytrain <- y[foldid!=k, ]
    xtest <- lapply(x, function(x) {x[foldid==k, ]})
    ytest <- y[foldid==k, ]
    
    # estimate the model on the training data
    fit <- enig(x=xtrain, y=ytrain, C=C, unpenalized=unpenalized, 
                standardize=standardize, intercept=intercept, 
                intercept.eb=intercept.eb, mult.lambda=mult.lambda, 
                fixed.eb=fixed.eb, full.post=full.post, init=init, 
                control=control)
    
    # estimate CPOi on the test data
    est <- .cpo(fit, xtest, ytest, n, p, D)
    val[foldid==k, ] <- est$value
    err[foldid==k, ] <- est$abs.error
    
  }
  
  out <- list(lpml=colSums(log(val))/n, cpo=val, abs.error=err)
  return(out)
  
}