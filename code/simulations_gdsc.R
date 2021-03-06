#!/usr/bin/env Rscript

### libraries
packages <- c("foreach", "doParallel", "NIG", "statmod", "glmnet", "rstan")
sapply(packages, library, character.only=TRUE)

### load and preprocess data
load(file="data/data_gdsc_dat1.Rdata")

################################################################################
################################# simulation 1 #################################
################################################################################
### simulation settings
D <- length(resp)
p <- 100
n <- floor(nrow(expr$expr))/2
ntest <- nrow(expr$expr) - n
alphaf <- c(1, 1, 3, 7)
G <- length(alphaf)
shape <- 3
rate <- 2
lambdaf <- 1
C <- replicate(D, list(unname(model.matrix(~ factor(rep(1:G, each=p/G))))))
x <- scale(expr$expr[, order(-apply(expr$expr, 2, sd))[1:p]])

# estimation settings
nfolds <- 10
control <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-5, 
                epsilon.vb=1e-3, maxit.eb=500, maxit.vb=1, 
                maxit.post=100, maxit.block=0)


# setup cluster
if(Sys.info()["sysname"]=="Windows") {
  cl <- makeCluster(min(ncores, nreps), type="PSOCK")
} else {
  cl <- makeCluster(min(ncores, nreps), type="FORK")
}
if(min(ncores, nreps) > 1) {
  registerDoParallel(cl)
} else {
  registerDoSEQ()
}
res <- foreach(r=1:nreps, .packages=packages) %dopar% {
  cat("\r", "rep", r)
  set.seed(2020 + r)

  ### simulate parameters
  gamma <- sapply(C, function(s) {
    sqrt(rinvgauss(p, 1/as.numeric(s %*% alphaf), lambdaf))})
  sigma <- 1/sqrt(rgamma(D, shape, rate))
  
  beta <- sapply(1:D, function(d) {rnorm(p, 0, gamma[, d]*sigma[d])})

  # draw and simulate data
  idtrain <- sample(1:nrow(x), n)
  xtrain <- scale(x[idtrain, ])
  rownames(xtrain) <- c(1:n)
  xtest <- scale(x[-idtrain, ])
  rownames(xtest) <- c(1:ntest)
  ytrain <- lapply(1:D, function(d) {
    s <- as.numeric(scale(rnorm(n, xtrain %*% beta[, d], sigma[d])))
    names(s) <- c(1:length(s))
    return(s)})
  ytest <- lapply(1:D, function(d) {
    s <- as.numeric(scale(rnorm(ntest, xtest %*% beta[, d], sigma[d])))
    names(s) <- c(1:length(s))
    return(s)})
  
  # creating folds
  foldid <- sample(c(rep(1:nfolds, each=n %/% nfolds),
                     rep(1:(n %% nfolds), (n %% nfolds)!=0)))

  # fitting models
  fit.nig1 <- nig(x=rep(list(xtrain), D), y=ytrain, 
                        C=replicate(D, list(matrix(rep(1, p)))), 
                        Z=NULL, full.post=TRUE, control=control)
  fit.nig2 <- nig(x=rep(list(xtrain), D), y=ytrain, C=C, Z=NULL, 
                        full.post=TRUE, control=control)

  # standard penalized methods
  fit.ridge1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain, ytrain[[d]], alpha=0, intercept=FALSE, standardize=TRUE,
              foldid=foldid)})
  fit.lasso1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain, ytrain[[d]], intercept=FALSE, standardize=TRUE,
              foldid=foldid)})
  
  # estimates
  best <- list(nig1=Reduce("cbind", fit.nig1$vb$mpost$beta), 
               nig2=Reduce("cbind", fit.nig2$vb$mpost$beta),
               ridge1=sapply(fit.ridge1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}),
               lasso1=sapply(fit.lasso1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}))
  
  # predictions
  pred <- lapply(best, function(b) {xtest %*% b})
  predt <- lapply(best, function(b) {xtrain %*% b})
  
  # performance measures
  emse <- sapply(best, function(b) {colMeans((b - beta)^2)})
  
  idh <- apply(-abs(beta), 2, order)[1:floor(0.1*p), ]
  emseh <- sapply(best, function(b) {
    sapply(1:D, function(d) {mean((b[idh[, d], d] - beta[idh[, d], d])^2)})})
    
  idl <- apply(abs(beta), 2, order)[1:floor(0.1*p), ]
  emsel <- sapply(best, function(b) {
    sapply(1:D, function(d) {mean((b[idl[, d], d] - beta[idl[, d], d])^2)})})
  
  pmse <- sapply(pred, function(s) {colMeans((s - Reduce("cbind", ytest))^2)})
  
  pmset <- sapply(predt, function(s) {
    colMeans((s - Reduce("cbind", ytrain))^2)})
  
  # save alpha estimates
  est <- cbind(c(fit.nig1$eb$alphaf, rep(NA, 3), fit.nig1$eb$lambdaf),
               c(fit.nig2$eb$alphaf, fit.nig2$eb$lambdaf),
               rep(NA, 5), rep(NA, 5))
  
  list(emse=emse, emsel=emsel, emseh=emseh, pmse=pmse, pmset=pmset, est=est)
}
stopCluster(cl=cl)

# prepare and save results table
emse <- Reduce("rbind", lapply(res, "[[", "emse"))
emsel <- Reduce("rbind", lapply(res, "[[", "emsel"))
emseh <- Reduce("rbind", lapply(res, "[[", "emseh"))
pmse <- Reduce("rbind", lapply(res, "[[", "pmse"))
pmset <- Reduce("rbind", lapply(res, "[[", "pmset"))
est <- Reduce("rbind", lapply(res, "[[", "est"))
  
res2 <- rbind(emse, emsel, emseh, pmse, pmset, est)
colnames(res2) <- c("nig1", "nig2", "ridge1", "lasso1")
rownames(res2) <- c(paste0(rep(c("emse", "emsel", "emseh", "pmse", "pmset"), 
                               each=D*nreps),
                           rep(rep(paste0(".drug", c(1:D)), nreps), 5)),
                    rep(c(paste0("alphaf", 0:3), "lambdaf"), times=nreps))
write.table(res2, file="results/simulations_gdsc_res1.txt")

################################################################################
################################# simulation 2 #################################
################################################################################
### simulation settings
D <- length(resp)
p <- 100
n <- floor(nrow(expr$expr))/2
ntest <- nrow(expr$expr) - n
alphad <- c(1, 1, 3, 7)
H <- length(alphad)
shape <- 3
rate <- 2
lambdad <- 10
Z <- unname(model.matrix(~ factor(sort(c(rep(1:H, each=D %/% H), 
                                         rep(1:(D %% H), (D %% H)!=0))))))
x <- scale(expr$expr[, order(-apply(expr$expr, 2, sd))[1:p]])

# estimation settings
nfolds <- 10
control <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-5, 
                epsilon.vb=1e-3, maxit.eb=500, maxit.vb=1, 
                maxit.post=100, maxit.block=0)

# setup cluster
if(Sys.info()["sysname"]=="Windows") {
  cl <- makeCluster(min(ncores, nreps), type="PSOCK")
} else {
  cl <- makeCluster(min(ncores, nreps), type="FORK")
}
if(min(ncores, nreps) > 1) {
  registerDoParallel(cl)
} else {
  registerDoSEQ()
}
res <- foreach(r=1:nreps, .packages=packages, .errorhandling="pass") %dopar% {
  cat("\r", "rep", r)
  set.seed(2020 + r)
  
  ### simulate parameters
  tau <- sqrt(rinvgauss(D, 1/as.numeric(Z %*% alphad), lambdad))
  sigma <- 1/sqrt(rgamma(D, shape, rate))
  beta <- sapply(1:D, function(d) {rnorm(p, 0, tau[d]*sigma[d])})
  
  # draw and simulate data
  idtrain <- sample(1:nrow(x), n)
  xtrain <- scale(x[idtrain, ])
  rownames(xtrain) <- c(1:n)
  xtest <- scale(x[-idtrain, ])
  rownames(xtest) <- c(1:ntest)
  ytrain <- lapply(1:D, function(d) {
    s <- as.numeric(scale(rnorm(n, xtrain %*% beta[, d], sigma[d])))
    names(s) <- c(1:length(s))
    return(s)})
  ytest <- lapply(1:D, function(d) {
    s <- as.numeric(scale(rnorm(ntest, xtest %*% beta[, d], sigma[d])))
    names(s) <- c(1:length(s))
    return(s)})
  
  # creating folds
  foldid <- sample(c(rep(1:nfolds, each=n %/% nfolds),
                     rep(1:(n %% nfolds), (n %% nfolds)!=0)))
  
  # fitting models
  fit.nig1 <- nig(x=rep(list(xtrain), D), y=ytrain,
                        C=NULL, Z=matrix(rep(1, D)), full.post=TRUE,
                        control=control)
  fit.nig2 <- nig(x=rep(list(xtrain), D), y=ytrain, C=NULL, Z=Z, 
                        full.post=TRUE, control=control)
  
  # standard penalized methods
  fit.ridge1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain, ytrain[[d]], alpha=0, intercept=FALSE, standardize=TRUE,
              foldid=foldid)})
  fit.lasso1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain, ytrain[[d]], intercept=FALSE, standardize=TRUE,
              foldid=foldid)})
  
  # estimates
  best <- list(nig1=Reduce("cbind", fit.nig1$vb$mpost$beta), 
               nig2=Reduce("cbind", fit.nig2$vb$mpost$beta),
               ridge1=sapply(fit.ridge1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}),
               lasso1=sapply(fit.lasso1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}))
  
  # predictions
  pred <- lapply(best, function(b) {xtest %*% b})
  predt <- lapply(best, function(b) {xtrain %*% b})
  
  # performance measures
  emse <- sapply(best, function(b) {colMeans((b - beta)^2)})
  
  idh <- apply(-abs(beta), 2, order)[1:floor(0.1*p), ]
  emseh <- sapply(best, function(b) {
    sapply(1:D, function(d) {mean((b[idh[, d], d] - beta[idh[, d], d])^2)})})
  
  idl <- apply(abs(beta), 2, order)[1:floor(0.1*p), ]
  emsel <- sapply(best, function(b) {
    sapply(1:D, function(d) {mean((b[idl[, d], d] - beta[idl[, d], d])^2)})})
  
  pmse <- sapply(pred, function(s) {colMeans((s - Reduce("cbind", ytest))^2)})
  
  pmset <- sapply(predt, function(s) {
    colMeans((s - Reduce("cbind", ytrain))^2)})
  
  # save alpha estimates
  est <- cbind(c(fit.nig1$eb$alphad, rep(NA, 3), fit.nig1$eb$lambdad),
               c(fit.nig2$eb$alphad, fit.nig2$eb$lambdad),
               rep(NA, 5), rep(NA, 5))
  
  list(emse=emse, emsel=emsel, emseh=emseh, pmse=pmse, pmset=pmset, est=est)
}
stopCluster(cl=cl)

# prepare and save results table
emse <- Reduce("rbind", lapply(res, "[[", "emse"))
emsel <- Reduce("rbind", lapply(res, "[[", "emsel"))
emseh <- Reduce("rbind", lapply(res, "[[", "emseh"))
pmse <- Reduce("rbind", lapply(res, "[[", "pmse"))
pmset <- Reduce("rbind", lapply(res, "[[", "pmset"))
est <- Reduce("rbind", lapply(res, "[[", "est"))

res2 <- rbind(emse, emsel, emseh, pmse, pmset, est)
colnames(res2) <- c("nig1", "nig2", "ridge1", "lasso1")
rownames(res2) <- c(paste0(rep(c("emse", "emsel", "emseh", "pmse", "pmset"), 
                               each=D*nreps),
                           rep(rep(paste0(".drug", c(1:D)), nreps), 5)),
                    rep(c(paste0("alphad", 0:3), "lambdad"), times=nreps))
write.table(res2, file="results/simulations_gdsc_res2.txt")

################################################################################
################################# simulation 3 #################################
################################################################################
### simulation settings
D <- length(resp)
p <- 100
n <- floor(nrow(expr$expr))/2
ntest <- nrow(expr$expr) - n
alphaf <- c(1, 1, 3, 7)
alphad <- c(1, 1, 3, 7)
G <- length(alphaf)
H <- length(alphad)
shape <- 3
rate <- 2
lambdaf <- 1
lambdad <- 1
C <- replicate(D, list(unname(model.matrix(~ factor(rep(1:G, each=p/G))))))
Z <- unname(model.matrix(~ factor(sort(c(rep(1:H, each=D %/% H), 
                                         rep(1:(D %% H), (D %% H)!=0))))))
x <- scale(expr$expr[, order(-apply(expr$expr, 2, sd))[1:p]])

# estimation settings
nfolds <- 10
control <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-5, 
                epsilon.vb=1e-3, maxit.eb=500, maxit.vb=1, 
                maxit.post=100, maxit.block=0)

# setup cluster
if(Sys.info()["sysname"]=="Windows") {
  cl <- makeCluster(min(ncores, nreps), type="PSOCK")
} else {
  cl <- makeCluster(min(ncores, nreps), type="FORK")
}
if(min(ncores, nreps) > 1) {
  registerDoParallel(cl)
} else {
  registerDoSEQ()
}
res <- foreach(r=1:nreps, .packages=packages) %dopar% {
  cat("\r", "rep", r)
  set.seed(2020 + r)

  ### simulate parameters
  gamma <- sapply(C, function(s) {
    sqrt(rinvgauss(p, 1/as.numeric(s %*% alphaf), lambdaf))})
  tau <- sqrt(rinvgauss(D, 1/as.numeric(Z %*% alphad), lambdad))
  sigma <- 1/sqrt(rgamma(D, shape, rate))
  beta <- sapply(1:D, function(d) {rnorm(p, 0, tau[d]*gamma[, d]*sigma[d])})
  
  # draw and simulate data
  idtrain <- sample(1:nrow(x), n)
  xtrain <- scale(x[idtrain, ])
  rownames(xtrain) <- c(1:n)
  xtest <- scale(x[-idtrain, ])
  rownames(xtest) <- c(1:ntest)
  ytrain <- lapply(1:D, function(d) {
    s <- as.numeric(scale(rnorm(n, xtrain %*% beta[, d], sigma[d])))
    names(s) <- c(1:length(s))
    return(s)})
  ytest <- lapply(1:D, function(d) {
    s <- as.numeric(scale(rnorm(ntest, xtest %*% beta[, d], sigma[d])))
    names(s) <- c(1:length(s))
    return(s)})
  
  # creating folds
  foldid <- sample(c(rep(1:nfolds, each=n %/% nfolds),
                     rep(1:(n %% nfolds), (n %% nfolds)!=0)))
  
  # fitting models
  fit.nig1 <- nig(x=rep(list(xtrain), D), y=ytrain, 
                        C=replicate(D, list(matrix(rep(1, p)))), 
                        Z=matrix(rep(1, D)), full.post=TRUE, 
                        control=control)
  fit.nig2 <- nig(x=rep(list(xtrain), D), y=ytrain, C=C, Z=Z, 
                        full.post=TRUE, control=control)
  
  # standard penalized methods
  fit.ridge1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain, ytrain[[d]], alpha=0, intercept=FALSE, standardize=TRUE,
              foldid=foldid)})
  fit.lasso1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain, ytrain[[d]], intercept=FALSE, standardize=TRUE,
              foldid=foldid)})
  
  # estimates
  best <- list(nig1=Reduce("cbind", fit.nig1$vb$mpost$beta), 
               nig2=Reduce("cbind", fit.nig2$vb$mpost$beta),
               ridge1=sapply(fit.ridge1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}),
               lasso1=sapply(fit.lasso1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}))
  
  # predictions
  pred <- lapply(best, function(b) {xtest %*% b})
  predt <- lapply(best, function(b) {xtrain %*% b})
  
  # performance measures
  emse <- sapply(best, function(b) {colMeans((b - beta)^2)})
  
  idh <- apply(-abs(beta), 2, order)[1:floor(0.1*p), ]
  emseh <- sapply(best, function(b) {
    sapply(1:D, function(d) {mean((b[idh[, d], d] - beta[idh[, d], d])^2)})})
  
  idl <- apply(abs(beta), 2, order)[1:floor(0.1*p), ]
  emsel <- sapply(best, function(b) {
    sapply(1:D, function(d) {mean((b[idl[, d], d] - beta[idl[, d], d])^2)})})
  
  pmse <- sapply(pred, function(s) {colMeans((s - Reduce("cbind", ytest))^2)})
  
  pmset <- sapply(predt, function(s) {
    colMeans((s - Reduce("cbind", ytrain))^2)})
  
  # save alpha estimates
  est <- cbind(c(fit.nig1$eb$alphaf, rep(NA, 3), fit.nig1$eb$alphad, 
                 rep(NA, 3), fit.nig1$eb$lambdaf, fit.nig1$eb$lambdad),
               c(fit.nig2$eb$alphaf, fit.nig2$eb$alphad, 
                 fit.nig2$eb$lambdaf, fit.nig2$eb$lambdad),
               rep(NA, 10), rep(NA, 10))
  
  list(emse=emse, emsel=emsel, emseh=emseh, pmse=pmse, pmset=pmset, est=est)
}
stopCluster(cl=cl)

# prepare and save results table
emse <- Reduce("rbind", lapply(res, "[[", "emse"))
emsel <- Reduce("rbind", lapply(res, "[[", "emsel"))
emseh <- Reduce("rbind", lapply(res, "[[", "emseh"))
pmse <- Reduce("rbind", lapply(res, "[[", "pmse"))
pmset <- Reduce("rbind", lapply(res, "[[", "pmset"))
est <- Reduce("rbind", lapply(res, "[[", "est"))

res2 <- rbind(emse, emsel, emseh, pmse, pmset, est)
colnames(res2) <- c("nig1", "nig2", "ridge1", "lasso1")
rownames(res2) <- c(paste0(rep(c("emse", "emsel", "emseh", "pmse", "pmset"), 
                               each=D*nreps),
                           rep(rep(paste0(".drug", c(1:D)), nreps), 5)),
                    rep(c(paste0("alphaf", 0:3), paste0("alphad", 0:3), 
                          "lambdaf", "lambdad"), times=nreps))
write.table(res2, file="results/simulations_gdsc_res3.txt")

################################################################################
################################# simulation 4 #################################
################################################################################
### simulation settings
D <- length(resp)
p <- 100
n <- floor(nrow(expr$expr))/2
ntest <- nrow(expr$expr) - n
alphaf <- c(1, 1, 3, 7)
alphad <- c(1, 1, 3, 7)
G <- length(alphaf)
H <- length(alphad)
shape <- 3
rate <- 2
lambdaf <- 1
lambdad <- 1
fracs <- c(1/10, 1/5, 1/3, 1/2, 2/3, 4/5, 1) # fraction of ext. data permuted
C <- replicate(D, list(unname(model.matrix(~ factor(rep(1:G, each=p/G))))))
Z <- unname(model.matrix(~ factor(sort(c(rep(1:H, each=D %/% H), 
                                         rep(1:(D %% H), (D %% H)!=0))))))
x <- scale(expr$expr[, order(-apply(expr$expr, 2, sd))[1:p]])

# estimation settings
nfolds <- 10
control <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-5, 
                epsilon.vb=1e-3, maxit.eb=500, maxit.vb=1, 
                maxit.post=100, maxit.block=0)

# setup cluster
if(Sys.info()["sysname"]=="Windows") {
  cl <- makeCluster(min(ncores, nreps), type="PSOCK")
} else {
  cl <- makeCluster(min(ncores, nreps), type="FORK")
}
if(min(ncores, nreps) > 1) {
  registerDoParallel(cl)
} else {
  registerDoSEQ()
}

res <- foreach(q=fracs) %:%
  foreach(r=1:nreps, .packages=packages) %dopar% {
  cat("\r", "rep", r)
  set.seed(2019 + r)
  
  ### simulate parameters
  gamma <- sapply(C, function(s) {
    sqrt(rinvgauss(p, 1/as.numeric(s %*% alphaf), lambdaf))})
  tau <- sqrt(rinvgauss(D, 1/as.numeric(Z %*% alphad), lambdad))
  sigma <- 1/sqrt(rgamma(D, shape, rate))
  beta <- sapply(1:D, function(d) {rnorm(p, 0, tau[d]*gamma[, d]*sigma[d])})
  
  # draw and simulate data
  idtrain <- sample(1:nrow(x), n)
  xtrain <- scale(x[idtrain, ])
  rownames(xtrain) <- c(1:n)
  xtest <- scale(x[-idtrain, ])
  rownames(xtest) <- c(1:ntest)
  ytrain <- lapply(1:D, function(d) {
    s <- as.numeric(scale(rnorm(n, xtrain %*% beta[, d], sigma[d])))
    names(s) <- c(1:length(s))
    return(s)})
  ytest <- lapply(1:D, function(d) {
    s <- as.numeric(scale(rnorm(ntest, xtest %*% beta[, d], sigma[d])))
    names(s) <- c(1:length(s))
    return(s)})
  
  # creating folds
  foldid <- sample(c(rep(1:nfolds, each=n %/% nfolds),
                     rep(1:(n %% nfolds), (n %% nfolds)!=0)))
  
  # randomly permuting q proportion of rows in external data
  Ctilde <- lapply(C, function(s) {
    ids <- sample(1:p, p*q)
    s[ids, ] <- s[ids, ][sample(1:(q*p)), ]
    s})
  Ztilde <- Z
  ids <- ids <- sample(1:D, D*q)
  Ztilde[ids, ] <- Ztilde[ids, ][sample(1:(q*D)), ]
  
  # fitting models
  fit.nig1 <- nig(x=rep(list(xtrain), D), y=ytrain,
                  C=replicate(D, list(matrix(rep(1, p)))),
                  Z=matrix(rep(1, D)), full.post=TRUE, control=control)
  fit.nig2 <- nig(x=rep(list(xtrain), D), y=ytrain, C=Ctilde, Z=Ztilde,
                  full.post=TRUE, control=control)

  # standard penalized methods
  fit.ridge1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain, ytrain[[d]], alpha=0, intercept=FALSE, standardize=TRUE,
              foldid=foldid)})
  fit.lasso1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain, ytrain[[d]], intercept=FALSE, standardize=TRUE,
              foldid=foldid)})

  # estimates
  best <- list(nig1=Reduce("cbind", fit.nig1$vb$mpost$beta),
               nig2=Reduce("cbind", fit.nig2$vb$mpost$beta),
               ridge1=sapply(fit.ridge1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}),
               lasso1=sapply(fit.lasso1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}))

  # predictions
  pred <- lapply(best, function(b) {xtest %*% b})
  predt <- lapply(best, function(b) {xtrain %*% b})

  # performance measures
  emse <- sapply(best, function(b) {colMeans((b - beta)^2)})

  idh <- apply(-abs(beta), 2, order)[1:floor(0.1*p), ]
  emseh <- sapply(best, function(b) {
    sapply(1:D, function(d) {mean((b[idh[, d], d] - beta[idh[, d], d])^2)})})

  idl <- apply(abs(beta), 2, order)[1:floor(0.1*p), ]
  emsel <- sapply(best, function(b) {
    sapply(1:D, function(d) {mean((b[idl[, d], d] - beta[idl[, d], d])^2)})})

  pmse <- sapply(pred, function(s) {colMeans((s - Reduce("cbind", ytest))^2)})

  pmset <- sapply(predt, function(s) {
    colMeans((s - Reduce("cbind", ytrain))^2)})

  # save alpha estimates
  est <- cbind(c(fit.nig1$eb$alphaf, rep(NA, 3), fit.nig1$eb$alphad,
                 rep(NA, 3), fit.nig1$eb$lambdaf, fit.nig1$eb$lambdad),
               c(fit.nig2$eb$alphaf, fit.nig2$eb$alphad,
                 fit.nig2$eb$lambdaf, fit.nig2$eb$lambdad),
               rep(NA, 10), rep(NA, 10))

  list(emse=emse, emsel=emsel, emseh=emseh, pmse=pmse, pmset=pmset, est=est,
       q=q)
}
stopCluster(cl=cl)

# prepare and save results table
res2.1 <- Reduce("rbind", lapply(res, function(r) {
  emse <- Reduce("rbind", lapply(r, "[[", "emse"))
  emsel <- Reduce("rbind", lapply(r, "[[", "emsel"))
  emseh <- Reduce("rbind", lapply(r, "[[", "emseh"))
  rbind(emse, emsel, emseh)}))
res2.2 <- Reduce("rbind", lapply(res, function(r) {
  pmse <- Reduce("rbind", lapply(r, "[[", "pmse"))
  pmset <- Reduce("rbind", lapply(r, "[[", "pmset"))
  est <- Reduce("rbind", lapply(r, "[[", "est"))
  rbind(pmse, pmset, est)}))
colnames(res2.1) <- colnames(res2.2) <- c("nig1", "nig2", "ridge1", "lasso1")

rownames(res2.1) <- 
  paste0(rep(paste0("frac", format(round(fracs, 2))), each=3*D*nreps), ".",
         paste0(rep(c("emse", "emsel", "emseh"), each=D*nreps),
                rep(rep(paste0(".drug", c(1:D)), nreps), 3)))
rownames(res2.2) <- 
  paste0(rep(paste0("frac", format(round(fracs, 2))), 
             each=2*D*nreps + 10*nreps), ".",
         c(paste0(rep(c("pmse", "pmset"), each=D*nreps),
                  rep(rep(paste0(".drug", c(1:D)), nreps), 2)),
           rep(c(paste0("alphaf", 0:3), paste0("alphad", 0:3), 
                 "lambdaf", "lambdad"), times=nreps)))
write.table(res2.1, file="results/simulations_gdsc_res4.1.txt")
write.table(res2.2, file="results/simulations_gdsc_res4.2.txt")

################################################################################
################################# simulation 5 #################################
################################################################################
### simulation settings
D <- length(resp)
p <- 100
n <- floor(nrow(expr$expr))/2
ntest <- nrow(expr$expr) - n
alphaf <- c(1, 1, 3, 7)
alphad <- c(1, 1, 3, 7)
G <- length(alphaf)
H <- length(alphad)
shape <- 3
rate <- 2
lambdaf <- 1
lambdad <- 1
C <- replicate(D, list(unname(model.matrix(~ factor(rep(1:G, each=p/G))))))
Z <- unname(model.matrix(~ factor(sort(c(rep(1:H, each=D %/% H), 
                                         rep(1:(D %% H), (D %% H)!=0))))))

# estimation settings
methods <- c("NIG", "MCMC")
control.nig <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-5, 
                       epsilon.vb=1e-3, maxit.eb=2000, maxit.vb=1, 
                       maxit.post=100, maxit.block=0)
control.stan <- list(verbose=FALSE, show_messages=FALSE, iter=3000, 
                     warmup=1000, chains=1)

### simulate parameters
set.seed(2020)
gamma <- sapply(C, function(s) {
  sqrt(rinvgauss(p, 1/as.numeric(s %*% alphaf), lambdaf))})
tau <- sqrt(rinvgauss(D, 1/as.numeric(Z %*% alphad), lambdad))
sigma <- 1/sqrt(rgamma(D, shape, rate))
beta <- sapply(1:D, function(d) {rnorm(p, 0, tau[d]*gamma[, d]*sigma[d])})

# simulate data
id <- sample(1:nrow(expr$expr), n)
x <- scale(expr$expr[, order(-apply(expr$expr, 2, sd))[1:p]])[id, ]
y <- lapply(1:D, function(d) {
  s <- as.numeric(scale(rnorm(n, x %*% beta[, d], sigma[d])))
  names(s) <- c(1:length(s))
  return(s)})

# fitting models
ct <- proc.time()[3]
fit.nig1 <- nig(x=rep(list(x), D), y=y, C=C, Z=Z, 
                      full.post=TRUE, control=control.nig)
time.nig1 <- proc.time()[3] - ct
time.nig1.post <- fit.nig1$time[2]
ct <- proc.time()[3]
fit.mcmc1 <- lapply(1:D, function(d) {
  sampling(stanmodels$nig, 
           data=list(p=p, n=n, x=x, y=y[[d]], 
                     phi=1/as.numeric(C[[d]] %*% fit.nig1$eb$alphaf),
                     chi=1/as.numeric(Z[d, ] %*% fit.nig1$eb$alphad),
                     lambdaf=fit.nig1$eb$lambdaf, 
                     lambdad=fit.nig1$eb$lambdad), 
           iter=control.stan$iter, warmup=control.stan$warmup, 
           chains=control.stan$chains,
           verbose=control.stan$verbose, 
           show_messages=control.stan$show_messages)})
time.mcmc1 <- proc.time()[3] - ct
ct <- proc.time()[3]
fit.mcmc2 <- sampling(stanmodels$nig_full,
                      data=list(D=D, p=rep(p, D), sump=p*D, n=n, G=G, H=H,
                                x=Reduce("cbind", replicate(D, list(x))),
                                y=Reduce("cbind", y), C=Reduce("rbind", C),
                                Z=Z, nuf=10, nud=10, kappaf=1, 
                                kappad=1, xif=1, xid=1),
                      iter=control.stan$iter, warmup=control.stan$warmup, 
                      chains=control.stan$chains,
                      verbose=control.stan$verbose, 
                      show_messages=control.stan$show_messages)
time.mcmc2 <- proc.time()[3] - ct
save(fit.nig1, fit.mcmc1, fit.mcmc2, 
     file="results/simulations_gdsc_fit5.Rdata")
hid <- sapply(1:H, function(h) {which((Z[, -1] %*% c(1:(H - 1)) + 1)==h)[1]})
gid <- p/G*(c(1:G) - 1) + 1
post.nig1 <- list(mu=fit.nig1$vb$mu[hid], Sigma=fit.nig1$vb$Sigma[hid])
save(post.nig1, file="results/simulations_gdsc_fit5_nig1.Rdata")
post.mcmc1 <- lapply(fit.mcmc1[hid], extract, pars=paste0("beta[", gid, "]"))
post.mcmc2 <- extract(fit.mcmc2, pars=c("beta", "alphaf", "alphad", 
                                        "lambdaf", "lambdad"))
post.mcmc2 <- list(beta=lapply(split(
  as.data.frame(t(post.mcmc2$beta)), rep(1:D, each=p))[hid],
  function(s) {unname(as.matrix(s)[gid, ])}),
  alphaf=t(post.mcmc2$alphaf), alphad=t(post.mcmc2$alphad),
  lambdaf=t(post.mcmc2$lambdaf), lambdad=t(post.mcmc2$lambdad))
summary.mcmc2 <- summary(fit.mcmc2)
save(time.nig1, time.nig1.post, time.mcmc1, time.mcmc2, post.mcmc1, post.mcmc2, 
     post.nig1, summary.mcmc2, file="results/simulations_gdsc_res5.Rdata")

