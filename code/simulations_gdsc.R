#!/usr/bin/env Rscript

# number of cores to use
ncores <- 100

### libraries
packages <- c("foreach", "doParallel", "cambridge", "statmod", "glmnet", 
              "xtune", "rstan")
sapply(packages, library, character.only=TRUE)

### load and preprocess data
load(file="data/data_gdsc_dat1.Rdata")

# number of reps
nreps <- 100

################################################################################
################################# simulation 1 #################################
################################################################################
### simulation settings
D <- 100
p <- 100
n <- 100
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
methods <- c("NIG-", "NIG", "ridge", "lasso", "xtune", "ebridge")
control.semnig <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-5, 
                       epsilon.vb=1e-3, maxit.eb=2000, maxit.vb=1, 
                       maxit.post=100, maxit.block=0)
control.ebridge <-list(epsilon=sqrt(.Machine$double.eps), maxit=500, 
                       trace=FALSE, glmnet.fit2=FALSE, beta2=FALSE)

# setup cluster
cl <- makeCluster(ncores)
if(ncores > 1) {
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
  fit.semnig1 <- semnig(x=rep(list(xtrain), D), y=ytrain, 
                        C=replicate(D, list(matrix(rep(1, p)))), 
                        Z=NULL, full.post=TRUE, control=control.semnig)
  fit.semnig2 <- semnig(x=rep(list(xtrain), D), y=ytrain, C=C, Z=NULL, 
                        full.post=TRUE, control=control.semnig)

  # standard penalized methods
  fit.lasso1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain, ytrain[[d]], intercept=FALSE, standardize=TRUE,
              foldid=foldid)})
  fit.ridge1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain, ytrain[[d]], alpha=0, intercept=FALSE, standardize=TRUE,
              foldid=foldid)})
  
  # empirical bayes penalized methods
  fit.xtune1 <- lapply(1:D, function(d) {
    xtune(xtrain, ytrain[[d]], C[[d]][, -1], family="linear", method="ridge", 
          message=FALSE, control=list(intercept=FALSE))})
  fit.ebridge1 <- ebridge(xtrain, ytrain, lapply(C, function(s) {s[, -1]}), 
                          NULL, foldid=rep(list(foldid), D), 
                          hyper=list(lambda=1/sqrt(n*sapply(
                            fit.ridge1, "[[", "lambda.min")), zeta=0, nu=0),
                          control=control.ebridge)
  
  # estimates
  best <- list(semnig1=Reduce("cbind", fit.semnig1$vb$mpost$beta), 
               semnig2=Reduce("cbind", fit.semnig2$vb$mpost$beta),
               ridge1=sapply(fit.ridge1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}),
               lasso1=sapply(fit.lasso1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}),
               xtune1=sapply(fit.xtune1, function(s) {
                 unname(s$beta.est[-1, ])}),
               ebridge1=Reduce("cbind", fit.ebridge1$beta1))
  
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
  est <- cbind(c(fit.semnig1$eb$alphaf, rep(NA, 3), fit.semnig1$eb$lambdaf),
               c(fit.semnig2$eb$alphaf, fit.semnig2$eb$lambdaf),
               rep(NA, 5), rep(NA, 5),
               c(colMeans(-t(sapply(fit.xtune1, "[[", "alpha.est"))), NA),
               c(NA, fit.ebridge1$alphaf, NA))

  # calculate ELBO and conditional predictive ordinate for semnig models
  elbot <- cbind(fit.semnig1$seq.elbo[nrow(fit.semnig1$seq.elbo), ],
                 fit.semnig2$seq.elbo[nrow(fit.semnig2$seq.elbo), ],
                 rep(NA, D), rep(NA, D), rep(NA, D), rep(NA, D))

  elbo <- cbind(new.elbo(fit.semnig1, replicate(D, list(xtest)), 
                         Reduce("cbind", ytest)),
                new.elbo(fit.semnig1, replicate(D, list(xtest)), 
                         Reduce("cbind", ytest)), rep(NA, D), rep(NA, D), 
                rep(NA, D), rep(NA, D))

  lpml <- cbind(colMeans(logcpo(replicate(D, list(xtest)), 
                                Reduce("cbind", ytest), n, fit.semnig1)),
                colMeans(logcpo(replicate(D, list(xtest)), 
                                Reduce("cbind", ytest), n, fit.semnig2)), 
                rep(NA, D), rep(NA, D), rep(NA, D), rep(NA, D))
  
  list(emse=emse, emsel=emsel, emseh=emseh, pmse=pmse, pmset=pmset, est=est, 
       elbo=elbo, elbot=elbot, lpml=lpml)
}
stopCluster(cl=cl)

# prepare and save results table
emse <- Reduce("rbind", lapply(res, "[[", "emse"))
emsel <- Reduce("rbind", lapply(res, "[[", "emsel"))
emseh <- Reduce("rbind", lapply(res, "[[", "emseh"))
pmse <- Reduce("rbind", lapply(res, "[[", "pmse"))
pmset <- Reduce("rbind", lapply(res, "[[", "pmset"))
elbo <- Reduce("rbind", lapply(res, "[[", "elbo"))
elbot <- Reduce("rbind", lapply(res, "[[", "elbot"))
lpml <- Reduce("rbind", lapply(res, "[[", "lpml"))
est <- Reduce("rbind", lapply(res, "[[", "est"))
  

res2 <- rbind(emse, emsel, emseh, pmse, pmset, elbo, elbot, lpml, est)
colnames(res2) <- c(methods)
rownames(res2) <- c(paste0(rep(c("emse", "emsel", "emseh", "pmse", "pmset", 
                                 "elbo", "elbot", "lpml"), each=D*nreps),
                           rep(rep(paste0(".drug", c(1:D)), nreps), 8)),
                    rep(c(paste0("alphaf", 0:3), "lambdaf"), times=nreps))
write.table(res2, file="results/simulations_gdsc_res1.txt")

################################################################################
################################# simulation 2 #################################
################################################################################
### simulation settings
D <- 100
p <- 100
n <- 100
ntest <- nrow(expr$expr) - n
alphad <- c(1, 1, 3, 7)
H <- length(alphad)
shape <- 3
rate <- 2
lambdad <- 10
Z <- unname(model.matrix(~ factor(rep(1:H, each=D/H))))
x <- scale(expr$expr[, order(-apply(expr$expr, 2, sd))[1:p]])

# estimation settings
nfolds <- 10
methods <- c("NIG-", "NIG", "ridge", "lasso", "xtune")
control.semnig <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-5, 
                       epsilon.vb=1e-3, maxit.eb=2000, maxit.vb=1, 
                       maxit.post=100, maxit.block=0)
control.ebridge <-list(epsilon=sqrt(.Machine$double.eps), maxit=500, 
                       trace=FALSE, glmnet.fit2=FALSE, beta2=FALSE)

# setup cluster
cl <- makeCluster(ncores)
if(ncores > 1) {
  registerDoParallel(cl)
} else {
  registerDoSEQ()
}
res <- foreach(r=1:nreps, .packages=packages) %dopar% {
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
  fit.semnig1 <- semnig(x=rep(list(xtrain), D), y=ytrain,
                        C=NULL, Z=matrix(rep(1, D)), full.post=TRUE,
                        control=control.semnig)
  fit.semnig2 <- semnig(x=rep(list(xtrain), D), y=ytrain, C=NULL, Z=Z, 
                        full.post=TRUE, control=control.semnig)
  
  # standard penalized methods
  fit.lasso1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain, ytrain[[d]], intercept=FALSE, standardize=TRUE,
              foldid=foldid)})
  fit.ridge1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain, ytrain[[d]], alpha=0, intercept=FALSE, standardize=TRUE,
              foldid=foldid)})
  
  # empirical bayes penalized methods
  fit.xtune1 <- lapply(1:D, function(d) {
    xtune(xtrain, ytrain[[d]], NULL, family="linear", method="ridge", 
          message=FALSE, control=list(intercept=FALSE))})
  
  # estimates
  best <- list(semnig1=Reduce("cbind", fit.semnig1$vb$mpost$beta), 
               semnig2=Reduce("cbind", fit.semnig2$vb$mpost$beta),
               ridge1=sapply(fit.ridge1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}),
               lasso1=sapply(fit.lasso1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}),
               xtune1=sapply(fit.xtune1, function(s) {
                 unname(s$beta.est[-1, ])}))
  
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
  est <- cbind(c(fit.semnig1$eb$alphad, rep(NA, 3), fit.semnig1$eb$lambdad),
               c(fit.semnig2$eb$alphad, fit.semnig2$eb$lambdad),
               rep(NA, 5), rep(NA, 5),
               c(mean(-t(sapply(fit.xtune1, "[[", "alpha.est"))), 
                 rep(NA, 4)))
  
  # calculate ELBO and conditional predictive ordinate for semnig models
  elbot <- cbind(fit.semnig1$seq.elbo[nrow(fit.semnig1$seq.elbo), ],
                 fit.semnig2$seq.elbo[nrow(fit.semnig2$seq.elbo), ],
                 rep(NA, D), rep(NA, D), rep(NA, D))
  
  elbo <- cbind(new.elbo(fit.semnig1, replicate(D, list(xtest)), 
                         Reduce("cbind", ytest)),
                new.elbo(fit.semnig1, replicate(D, list(xtest)), 
                         Reduce("cbind", ytest)), rep(NA, D), rep(NA, D), 
                rep(NA, D))
  
  lpml <- cbind(colMeans(logcpo(replicate(D, list(xtest)), 
                                Reduce("cbind", ytest), n, fit.semnig1)),
                colMeans(logcpo(replicate(D, list(xtest)), 
                                Reduce("cbind", ytest), n, fit.semnig2)), 
                rep(NA, D), rep(NA, D), rep(NA, D))
  
  list(emse=emse, emsel=emsel, emseh=emseh, pmse=pmse, pmset=pmset, est=est, 
       elbo=elbo, elbot=elbot, lpml=lpml)
}
stopCluster(cl=cl)

# prepare and save results table
emse <- Reduce("rbind", lapply(res, "[[", "emse"))
emsel <- Reduce("rbind", lapply(res, "[[", "emsel"))
emseh <- Reduce("rbind", lapply(res, "[[", "emseh"))
pmse <- Reduce("rbind", lapply(res, "[[", "pmse"))
pmset <- Reduce("rbind", lapply(res, "[[", "pmset"))
elbo <- Reduce("rbind", lapply(res, "[[", "elbo"))
elbot <- Reduce("rbind", lapply(res, "[[", "elbot"))
lpml <- Reduce("rbind", lapply(res, "[[", "lpml"))
est <- Reduce("rbind", lapply(res, "[[", "est"))

res2 <- rbind(emse, emsel, emseh, pmse, pmset, elbo, elbot, lpml, est)
colnames(res2) <- c(methods)
rownames(res2) <- c(paste0(rep(c("emse", "emsel", "emseh", "pmse", "pmset", 
                                 "elbo", "elbot", "lpml"), each=D*nreps),
                           rep(rep(paste0(".drug", c(1:D)), nreps), 8)),
                    rep(c(paste0("alphad", 0:3), "lambdad"), times=nreps))
write.table(res2, file="results/simulations_gdsc_res2.txt")

################################################################################
################################# simulation 3 #################################
################################################################################
### simulation settings
D <- 100
p <- 100
n <- 100
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
Z <- unname(model.matrix(~ factor(rep(1:H, each=D/H))))
x <- scale(expr$expr[, order(-apply(expr$expr, 2, sd))[1:p]])

# estimation settings
nfolds <- 10
methods <- c("NIG-", "NIG", "ridge", "lasso", "xtune", "ebridge")
control.semnig <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-5, 
                       epsilon.vb=1e-3, maxit.eb=1, maxit.vb=1, 
                       maxit.post=100, maxit.block=0)
control.ebridge <-list(epsilon=sqrt(.Machine$double.eps), maxit=500, 
                       trace=FALSE, glmnet.fit2=FALSE, beta2=FALSE)
control.stan <- list(verbose=FALSE, show_messages=FALSE)

# setup cluster
cl <- makeCluster(ncores)
if(ncores > 1) {
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
  fit.semnig1 <- semnig(x=rep(list(xtrain), D), y=ytrain, 
                        C=replicate(D, list(matrix(rep(1, p)))), 
                        Z=matrix(rep(1, D)), full.post=TRUE, 
                        control=control.semnig)
  fit.semnig2 <- semnig(x=rep(list(xtrain), D), y=ytrain, C=C, Z=Z, 
                        full.post=TRUE, control=control.semnig)
  
  # standard penalized methods
  fit.lasso1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain, ytrain[[d]], intercept=FALSE, standardize=TRUE,
              foldid=foldid)})
  fit.ridge1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain, ytrain[[d]], alpha=0, intercept=FALSE, standardize=TRUE,
              foldid=foldid)})
  
  # empirical bayes penalized methods
  fit.xtune1 <- lapply(1:D, function(d) {
    xtune(xtrain, ytrain[[d]], C[[d]][, -1], family="linear", method="ridge", 
          message=FALSE, control=list(intercept=FALSE))})
  fit.ebridge1 <- ebridge(xtrain, ytrain, lapply(C, function(s) {s[, -1]}), 
                          Z[, -1], foldid=rep(list(foldid), D), 
                          hyper=list(lambda=1/sqrt(n*sapply(
                            fit.ridge1, "[[", "lambda.min")), zeta=0, nu=0),
                          control=control.ebridge)
  
  # estimates
  best <- list(semnig1=Reduce("cbind", fit.semnig1$vb$mpost$beta), 
               semnig2=Reduce("cbind", fit.semnig2$vb$mpost$beta),
               ridge1=sapply(fit.ridge1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}),
               lasso1=sapply(fit.lasso1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}),
               xtune1=sapply(fit.xtune1, function(s) {
                 unname(s$beta.est[-1, ])}),
               ebridge1=Reduce("cbind", fit.ebridge1$beta1))
  
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
  est <- cbind(c(fit.semnig1$eb$alphaf, rep(NA, 3), fit.semnig1$eb$alphad, 
                 rep(NA, 3), fit.semnig1$eb$lambdaf, fit.semnig1$eb$lambdad),
               c(fit.semnig2$eb$alphaf, fit.semnig2$eb$alphad, 
                 fit.semnig2$eb$lambdaf, fit.semnig2$eb$lambdad),
               rep(NA, 10), rep(NA, 10),
               c(colMeans(-t(sapply(fit.xtune1, "[[", "alpha.est"))), 
                 rep(NA, 6)),
               c(NA, fit.ebridge1$alphaf, rep(NA, 6)))
  
  # calculate ELBO and conditional predictive ordinate for semnig models
  elbot <- cbind(fit.semnig1$seq.elbo[nrow(fit.semnig1$seq.elbo), ],
                 fit.semnig2$seq.elbo[nrow(fit.semnig2$seq.elbo), ],
                 rep(NA, D), rep(NA, D), rep(NA, D), rep(NA, D))
  
  elbo <- cbind(new.elbo(fit.semnig1, replicate(D, list(xtest)), 
                         Reduce("cbind", ytest)),
                new.elbo(fit.semnig1, replicate(D, list(xtest)), 
                         Reduce("cbind", ytest)), rep(NA, D), rep(NA, D), 
                rep(NA, D), rep(NA, D))
  
  lpml <- cbind(colMeans(logcpo(replicate(D, list(xtest)), 
                                Reduce("cbind", ytest), n, fit.semnig1)),
                colMeans(logcpo(replicate(D, list(xtest)), 
                                Reduce("cbind", ytest), n, fit.semnig2)), 
                rep(NA, D), rep(NA, D), rep(NA, D), rep(NA, D))
  
  list(emse=emse, emsel=emsel, emseh=emseh, pmse=pmse, pmset=pmset, est=est, 
       elbo=elbo, elbot=elbot, lpml=lpml)
}
stopCluster(cl=cl)

# prepare and save results table
emse <- Reduce("rbind", lapply(res, "[[", "emse"))
emsel <- Reduce("rbind", lapply(res, "[[", "emsel"))
emseh <- Reduce("rbind", lapply(res, "[[", "emseh"))
pmse <- Reduce("rbind", lapply(res, "[[", "pmse"))
pmset <- Reduce("rbind", lapply(res, "[[", "pmset"))
elbo <- Reduce("rbind", lapply(res, "[[", "elbo"))
elbot <- Reduce("rbind", lapply(res, "[[", "elbot"))
lpml <- Reduce("rbind", lapply(res, "[[", "lpml"))
est <- Reduce("rbind", lapply(res, "[[", "est"))


res2 <- rbind(emse, emsel, emseh, pmse, pmset, elbo, elbot, lpml, est)
colnames(res2) <- c(methods)
rownames(res2) <- c(paste0(rep(c("emse", "emsel", "emseh", "pmse", "pmset", 
                                 "elbo", "elbot", "lpml"), each=D*nreps),
                           rep(rep(paste0(".drug", c(1:D)), nreps), 8)),
                    rep(c(paste0("alphaf", 0:3), paste0("alphad", 0:3), 
                          "lambdaf", "lambdad"), times=nreps))
write.table(res2, file="results/simulations_gdsc_res3.txt")

################################################################################
################################# simulation 4 #################################
################################################################################
### simulation settings
D <- 100
p <- 100
n <- 100
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
Z <- unname(model.matrix(~ factor(rep(1:H, each=D/H))))
x <- scale(expr$expr[, order(-apply(expr$expr, 2, sd))[1:p]])

# estimation settings
nfolds <- 10
methods <- c("NIG-", "NIG", "ridge", "lasso", "xtune", "ebridge")
control.semnig <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-5, 
                       epsilon.vb=1e-3, maxit.eb=2000, maxit.vb=1, 
                       maxit.post=100, maxit.block=0)
control.ebridge <-list(epsilon=sqrt(.Machine$double.eps), maxit=500, 
                       trace=FALSE, glmnet.fit2=FALSE, beta2=FALSE)

# setup cluster
cl <- makeCluster(ncores)
if(ncores > 1) {
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
  fit.semnig1 <- semnig(x=rep(list(xtrain), D), y=ytrain,
                        C=replicate(D, list(matrix(rep(1, p)))),
                        Z=matrix(rep(1, D)), full.post=TRUE,
                        control=control.semnig)
  fit.semnig2 <- semnig(x=rep(list(xtrain), D), y=ytrain, C=Ctilde, Z=Ztilde,
                        full.post=TRUE, control=control.semnig)

  # standard penalized methods
  fit.lasso1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain, ytrain[[d]], intercept=FALSE, standardize=TRUE,
              foldid=foldid)})
  fit.ridge1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain, ytrain[[d]], alpha=0, intercept=FALSE, standardize=TRUE,
              foldid=foldid)})

  # empirical bayes penalized methods
  fit.xtune1 <- lapply(1:D, function(d) {
    xtune(xtrain, ytrain[[d]], C[[d]][, -1], family="linear", method="ridge",
          message=FALSE, control=list(intercept=FALSE))})
  fit.ebridge1 <- ebridge(xtrain, ytrain, lapply(C, function(s) {s[, -1]}),
                          Z[, -1], foldid=rep(list(foldid), D),
                          hyper=list(lambda=1/sqrt(n*sapply(
                            fit.ridge1, "[[", "lambda.min")), zeta=0, nu=0),
                          control=control.ebridge)

  # estimates
  best <- list(semnig1=Reduce("cbind", fit.semnig1$vb$mpost$beta),
               semnig2=Reduce("cbind", fit.semnig2$vb$mpost$beta),
               ridge1=sapply(fit.ridge1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}),
               lasso1=sapply(fit.lasso1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}),
               xtune1=sapply(fit.xtune1, function(s) {
                 unname(s$beta.est[-1, ])}),
               ebridge1=Reduce("cbind", fit.ebridge1$beta1))

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
  est <- cbind(c(fit.semnig1$eb$alphaf, rep(NA, 3), fit.semnig1$eb$alphad,
                 rep(NA, 3), fit.semnig1$eb$lambdaf, fit.semnig1$eb$lambdad),
               c(fit.semnig2$eb$alphaf, fit.semnig2$eb$alphad,
                 fit.semnig2$eb$lambdaf, fit.semnig2$eb$lambdad),
               rep(NA, 10), rep(NA, 10),
               c(colMeans(-t(sapply(fit.xtune1, "[[", "alpha.est"))),
                 rep(NA, 6)),
               c(NA, fit.ebridge1$alphaf, rep(NA, 6)))

  # calculate ELBO and conditional predictive ordinate for semnig models
  elbot <- cbind(fit.semnig1$seq.elbo[nrow(fit.semnig1$seq.elbo), ],
                 fit.semnig2$seq.elbo[nrow(fit.semnig2$seq.elbo), ],
                 rep(NA, D), rep(NA, D), rep(NA, D), rep(NA, D))

  elbo <- cbind(new.elbo(fit.semnig1, replicate(D, list(xtest)),
                         Reduce("cbind", ytest)),
                new.elbo(fit.semnig1, replicate(D, list(xtest)),
                         Reduce("cbind", ytest)), rep(NA, D), rep(NA, D),
                rep(NA, D), rep(NA, D))

  lpml <- cbind(colMeans(logcpo(replicate(D, list(xtest)),
                                Reduce("cbind", ytest), n, fit.semnig1)),
                colMeans(logcpo(replicate(D, list(xtest)),
                                Reduce("cbind", ytest), n, fit.semnig2)),
                rep(NA, D), rep(NA, D), rep(NA, D), rep(NA, D))

  list(emse=emse, emsel=emsel, emseh=emseh, pmse=pmse, pmset=pmset, est=est,
       elbo=elbo, elbot=elbot, lpml=lpml, q=q)
}
stopCluster(cl=cl)

# prepare and save results table
res2 <- Reduce("rbind", lapply(res, function(r) {
  emse <- Reduce("rbind", lapply(r, "[[", "emse"))
  emsel <- Reduce("rbind", lapply(r, "[[", "emsel"))
  emseh <- Reduce("rbind", lapply(r, "[[", "emseh"))
  pmse <- Reduce("rbind", lapply(r, "[[", "pmse"))
  pmset <- Reduce("rbind", lapply(r, "[[", "pmset"))
  elbo <- Reduce("rbind", lapply(r, "[[", "elbo"))
  elbot <- Reduce("rbind", lapply(r, "[[", "elbot"))
  lpml <- Reduce("rbind", lapply(r, "[[", "lpml"))
  est <- Reduce("rbind", lapply(r, "[[", "est"))
  rbind(emse, emsel, emseh, pmse, pmset, elbo, elbot, lpml, est)}))

colnames(res2) <- c(methods)
rownames(res2) <- 
  paste0(rep(paste0("frac", fracs), each=8*D*nreps + 10*nreps), ".",
         c(paste0(rep(c("emse", "emsel", "emseh", "pmse", "pmset", 
                        "elbo", "elbot", "lpml"), each=D*nreps),
                  rep(rep(paste0(".drug", c(1:D)), nreps), 8)),
           rep(c(paste0("alphaf", 0:3), paste0("alphad", 0:3), 
                 "lambdaf", "lambdad"), times=nreps)))
write.table(res2, file="results/simulations_gdsc_res4.txt")

################################################################################
################################# simulation 5 #################################
################################################################################
### simulation settings
D <- 100
p <- 100
n <- 100
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
Z <- unname(model.matrix(~ factor(rep(1:H, each=D/H))))

# estimation settings
methods <- c("NIG", "MCMC")
control.semnig <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-5, 
                       epsilon.vb=1e-3, maxit.eb=2000, maxit.vb=1, 
                       maxit.post=100, maxit.block=0)
control.stan <- list(verbose=FALSE, show_messages=FALSE, iter=3000, 
                     warmup=1000)

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
fit.semnig1 <- semnig(x=rep(list(x), D), y=y, C=C, Z=Z, 
                      full.post=TRUE, control=control.semnig)
time.semnig1 <- proc.time()[3] - ct
ct <- proc.time()[3]
fit.mcmc1 <- lapply(1:D, function(d) {
  sampling(stanmodels$nig, 
           data=list(p=p, n=n, x=x, y=y[[d]], 
                     phi=1/as.numeric(C[[d]] %*% fit.semnig1$eb$alphaf),
                     chi=1/as.numeric(Z[d, ] %*% fit.semnig1$eb$alphad),
                     lambdaf=fit.semnig1$eb$lambdaf, 
                     lambdad=fit.semnig1$eb$lambdad), 
           iter=control.stan$iter, warmup=control.stan$warmup,
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
                      chains=1, iter=2000, warmup=1000)
time.mcmc2 <- proc.time()[3] - ct

save(fit.semnig1, fit.mcmc1, fit.mcmc2, 
     file="results/simulations_gdsc_fit5.Rdata")
load(file="results/simulations_gdsc_fit5.Rdata")
library(rstan)
get_stancode(fit.mcmc2)

post.mcmc1 <- lapply(fit.mcmc1[D/H*(c(1:H) - 1) + 1], extract, 
                     pars=paste0("beta[", p/G*(c(1:G) - 1) + 1, "]"))
post.mcmc2 <- extract(fit.mcmc2, pars=c("beta", "alphaf", "alphad", 
                                        "lambdaf", "lambdad"))
post.mcmc2 <- list(beta=lapply(split(as.data.frame(t(post.mcmc2$beta)), 
                                     rep(1:D, each=p))[D/H*(c(1:H) - 1) + 1],
                               function(s) {
                                 unname(as.matrix(s)[p/G*(c(1:G) - 1) + 1, ])}),
                   alphaf=t(post.mcmc2$alphaf), alphad=t(post.mcmc2$alphad),
                   lambdaf=t(post.mcmc2$lambdaf), lambdad=t(post.mcmc2$lambdad))
save(time.semnig1, time.mcmc1, time.mcmc2, post.mcmc1, post.mcmc2, post.semnig1, 
     file="results/simulations_gdsc_res5.Rdata")

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(c(rep(rep(c(1:2), each=2), 2), rep(rep(c(3:4), each=2), 2)), 
              nrow=4, ncol=4, byrow=TRUE))
hist(post.mcmc2$alphaf[1, ])
hist(post.mcmc2$alphaf[2, ])
hist(post.mcmc2$alphaf[3, ])
hist(post.mcmc2$alphaf[4, ])
par(opar)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(c(rep(rep(c(1:2), each=2), 2), rep(rep(c(3:4), each=2), 2)), 
              nrow=4, ncol=4, byrow=TRUE))
hist(1/(t(post.mcmc2$alphaf) %*% t(cbind(1, rbind(0, diag(3)))))[, 1])
hist(1/(t(post.mcmc2$alphaf) %*% t(cbind(1, rbind(0, diag(3)))))[, 2])
hist(1/(t(post.mcmc2$alphaf) %*% t(cbind(1, rbind(0, diag(3)))))[, 3])
hist(1/(t(post.mcmc2$alphaf) %*% t(cbind(1, rbind(0, diag(3)))))[, 4])
par(opar)

rowMeans(post.mcmc2$alphaf)
alphaf
test <- stan_model("rpackage/inst/stan/nig_full.stan")
1/as.numeric(alphaf %*% t(cbind(1, rbind(0, diag(3)))))

stan_ess(fit.mcmc2)
stan_rhat
stan_diag
stan_mcse
stan_ac
stan_trace(fit.mcmc2, pars=c("beta[2]"))
