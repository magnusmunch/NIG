#!/usr/bin/env Rscript

# number of cores to use
ncores <- 100

### libraries
packages <- c("foreach", "doParallel", "cambridge", "statmod", "glmnet", 
              "xtune")
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
control.semnig <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-3, 
                       epsilon.vb=1e-3, maxit.eb=500, maxit.vb=1, 
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
  set.seed(2019 + r)

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
lambdad <- 1
Z <- unname(model.matrix(~ factor(rep(1:H, each=D/H))))
x <- scale(expr$expr[, order(-apply(expr$expr, 2, sd))[1:p]])

# estimation settings
nfolds <- 10
methods <- c("NIG-", "NIG", "ridge", "lasso", "xtune")
control.semnig <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-3, 
                       epsilon.vb=1e-3, maxit.eb=500, maxit.vb=1, 
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
  set.seed(2019 + r)
  
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
control.semnig <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-3, 
                       epsilon.vb=1e-3, maxit.eb=500, maxit.vb=1, 
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
fracs <- c(1/10, 1/5, 1/3) # fraction of permuted rows in external covariates
C <- replicate(D, list(unname(model.matrix(~ factor(rep(1:G, each=p/G))))))
Z <- unname(model.matrix(~ factor(rep(1:H, each=D/H))))
x <- scale(expr$expr[, order(-apply(expr$expr, 2, sd))[1:p]])

# estimation settings
nfolds <- 10
methods <- c("NIG-", "NIG", "ridge", "lasso", "xtune", "ebridge")
control.semnig <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, 
                       epsilon.vb=1e-3, maxit.eb=500, maxit.vb=1, 
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
write.table(res2, file="results/simulations_gdsc_res4.txt")
  
