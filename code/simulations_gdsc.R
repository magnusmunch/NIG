#!/usr/bin/env Rscript

### installation of packages
# if(!("devtools" %in% installed.packages())) {
#   install.packages("devtools")
# }
# library(devtools)
# install_github("magnusmunch/cambridge/rpackage", local=FALSE,
#                auth_token=Sys.getenv("GITHUB_PAT"))

### libraries
detach("package:cambridge", unload=TRUE)
library(cambridge)
library(statmod)
library(glmnet)

### load data
load(file="data/data_gdsc_dat1.Rdata")
load(file="results/data_gdsc_fit1.Rdata")
# expr.sel <- sapply(idsel, function(id) {expr.prep[, id]})
psel <- 300
idsel <- order(apply(expr.prep, 2, sd), decreasing=TRUE)[c(1:psel)]
expr.sel <- expr.prep[, idsel]
rm(drug.prep, expr.prep, fit.enet, idsel, resp.prep)

### data preparation
# x <- lapply(expr.sel, scale)
x <- scale(expr.sel)
# D <- length(x)
D <- 500
# p <- sapply(x, ncol)
p <- ncol(x)
# n <- nrow(x[[1]])
n <- nrow(x)
ntrain <- floor(n/2)

################################# simulation 1 #################################
### simulation settings
alphaf <- c(1, 1, 3, 7)
lambdaf <- 1
shape <- 3
rate <- 2
set <- list(alphaf=alphaf, lambdaf=lambdaf, shape=shape, rate=rate)

nreps <- 50
methods <- c(paste0("semnig", 1:4), "lasso", "ridge")
res <- matrix(NA, nrow=nreps, ncol=3*length(methods))
colnames(res) <- c(paste0("aemse.", methods), paste0("apmse.", methods),
                   paste0("apmset.", methods))
est <- matrix(NA, nrow=nreps, ncol=14)
colnames(est) <- paste0(c(rep(methods[c(1, 2)], each=5),
                          rep(methods[c(3, 4)], each=2)), ".",
                        c(rep(c(paste0("alphaf.", 0:3), "lambdaf"), 2),
                          rep(c("alphaf.0", "lambdaf"), 2)))
elbo <- matrix(NA, nrow=nreps, ncol=4)
colnames(elbo) <- methods[c(1:4)]

for(r in 1:nreps) {
  cat("\r", "replication", r)
  set.seed(2019 + r)
  
  ### simulate parameters
  C <- lapply(1:D, function(d) {
    cbind(1, t(replicate(p[d], sample(c(1, 0, 0, 0), 4)[-1])))})
  gamma <- lapply(1:D, function(d) {
    sqrt(rinvgauss(p[d], 1/as.numeric(C[[d]] %*% alphaf), lambdaf))})
  sigma <- 1/sqrt(rgamma(D, shape, rate))
  beta <- lapply(1:D, function(d) {rnorm(p[d], 0, gamma[[d]]*sigma[d])})
  
  # simulate data
  idtrain <- sample(1:n, ntrain)
  xtrain <- lapply(x, function(s) {scale(s[idtrain, ])})
  xtest <- lapply(x, function(s) {scale(s[-idtrain, ])})
  ytrain <- scale(sapply(1:D, function(d) {
    rnorm(ntrain, as.numeric(xtrain[[d]] %*% beta[[d]]), sigma[d])}))
  ytest <- scale(sapply(1:D, function(d) {
    rnorm(n - ntrain, as.numeric(xtest[[d]] %*% beta[[d]]), sigma[d])}))
  
  ### fitting models
  control <- list(conv.post=TRUE, trace=FALSE,
                  epsilon.eb=1e-3, epsilon.vb=1e-3,
                  maxit.eb=100, maxit.vb=1, maxit.post=200)
  fit.semnig1 <- semnig(x=xtrain, y=ytrain, C=C, Z=NULL, unpenalized=NULL,
                        standardize=FALSE, intercept=FALSE, fixed.eb="none", 
                        full.post=FALSE, init=NULL, control=control)
  init <- list(lambdaf=lambdaf, lambdad=1000)
  fit.semnig2 <- semnig(x=xtrain, y=ytrain, C=C, Z=NULL, unpenalized=NULL,
                        standardize=FALSE, intercept=FALSE, fixed.eb="lambda", 
                        full.post=FALSE, init=init, control=control)
  C <- lapply(1:D, function(d) {matrix(1, nrow=p[d], ncol=1)})
  fit.semnig3 <- semnig(x=xtrain, y=ytrain, C=C, Z=NULL, unpenalized=NULL,
                        standardize=FALSE, intercept=FALSE, fixed.eb="none", 
                        full.post=FALSE, init=init, control=control)
  fit.semnig4 <- semnig(x=xtrain, y=ytrain, C=C, Z=NULL, unpenalized=NULL,
                        standardize=FALSE, intercept=FALSE, fixed.eb="lambda", 
                        full.post=FALSE, init=init, control=control)
  
  fit.lasso <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[, d], intercept=FALSE, standardize=FALSE)})
  fit.ridge <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[, d], alpha=0, intercept=FALSE, 
              standardize=FALSE)})
  
  best.semnig1 <- fit.semnig1$vb$mpost$beta
  best.semnig2 <- fit.semnig2$vb$mpost$beta
  best.semnig3 <- fit.semnig3$vb$mpost$beta
  best.semnig4 <- fit.semnig4$vb$mpost$beta
  best.lasso <- lapply(fit.lasso, function(s) {coef(s, s="lambda.min")[-1, 1]})
  best.ridge <- lapply(fit.ridge, function(s) {coef(s, s="lambda.min")[-1, 1]})
  
  pred.semnig1 <- sapply(1:D, function(d) {xtest[[d]] %*% best.semnig1[[d]]})
  pred.semnig2 <- sapply(1:D, function(d) {xtest[[d]] %*% best.semnig2[[d]]})
  pred.semnig3 <- sapply(1:D, function(d) {xtest[[d]] %*% best.semnig3[[d]]})
  pred.semnig4 <- sapply(1:D, function(d) {xtest[[d]] %*% best.semnig4[[d]]})
  pred.lasso <- sapply(1:D, function(d) {xtest[[d]] %*% best.lasso[[d]]})
  pred.ridge <- sapply(1:D, function(d) {xtest[[d]] %*% best.ridge[[d]]})
  
  psel.lasso <- sapply(fit.lasso, function(s) {s$nzero[which.min(s$cvm)]})
  
  emse.semnig1 <- sapply(1:D, function(d) {
    mean((best.semnig1[[d]] - beta[[d]])^2)})
  emse.semnig2 <- sapply(1:D, function(d) {
    mean((best.semnig2[[d]] - beta[[d]])^2)})
  emse.semnig3 <- sapply(1:D, function(d) {
    mean((best.semnig3[[d]] - beta[[d]])^2)})
  emse.semnig4 <- sapply(1:D, function(d) {
    mean((best.semnig4[[d]] - beta[[d]])^2)})
  emse.lasso <- sapply(1:D, function(d) {mean((best.lasso[[d]] - beta[[d]])^2)})
  emse.ridge <- sapply(1:D, function(d) {mean((best.ridge[[d]] - beta[[d]])^2)})
  
  pmse.semnig1 <- colMeans((pred.semnig1 - ytest)^2)
  pmse.semnig2 <- colMeans((pred.semnig2 - ytest)^2)
  pmse.semnig3 <- colMeans((pred.semnig3 - ytest)^2)
  pmse.semnig4 <- colMeans((pred.semnig4 - ytest)^2)
  pmse.lasso <- colMeans((pred.lasso - ytest)^2)
  pmse.ridge <- colMeans((pred.ridge - ytest)^2)
  
  predt.semnig1 <- sapply(1:D, function(d) {xtrain[[d]] %*% best.semnig1[[d]]})
  predt.semnig2 <- sapply(1:D, function(d) {xtrain[[d]] %*% best.semnig2[[d]]})
  predt.semnig3 <- sapply(1:D, function(d) {xtrain[[d]] %*% best.semnig3[[d]]})
  predt.semnig4 <- sapply(1:D, function(d) {xtrain[[d]] %*% best.semnig4[[d]]})
  predt.lasso <- sapply(1:D, function(d) {xtrain[[d]] %*% best.lasso[[d]]})
  predt.ridge <- sapply(1:D, function(d) {xtrain[[d]] %*% best.ridge[[d]]})
  
  pmset.semnig1 <- colMeans((predt.semnig1 - ytrain)^2)
  pmset.semnig2 <- colMeans((predt.semnig2 - ytrain)^2)
  pmset.semnig3 <- colMeans((predt.semnig3 - ytrain)^2)
  pmset.semnig4 <- colMeans((predt.semnig4 - ytrain)^2)
  pmset.lasso <- colMeans((predt.lasso - ytrain)^2)
  pmset.ridge <- colMeans((predt.ridge - ytrain)^2)
  
  res[r, 1:6] <- c(mean(emse.semnig1), mean(emse.semnig2), mean(emse.semnig3),
                   mean(emse.semnig4), mean(emse.lasso), mean(emse.ridge))
  res[r, 7:12] <- c(mean(pmse.semnig1), mean(pmse.semnig2), mean(pmse.semnig3),
                   mean(pmse.semnig4), mean(pmse.lasso), mean(pmse.ridge))
  res[r, 13:18] <- c(mean(pmset.semnig1), mean(pmset.semnig2), 
                     mean(pmset.semnig3), mean(pmset.semnig4), 
                     mean(pmset.lasso), mean(pmse.ridge))
  est[r, ] <- c(fit.semnig1$eb$alphaf, fit.semnig1$eb$lambdaf,
                fit.semnig2$eb$alphaf, fit.semnig2$eb$lambdaf,
                fit.semnig3$eb$alphaf, fit.semnig3$eb$lambdaf,
                fit.semnig4$eb$alphaf, fit.semnig4$eb$lambdaf)
  elbo[r, ] <- c(tail(rowMeans(fit.semnig1$seq.elbo), n=1L), 
                 tail(rowMeans(fit.semnig2$seq.elbo), n=1L),
                 tail(rowMeans(fit.semnig3$seq.elbo), n=1L), 
                 tail(rowMeans(fit.semnig4$seq.elbo), n=1L))
  save(res, est, set, elbo, file="results/simulations_gdsc_set1.Rdata")
}  


################################# simulation 2 #################################
### simulation settings
alphaf <- c(1, 1, 3, 7)
alphad <- c(1, 1, 3, 7)
lambdaf <- 1
lambdad <- 1
shape <- 3
rate <- 2
set <- list(alphaf=alphaf, alphad=alphad, lambdaf=lambdaf, lambdad=lambdad,
            shape=shape, rate=rate)

nreps <- 50
methods <- c(paste0("semnig", 1:4), "lasso", "ridge")
res <- matrix(NA, nrow=nreps, ncol=3*length(methods))
colnames(res) <- c(paste0("aemse.", methods), paste0("apmse.", methods),
                   paste0("apmset.", methods))
est <- matrix(NA, nrow=nreps, ncol=28)
colnames(est) <- paste0(c(rep(methods[c(1, 2)], each=10),
                          rep(methods[c(3, 4)], each=4)), ".",
                        c(rep(c(paste0("alphaf.", 0:3), paste0("alphad.", 0:3), 
                                "lambdaf", "lambdad"), 2),
                          rep(c("alphaf.0", "alphad.0", "lambdaf", "lambdad"), 
                              2)))
elbo <- matrix(NA, nrow=nreps, ncol=4)
mprod <- matrix(NA, nrow=nreps, ncol=34)
colnames(mprod) <- c(rep(c("semnig1", "semnig2"), each=16), "semnig3", 
                     "semnig4")
colnames(elbo) <- methods[c(1:4)]

for(r in 1:nreps) {
  cat("\r", "replication", r)
  set.seed(2019 + r)
  
  ### simulate parameters
  C <- lapply(1:D, function(d) {
    cbind(1, t(replicate(p, sample(c(1, 0, 0, 0), 4)[-1])))})
  Z <- cbind(1, t(replicate(D, sample(c(1, 0, 0, 0), 4)[-1])))
  gamma <- lapply(1:D, function(d) {
    sqrt(rinvgauss(p, 1/as.numeric(C[[d]] %*% alphaf), lambdaf))})
  tau <- sqrt(rinvgauss(D, 1/as.numeric(Z %*% alphad), lambdad))
  sigma <- 1/sqrt(rgamma(D, shape, rate))
  beta <- lapply(1:D, function(d) {rnorm(p, 0, gamma[[d]]*tau[d]*sigma[d])})
  
  # simulate data
  idtrain <- sample(1:n, ntrain)
  xtrain <- scale(x[idtrain, ])
  xtest <- scale(x[-idtrain, ])
  ytrain <- scale(sapply(1:D, function(d) {
    rnorm(ntrain, as.numeric(xtrain %*% beta[[d]]), sigma[d])}))
  ytest <- scale(sapply(1:D, function(d) {
    rnorm(n - ntrain, as.numeric(xtest %*% beta[[d]]), sigma[d])}))
  
  ### fitting models
  control <- list(conv.post=TRUE, trace=TRUE,
                  epsilon.eb=1e-3, epsilon.vb=1e-3,
                  maxit.eb=100, maxit.vb=1, maxit.post=200)
  fit.semnig1 <- semnig(x=rep(list(xtrain), D), y=ytrain, C=C, Z=Z, unpenalized=NULL,
                        standardize=FALSE, intercept=FALSE, fixed.eb="none", 
                        full.post=FALSE, init=NULL, control=control)
  
  plot(1/as.numeric(set$alphad %*% t(unname(model.matrix(~ factor(1:4))))),
       1/as.numeric(fit.semnig1$eb$alphad %*% t(unname(model.matrix(~ factor(1:4))))))
  plot(1/as.numeric(set$alphaf %*% t(unname(model.matrix(~ factor(1:4))))),
       1/as.numeric(fit.semnig1$eb$alphaf %*% t(unname(model.matrix(~ factor(1:4))))))
  
  init <- list(lambdaf=lambdaf, lambdad=lambdad)
  fit.semnig2 <- semnig(x=xtrain, y=ytrain, C=C, Z=Z, unpenalized=NULL,
                        standardize=FALSE, intercept=FALSE, fixed.eb="lambda", 
                        full.post=FALSE, init=init, control=control)
  C <- lapply(1:D, function(d) {matrix(1, nrow=p[d], ncol=1)})
  Z <- matrix(1, nrow=D, ncol=1)
  fit.semnig3 <- semnig(x=xtrain, y=ytrain, C=C, Z=Z, unpenalized=NULL,
                        standardize=FALSE, intercept=FALSE, fixed.eb="none", 
                        full.post=FALSE, init=NULL, control=control)
  fit.semnig4 <- semnig(x=xtrain, y=ytrain, C=C, Z=Z, unpenalized=NULL,
                        standardize=FALSE, intercept=FALSE, fixed.eb="lambda", 
                        full.post=FALSE, init=init, control=control)
  
  fit.lasso <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[, d], intercept=FALSE, standardize=FALSE)})
  fit.ridge <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[, d], alpha=0, intercept=FALSE, 
              standardize=FALSE)})
  
  best.semnig1 <- fit.semnig1$vb$mpost$beta
  best.semnig2 <- fit.semnig2$vb$mpost$beta
  best.semnig3 <- fit.semnig3$vb$mpost$beta
  best.semnig4 <- fit.semnig4$vb$mpost$beta
  best.lasso <- lapply(fit.lasso, function(s) {coef(s, s="lambda.min")[-1, 1]})
  best.ridge <- lapply(fit.ridge, function(s) {coef(s, s="lambda.min")[-1, 1]})
  
  pred.semnig1 <- sapply(1:D, function(d) {xtest[[d]] %*% best.semnig1[[d]]})
  pred.semnig2 <- sapply(1:D, function(d) {xtest[[d]] %*% best.semnig2[[d]]})
  pred.semnig3 <- sapply(1:D, function(d) {xtest[[d]] %*% best.semnig3[[d]]})
  pred.semnig4 <- sapply(1:D, function(d) {xtest[[d]] %*% best.semnig4[[d]]})
  pred.lasso <- sapply(1:D, function(d) {xtest[[d]] %*% best.lasso[[d]]})
  pred.ridge <- sapply(1:D, function(d) {xtest[[d]] %*% best.ridge[[d]]})
  
  psel.lasso <- sapply(fit.lasso, function(s) {s$nzero[which.min(s$cvm)]})
  
  emse.semnig1 <- sapply(1:D, function(d) {
    mean((best.semnig1[[d]] - beta[[d]])^2)})
  emse.semnig2 <- sapply(1:D, function(d) {
    mean((best.semnig2[[d]] - beta[[d]])^2)})
  emse.semnig3 <- sapply(1:D, function(d) {
    mean((best.semnig3[[d]] - beta[[d]])^2)})
  emse.semnig4 <- sapply(1:D, function(d) {
    mean((best.semnig4[[d]] - beta[[d]])^2)})
  emse.lasso <- sapply(1:D, function(d) {mean((best.lasso[[d]] - beta[[d]])^2)})
  emse.ridge <- sapply(1:D, function(d) {mean((best.ridge[[d]] - beta[[d]])^2)})
  
  pmse.semnig1 <- colMeans((pred.semnig1 - ytest)^2)
  pmse.semnig2 <- colMeans((pred.semnig2 - ytest)^2)
  pmse.semnig3 <- colMeans((pred.semnig3 - ytest)^2)
  pmse.semnig4 <- colMeans((pred.semnig4 - ytest)^2)
  pmse.lasso <- colMeans((pred.lasso - ytest)^2)
  pmse.ridge <- colMeans((pred.ridge - ytest)^2)
  
  predt.semnig1 <- sapply(1:D, function(d) {xtrain[[d]] %*% best.semnig1[[d]]})
  predt.semnig2 <- sapply(1:D, function(d) {xtrain[[d]] %*% best.semnig2[[d]]})
  predt.semnig3 <- sapply(1:D, function(d) {xtrain[[d]] %*% best.semnig3[[d]]})
  predt.semnig4 <- sapply(1:D, function(d) {xtrain[[d]] %*% best.semnig4[[d]]})
  predt.lasso <- sapply(1:D, function(d) {xtrain[[d]] %*% best.lasso[[d]]})
  predt.ridge <- sapply(1:D, function(d) {xtrain[[d]] %*% best.ridge[[d]]})
  
  pmset.semnig1 <- colMeans((predt.semnig1 - ytrain)^2)
  pmset.semnig2 <- colMeans((predt.semnig2 - ytrain)^2)
  pmset.semnig3 <- colMeans((predt.semnig3 - ytrain)^2)
  pmset.semnig4 <- colMeans((predt.semnig4 - ytrain)^2)
  pmset.lasso <- colMeans((predt.lasso - ytrain)^2)
  pmset.ridge <- colMeans((predt.ridge - ytrain)^2)
  
  res[r, 1:6] <- c(mean(emse.semnig1), mean(emse.semnig2), mean(emse.semnig3),
                   mean(emse.semnig4), mean(emse.lasso), mean(emse.ridge))
  res[r, 7:12] <- c(mean(pmse.semnig1), mean(pmse.semnig2), mean(pmse.semnig3),
                    mean(pmse.semnig4), mean(pmse.lasso), mean(pmse.ridge))
  res[r, 13:18] <- c(mean(pmset.semnig1), mean(pmset.semnig2), 
                     mean(pmset.semnig3), mean(pmset.semnig4), 
                     mean(pmset.lasso), mean(pmse.ridge))
  est[r, ] <- c(fit.semnig1$eb$alphaf, fit.semnig1$eb$alphad,
                fit.semnig1$eb$lambdaf, fit.semnig1$eb$lambdad,
                fit.semnig2$eb$alphaf, fit.semnig2$eb$alphad,
                fit.semnig2$eb$lambdaf, fit.semnig2$eb$lambdad,
                fit.semnig3$eb$alphaf, fit.semnig3$eb$alphad,
                fit.semnig3$eb$lambdaf, fit.semnig3$eb$lambdad,
                fit.semnig4$eb$alphaf, fit.semnig4$eb$alphad, 
                fit.semnig4$eb$lambdaf, fit.semnig4$eb$lambdad)
  elbo[r, ] <- c(tail(rowMeans(fit.semnig1$seq.elbo), n=1L), 
                 tail(rowMeans(fit.semnig2$seq.elbo), n=1L),
                 tail(rowMeans(fit.semnig3$seq.elbo), n=1L), 
                 tail(rowMeans(fit.semnig4$seq.elbo), n=1L))
  mprod[r, ] <- c(outer(1/as.numeric(fit.semnig1$eb$alphad %*% t(
    unname(model.matrix(~ factor(1:4))))), 1/as.numeric(
      fit.semnig1$eb$alphaf %*% t(unname(model.matrix(~ factor(1:4))))), "*"),
    outer(1/as.numeric(fit.semnig2$eb$alphad %*% t(unname(
      model.matrix(~ factor(1:4))))),
      1/as.numeric(fit.semnig2$eb$alphaf %*% t(unname(model.matrix(~ factor(
        1:4))))), "*"),
    1/(fit.semnig3$eb$alphad*fit.semnig3$eb$alphaf),
    1/(fit.semnig4$eb$alphad*fit.semnig4$eb$alphaf))
  save(res, est, set, elbo, mprod, file="results/simulations_gdsc_set2.Rdata")
}  



################################# simulation 3 #################################
### simulation settings
alphad <- c(1, 1, 3, 7)
lambdad <- 1
shape <- 3
rate <- 2
set <- list(alphad=alphad, lambdad=lambdad, shape=shape, rate=rate)

nreps <- 50
methods <- c(paste0("semnig", 1:4), "lasso", "ridge")
res <- matrix(NA, nrow=nreps, ncol=3*length(methods))
colnames(res) <- c(paste0("aemse.", methods), paste0("apmse.", methods),
                   paste0("apmset.", methods))
est <- matrix(NA, nrow=nreps, ncol=14)
colnames(est) <- paste0(c(rep(methods[c(1, 2)], each=5),
                          rep(methods[c(3, 4)], each=2)), ".",
                        c(rep(c(paste0("alphad.", 0:3), "lambdad"), 2),
                          rep(c("alphad.0", "lambdad"), 2)))
elbo <- matrix(NA, nrow=nreps, ncol=4)
colnames(elbo) <- methods[c(1:4)]

for(r in 1:nreps) {
  cat("\r", "replication", r)
  set.seed(2019 + r)
  
  ### simulate parameters
  Z <- cbind(1, t(replicate(D, sample(c(1, 0, 0, 0), 4)[-1])))
  tau <- sqrt(rinvgauss(D, 1/as.numeric(Z %*% alphad), lambdad))
  sigma <- 1/sqrt(rgamma(D, shape, rate))
  beta <- lapply(1:D, function(d) {rnorm(p[d], 0, tau[d]*sigma[d])})
  
  # simulate data
  idtrain <- sample(1:n, ntrain)
  xtrain <- lapply(x, function(s) {scale(s[idtrain, ])})
  xtest <- lapply(x, function(s) {scale(s[-idtrain, ])})
  ytrain <- scale(sapply(1:D, function(d) {
    rnorm(ntrain, as.numeric(xtrain[[d]] %*% beta[[d]]), sigma[d])}))
  ytest <- scale(sapply(1:D, function(d) {
    rnorm(n - ntrain, as.numeric(xtest[[d]] %*% beta[[d]]), sigma[d])}))
  
  ### fitting models
  control <- list(conv.post=TRUE, trace=TRUE,
                  epsilon.eb=1e-3, epsilon.vb=1e-3,
                  maxit.eb=50, maxit.vb=1, maxit.post=200)
  fit.semnig1 <- semnig(x=xtrain, y=ytrain, C=NULL, Z=Z, unpenalized=NULL,
                        standardize=FALSE, intercept=FALSE, fixed.eb="none", 
                        full.post=FALSE, init=NULL, control=control)
  plot(1/as.numeric(set$alphad %*% t(unname(model.matrix(~ factor(1:4))))),
       1/as.numeric(fit.semnig1$eb$alphad %*% t(unname(model.matrix(~ factor(1:4))))))
       
  init <- list(lambdad=lambdad)
  fit.semnig2 <- semnig(x=xtrain, y=ytrain, C=NULL, Z=Z, unpenalized=NULL,
                        standardize=FALSE, intercept=FALSE, fixed.eb="lambda", 
                        full.post=FALSE, init=init, control=control)
  C <- lapply(1:D, function(d) {matrix(1, nrow=p[d], ncol=1)})
  Z <- matrix(1, nrow=D, ncol=1)
  fit.semnig3 <- semnig(x=xtrain, y=ytrain, C=NULL, Z=Z, unpenalized=NULL,
                        standardize=FALSE, intercept=FALSE, fixed.eb="none", 
                        full.post=FALSE, init=NULL, control=control)
  fit.semnig4 <- semnig(x=xtrain, y=ytrain, C=NULL, Z=Z, unpenalized=NULL,
                        standardize=FALSE, intercept=FALSE, fixed.eb="lambda", 
                        full.post=FALSE, init=init, control=control)
  
  fit.lasso <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[, d], intercept=FALSE, standardize=FALSE)})
  fit.ridge <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[, d], alpha=0, intercept=FALSE, 
              standardize=FALSE)})
  
  best.semnig1 <- fit.semnig1$vb$mpost$beta
  best.semnig2 <- fit.semnig2$vb$mpost$beta
  best.semnig3 <- fit.semnig3$vb$mpost$beta
  best.semnig4 <- fit.semnig4$vb$mpost$beta
  best.lasso <- lapply(fit.lasso, function(s) {coef(s, s="lambda.min")[-1, 1]})
  best.ridge <- lapply(fit.ridge, function(s) {coef(s, s="lambda.min")[-1, 1]})
  
  pred.semnig1 <- sapply(1:D, function(d) {xtest[[d]] %*% best.semnig1[[d]]})
  pred.semnig2 <- sapply(1:D, function(d) {xtest[[d]] %*% best.semnig2[[d]]})
  pred.semnig3 <- sapply(1:D, function(d) {xtest[[d]] %*% best.semnig3[[d]]})
  pred.semnig4 <- sapply(1:D, function(d) {xtest[[d]] %*% best.semnig4[[d]]})
  pred.lasso <- sapply(1:D, function(d) {xtest[[d]] %*% best.lasso[[d]]})
  pred.ridge <- sapply(1:D, function(d) {xtest[[d]] %*% best.ridge[[d]]})
  
  psel.lasso <- sapply(fit.lasso, function(s) {s$nzero[which.min(s$cvm)]})
  
  emse.semnig1 <- sapply(1:D, function(d) {
    mean((best.semnig1[[d]] - beta[[d]])^2)})
  emse.semnig2 <- sapply(1:D, function(d) {
    mean((best.semnig2[[d]] - beta[[d]])^2)})
  emse.semnig3 <- sapply(1:D, function(d) {
    mean((best.semnig3[[d]] - beta[[d]])^2)})
  emse.semnig4 <- sapply(1:D, function(d) {
    mean((best.semnig4[[d]] - beta[[d]])^2)})
  emse.lasso <- sapply(1:D, function(d) {mean((best.lasso[[d]] - beta[[d]])^2)})
  emse.ridge <- sapply(1:D, function(d) {mean((best.ridge[[d]] - beta[[d]])^2)})
  
  pmse.semnig1 <- colMeans((pred.semnig1 - ytest)^2)
  pmse.semnig2 <- colMeans((pred.semnig2 - ytest)^2)
  pmse.semnig3 <- colMeans((pred.semnig3 - ytest)^2)
  pmse.semnig4 <- colMeans((pred.semnig4 - ytest)^2)
  pmse.lasso <- colMeans((pred.lasso - ytest)^2)
  pmse.ridge <- colMeans((pred.ridge - ytest)^2)
  
  predt.semnig1 <- sapply(1:D, function(d) {xtrain[[d]] %*% best.semnig1[[d]]})
  predt.semnig2 <- sapply(1:D, function(d) {xtrain[[d]] %*% best.semnig2[[d]]})
  predt.semnig3 <- sapply(1:D, function(d) {xtrain[[d]] %*% best.semnig3[[d]]})
  predt.semnig4 <- sapply(1:D, function(d) {xtrain[[d]] %*% best.semnig4[[d]]})
  predt.lasso <- sapply(1:D, function(d) {xtrain[[d]] %*% best.lasso[[d]]})
  predt.ridge <- sapply(1:D, function(d) {xtrain[[d]] %*% best.ridge[[d]]})
  
  pmset.semnig1 <- colMeans((predt.semnig1 - ytrain)^2)
  pmset.semnig2 <- colMeans((predt.semnig2 - ytrain)^2)
  pmset.semnig3 <- colMeans((predt.semnig3 - ytrain)^2)
  pmset.semnig4 <- colMeans((predt.semnig4 - ytrain)^2)
  pmset.lasso <- colMeans((predt.lasso - ytrain)^2)
  pmset.ridge <- colMeans((predt.ridge - ytrain)^2)
  
  res[r, 1:6] <- c(mean(emse.semnig1), mean(emse.semnig2), mean(emse.semnig3),
                   mean(emse.semnig4), mean(emse.lasso), mean(emse.ridge))
  res[r, 7:12] <- c(mean(pmse.semnig1), mean(pmse.semnig2), mean(pmse.semnig3),
                    mean(pmse.semnig4), mean(pmse.lasso), mean(pmse.ridge))
  res[r, 13:18] <- c(mean(pmset.semnig1), mean(pmset.semnig2), 
                     mean(pmset.semnig3), mean(pmset.semnig4), 
                     mean(pmset.lasso), mean(pmse.ridge))
  est[r, ] <- c(fit.semnig1$eb$alphaf, fit.semnig1$eb$alphad,
                fit.semnig1$eb$lambdaf, fit.semnig1$eb$lambdad,
                fit.semnig2$eb$alphaf, fit.semnig2$eb$alphad,
                fit.semnig2$eb$lambdaf, fit.semnig2$eb$lambdad,
                fit.semnig3$eb$alphaf, fit.semnig3$eb$alphad,
                fit.semnig3$eb$lambdaf, fit.semnig3$eb$lambdad,
                fit.semnig4$eb$alphaf, fit.semnig4$eb$alphad, 
                fit.semnig4$eb$lambdaf, fit.semnig4$eb$lambdad)
  elbo[r, ] <- c(tail(rowMeans(fit.semnig1$seq.elbo), n=1L), 
                 tail(rowMeans(fit.semnig2$seq.elbo), n=1L),
                 tail(rowMeans(fit.semnig3$seq.elbo), n=1L), 
                 tail(rowMeans(fit.semnig4$seq.elbo), n=1L))
  mprod[r, ] <- c(outer(1/as.numeric(fit.semnig1$eb$alphad %*% t(
    unname(model.matrix(~ factor(1:4))))), 1/as.numeric(
      fit.semnig1$eb$alphaf %*% t(unname(model.matrix(~ factor(1:4))))), "*"),
    outer(1/as.numeric(fit.semnig2$eb$alphad %*% t(unname(
      model.matrix(~ factor(1:4))))),
      1/as.numeric(fit.semnig2$eb$alphaf %*% t(unname(model.matrix(~ factor(
        1:4))))), "*"),
    1/(fit.semnig3$eb$alphad*fit.semnig3$eb$alphaf),
    1/(fit.semnig4$eb$alphad*fit.semnig4$eb$alphaf))
  save(res, est, set, elbo, mprod, file="results/simulations_gdsc_set2.Rdata")
}  

