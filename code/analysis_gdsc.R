#!/usr/bin/env Rscript

### libraries
library(cambridge)
library(pInc)
library(glmnet)

################################## analysis 1 ##################################
### data preparation
# load prepared data
load(file="data/data_gdsc_dat1.Rdata")

# select features
D <- ncol(resp.prep)
psel <- 500
o <- order(-apply(expr.prep, 2, sd))
idsel <- lapply(feat.prep$inpathway, function(s) {
  m1 <- o[o %in% which(s==1)]; id1 <- head(m1, n=min(length(m1), psel/2))
  m0 <- o[o %in% which(s==0)]; id0 <- head(m0, n=psel - length(id1))
  c(id1, id0)})
expr.sel <- lapply(idsel, function(s) {expr.prep[, s]})
inpathway <- lapply(1:D, function(d) {feat.prep$inpathway[[d]][idsel[[d]]]})
# idsel <- o[1:psel]
# expr.sel <- lapply(1:D, function(d) {expr.prep[, idsel]})
# inpathway <- lapply(feat.prep$inpathway, function(s) {s[idsel]})
x <- lapply(expr.sel, function(s) {scale(s)})
y <- scale(resp.prep)

# number of equations, features, and observations
p <- sapply(expr.sel, ncol)
n <- nrow(resp.prep)

### model fitting
# setting seed
set.seed(2019)

# estimation settings
control <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, epsilon.vb=1e-3,
                maxit.eb=200, maxit.vb=1, maxit.post=100)
methods <- c("NIG1", "NIG2", "NIG3", "NIG4", "lasso", "ridge")

# model without external covariates
Z <- matrix(1, nrow=D)
colnames(Z) <- c("intercept")
C <- lapply(inpathway, function(s) {
  s <- matrix(1, nrow=length(s)); colnames(s) <- c("intercept"); return(s)})
fit1.semnig <- semnig(x=x, y=y, C=C, Z=Z, unpenalized=NULL,
                      standardize=FALSE, intercept=FALSE, fixed.eb="none",
                      full.post=TRUE, init=NULL, control=control)

# models with external covariates
Z <- model.matrix(~ replace(drug.prep$stage, is.na(drug.prep$stage),
                            "experimental") + replace(drug.prep$action, is.na(
                              drug.prep$action), "unknown"))
colnames(Z) <- c("intercept", "experimental", "in clinical development",
                 "targeted", "unknown")
C <- lapply(inpathway, function(s) {
  s <- cbind(1, s); colnames(s) <- c("intercept", "in pathway"); return(s)})
fit2.semnig <- semnig(x=x, y=y, C=C, Z=Z, unpenalized=NULL,
                      standardize=FALSE, intercept=FALSE, fixed.eb=FALSE,
                      full.post=TRUE, init=NULL, control=control)
fit3.semnig <- semnig(x=x, y=y, C=C, Z=NULL, unpenalized=NULL,
                      standardize=FALSE, intercept=FALSE, fixed.eb=FALSE,
                      full.post=TRUE, init=NULL, control=control)
fit4.semnig <- semnig(x=x, y=y, C=NULL, Z=Z, unpenalized=NULL,
                      standardize=FALSE, intercept=FALSE, fixed.eb=FALSE,
                      full.post=TRUE, init=NULL, control=control)

# bSEM model
# fit1.bSEM <- bSEM(rep(list(x), D), split(y, 1:ncol(y)),
#                   lapply(inpathway, "+", 1),
#                   control=list(maxit=200, trace=TRUE, epsilon=1e-3))

fit1.lasso <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[, d], intercept=FALSE, standardize=FALSE)})
fit1.ridge <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[, d], alpha=0, intercept=FALSE, standardize=FALSE)})

save(fit1.semnig, fit2.semnig, fit3.semnig, fit4.semnig, fit1.lasso, fit1.ridge,
     file="results/data_gdsc_fit1.Rdata")

# EB estimates
tab <- rbind(c(fit1.semnig$eb$alphaf, NA, fit1.semnig$eb$alphad, rep(NA, 4), 
               fit1.semnig$eb$lambdaf, fit1.semnig$eb$lambdad),
             c(fit2.semnig$eb$alphaf, fit2.semnig$eb$alphad, 
               fit2.semnig$eb$lambdaf, fit2.semnig$eb$lambdad),
             c(fit3.semnig$eb$alphaf, rep(NA, 5), fit3.semnig$eb$lambdaf, NA),
             c(rep(NA, 2), fit4.semnig$eb$alphad, NA, fit2.semnig$eb$lambdad))
             # c(rep(NA, 5), fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 1]/
             #     fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 1],
             #   fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 2]/
             #     fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 2] -
             #     fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 1]/
             #     fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 1], rep(NA, 7))
colnames(tab) <- c(colnames(C[[1]]), colnames(Z), "lambdaf", "lambdad")
rownames(tab) <- methods[c(1:4)]
write.table(tab, file="results/data_gdsc_fit1.txt")

### cross-validation
# estimation settings
control <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-3, epsilon.vb=1e-3, 
                maxit.eb=200, maxit.vb=1, maxit.post=100)

nreps <- 50
ntrain <- floor(nrow(resp.prep)/2)

# results matrices
methods <- c("NIG1", "NIG2", "NIG3", "NIG4", "lasso", "ridge")
pmse <- matrix(NA, nrow=nreps, ncol=length(methods))
colnames(pmse) <- paste0("pmse.", methods)
brank <- matrix(NA, nrow=nreps, ncol=sum(p)*length(methods))
colnames(brank) <- paste0(
  "brank.", paste0(rep(methods, each=sum(p)), ".", 
                   rep(paste0(rep(1:D, times=p), ".", 
                              unlist(sapply(p, function(s) {
                                return(1:s)}))), length(methods))))
elbo <- matrix(NA, ncol=length(methods) - 2, nrow=nreps)
colnames(elbo) <- c("NIG1", "NIG2", "NIG3", "NIG4")

for(r in 1:nreps) {
  cat("\r", "rep", r)
  set.seed(2019 + r)
  
  # splitting data
  idtrain <- sample(1:nrow(y), ntrain)
  xtrain <- lapply(x, function(s) {scale(s[idtrain, ])})
  ytrain <- scale(y[idtrain, ])
  xtest <- lapply(x, function(s) {scale(s[-idtrain, ])})
  ytest <- scale(y[-idtrain, ])
  
  # model without external covariates
  Z <- matrix(1, nrow=D)
  colnames(Z) <- c("intercept")
  C <- lapply(inpathway, function(s) {
    s <- matrix(1, nrow=length(s)); colnames(s) <- c("intercept"); return(s)})
  cv1.semnig <- semnig(x=xtrain, y=ytrain, C=C, Z=Z, unpenalized=NULL, 
                       standardize=FALSE, intercept=FALSE, fixed.eb="none", 
                       full.post=TRUE, init=NULL, control=control)
  
  # models with external covariates
  Z <- model.matrix(~ replace(drug.prep$stage, is.na(drug.prep$stage),
                              "experimental") + replace(drug.prep$action, is.na(
                                drug.prep$action), "unknown"))
  colnames(Z) <- c("intercept", "experimental", "in clinical development",
                   "targeted", "unknown")
  C <- lapply(inpathway, function(s) {
    s <- cbind(1, s); colnames(s) <- c("intercept", "in pathway"); return(s)})
  cv2.semnig <- semnig(x=xtrain, y=ytrain, C=C, Z=Z, unpenalized=NULL, 
                       standardize=FALSE, intercept=FALSE, fixed.eb="none", 
                       full.post=TRUE, init=NULL, control=control)
  cv3.semnig <- semnig(x=xtrain, y=ytrain, C=C, Z=NULL, unpenalized=NULL, 
                       standardize=FALSE, intercept=FALSE, fixed.eb="none", 
                       full.post=TRUE, init=NULL, control=control)
  cv4.semnig <- semnig(x=xtrain, y=ytrain, C=NULL, Z=Z, unpenalized=NULL, 
                       standardize=FALSE, intercept=FALSE, fixed.eb="none", 
                       full.post=TRUE, init=NULL, control=control)
  
  # bSEM model
  # cv1.bSEM <- bSEM(xtrain, split(ytrain, 1:ncol(ytrain)), 
  #                  lapply(inpathway, "+", 1), 
  #                   control=list(maxit=200, trace=FALSE, epsilon=1e-3))
  
  # lasso and ridge
  cv1.lasso <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[, d], intercept=FALSE, standardize=FALSE)})
  cv1.ridge <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[, d], alpha=0, intercept=FALSE,
              standardize=FALSE)})
  
  # calculating average prediction mean squared error
  pmse[r, 1] <- mean(sapply(1:D, function(d) {
    mean((ytest[, d] - xtest[[d]] %*% cv1.semnig$vb$mpost$beta[[d]])^2)}))
  pmse[r, 2] <- mean(sapply(1:D, function(d) {
    mean((ytest[, d] - xtest[[d]] %*% cv2.semnig$vb$mpost$beta[[d]])^2)}))
  pmse[r, 3] <- mean(sapply(1:D, function(d) {
    mean((ytest[, d] - xtest[[d]] %*% cv3.semnig$vb$mpost$beta[[d]])^2)}))
  pmse[r, 4] <- mean(sapply(1:D, function(d) {
    mean((ytest[, d] - xtest[[d]] %*% cv4.semnig$vb$mpost$beta[[d]])^2)}))
  pmse[r, 5] <- mean((ytest - sapply(1:D, function(d) {
    predict(cv1.lasso[[d]], xtest[[d]], s="lambda.min")}))^2)
  pmse[r, 6] <- mean((ytest - sapply(1:D, function(d) {
    predict(cv1.ridge[[d]], xtest[[d]], s="lambda.min")}))^2)
  # pmse[r, 5] <- mean(sapply(1:D, function(d) {
  #   mean((ytest[, d] - xtest[[d]] %*% cv1.bSEM$vb$beta[[d]][, 1])^2)}))
  
  # calculate ELBO for semnig models
  elbo[r, ] <- c(mean(cv1.semnig$seq.elbo[nrow(cv1.semnig$seq.elbo), ]),
                 mean(cv2.semnig$seq.elbo[nrow(cv2.semnig$seq.elbo), ]),
                 mean(cv3.semnig$seq.elbo[nrow(cv3.semnig$seq.elbo), ]),
                 mean(cv4.semnig$seq.elbo[nrow(cv4.semnig$seq.elbo), ]))
                 # mean(cv1.bSEM$mll[, ncol(cv1.bSEM$mll)]))
  
  
  # determine the ranks of the model parameter point estimates
  brank[r, 1:sum(p)] <- unlist(lapply(cv1.semnig$vb$mpost$beta, function(s) {
    rank(abs(s), ties.method="average")}))
  brank[r, (sum(p) + 1):(2*sum(p))] <- 
    unlist(lapply(cv2.semnig$vb$mpost$beta, function(s) {
      rank(abs(s), ties.method="average")}))
  brank[r, (2*sum(p) + 1):(3*sum(p))] <-
    unlist(lapply(cv3.semnig$vb$mpost$beta, function(s) {
      rank(abs(s), ties.method="average")}))
  brank[r, (3*sum(p) + 1):(4*sum(p))] <- 
    unlist(lapply(cv4.semnig$vb$mpost$beta, function(s) {
      rank(abs(s), ties.method="average")}))
  brank[r, (4*sum(p) + 1):(5*sum(p))] <- unlist(sapply(cv1.lasso, function(s) {
    rank(abs(as.numeric(coef(s, s="lambda.min"))[-1]), ties.method="average")}))
  brank[r, (5*sum(p) + 1):(6*sum(p))] <- unlist(sapply(cv1.ridge, function(s) {
    rank(abs(as.numeric(coef(s, s="lambda.min"))[-1]), ties.method="average")}))
  # brank[r, (4*sum(p) + 1):(5*sum(p))] <-
  #   unlist(lapply(cv1.bSEM$vb$beta, function(s) {
  #     rank(abs(s[, 1]), ties.method="average")}))

}

# Euclidian distance between rank vectors
brankdist <- sapply(methods, function(s) {
  dist(brank[, substr(colnames(brank), 7, nchar(s) + 6)==s])})

# combine tables and save
res <- rbind(pmse, cbind(elbo, NA, NA), brankdist)
colnames(res) <- methods
rownames(res) <- c(rep("pmse", nreps), rep("elbo", nreps), 
                   rep("brankdist", nrow(brankdist)))
write.table(res, file="results/data_gdsc_cv1.txt")
