#!/usr/bin/env Rscript

################################################################################
# regression GDSC log IC50 values on gene expression/mutations of cell lines
################################################################################

### libraries
packages <- c("foreach", "doParallel", "rstan", "glmnet", "pInc", "Hyperpar", 
              "cambridge")
sapply(packages, library, character.only=TRUE)

# number of CV folds
nfolds <- 10

# number of cores to use
ncores <- min(nfolds, 100)
# options(mc.cores=max(floor(100/ncores), 1))

### load data
load(file="data/data_gdsc_dat1.Rdata")

# reassign data
resp <- resp.prep
expr <- expr.prep
meth <- meth.prep
mut <- mut.prep

################################## analysis 1 ##################################
### data preparation
# only retain the cell lines that are available in both response and expression
expr.prep <- expr[rownames(expr) %in% rownames(resp), ]
resp.prep <- resp[rownames(resp) %in% rownames(expr), ]

# # select features with largest variance separately for in and out pathway genes
# D <- ncol(resp.prep)
# psel <- 100
# o <- order(-apply(expr.prep, 2, sd))
# idsel <- lapply(feat.prep$inpathway, function(s) {
#   m1 <- o[o %in% which(s==1)]; id1 <- head(m1, n=min(length(m1), psel/2))
#   m0 <- o[o %in% which(s==0)]; id0 <- head(m0, n=psel - length(id1))
#   c(id1, id0)})
# expr.sel <- lapply(idsel, function(s) {expr.prep[, s]})
# inpathway <- lapply(1:D, function(d) {feat.prep$inpathway[[d]][idsel[[d]]]})
#
# # create data objects used in fitting
# x <- lapply(expr.sel, function(s) {scale(s)})
# y <- scale(resp.prep)
# p <- sapply(x, ncol)
# n <- nrow(y)

# create data objects used in fitting
set.seed(2020)
D <- ncol(resp.prep)
psel <- ncol(expr.prep)
o <- order(-apply(expr.prep, 2, sd))
idsel <- which(o %in% c(1:psel))
expr.sel <- expr.prep[, idsel]
n <- nrow(resp.prep)
ntrain <- floor(n/2)
idtrain <- sample(1:n, ntrain)
xtrain <- scale(expr.sel[idtrain, ])
xtest <- scale(expr.sel[-idtrain, ])
ytrain <- scale(resp.prep[idtrain, ])
ytest <- scale(resp.prep[-idtrain, ])
if(is.list(xtrain)) {
  p <- sapply(xtrain, ncol)
} else {
  p <- ncol(xtrain)
}
target <- lapply(1:D, function(d) {
  as.numeric((feat.prep$inpathway[[d]] + feat.prep$inpathway[[d]])!=0)[idsel]})

### model fitting
# setting seed
set.seed(2020)
foldid <- sample(c(rep(1:nfolds, each=ntrain %/% nfolds),
                   rep(1:(ntrain %% nfolds), (ntrain %% nfolds)!=0)))
# foldid <- sample(c(rep(1:nfolds, each=n %/% nfolds),
#                    rep(1:(n %% nfolds), (n %% nfolds)!=0)))
# 
# # # estimation settings
# control <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, epsilon.vb=1e-3,
#                 maxit.eb=200, maxit.vb=1, maxit.post=100, maxit.block=0)
# methods <- c("NIG1", "NIG2", "NIG3", "NIG4", "NIG5", "lasso", "ridge", "bSEM")
# 
# # model without external covariates
# Z <- matrix(1, nrow=D)
# colnames(Z) <- c("intercept")
# C <- lapply(inpathway, function(s) {
#   s <- matrix(1, nrow=length(s)); colnames(s) <- c("intercept"); return(s)})
# fit1.semnig <- semnig(x=x, y=y, C=C, Z=Z, unpenalized=NULL,
#                       standardize=FALSE, intercept=FALSE, fixed.eb="none",
#                       full.post=TRUE, init=NULL, control=control)
# 
# # models with external covariates
# Z <- model.matrix(~ replace(drug.prep$stage, is.na(drug.prep$stage),
#                             "experimental") + replace(drug.prep$action, is.na(
#                               drug.prep$action), "unknown"))
# colnames(Z) <- c("intercept", "experimental", "in clinical development",
#                  "targeted", "unknown")
# C <- lapply(inpathway, function(s) {
#   s <- cbind(1, s); colnames(s) <- c("intercept", "in pathway"); return(s)})
# fit2.semnig <- semnig(x=x, y=y, C=C, Z=Z, unpenalized=NULL,
#                       standardize=FALSE, intercept=FALSE, fixed.eb=FALSE,
#                       full.post=TRUE, init=NULL, control=control)
# fit3.semnig <- semnig(x=x, y=y, C=C, Z=NULL, unpenalized=NULL,
#                       standardize=FALSE, intercept=FALSE, fixed.eb=FALSE,
#                       full.post=TRUE, init=NULL, control=control)
# fit4.semnig <- semnig(x=x, y=y, C=NULL, Z=Z, unpenalized=NULL,
#                       standardize=FALSE, intercept=FALSE, fixed.eb=FALSE,
#                       full.post=TRUE, init=NULL, control=control)
# 
# # blockwise updating of hyperparameters
# control <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, epsilon.vb=1e-3,
#                 maxit.eb=200, maxit.vb=1, maxit.post=100, maxit.block=10)
# fit5.semnig <- semnig(x=x, y=y, C=C, Z=Z, unpenalized=NULL,
#                       standardize=FALSE, intercept=FALSE, fixed.eb=FALSE,
#                       full.post=TRUE, init=NULL, control=control)

# EB ridge
Z <- cbind(model.matrix(~ 0 + replace(drug.prep$stage, is.na(drug.prep$stage), 
                                      "experimental")),
           model.matrix(~ 0 + replace(drug.prep$action, 
                                      is.na(drug.prep$action), "unknown")))
colnames(Z) <- c("clinically approved", "experimental", 
                 "in clinical development", "cytotoxic", "targeted", 
                 "unknown")
Z <- scale(Z, scale=FALSE)
C <- lapply(target, function(s) {
  s <- scale(as.matrix(s), scale=FALSE); colnames(s) <- "target"; s})
fit1.ebridge <- ebridge(xtrain, ytrain, C, Z, mult.lambda=TRUE, foldid=foldid,
                        hyper=list(lambda=NULL, zeta=0, nu=0),
                        control=list(epsilon=sqrt(.Machine$double.eps), 
                                     maxit=500, trace=TRUE, glmnet.fit2=FALSE))
Z <- cbind(model.matrix(~ 0 + replace(drug.prep$stage, is.na(drug.prep$stage), 
                                      "experimental")),
           model.matrix(~ 0 + replace(drug.prep$action, 
                                      is.na(drug.prep$action), "unknown")))
colnames(Z) <- c("clinically approved", "experimental", 
                 "in clinical development", "cytotoxic", "targeted", 
                 "unknown")
C <- lapply(target, function(s) {
  s <- as.matrix(s); colnames(s) <- "target"; s})
fit2.ebridge <- ebridge(xtrain, ytrain, C, Z, mult.lambda=TRUE, foldid=foldid,
                        hyper=list(lambda=fit1.ebridge$lambda, zeta=0, nu=0),
                        control=list(epsilon=sqrt(.Machine$double.eps), 
                                     maxit=500, trace=TRUE, glmnet.fit2=FALSE))
save(fit1.ebridge, fit2.ebridge, file="results/analysis_gdsc_fit1.2.Rdata")
# load(file="results/analysis_gdsc_fit1.2.Rdata")


colnames(Z)
fit1.ebridge$alphad
fit2.ebridge$alphad
fit1.ebridge$alphaf
fit2.ebridge$alphaf

hist(fit1.ebridge$tau^2)
hist(unlist(fit1.ebridge$gamma)^2, breaks=40)
str(fit1.ebridge$gamma)

pred1.ebridge <- sapply(fit1.ebridge$beta1, function(b) {xtest %*% b})
pred2.ebridge <- sapply(fit1.ebridge$beta2, function(b) {xtest %*% b})
pred3.ebridge <- sapply(fit2.ebridge$beta1, function(b) {xtest %*% b})
pred4.ebridge <- sapply(fit2.ebridge$beta2, function(b) {xtest %*% b})
pred1.ridge <- sapply(fit1.ebridge$glmnet.fit1, function(s) {
  predict(s, xtest, s="lambda.min")})

pmse1.ebridge <- colMeans((pred1.ebridge - ytest)^2)
pmse2.ebridge <- colMeans((pred2.ebridge - ytest)^2)
pmse3.ebridge <- colMeans((pred3.ebridge - ytest)^2)
pmse4.ebridge <- colMeans((pred4.ebridge - ytest)^2)
pmse1.ridge <- colMeans((pred1.ridge - ytest)^2)
pmse.null <- colMeans(ytest^2)
sort(c(ebridge1=mean(pmse1.ebridge), ebridge2=mean(pmse2.ebridge), 
  ebridge3=mean(pmse3.ebridge), ebridge4=mean(pmse4.ebridge), 
  ridge=mean(pmse1.ridge), null=mean(pmse.null)))

mean(pmse1.ridge)
mean(pmse1.ebridge)
mean(pmse2.ebridge)
plot((pmse.null - pmse1.ridge)/pmse.null, cex=0.5, pch=16, 
     ylim=range((pmse.null - c(pmse1.ebridge, pmse1.ridge))/pmse.null),
     xlab="Drug", ylab="MSE reduction as fraction of null")
points((pmse.null - pmse1.ebridge)/pmse.null, cex=0.5, pch=16, col=2)
abline(v=order(abs(pmse1.ridge - pmse1.ebridge)/pmse.null, 
               decreasing=TRUE)[c(1:5)], lty=2)
legend("topleft", c("CV", "CV + EB"), pch=16, col=c(1, 2))

id.max <- which.max(abs(pmse1.ridge - pmse2.ebridge)/pmse.null)
plot(fit1.ebridge$beta2[[id.max]], 
     ylim=range(c(fit1.ebridge$beta2[[id.max]],
                  coef(fit1.ebridge$glmnet.fit1[[id.max]], 
                       s="lambda.min")[-1, ])), cex=0.5, pch=16)
points(coef(fit1.ebridge$glmnet.fit1[[id.max]], s="lambda.min")[-1, ],
       cex=0.5, pch=16, col=2)

plot(fit1.ebridge$beta2[[id.max]], ylim=range(fit1.ebridge$beta2[[id.max]]), 
     cex=0.5, pch=16, col=c(1:2)[as.numeric(target[[id.max]]==1) + 1])

pmse1.ridge[id.max]
pmse1.ebridge[id.max]

# bSEM model
fit1.bSEM <- bSEM(x, split(y, 1:ncol(y)), lapply(inpathway, "+", 1),
                  control=list(maxit=200, trace=TRUE, epsilon=1e-3))

# penalized regression models
fit1.lasso <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[, d], foldid=foldid, intercept=FALSE)})
fit1.ridge <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[, d], alpha=0, foldid=foldid, intercept=FALSE)})

# models with external covariates
Z <- model.matrix(~ replace(drug.prep$stage, is.na(drug.prep$stage),
                            "experimental") + replace(drug.prep$action, is.na(
                              drug.prep$action), "unknown"))[, -1]
colnames(Z) <- c("experimental", "in clinical development", "targeted", 
                 "unknown")
C <- lapply(inpathway, function(s) {
  s <- as.matrix(s); colnames(s) <- c("in pathway"); return(s)})
var.scale <- 1/(n*sapply(fit1.ridge, "[[", "lambda.min"))
fit1.semnig <- semnig(x=x, y=y, C=C, Z=Z, unpenalized=NULL,
                      standardize=FALSE, intercept=FALSE, fixed.eb="none",
                      full.post=TRUE, init=NULL, var.scale=var.scale, 
                      control=control)

# # saving fitted model objects
# save(fit1.semnig, fit2.semnig, fit3.semnig, fit4.semnig, fit5.semnig, 
#      fit1.lasso, fit1.ridge,
#      fit1.bSEM, file="results/analysis_gdsc_fit1.Rdata")

test1 <- sapply(fit1.ebridge$glmnet.fit1, function(s) {coef(s, s="lambda.min")[-1, ]})
test2 <- sapply(fit1.ridge, function(s) {coef(s, s="lambda.min")[-1, ]})
test3 <- sapply(fit1.ebridge$glmnet.fit2, function(s) {coef(s, s="lambda.min")[-1, ]})

pairs(cbind(c(test1), c(test2), c(test3)), panel=function(x, y) {
  points(x, y); abline(a=0, b=1)})


# saving fitted model objects
save(fit1.ebridge,
     fit1.lasso, fit1.ridge,
     fit1.bSEM, file="results/analysis_gdsc_fit1.2.Rdata")

# EB estimates
tab <- rbind(c(fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 1]/
                 fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 1],
               fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 2]/
                 fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 2] -
                 fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 1]/
                 fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 1]))
colnames(tab) <- c(colnames(C[[1]]), colnames(Z), "lambdaf", "lambdad")
rownames(tab) <- methods[c(1:5, 7)]
write.table(tab, file="results/analysis_gdsc_fit1.txt")

### cross-validation of performance measures
# settings seed
set.seed(2020)

# estimation settings
# lasso <- stan_model("code/lasso.stan", auto_write=TRUE)
# ridge <- stan_model("code/ridge.stan", auto_write=TRUE)
# methods <- c(paste0("NIG", 1:9), paste0("lasso", 1:4), paste0("ridge", 1:5), 
#              "bSEM")
methods <- c(paste0("EBridge", 1:2), paste0("lasso", 1), paste0("ridge", 1), 
             "bSEM")
foldid <- sample(c(rep(1:nfolds, each=n %/% nfolds), 
                   rep(1:(n %% nfolds), (n %% nfolds)!=0)))

# setup cluster
if(Sys.info()["sysname"]=="Windows") {
  cl <- makeCluster(ncores, type="PSOCK")
} else {
  cl <- makeCluster(ncores, type="FORK")
}
if(ncores > 1) {
  registerDoParallel(cl)
} else {
  registerDoSEQ()
}
res <- foreach(r=1:nfolds, .packages=packages, .errorhandling="pass", 
               .verbose=TRUE) %dopar% {
  cat("\r", "fold", r)
  set.seed(2020 + r)
  
  # splitting data
  ntrain <- sum(foldid!=r)
  xtrain <- lapply(x, function(s) {
    scale(as.matrix(s[foldid!=r, ], nrow=ntrain))})
  ytrain <- scale(as.matrix(y[foldid!=r, ], nrow=ntrain))
  xtest <- lapply(x, function(s) {
    scale(as.matrix(s[foldid==r, ], nrow=n - ntrain))})
  ytest <- scale(as.matrix(y[foldid==r, ], nrow=n - ntrain))
  
  # # model without external covariates
  # control <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-3, epsilon.vb=1e-3, 
  #                 maxit.eb=200, maxit.vb=1, maxit.post=100, maxit.block=0)
  # Z <- matrix(1, nrow=D)
  # colnames(Z) <- c("intercept")
  # C <- lapply(inpathway, function(s) {
  #   s <- matrix(1, nrow=length(s)); colnames(s) <- c("intercept"); return(s)})
  # cv1.semnig <- semnig(x=xtrain, y=ytrain, C=C, Z=Z, unpenalized=NULL, 
  #                      standardize=FALSE, intercept=FALSE, fixed.eb="none", 
  #                      full.post=TRUE, init=NULL, control=control)
  #
  # # models with external covariates
  # Z <- model.matrix(~ replace(drug.prep$stage, is.na(drug.prep$stage),
  #                             "experimental") + replace(drug.prep$action, is.na(
  #                               drug.prep$action), "unknown"))
  # colnames(Z) <- c("intercept", "experimental", "in clinical development",
  #                  "targeted", "unknown")
  # C <- lapply(inpathway, function(s) {
  #   s <- cbind(1, s); colnames(s) <- c("intercept", "in pathway"); return(s)})
  # cv2.semnig <- semnig(x=xtrain, y=ytrain, C=C, Z=Z, unpenalized=NULL,
  #                      standardize=FALSE, intercept=FALSE, fixed.eb="none",
  #                      full.post=TRUE, init=NULL, control=control)
  # cv3.semnig <- semnig(x=xtrain, y=ytrain, C=C, Z=NULL, unpenalized=NULL, 
  #                      standardize=FALSE, intercept=FALSE, fixed.eb="none", 
  #                      full.post=TRUE, init=NULL, control=control)
  # cv4.semnig <- semnig(x=xtrain, y=ytrain, C=NULL, Z=Z, unpenalized=NULL, 
  #                      standardize=FALSE, intercept=FALSE, fixed.eb="none", 
  #                      full.post=TRUE, init=NULL, control=control)
  # 
  # 
  # 
  # # blockwise updating of hyperparameters
  # control <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, epsilon.vb=1e-3, 
  #                 maxit.eb=1, maxit.vb=1, maxit.post=100, maxit.block=10)
  # cv5.semnig <- semnig(x=xtrain, y=ytrain, C=C, Z=Z, unpenalized=NULL,
  #                      standardize=FALSE, intercept=FALSE, fixed.eb=FALSE,
  #                      full.post=TRUE, init=NULL, control=control)
  # 
  # # optimisation of MAP
  # cv6.semnig <- lapply(1:D, function(d) {
  #   optimizing(stanmodels$nig, 
  #              data=list(p=p[d], n=ntrain, x=xtrain[[d]], y=ytrain[, d], 
  #                        phi=cv2.semnig$eb$mpriorf[[d]], 
  #                        lambdaf=cv2.semnig$eb$lambdaf, 
  #                        chi=cv2.semnig$eb$mpriord[d], 
  #                        lambdad=cv2.semnig$eb$lambdad))})
  # 
  # phi <- seq(0.01, 10, length.out=10)
  # chi <- seq(0.01, 10, length.out=10)
  # lambdaf <- cv2.semnig$eb$lambdaf
  # lambdad <- cv2.semnig$eb$lambdad
  # cv7.semnig <- cv.semnig(x=xtrain, y=ytrain, nfolds=5, foldid=NULL, 
  #                         seed=NULL, phi=phi, chi=chi, lambdaf=lambdaf, 
  #                         lambdad=lambdad, type.measure="mse", 
  #                         control=list(trace=FALSE))
  # cv8.semnig <- lapply(1:D, function(d) {
  #   sampling(stanmodels$nig, 
  #            data=list(p=p[d], n=ntrain, x=xtrain[[d]], y=ytrain[, d], 
  #                      phi=cv2.semnig$eb$mpriorf[[d]], 
  #                      lambdaf=cv2.semnig$eb$lambdaf, 
  #                      chi=cv2.semnig$eb$mpriord[d], 
  #                      lambdad=cv2.semnig$eb$lambdad),
  #            chains=1, iter=1000, warmup=500, verbose=FALSE, refresh=0)})
  
  # bSEM model
  cv1.bSEM <- bSEM(xtrain, split(ytrain, 1:ncol(ytrain)),
                   lapply(inpathway, "+", 1),
                   control=list(maxit=200, trace=FALSE, epsilon=1e-3))
  
  # EB ridge
  replace(drug.prep$stage, is.na(drug.prep$stage), "experimental")
  Z <- cbind(model.matrix(~ 0 + replace(drug.prep$stage, is.na(drug.prep$stage), 
                                        "experimental")),
             model.matrix(~ 0 + replace(drug.prep$action, 
                                        is.na(drug.prep$action), "unknown")))
  colnames(Z) <- c("clinically approved", "experimental", 
                   "in clinical development", "cytotoxic", "targeted", 
                   "unknown")
  Z <- scale(Z, scale=FALSE)
  C <- lapply(inpathway, function(s) {
    s <- scale(as.matrix(s), scale=FALSE); colnames(s) <- "inpathway"; s})
  cv1.ebridge <- ebridge(xtrain, ytrain, C, Z, mult.lambda=TRUE,
                         hyper=list(lambda=NULL, zeta=0, nu=0),
                         control=list(epsilon=sqrt(.Machine$double.eps), 
                                      maxit=500, trace=FALSE))

  # penalized regression models
  cv1.lasso <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[, d], intercept=FALSE)})
  # cv2.lasso <- lapply(1:D, function(d) {
  #   sampling(lasso, data=list(p=p[d], n=ntrain, x=xtrain[[d]], y=ytrain[, d], 
  #                             a0=0.001, b0=0.001),
  #            chains=1, iter=2000, warmup=1000, verbose=FALSE, refresh=0)})
  # cv3.lasso <- cv.glmnet(xtrain[[1]], ytrain, family="mgaussian", 
  #                        intercept=FALSE, standardize=FALSE)
  
  cv1.ridge <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[, d], alpha=0, intercept=FALSE)})
  # cv2.ridge <- lapply(1:D, function(d) {
  #   sampling(ridge, data=list(p=p[d], n=ntrain, x=xtrain[[d]], y=ytrain[, d], 
  #                             a0=0.001, b0=0.001),
  #            chains=1, iter=1000, warmup=500, verbose=FALSE, refresh=0)})
  # cv3.ridge <- cv.glmnet(xtrain[[1]], ytrain, family="mgaussian", alpha=0, 
  #                        intercept=FALSE, standardize=FALSE)
  # cv4.ridge <- lapply(1:D, function(d) {
  #   BayesRidge(ytrain[, d], xtrain[[d]], plotML=FALSE, a=0.001, b=0.001, 
  #              SVD=NULL)})
  
  # models with external covariates
  Z <- model.matrix(~ replace(drug.prep$stage, is.na(drug.prep$stage),
                              "experimental") + replace(drug.prep$action, is.na(
                                drug.prep$action), "unknown"))[, -1]
  colnames(Z) <- c("experimental", "in clinical development", "targeted", 
                   "unknown")
  C <- lapply(inpathway, function(s) {
    s <- as.matrix; colnames(s) <- c("in pathway"); return(s)})
  cv1.semnig <- semnig(x=xtrain, y=ytrain, C=C, Z=Z, unpenalized=NULL,
                       standardize=FALSE, intercept=FALSE, fixed.eb="none",
                       full.post=TRUE, init=NULL, control=control)
  
  # # extracting point estimates
  # best <- list(Reduce("cbind", cv1.semnig$vb$mpost$beta),
  #              Reduce("cbind", cv2.semnig$vb$mpost$beta),
  #              Reduce("cbind", cv3.semnig$vb$mpost$beta),
  #              Reduce("cbind", cv4.semnig$vb$mpost$beta),
  #              Reduce("cbind", cv5.semnig$vb$mpost$beta),
  #              sapply(cv6.semnig, function(s) {
  #                s$par[names(s$par) %in% paste0("beta[", 1:p[1], "]")]}),
  #              sapply(cv7.semnig, function(s) {
  #                s$fit$par[names(s$fit$par) %in% 
  #                            paste0("beta[", 1:p[1], "]")]}),
  #              sapply(cv8.semnig, get_posterior_mean, 
  #                     pars=paste0("beta[", 1:p[1], "]")),
  #              sapply(cv8.semnig, function(s) {
  #                sapply(extract(s, pars=paste0("beta[", 1:p[1], "]")), 
  #                       function(b) {d <- density(b); d$x[which.max(d$y)]})}),
  #              sapply(cv1.lasso, function(s) {
  #                as.numeric(coef(s, s="lambda.min"))[-1]}),
  #              sapply(cv2.lasso, get_posterior_mean, 
  #                     pars=paste0("beta[", 1:p[1], "]")),
  #              sapply(cv2.lasso, function(s) {
  #                sapply(extract(s, pars=paste0("beta[", 1:p[1], "]")), 
  #                       function(b) {d <- density(b); d$x[which.max(d$y)]})}),
  #              as.matrix(Reduce("cbind", 
  #                               coef(cv3.lasso, s="lambda.min")))[-1, ],
  #              sapply(cv1.ridge, function(s) {
  #                as.numeric(coef(s, s="lambda.min"))[-1]}),
  #              sapply(cv2.ridge, get_posterior_mean, 
  #                     pars=paste0("beta[", 1:p[1], "]")),
  #              sapply(cv2.ridge, function(s) {
  #                sapply(extract(s, pars=paste0("beta[", 1:p[1], "]")), 
  #                       function(b) {d <- density(b); d$x[which.max(d$y)]})}),
  #              as.matrix(Reduce("cbind", 
  #                               coef(cv3.ridge, s="lambda.min")))[-1, ],
  #              sapply(cv4.ridge, "[[", "postmeanbeta"),
  #              sapply(cv1.bSEM$vb$beta, "[", 1:p[1]))
  # 
  
  best <- list(cv1.ebridge$beta1,
               cv1.ebridge$beta2,
               lapply(cv1.lasso, function(s) {
                 as.numeric(coef(s, s="lambda.min"))[-1]}),
               lapply(cv1.ridge, function(s) {
                 as.numeric(coef(s, s="lambda.min"))[-1]}),
               lapply(cv1.bSEM$vb$beta, "[", 1:p[1]))
  
  # calculating average prediction mean squared error
  pmse <- cbind(sapply(best, function(b) {
    sapply(1:D, function(d) {mean((ytest[, d] - xtest[[d]] %*% b[[d]])^2)})}), 
    colMeans((ytest - colMeans(ytest))^2))
  
  # # calculate ELBO for semnig models on training and test data
  # elbot <- c(mean(cv1.semnig$seq.elbo[nrow(cv1.semnig$seq.elbo), ]),
  #            mean(cv2.semnig$seq.elbo[nrow(cv2.semnig$seq.elbo), ]),
  #            mean(cv3.semnig$seq.elbo[nrow(cv3.semnig$seq.elbo), ]),
  #            mean(cv4.semnig$seq.elbo[nrow(cv4.semnig$seq.elbo), ]),
  #            mean(cv5.semnig$seq.elbo[nrow(cv5.semnig$seq.elbo), ]),
  #            mean(cv1.bSEM$mll[, ncol(cv1.bSEM$mll)]))
  # elbo <- c(mean(new.elbo(cv1.semnig, xtest, ytest)),
  #           mean(new.elbo(cv2.semnig, xtest, ytest)),
  #           mean(new.elbo(cv3.semnig, xtest, ytest)),
  #           mean(new.elbo(cv4.semnig, xtest, ytest)),
  #           mean(new.elbo(cv5.semnig, xtest, ytest)))
  # 
  # # log pseudo marginal log likelihood
  # lpml <- c(sum(logcpo(xtest, ytest, ntrain, cv1.semnig), na.rm=TRUE),
  #           sum(logcpo(xtest, ytest, ntrain, cv2.semnig), na.rm=TRUE),
  #           sum(logcpo(xtest, ytest, ntrain, cv3.semnig), na.rm=TRUE),
  #           sum(logcpo(xtest, ytest, ntrain, cv4.semnig), na.rm=TRUE),
  #           sum(logcpo(xtest, ytest, ntrain, cv5.semnig), na.rm=TRUE))
  
  # determine the ranks of the model parameter point estimates
  brank <- c(sapply(best, function(s) {rank(unlist(s), ties.method="average")}))
  
  names(brank) <- paste0(
    "brank.", paste0(rep(methods, each=sum(p)), ".",
                     rep(paste0(rep(1:D, times=p), ".",
                                unlist(sapply(p, function(s) {
                                  return(1:s)}))), length(methods))))

  list(pmse=pmse, 
       # elbo=elbo, elbot=elbot, lpml=lpml, 
       brank=brank)
}
stopCluster(cl=cl)

errorid <- sapply(res, function(s) {is.numeric(s[[1]])})
res2 <- res[errorid]

# prepare and save results table
pmse <- Reduce("rbind", sapply(res2, function(s) {as.list(s["pmse"])}))
# elbo <- t(sapply(res2, "[[", "elbo"))
# elbot <- t(sapply(res2, "[[", "elbot"))
# lpml <- t(sapply(res2, "[[", "lpml"))
brank <- t(sapply(res2, "[[", "brank"))

# Euclidian distance between rank vectors
brankdist <- sapply(methods, function(s) {
  dist(brank[, substr(colnames(brank), 7, nchar(s) + 6)==s])})

# combine tables and save
res <- rbind(pmse, 
             # cbind(elbo, matrix(NA, nrow=nfolds, ncol=15)), 
             # cbind(elbot, matrix(NA, nrow=nfolds, ncol=14)), 
             # cbind(lpml, matrix(NA, nrow=nfolds, ncol=15)), 
             cbind(brankdist, NA))
colnames(res) <- c(methods, "null")
rownames(res) <- c(paste0(rep("pmse", nfolds*D), ".drug", 
                          rep(1:D, times=nfolds)), 
                   # rep("elbo", nfolds), rep("elbot", nfolds),
                   # rep("lpml", nfolds), 
                   rep("brankdist", nrow(brankdist)))
write.table(res, file="results/analysis_gdsc_res1.2.txt")


################################## analysis 2 ##################################
### data preparation
# only retain the cell lines that are available in both response and expression
expr.prep <- expr[rownames(expr) %in% rownames(resp), ]
resp.prep <- resp[rownames(resp) %in% rownames(expr), ]

D <- ncol(resp.prep)
psel <- 100
# o <- order(-colSums(Reduce("rbind", feat.prep$inpathway)))
o <- order(-apply(expr.prep, 2, sd))
idsel <- rep(list(o[1:psel]), D)
expr.sel <- lapply(idsel, function(s) {expr.prep[, s]})
inpathway <- lapply(1:D, function(d) {feat.prep$inpathway[[d]][idsel[[d]]]})

# create data objects used in fitting
x <- lapply(expr.sel, function(s) {scale(s)})
y <- scale(resp.prep)
p <- sapply(x, ncol)
n <- nrow(y)

### model fitting
# setting seed
set.seed(2019)

# estimation settings
control <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, epsilon.vb=1e-3,
                maxit.eb=200, maxit.vb=1, maxit.post=100, maxit.block=0)
methods <- c("NIG1", "NIG2", "NIG3", "NIG4", "NIG5", "lasso", "ridge", "bSEM")

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

# blockwise updating of hyperparameters
control <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, epsilon.vb=1e-3,
                maxit.eb=200, maxit.vb=1, maxit.post=100, maxit.block=10)
fit5.semnig <- semnig(x=x, y=y, C=C, Z=Z, unpenalized=NULL,
                      standardize=FALSE, intercept=FALSE, fixed.eb=FALSE,
                      full.post=TRUE, init=NULL, control=control)

# bSEM model
fit1.bSEM <- bSEM(x, split(y, 1:ncol(y)), lapply(inpathway, "+", 1),
                  control=list(maxit=200, trace=TRUE, epsilon=1e-3))

# penalized regression models
fit1.lasso <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[, d], intercept=FALSE, standardize=FALSE)})
fit1.ridge <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[, d], alpha=0, intercept=FALSE, standardize=FALSE)})

# saving fitted model objects
save(fit1.semnig, fit2.semnig, fit3.semnig, fit4.semnig, fit5.semnig, 
     fit1.lasso, fit1.ridge,
     fit1.bSEM, file="results/analysis_gdsc_fit2.Rdata")

# EB estimates
tab <- rbind(c(fit1.semnig$eb$alphaf, NA, fit1.semnig$eb$alphad, rep(NA, 4), 
               fit1.semnig$eb$lambdaf, fit1.semnig$eb$lambdad),
             c(fit2.semnig$eb$alphaf, fit2.semnig$eb$alphad, 
               fit2.semnig$eb$lambdaf, fit2.semnig$eb$lambdad),
             c(fit3.semnig$eb$alphaf, rep(NA, 5), fit3.semnig$eb$lambdaf, NA),
             c(rep(NA, 2), fit4.semnig$eb$alphad, NA, fit2.semnig$eb$lambdad),
             c(fit5.semnig$eb$alphaf, fit5.semnig$eb$alphad, 
               fit5.semnig$eb$lambdaf, fit5.semnig$eb$lambdad),
             c(fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 1]/
                 fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 1],
               fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 2]/
                 fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 2] -
                 fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 1]/
                 fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 1], rep(NA, 7)))
colnames(tab) <- c(colnames(C[[1]]), colnames(Z), "lambdaf", "lambdad")
rownames(tab) <- methods[c(1:5, 7)]
write.table(tab, file="results/analysis_gdsc_fit2.txt")

### cross-validation of performance measures
# settings seed
set.seed(2020)

# estimation settings
lasso <- stan_model("code/lasso.stan", auto_write=TRUE)
ridge <- stan_model("code/ridge.stan", auto_write=TRUE)
methods <- c(paste0("NIG", 1:9), paste0("lasso", 1:4), paste0("ridge", 1:5), 
             "bSEM")
foldid <- sample(c(rep(1:nfolds, each=n %/% nfolds), 
                   rep(1:(n %% nfolds), (n %% nfolds)!=0)))

# setup cluster
if(Sys.info()["sysname"]=="Windows") {
  cl <- makeCluster(ncores, type="PSOCK")
} else {
  cl <- makeCluster(ncores, type="FORK")
}
if(ncores > 1) {
  registerDoParallel(cl)
} else {
  registerDoSEQ()
}
res <- foreach(r=1:nfolds, .packages=packages, .errorhandling="pass", 
               .verbose=TRUE) %dopar% {
  cat("\r", "fold", r)
  set.seed(2020 + r)
  
  # splitting data
  ntrain <- sum(foldid!=r)
  xtrain <- lapply(x, function(s) {
   scale(as.matrix(s[foldid!=r, ], nrow=ntrain))})
  ytrain <- scale(as.matrix(y[foldid!=r, ], nrow=ntrain))
  xtest <- lapply(x, function(s) {
   scale(as.matrix(s[foldid==r, ], nrow=n - ntrain))})
  ytest <- scale(as.matrix(y[foldid==r, ], nrow=n - ntrain))
  
  # model without external covariates
  control <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-3, epsilon.vb=1e-3, 
                 maxit.eb=200, maxit.vb=1, maxit.post=100, maxit.block=0)
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
  
  # blockwise updating of hyperparameters
  control <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, epsilon.vb=1e-3, 
                 maxit.eb=1, maxit.vb=1, maxit.post=100, maxit.block=10)
  cv5.semnig <- semnig(x=xtrain, y=ytrain, C=C, Z=Z, unpenalized=NULL,
                      standardize=FALSE, intercept=FALSE, fixed.eb=FALSE,
                      full.post=TRUE, init=NULL, control=control)
  
  # optimisation of MAP
  cv6.semnig <- lapply(1:D, function(d) {
   optimizing(stanmodels$nig, 
              data=list(p=p[d], n=ntrain, x=xtrain[[d]], y=ytrain[, d], 
                        phi=cv2.semnig$eb$mpriorf[[d]], 
                        lambdaf=cv2.semnig$eb$lambdaf, 
                        chi=cv2.semnig$eb$mpriord[d], 
                        lambdad=cv2.semnig$eb$lambdad))})
  
  phi <- seq(0.01, 10, length.out=10)
  chi <- seq(0.01, 10, length.out=10)
  lambdaf <- cv2.semnig$eb$lambdaf
  lambdad <- cv2.semnig$eb$lambdad
  cv7.semnig <- cv.semnig(x=xtrain, y=ytrain, nfolds=5, foldid=NULL, 
                         seed=NULL, phi=phi, chi=chi, lambdaf=lambdaf, 
                         lambdad=lambdad, type.measure="mse", 
                         control=list(trace=FALSE))
  cv8.semnig <- lapply(1:D, function(d) {
   sampling(stanmodels$nig, 
            data=list(p=p[d], n=ntrain, x=xtrain[[d]], y=ytrain[, d], 
                      phi=cv2.semnig$eb$mpriorf[[d]], 
                      lambdaf=cv2.semnig$eb$lambdaf, 
                      chi=cv2.semnig$eb$mpriord[d], 
                      lambdad=cv2.semnig$eb$lambdad),
            chains=1, iter=1000, warmup=500, verbose=FALSE, refresh=0)})
  
  # bSEM model
  cv1.bSEM <- bSEM(xtrain, split(ytrain, 1:ncol(ytrain)),
                  lapply(inpathway, "+", 1),
                  control=list(maxit=200, trace=FALSE, epsilon=1e-3))
  
  # penalized regression models
  cv1.lasso <- lapply(1:D, function(d) {
   cv.glmnet(xtrain[[d]], ytrain[, d], intercept=FALSE, standardize=FALSE)})
  cv2.lasso <- lapply(1:D, function(d) {
   sampling(lasso, data=list(p=p[d], n=ntrain, x=xtrain[[d]], y=ytrain[, d], 
                             a0=0.001, b0=0.001),
            chains=1, iter=2000, warmup=1000, verbose=FALSE, refresh=0)})
  cv3.lasso <- cv.glmnet(xtrain[[1]], ytrain, family="mgaussian", 
                        intercept=FALSE, standardize=FALSE)
  
  cv1.ridge <- lapply(1:D, function(d) {
   cv.glmnet(xtrain[[d]], ytrain[, d], alpha=0, intercept=FALSE,
             standardize=FALSE)})
  cv2.ridge <- lapply(1:D, function(d) {
   sampling(ridge, data=list(p=p[d], n=ntrain, x=xtrain[[d]], y=ytrain[, d], 
                             a0=0.001, b0=0.001),
            chains=1, iter=1000, warmup=500, verbose=FALSE, refresh=0)})
  cv3.ridge <- cv.glmnet(xtrain[[1]], ytrain, family="mgaussian", alpha=0, 
                        intercept=FALSE, standardize=FALSE)
  cv4.ridge <- lapply(1:D, function(d) {
   BayesRidge(ytrain[, d], xtrain[[d]], plotML=FALSE, a=0.001, b=0.001, 
              SVD=NULL)})
  
  # extracting point estimates
  best <- list(Reduce("cbind", cv1.semnig$vb$mpost$beta),
              Reduce("cbind", cv2.semnig$vb$mpost$beta),
              Reduce("cbind", cv3.semnig$vb$mpost$beta),
              Reduce("cbind", cv4.semnig$vb$mpost$beta),
              Reduce("cbind", cv5.semnig$vb$mpost$beta),
              sapply(cv6.semnig, function(s) {
                s$par[names(s$par) %in% paste0("beta[", 1:p[1], "]")]}),
              sapply(cv7.semnig, function(s) {
                s$fit$par[names(s$fit$par) %in% 
                            paste0("beta[", 1:p[1], "]")]}),
              sapply(cv8.semnig, get_posterior_mean, 
                     pars=paste0("beta[", 1:p[1], "]")),
              sapply(cv8.semnig, function(s) {
                sapply(extract(s, pars=paste0("beta[", 1:p[1], "]")), 
                       function(b) {d <- density(b); d$x[which.max(d$y)]})}),
              sapply(cv1.lasso, function(s) {
                as.numeric(coef(s, s="lambda.min"))[-1]}),
              sapply(cv2.lasso, get_posterior_mean, 
                     pars=paste0("beta[", 1:p[1], "]")),
              sapply(cv2.lasso, function(s) {
                sapply(extract(s, pars=paste0("beta[", 1:p[1], "]")), 
                       function(b) {d <- density(b); d$x[which.max(d$y)]})}),
              as.matrix(Reduce("cbind", 
                               coef(cv3.lasso, s="lambda.min")))[-1, ],
              sapply(cv1.ridge, function(s) {
                as.numeric(coef(s, s="lambda.min"))[-1]}),
              sapply(cv2.ridge, get_posterior_mean, 
                     pars=paste0("beta[", 1:p[1], "]")),
              sapply(cv2.ridge, function(s) {
                sapply(extract(s, pars=paste0("beta[", 1:p[1], "]")), 
                       function(b) {d <- density(b); d$x[which.max(d$y)]})}),
              as.matrix(Reduce("cbind", 
                               coef(cv3.ridge, s="lambda.min")))[-1, ],
              sapply(cv4.ridge, "[[", "postmeanbeta"),
              sapply(cv1.bSEM$vb$beta, "[", 1:p[1]))
  
  # calculating average prediction mean squared error
  pmse <- cbind(sapply(best, function(b) {
   colMeans((ytest - xtest[[1]] %*% b)^2)}), 
   sapply(1:D, function(d) {mean((ytest[, d] - colMeans(ytest)[d])^2)}))
  
  # calculate ELBO for semnig models on training and test data
  elbot <- c(mean(cv1.semnig$seq.elbo[nrow(cv1.semnig$seq.elbo), ]),
            mean(cv2.semnig$seq.elbo[nrow(cv2.semnig$seq.elbo), ]),
            mean(cv3.semnig$seq.elbo[nrow(cv3.semnig$seq.elbo), ]),
            mean(cv4.semnig$seq.elbo[nrow(cv4.semnig$seq.elbo), ]),
            mean(cv5.semnig$seq.elbo[nrow(cv5.semnig$seq.elbo), ]),
            mean(cv1.bSEM$mll[, ncol(cv1.bSEM$mll)]))
  elbo <- c(mean(new.elbo(cv1.semnig, xtest, ytest)),
           mean(new.elbo(cv2.semnig, xtest, ytest)),
           mean(new.elbo(cv3.semnig, xtest, ytest)),
           mean(new.elbo(cv4.semnig, xtest, ytest)),
           mean(new.elbo(cv5.semnig, xtest, ytest)))
  
  # log pseudo marginal log likelihood
  lpml <- c(sum(logcpo(xtest, ytest, ntrain, cv1.semnig), na.rm=TRUE),
           sum(logcpo(xtest, ytest, ntrain, cv2.semnig), na.rm=TRUE),
           sum(logcpo(xtest, ytest, ntrain, cv3.semnig), na.rm=TRUE),
           sum(logcpo(xtest, ytest, ntrain, cv4.semnig), na.rm=TRUE),
           sum(logcpo(xtest, ytest, ntrain, cv5.semnig), na.rm=TRUE))
  
  # determine the ranks of the model parameter point estimates
  brank <- c(sapply(best, rank, ties.method="average"))
  
  names(brank) <- paste0(
   "brank.", paste0(rep(methods, each=sum(p)), ".",
                    rep(paste0(rep(1:D, times=p), ".",
                               unlist(sapply(p, function(s) {
                                 return(1:s)}))), length(methods))))
  
  list(pmse=pmse, elbo=elbo, elbot=elbot, lpml=lpml, brank=brank)
}
stopCluster(cl=cl)

errorid <- sapply(res, function(s) {is.numeric(s[[1]])})
res2 <- res[errorid]

# prepare and save results table
pmse <- Reduce("rbind", sapply(res2, function(s) {as.list(s["pmse"])}))
elbo <- t(sapply(res2, "[[", "elbo"))
elbot <- t(sapply(res2, "[[", "elbot"))
lpml <- t(sapply(res2, "[[", "lpml"))
brank <- t(sapply(res2, "[[", "brank"))

# Euclidian distance between rank vectors
brankdist <- sapply(methods, function(s) {
  dist(brank[, substr(colnames(brank), 7, nchar(s) + 6)==s])})

# combine tables and save
res <- rbind(pmse, cbind(elbo, matrix(NA, nrow=nfolds, ncol=15)), 
             cbind(elbot, matrix(NA, nrow=nfolds, ncol=14)), 
             cbind(lpml, matrix(NA, nrow=nfolds, ncol=15)), 
             cbind(brankdist, NA))
colnames(res) <- c(methods, "null")
rownames(res) <- c(paste0(rep("pmse", nfolds*D), ".drug", 
                          rep(1:D, times=nfolds)), 
                   rep("elbo", nfolds), rep("elbot", nfolds),
                   rep("lpml", nfolds), rep("brankdist", nrow(brankdist)))
write.table(res, file="results/analysis_gdsc_res2.txt")



################################## analysis 3 ##################################
### data preparation
# only retain the cell lines that are available in both response and expression
expr.prep <- expr[rownames(expr) %in% rownames(resp), ]
resp.prep <- resp[rownames(resp) %in% rownames(expr), ]

# select features with maximum variance
D <- ncol(resp.prep)
psel <- 1000
# o <- order(-colSums(Reduce("rbind", feat.prep$inpathway)))
o <- order(-apply(expr.prep, 2, sd))
idsel <- rep(list(o[1:psel]), D)
expr.sel <- lapply(idsel, function(s) {expr.prep[, s]})
inpathway <- lapply(1:D, function(d) {feat.prep$inpathway[[d]][idsel[[d]]]})

# create data objects used in fitting
x <- lapply(expr.sel, function(s) {scale(s)})
y <- scale(resp.prep)
p <- sapply(x, ncol)
n <- nrow(y)

### model fitting
# setting seed
set.seed(2019)

# estimation settings
control <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, epsilon.vb=1e-3,
                maxit.eb=200, maxit.vb=1, maxit.post=100, maxit.block=0)
methods <- c("NIG1", "NIG2", "NIG3", "NIG4", "NIG5", "lasso", "ridge", "bSEM")

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

# blockwise updating of hyperparameters
control <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, epsilon.vb=1e-3,
                maxit.eb=200, maxit.vb=1, maxit.post=100, maxit.block=10)
fit5.semnig <- semnig(x=x, y=y, C=C, Z=Z, unpenalized=NULL,
                      standardize=FALSE, intercept=FALSE, fixed.eb=FALSE,
                      full.post=TRUE, init=NULL, control=control)

# bSEM model
fit1.bSEM <- bSEM(x, split(y, 1:ncol(y)), lapply(inpathway, "+", 1),
                  control=list(maxit=200, trace=TRUE, epsilon=1e-3))

# penalized regression models
fit1.lasso <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[, d], intercept=FALSE, standardize=FALSE)})
fit1.ridge <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[, d], alpha=0, intercept=FALSE, standardize=FALSE)})

# saving fitted model objects
save(fit1.semnig, fit2.semnig, fit3.semnig, fit4.semnig, fit5.semnig, 
     fit1.lasso, fit1.ridge,
     fit1.bSEM, file="results/analysis_gdsc_fit3.Rdata")

# EB estimates
tab <- rbind(c(fit1.semnig$eb$alphaf, NA, fit1.semnig$eb$alphad, rep(NA, 4), 
               fit1.semnig$eb$lambdaf, fit1.semnig$eb$lambdad),
             c(fit2.semnig$eb$alphaf, fit2.semnig$eb$alphad, 
               fit2.semnig$eb$lambdaf, fit2.semnig$eb$lambdad),
             c(fit3.semnig$eb$alphaf, rep(NA, 5), fit3.semnig$eb$lambdaf, NA),
             c(rep(NA, 2), fit4.semnig$eb$alphad, NA, fit2.semnig$eb$lambdad),
             c(fit5.semnig$eb$alphaf, fit5.semnig$eb$alphad, 
               fit5.semnig$eb$lambdaf, fit5.semnig$eb$lambdad),
             c(fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 1]/
                 fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 1],
               fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 2]/
                 fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 2] -
                 fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 1]/
                 fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 1], rep(NA, 7)))
colnames(tab) <- c(colnames(C[[1]]), colnames(Z), "lambdaf", "lambdad")
rownames(tab) <- methods[c(1:5, 7)]
write.table(tab, file="results/analysis_gdsc_fit3.txt")

### cross-validation of performance measures
# settings seed
set.seed(2020)

# estimation settings
lasso <- stan_model("code/lasso.stan", auto_write=TRUE)
ridge <- stan_model("code/ridge.stan", auto_write=TRUE)
methods <- c(paste0("NIG", 1:9), paste0("lasso", 1:4), paste0("ridge", 1:5), 
             "bSEM")
foldid <- sample(c(rep(1:nfolds, each=n %/% nfolds), 
                   rep(1:(n %% nfolds), (n %% nfolds)!=0)))

# setup cluster
if(Sys.info()["sysname"]=="Windows") {
  cl <- makeCluster(ncores, type="PSOCK")
} else {
  cl <- makeCluster(ncores, type="FORK")
}
if(ncores > 1) {
  registerDoParallel(cl)
} else {
  registerDoSEQ()
}
res <- foreach(r=1:nfolds, .packages=packages, .errorhandling="pass", 
               .verbose=TRUE) %dopar% {
  cat("\r", "fold", r)
  set.seed(2020 + r)
  
  # splitting data
  ntrain <- sum(foldid!=r)
  xtrain <- lapply(x, function(s) {
   scale(as.matrix(s[foldid!=r, ], nrow=ntrain))})
  ytrain <- scale(as.matrix(y[foldid!=r, ], nrow=ntrain))
  xtest <- lapply(x, function(s) {
   scale(as.matrix(s[foldid==r, ], nrow=n - ntrain))})
  ytest <- scale(as.matrix(y[foldid==r, ], nrow=n - ntrain))
  
  # model without external covariates
  control <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-3, epsilon.vb=1e-3, 
                 maxit.eb=200, maxit.vb=1, maxit.post=100, maxit.block=0)
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
  
  # blockwise updating of hyperparameters
  control <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, epsilon.vb=1e-3, 
                 maxit.eb=1, maxit.vb=1, maxit.post=100, maxit.block=10)
  cv5.semnig <- semnig(x=xtrain, y=ytrain, C=C, Z=Z, unpenalized=NULL,
                      standardize=FALSE, intercept=FALSE, fixed.eb=FALSE,
                      full.post=TRUE, init=NULL, control=control)
  
  # optimisation of MAP
  cv6.semnig <- lapply(1:D, function(d) {
   optimizing(stanmodels$nig, 
              data=list(p=p[d], n=ntrain, x=xtrain[[d]], y=ytrain[, d], 
                        phi=cv2.semnig$eb$mpriorf[[d]], 
                        lambdaf=cv2.semnig$eb$lambdaf, 
                        chi=cv2.semnig$eb$mpriord[d], 
                        lambdad=cv2.semnig$eb$lambdad))})
  
  phi <- seq(0.01, 10, length.out=10)
  chi <- seq(0.01, 10, length.out=10)
  lambdaf <- cv2.semnig$eb$lambdaf
  lambdad <- cv2.semnig$eb$lambdad
  cv7.semnig <- cv.semnig(x=xtrain, y=ytrain, nfolds=5, foldid=NULL, 
                         seed=NULL, phi=phi, chi=chi, lambdaf=lambdaf, 
                         lambdad=lambdad, type.measure="mse", 
                         control=list(trace=FALSE))
  cv8.semnig <- lapply(1:D, function(d) {
   sampling(stanmodels$nig, 
            data=list(p=p[d], n=ntrain, x=xtrain[[d]], y=ytrain[, d], 
                      phi=cv2.semnig$eb$mpriorf[[d]], 
                      lambdaf=cv2.semnig$eb$lambdaf, 
                      chi=cv2.semnig$eb$mpriord[d], 
                      lambdad=cv2.semnig$eb$lambdad),
            chains=1, iter=1000, warmup=500, verbose=FALSE, refresh=0)})
  
  # bSEM model
  cv1.bSEM <- bSEM(xtrain, split(ytrain, 1:ncol(ytrain)),
                  lapply(inpathway, "+", 1),
                  control=list(maxit=200, trace=FALSE, epsilon=1e-3))
  
  # penalized regression models
  cv1.lasso <- lapply(1:D, function(d) {
   cv.glmnet(xtrain[[d]], ytrain[, d], intercept=FALSE, standardize=FALSE)})
  cv2.lasso <- lapply(1:D, function(d) {
   sampling(lasso, data=list(p=p[d], n=ntrain, x=xtrain[[d]], y=ytrain[, d], 
                             a0=0.001, b0=0.001),
            chains=1, iter=2000, warmup=1000, verbose=FALSE, refresh=0)})
  cv3.lasso <- cv.glmnet(xtrain[[1]], ytrain, family="mgaussian", 
                        intercept=FALSE, standardize=FALSE)
  
  cv1.ridge <- lapply(1:D, function(d) {
   cv.glmnet(xtrain[[d]], ytrain[, d], alpha=0, intercept=FALSE,
             standardize=FALSE)})
  cv2.ridge <- lapply(1:D, function(d) {
   sampling(ridge, data=list(p=p[d], n=ntrain, x=xtrain[[d]], y=ytrain[, d], 
                             a0=0.001, b0=0.001),
            chains=1, iter=1000, warmup=500, verbose=FALSE, refresh=0)})
  cv3.ridge <- cv.glmnet(xtrain[[1]], ytrain, family="mgaussian", alpha=0, 
                        intercept=FALSE, standardize=FALSE)
  cv4.ridge <- lapply(1:D, function(d) {
   BayesRidge(ytrain[, d], xtrain[[d]], plotML=FALSE, a=0.001, b=0.001, 
              SVD=NULL)})
  
  # extracting point estimates
  best <- list(Reduce("cbind", cv1.semnig$vb$mpost$beta),
              Reduce("cbind", cv2.semnig$vb$mpost$beta),
              Reduce("cbind", cv3.semnig$vb$mpost$beta),
              Reduce("cbind", cv4.semnig$vb$mpost$beta),
              Reduce("cbind", cv5.semnig$vb$mpost$beta),
              sapply(cv6.semnig, function(s) {
                s$par[names(s$par) %in% paste0("beta[", 1:p[1], "]")]}),
              sapply(cv7.semnig, function(s) {
                s$fit$par[names(s$fit$par) %in% 
                            paste0("beta[", 1:p[1], "]")]}),
              sapply(cv8.semnig, get_posterior_mean, 
                     pars=paste0("beta[", 1:p[1], "]")),
              sapply(cv8.semnig, function(s) {
                sapply(extract(s, pars=paste0("beta[", 1:p[1], "]")), 
                       function(b) {d <- density(b); d$x[which.max(d$y)]})}),
              sapply(cv1.lasso, function(s) {
                as.numeric(coef(s, s="lambda.min"))[-1]}),
              sapply(cv2.lasso, get_posterior_mean, 
                     pars=paste0("beta[", 1:p[1], "]")),
              sapply(cv2.lasso, function(s) {
                sapply(extract(s, pars=paste0("beta[", 1:p[1], "]")), 
                       function(b) {d <- density(b); d$x[which.max(d$y)]})}),
              as.matrix(Reduce("cbind", 
                               coef(cv3.lasso, s="lambda.min")))[-1, ],
              sapply(cv1.ridge, function(s) {
                as.numeric(coef(s, s="lambda.min"))[-1]}),
              sapply(cv2.ridge, get_posterior_mean, 
                     pars=paste0("beta[", 1:p[1], "]")),
              sapply(cv2.ridge, function(s) {
                sapply(extract(s, pars=paste0("beta[", 1:p[1], "]")), 
                       function(b) {d <- density(b); d$x[which.max(d$y)]})}),
              as.matrix(Reduce("cbind", 
                               coef(cv3.ridge, s="lambda.min")))[-1, ],
              sapply(cv4.ridge, "[[", "postmeanbeta"),
              sapply(cv1.bSEM$vb$beta, "[", 1:p[1]))
  
  # calculating average prediction mean squared error
  pmse <- cbind(sapply(best, function(b) {
   colMeans((ytest - xtest[[1]] %*% b)^2)}), 
   sapply(1:D, function(d) {mean((ytest[, d] - colMeans(ytest)[d])^2)}))
  
  # calculate ELBO for semnig models on training and test data
  elbot <- c(mean(cv1.semnig$seq.elbo[nrow(cv1.semnig$seq.elbo), ]),
            mean(cv2.semnig$seq.elbo[nrow(cv2.semnig$seq.elbo), ]),
            mean(cv3.semnig$seq.elbo[nrow(cv3.semnig$seq.elbo), ]),
            mean(cv4.semnig$seq.elbo[nrow(cv4.semnig$seq.elbo), ]),
            mean(cv5.semnig$seq.elbo[nrow(cv5.semnig$seq.elbo), ]),
            mean(cv1.bSEM$mll[, ncol(cv1.bSEM$mll)]))
  elbo <- c(mean(new.elbo(cv1.semnig, xtest, ytest)),
           mean(new.elbo(cv2.semnig, xtest, ytest)),
           mean(new.elbo(cv3.semnig, xtest, ytest)),
           mean(new.elbo(cv4.semnig, xtest, ytest)),
           mean(new.elbo(cv5.semnig, xtest, ytest)))
  
  # log pseudo marginal log likelihood
  lpml <- c(sum(logcpo(xtest, ytest, ntrain, cv1.semnig), na.rm=TRUE),
           sum(logcpo(xtest, ytest, ntrain, cv2.semnig), na.rm=TRUE),
           sum(logcpo(xtest, ytest, ntrain, cv3.semnig), na.rm=TRUE),
           sum(logcpo(xtest, ytest, ntrain, cv4.semnig), na.rm=TRUE),
           sum(logcpo(xtest, ytest, ntrain, cv5.semnig), na.rm=TRUE))
  
  # determine the ranks of the model parameter point estimates
  brank <- c(sapply(best, rank, ties.method="average"))
  
  names(brank) <- paste0(
   "brank.", paste0(rep(methods, each=sum(p)), ".",
                    rep(paste0(rep(1:D, times=p), ".",
                               unlist(sapply(p, function(s) {
                                 return(1:s)}))), length(methods))))
  
  list(pmse=pmse, elbo=elbo, elbot=elbot, lpml=lpml, brank=brank)
}
stopCluster(cl=cl)
str(res)
errorid <- sapply(res, function(s) {is.numeric(s[[1]])})
res2 <- res[errorid]

# prepare and save results table
pmse <- Reduce("rbind", sapply(res2, function(s) {as.list(s["pmse"])}))
elbo <- t(sapply(res2, "[[", "elbo"))
elbot <- t(sapply(res2, "[[", "elbot"))
lpml <- t(sapply(res2, "[[", "lpml"))
brank <- t(sapply(res2, "[[", "brank"))

# Euclidian distance between rank vectors
brankdist <- sapply(methods, function(s) {
  dist(brank[, substr(colnames(brank), 7, nchar(s) + 6)==s])})

# combine tables and save
res <- rbind(pmse, cbind(elbo, matrix(NA, nrow=nfolds, ncol=15)), 
             cbind(elbot, matrix(NA, nrow=nfolds, ncol=14)), 
             cbind(lpml, matrix(NA, nrow=nfolds, ncol=15)), 
             cbind(brankdist, NA))
colnames(res) <- c(methods, "null")
rownames(res) <- c(paste0(rep("pmse", nfolds*D), ".drug", 
                          rep(1:D, times=nfolds)), 
                   rep("elbo", nfolds), rep("elbot", nfolds),
                   rep("lpml", nfolds), rep("brankdist", nrow(brankdist)))
write.table(res, file="results/analysis_gdsc_res3.txt")



################################## analysis 4 ##################################
### data preparation
# only retain the cell lines that are available in both response and methylation
expr.prep <- expr[(rownames(expr) %in% rownames(resp)) & 
                    (rownames(expr) %in% rownames(mut)), ]
mut.prep <- mut[(rownames(mut) %in% rownames(resp)) & 
                    (rownames(mut) %in% rownames(expr)), ]
resp.prep <- resp[rownames(resp) %in% rownames(expr) &
                    rownames(resp) %in% rownames(mut), ]

# selecting expression values and mutations with largest variance
psel <- 100
o <- order(-apply(expr.prep, 2, sd))
idsel <- o[1:psel]
expr.sel <- expr.prep[, idsel]
o <- order(-apply(mut.prep, 2, sd))
idsel <- o[1:psel]
mut.sel <- mut.prep[, idsel]
mutation <- replicate(ncol(resp.prep), rep(c(0, 1), each=psel), simplify=FALSE)

# transform methylation values
x <- lapply(1:ncol(resp.prep), function(s) {
  cbind(scale(log(expr.sel)), scale(mut.sel))})
y <- scale(resp.prep)
D <- ncol(y)
p <- sapply(x, ncol)
n <- nrow(y)

### model fitting
# setting seed
set.seed(2020)

# estimation settings
control <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, epsilon.vb=1e-3,
                maxit.eb=200, maxit.vb=1, maxit.post=100, maxit.block=0)
methods <- c("NIG1", "NIG2", "NIG3", "NIG4", "NIG5", "lasso", "ridge", "bSEM")

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

# blockwise updating of hyperparameters
control <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, epsilon.vb=1e-3,
                maxit.eb=200, maxit.vb=1, maxit.post=100, maxit.block=10)
fit5.semnig <- semnig(x=x, y=y, C=C, Z=Z, unpenalized=NULL,
                      standardize=FALSE, intercept=FALSE, fixed.eb=FALSE,
                      full.post=TRUE, init=NULL, control=control)

# bSEM model
fit1.bSEM <- bSEM(x, split(y, 1:ncol(y)), lapply(inpathway, "+", 1),
                  control=list(maxit=200, trace=TRUE, epsilon=1e-3))

# penalized regression models
fit1.lasso <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[, d], intercept=FALSE, standardize=FALSE)})
fit1.ridge <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[, d], alpha=0, intercept=FALSE, standardize=FALSE)})

# saving fitted model objects
save(fit1.semnig, fit2.semnig, fit3.semnig, fit4.semnig, fit5.semnig, 
     fit1.lasso, fit1.ridge,
     fit1.bSEM, file="results/analysis_gdsc_fit4.Rdata")

# EB estimates
tab <- rbind(c(fit1.semnig$eb$alphaf, NA, fit1.semnig$eb$alphad, rep(NA, 4), 
               fit1.semnig$eb$lambdaf, fit1.semnig$eb$lambdad),
             c(fit2.semnig$eb$alphaf, fit2.semnig$eb$alphad, 
               fit2.semnig$eb$lambdaf, fit2.semnig$eb$lambdad),
             c(fit3.semnig$eb$alphaf, rep(NA, 5), fit3.semnig$eb$lambdaf, NA),
             c(rep(NA, 2), fit4.semnig$eb$alphad, NA, fit2.semnig$eb$lambdad),
             c(fit5.semnig$eb$alphaf, fit5.semnig$eb$alphad, 
               fit5.semnig$eb$lambdaf, fit5.semnig$eb$lambdad),
             c(fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 1]/
                 fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 1],
               fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 2]/
                 fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 2] -
                 fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 1]/
                 fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 1], rep(NA, 7)))
colnames(tab) <- c(colnames(C[[1]]), colnames(Z), "lambdaf", "lambdad")
rownames(tab) <- methods[c(1:5, 7)]
write.table(tab, file="results/analysis_gdsc_fit4.txt")

### cross-validation of performance measures
# settings seed
set.seed(2020)

# estimation settings
lasso <- stan_model("code/lasso.stan", auto_write=TRUE)
ridge <- stan_model("code/ridge.stan", auto_write=TRUE)
methods <- c(paste0("NIG", c(1:5, 7:8)),
             # paste0("NIG", 1:8), 
             paste0("lasso", 1:4), paste0("ridge", 1:5), 
             "bSEM")
foldid <- sample(c(rep(1:nfolds, each=n %/% nfolds), 
                   rep(1:(n %% nfolds), (n %% nfolds)!=0)))

# setup cluster
if(Sys.info()["sysname"]=="Windows") {
  cl <- makeCluster(ncores, type="PSOCK")
} else {
  cl <- makeCluster(ncores, type="FORK")
}
if(ncores > 1) {
  registerDoParallel(cl)
} else {
  registerDoSEQ()
}
res <- foreach(r=1:nfolds, .packages=packages, .errorhandling="pass", 
               .verbose=TRUE) %dopar% {
  cat("\r", "fold", r)
  set.seed(2020 + r)
 
  # splitting data
  ntrain <- sum(foldid!=r)
  xtrain <- lapply(x, function(s) {
   scale(as.matrix(s[foldid!=r, ], nrow=ntrain))})
  ytrain <- scale(as.matrix(y[foldid!=r, ], nrow=ntrain))
  xtest <- lapply(x, function(s) {
   scale(as.matrix(s[foldid==r, ], nrow=n - ntrain))})
  ytest <- scale(as.matrix(y[foldid==r, ], nrow=n - ntrain))
  
  # model without external covariates
  control <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-3, epsilon.vb=1e-3, 
                 maxit.eb=200, maxit.vb=1, maxit.post=100, maxit.block=0)
  Z <- matrix(1, nrow=D)
  colnames(Z) <- c("intercept")
  C <- lapply(mutation, function(s) {
   s <- matrix(1, nrow=length(s)); colnames(s) <- c("intercept"); return(s)})
  cv1.semnig <- semnig(x=xtrain, y=ytrain, C=C, Z=Z, unpenalized=NULL, 
                      standardize=FALSE, intercept=FALSE, fixed.eb="none", 
                      full.post=TRUE, init=NULL, control=control)
  C <- lapply(mutation, function(s) {
   s <- cbind(1, s); colnames(s) <- c("intercept", "mutation"); return(s)})
  cv2.semnig <- semnig(x=xtrain, y=ytrain, C=C, Z=Z, unpenalized=NULL, 
                      standardize=FALSE, intercept=FALSE, fixed.eb="none", 
                      full.post=TRUE, init=NULL, control=control)
  
  # models with external covariates
  Z <- model.matrix(~ replace(drug.prep$stage, is.na(drug.prep$stage),
                             "experimental") + replace(drug.prep$action, is.na(
                               drug.prep$action), "unknown"))
  colnames(Z) <- c("intercept", "experimental", "in clinical development",
                  "targeted", "unknown")
  cv3.semnig <- semnig(x=xtrain, y=ytrain, C=C, Z=Z, unpenalized=NULL, 
                      standardize=FALSE, intercept=FALSE, fixed.eb="none", 
                      full.post=TRUE, init=NULL, control=control)
  
  # blockwise updating of hyperparameters
  control <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-3, epsilon.vb=1e-3, 
                 maxit.eb=1, maxit.vb=1, maxit.post=100, maxit.block=10)
  cv4.semnig <- semnig(x=xtrain, y=ytrain, C=C, Z=Z, unpenalized=NULL,
                      standardize=FALSE, intercept=FALSE, fixed.eb=FALSE,
                      full.post=TRUE, init=NULL, control=control)
  
  # optimisation of MAP
  cv5.semnig <- lapply(1:D, function(d) {
   optimizing(stanmodels$nig, 
              data=list(p=p[d], n=ntrain, x=xtrain[[d]], y=ytrain[, d], 
                        phi=cv3.semnig$eb$mpriorf[[d]], 
                        lambdaf=cv3.semnig$eb$lambdaf, 
                        chi=cv3.semnig$eb$mpriord[d], 
                        lambdad=cv3.semnig$eb$lambdad))})
  
  # phi <- seq(0.01, 10, length.out=10)
  # chi <- seq(0.01, 10, length.out=10)
  # lambdaf <- cv3.semnig$eb$lambdaf
  # lambdad <- cv3.semnig$eb$lambdad
  # cv6.semnig <- cv.semnig(x=xtrain, y=ytrain, nfolds=5, foldid=NULL, 
  #                        seed=NULL, phi=phi, chi=chi, lambdaf=lambdaf, 
  #                        lambdad=lambdad, type.measure="mse", 
  #                        control=list(trace=FALSE))
  cv7.semnig <- lapply(1:D, function(d) {
   sampling(stanmodels$nig, 
            data=list(p=p[d], n=ntrain, x=xtrain[[d]], y=ytrain[, d], 
                      phi=cv3.semnig$eb$mpriorf[[d]], 
                      lambdaf=cv3.semnig$eb$lambdaf, 
                      chi=cv3.semnig$eb$mpriord[d], 
                      lambdad=cv3.semnig$eb$lambdad),
            chains=1, iter=1000, warmup=500, verbose=FALSE, refresh=0)})
  
  # bSEM model
  cv1.bSEM <- bSEM(xtrain, split(ytrain, 1:ncol(ytrain)),
                  lapply(mutation, "+", 1),
                  control=list(maxit=200, trace=FALSE, epsilon=1e-3))
  
  # penalized regression models
  cv1.lasso <- lapply(1:D, function(d) {
   cv.glmnet(xtrain[[d]], ytrain[, d], intercept=FALSE, standardize=FALSE)})
  cv2.lasso <- lapply(1:D, function(d) {
   sampling(lasso, data=list(p=p[d], n=ntrain, x=xtrain[[d]], y=ytrain[, d], 
                             a0=0.001, b0=0.001),
            chains=1, iter=2000, warmup=1000, verbose=FALSE, refresh=0)})
  cv3.lasso <- cv.glmnet(xtrain[[1]], ytrain, family="mgaussian", 
                        intercept=FALSE, standardize=FALSE)
  
  cv1.ridge <- lapply(1:D, function(d) {
   cv.glmnet(xtrain[[d]], ytrain[, d], alpha=0, intercept=FALSE,
             standardize=FALSE)})
  cv2.ridge <- lapply(1:D, function(d) {
   sampling(ridge, data=list(p=p[d], n=ntrain, x=xtrain[[d]], y=ytrain[, d], 
                             a0=0.001, b0=0.001),
            chains=1, iter=1000, warmup=500, verbose=FALSE, refresh=0)})
  cv3.ridge <- cv.glmnet(xtrain[[1]], ytrain, family="mgaussian", alpha=0, 
                        intercept=FALSE, standardize=FALSE)
  cv4.ridge <- lapply(1:D, function(d) {
   BayesRidge(ytrain[, d], xtrain[[d]], plotML=FALSE, a=0.001, b=0.001, 
              SVD=NULL)})
  
  # extracting point estimates
  best <- list(Reduce("cbind", cv1.semnig$vb$mpost$beta),
              Reduce("cbind", cv2.semnig$vb$mpost$beta),
              Reduce("cbind", cv3.semnig$vb$mpost$beta),
              Reduce("cbind", cv4.semnig$vb$mpost$beta),
              sapply(cv5.semnig, function(s) {
                s$par[names(s$par) %in% paste0("beta[", 1:p[1], "]")]}),
              # sapply(cv6.semnig, function(s) {
              #   s$fit$par[names(s$fit$par) %in% 
              #               paste0("beta[", 1:p[1], "]")]}),
              sapply(cv7.semnig, get_posterior_mean, 
                     pars=paste0("beta[", 1:p[1], "]")),
              sapply(cv7.semnig, function(s) {
                sapply(extract(s, pars=paste0("beta[", 1:p[1], "]")), 
                       function(b) {d <- density(b); d$x[which.max(d$y)]})}),
              sapply(cv1.lasso, function(s) {
                as.numeric(coef(s, s="lambda.min"))[-1]}),
              sapply(cv2.lasso, get_posterior_mean, 
                     pars=paste0("beta[", 1:p[1], "]")),
              sapply(cv2.lasso, function(s) {
                sapply(extract(s, pars=paste0("beta[", 1:p[1], "]")), 
                       function(b) {d <- density(b); d$x[which.max(d$y)]})}),
              as.matrix(Reduce("cbind", 
                               coef(cv3.lasso, s="lambda.min")))[-1, ],
              sapply(cv1.ridge, function(s) {
                as.numeric(coef(s, s="lambda.min"))[-1]}),
              sapply(cv2.ridge, get_posterior_mean, 
                     pars=paste0("beta[", 1:p[1], "]")),
              sapply(cv2.ridge, function(s) {
                sapply(extract(s, pars=paste0("beta[", 1:p[1], "]")), 
                       function(b) {d <- density(b); d$x[which.max(d$y)]})}),
              as.matrix(Reduce("cbind", 
                               coef(cv3.ridge, s="lambda.min")))[-1, ],
              sapply(cv4.ridge, "[[", "postmeanbeta"),
              sapply(cv1.bSEM$vb$beta, "[", 1:p[1]))
  
  # calculating average prediction mean squared error
  pmse <- cbind(sapply(best, function(b) {
   colMeans((ytest - xtest[[1]] %*% b)^2)}), 
   sapply(1:D, function(d) {mean((ytest[, d] - colMeans(ytest)[d])^2)}))
  
  # calculate ELBO for semnig models on training and test data
  elbot <- c(mean(cv1.semnig$seq.elbo[nrow(cv1.semnig$seq.elbo), ]),
            mean(cv2.semnig$seq.elbo[nrow(cv2.semnig$seq.elbo), ]),
            mean(cv3.semnig$seq.elbo[nrow(cv3.semnig$seq.elbo), ]),
            mean(cv4.semnig$seq.elbo[nrow(cv4.semnig$seq.elbo), ]),
            mean(cv1.bSEM$mll[, ncol(cv1.bSEM$mll)]))
  elbo <- c(mean(new.elbo(cv1.semnig, xtest, ytest)),
           mean(new.elbo(cv2.semnig, xtest, ytest)),
           mean(new.elbo(cv3.semnig, xtest, ytest)),
           mean(new.elbo(cv4.semnig, xtest, ytest)))
  
  # log pseudo marginal log likelihood
  # lpml <- c(sum(logcpo(xtest, ytest, ntrain, cv1.semnig), na.rm=TRUE),
  #          sum(logcpo(xtest, ytest, ntrain, cv2.semnig), na.rm=TRUE),
  #          sum(logcpo(xtest, ytest, ntrain, cv3.semnig), na.rm=TRUE),
  #          sum(logcpo(xtest, ytest, ntrain, cv4.semnig), na.rm=TRUE))
  
  # determine the ranks of the model parameter point estimates
  brank <- c(sapply(best, rank, ties.method="average"))
  
  names(brank) <- paste0(
   "brank.", paste0(rep(methods, each=sum(p)), ".",
                    rep(paste0(rep(1:D, times=p), ".",
                               unlist(sapply(p, function(s) {
                                 return(1:s)}))), length(methods))))
  
  list(pmse=pmse, elbo=elbo, elbot=elbot, 
       # lpml=lpml, 
       brank=brank)
}
stopCluster(cl=cl)
str(res)
errorid <- sapply(res, function(s) {is.numeric(s[[1]])})
res2 <- res[errorid]

# prepare and save results table
pmse <- Reduce("rbind", sapply(res2, function(s) {as.list(s["pmse"])}))
elbo <- t(sapply(res2, "[[", "elbo"))
elbot <- t(sapply(res2, "[[", "elbot"))
# lpml <- t(sapply(res2, "[[", "lpml"))
brank <- t(sapply(res2, "[[", "brank"))

# Euclidian distance between rank vectors
brankdist <- sapply(methods, function(s) {
  dist(brank[, substr(colnames(brank), 7, nchar(s) + 6)==s])})

# combine tables and save
res <- rbind(pmse, cbind(elbo, matrix(NA, nrow=nfolds, ncol=14)), 
             cbind(elbot, matrix(NA, nrow=nfolds, ncol=13)), 
             # cbind(lpml, matrix(NA, nrow=nfolds, ncol=14)), 
             cbind(brankdist, NA))
colnames(res) <- c(methods, "null")
rownames(res) <- c(paste0(rep("pmse", nfolds*D), ".drug", 
                          rep(1:D, times=nfolds)), 
                   rep("elbo", nfolds), rep("elbot", nfolds),
                   # rep("lpml", nfolds), 
                   rep("brankdist", nrow(brankdist)))
write.table(res, file="results/analysis_gdsc_res4.txt")



################################## analysis 5 ##################################
### data preparation
# only retain the cell lines that are available in both response and expression
expr.prep <- expr[rownames(expr) %in% rownames(resp), ]
resp.prep <- resp[rownames(resp) %in% rownames(expr), ]

# select features with largest variance
psel <- 100
o <- order(-apply(expr.prep, 2, sd))
idsel <- o[1:psel]
expr.sel <- expr.prep[, idsel]

# select half of the drugs
set.seed(2020)
Dext <- floor(ncol(resp.prep)/2)
D <- ncol(resp.prep) - Dext
dsel <- sample(1:(Dext + D), Dext)

# create data objects used in fitting
xext <- lapply(1:Dext, function(d) {scale(expr.sel)})
x <- lapply(1:D, function(d) {scale(expr.sel)})
yext <- scale(resp.prep[, dsel])
y <- scale(resp.prep[, -dsel])
p <- sapply(x, ncol)
n <- nrow(y)

# estimation settings
control <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, epsilon.vb=1e-3,
                maxit.eb=200, maxit.vb=1, maxit.post=100, maxit.block=0)
methods <- c("NIG1", "NIG2", "NIG3", "NIG4", "NIG5", "lasso", "ridge")

# fit external data
C <- lapply(1:D, function(d) {
  s <- matrix(1, nrow=p[d]); colnames(s) <- c("intercept"); return(s)})
Z <- matrix(1, nrow=D)
colnames(Z) <- c("intercept")
pvalues <- sapply(1:Dext, function(d) {
  apply(xext[[d]], 2, function(s) {cor.test(s, yext[, d])$p.value})})
C1 <- lapply(1:D, function(d) {
  s <- cbind(1, -rowSums(log(pvalues))/Dext); 
  colnames(s) <- c("intercept", "extern"); return(s)})
C2 <- lapply(C1, function(s) {
  s <- cbind(1, -pchisq(2*Dext*s[, 2], 2*Dext, lower.tail=FALSE, log.p=TRUE)); 
  colnames(s) <- c("intercept", "extern"); return(s)})
C3 <- lapply(1:D, function(d) {
  s <- cbind(1, log(rowSums(1/(pvalues*Dext)))); 
  colnames(s) <- c("intercept", "extern"); return(s)})

### model fitting
# setting seed
set.seed(2020)

# model without external covariates
fit1.semnig <- semnig(x=x, y=y, C=C, Z=Z, unpenalized=NULL,
                      standardize=FALSE, intercept=FALSE, fixed.eb="none",
                      full.post=TRUE, init=NULL, control=control)

# models with external covariates
fit2.semnig <- semnig(x=x, y=y, C=C1, Z=Z, unpenalized=NULL,
                      standardize=FALSE, intercept=FALSE, fixed.eb=FALSE,
                      full.post=TRUE, init=NULL, control=control)
fit3.semnig <- semnig(x=x, y=y, C=C2, Z=Z, unpenalized=NULL,
                      standardize=FALSE, intercept=FALSE, fixed.eb=FALSE,
                      full.post=TRUE, init=NULL, control=control)
fit4.semnig <- semnig(x=x, y=y, C=C3, Z=Z, unpenalized=NULL,
                      standardize=FALSE, intercept=FALSE, fixed.eb=FALSE,
                      full.post=TRUE, init=NULL, control=control)
phi <- seq(0.01, 10, length.out=10)
chi <- seq(0.01, 10, length.out=10)
lambdaf <- fit2.semnig$eb$lambdaf
lambdad <- fit2.semnig$eb$lambdad
fit5.semnig <- cv.semnig(x=x, y=y, nfolds=5, foldid=NULL, 
                         seed=NULL, phi=phi, chi=chi, lambdaf=lambdaf, 
                         lambdad=lambdad, type.measure="mse", 
                         control=list(trace=TRUE))

# penalized regression models
fit1.lasso <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[, d], intercept=FALSE, standardize=FALSE)})
fit1.ridge <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[, d], alpha=0, intercept=FALSE, standardize=FALSE)})

save(fit1.semnig, fit2.semnig, fit3.semnig, fit4.semnig, fit5.semnig, 
     fit1.lasso, fit1.ridge, file="results/analysis_gdsc_fit5.Rdata")


# EB estimates
par.min <- sapply(fit5.semnig, "[[", "par.min")
tab <- rbind(c(fit1.semnig$eb$alphaf, NA, fit1.semnig$eb$alphad,
               fit1.semnig$eb$lambdaf, fit1.semnig$eb$lambdad),
             c(fit2.semnig$eb$alphaf, fit2.semnig$eb$alphad, 
               fit2.semnig$eb$lambdaf, fit2.semnig$eb$lambdad),
             c(fit3.semnig$eb$alphaf, fit3.semnig$eb$alphad, 
               fit3.semnig$eb$lambdaf, fit3.semnig$eb$lambdad),
             c(fit4.semnig$eb$alphaf, fit4.semnig$eb$alphad, 
               fit4.semnig$eb$lambdaf, fit4.semnig$eb$lambdad),
             c(mean(1/par.min["phi", ]), NA, mean(1/par.min["chi", ]),
               mean(par.min["lambdaf", ]), mean(par.min["lambdad", ])))
colnames(tab) <- c(colnames(C1[[1]]), colnames(Z), "lambdaf", "lambdad")
rownames(tab) <- methods[c(1:5)]
write.table(tab, file="results/analysis_gdsc_fit5.txt")

### cross-validation of performance measures
# settings seed
set.seed(2020)

# estimation settings
foldid <- sample(c(rep(1:nfolds, each=n %/% nfolds), 
                   rep(1:(n %% nfolds), (n %% nfolds)!=0)))
control$trace <- FALSE

# setup cluster
cl <- makeCluster(ncores) 
if(ncores > 1) {
  registerDoParallel(cl)
} else {
  registerDoSEQ()
}
res <- foreach(r=1:nfolds, .packages=packages, .errorhandling="pass",
               .verbose=TRUE) %dopar% {
  cat("\r", "fold", r)
  
  # splitting data
  ntrain <- sum(foldid!=r)
  xtrain <- lapply(x, function(s) {
    scale(as.matrix(s[foldid!=r, ], nrow=ntrain))})
  ytrain <- scale(as.matrix(y[foldid!=r, ], nrow=ntrain))
  xtest <- lapply(x, function(s) {
    scale(as.matrix(s[foldid==r, ], nrow=sum(foldid==r)))})
  ytest <- scale(as.matrix(y[foldid==r, ], nrow=sum(foldid==r)))
  
  # model without external covariates
  cv1.semnig <- semnig(x=xtrain, y=ytrain, C=C, Z=Z, unpenalized=NULL,
                       standardize=FALSE, intercept=FALSE, fixed.eb="none",
                       full.post=TRUE, init=NULL, control=control)
  
  # models with external covariates
  cv2.semnig <- semnig(x=xtrain, y=ytrain, C=C1, Z=Z, unpenalized=NULL,
                       standardize=FALSE, intercept=FALSE, fixed.eb=FALSE,
                       full.post=TRUE, init=NULL, control=control)
  cv3.semnig <- semnig(x=xtrain, y=ytrain, C=C2, Z=Z, unpenalized=NULL,
                       standardize=FALSE, intercept=FALSE, fixed.eb=FALSE,
                       full.post=TRUE, init=NULL, control=control)
  cv4.semnig <- semnig(x=xtrain, y=ytrain, C=C3, Z=Z, unpenalized=NULL,
                       standardize=FALSE, intercept=FALSE, fixed.eb=FALSE,
                       full.post=TRUE, init=NULL, control=control)
  phi <- seq(0.01, 10, length.out=10)
  chi <- seq(0.01, 10, length.out=10)
  lambdaf <- cv2.semnig$eb$lambdaf
  lambdad <- cv2.semnig$eb$lambdad
  cv5.semnig <- cv.semnig(x=xtrain, y=ytrain, nfolds=5, foldid=NULL, 
                          seed=NULL, phi=phi, chi=chi, lambdaf=lambdaf, 
                          lambdad=lambdad, type.measure="mse", 
                          control=list(trace=TRUE))
  
  # penalized regression models
  cv1.lasso <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[, d], intercept=FALSE, standardize=FALSE)})
  cv1.ridge <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[, d], alpha=0, intercept=FALSE, 
              standardize=FALSE)})
  
  # calculating average prediction mean squared error
  pmse <- c(mean(sapply(1:D, function(d) {
    mean((ytest[, d] - xtest[[d]] %*% cv1.semnig$vb$mpost$beta[[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((ytest[, d] - xtest[[d]] %*% cv2.semnig$vb$mpost$beta[[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((ytest[, d] - xtest[[d]] %*% cv3.semnig$vb$mpost$beta[[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((ytest[, d] - xtest[[d]] %*% cv4.semnig$vb$mpost$beta[[d]])^2)})),
    mean(sapply(1:D, function(d) {
      best <- cv5.semnig[[d]]$fit$par[
        names(cv5.semnig[[d]]$fit$par) %in% paste0("beta[", 1:p[d], "]")];
      mean((ytest[, d] - xtest[[d]] %*% best)^2)})),
    mean((ytest - sapply(1:D, function(d) {
      predict(cv1.lasso[[d]], xtest[[d]], s="lambda.min")}))^2),
    mean((ytest - sapply(1:D, function(d) {
      predict(cv1.ridge[[d]], xtest[[d]], s="lambda.min")}))^2),
    mean(apply(ytest, 1, "-", colMeans(ytrain))^2))
  
  # calculate ELBO for semnig models on training and test data
  elbot <- c(mean(cv1.semnig$seq.elbo[nrow(cv1.semnig$seq.elbo), ]),
             mean(cv2.semnig$seq.elbo[nrow(cv2.semnig$seq.elbo), ]),
             mean(cv3.semnig$seq.elbo[nrow(cv3.semnig$seq.elbo), ]),
             mean(cv4.semnig$seq.elbo[nrow(cv4.semnig$seq.elbo), ]))
  elbo <- c(mean(new.elbo(cv1.semnig, xtest, ytest)),
            mean(new.elbo(cv2.semnig, xtest, ytest)),
            mean(new.elbo(cv3.semnig, xtest, ytest)),
            mean(new.elbo(cv4.semnig, xtest, ytest)))
  lpml <- c(sum(logcpo(xtest, ytest, ntrain, cv1.semnig), na.rm=TRUE),
            sum(logcpo(xtest, ytest, ntrain, cv2.semnig), na.rm=TRUE),
            sum(logcpo(xtest, ytest, ntrain, cv3.semnig), na.rm=TRUE),
            sum(logcpo(xtest, ytest, ntrain, cv4.semnig), na.rm=TRUE))
  
  # determine the ranks of the model parameter point estimates
  brank <- c(unlist(lapply(cv1.semnig$vb$mpost$beta, function(s) {
    rank(abs(s), ties.method="average")})),
    unlist(lapply(cv2.semnig$vb$mpost$beta, function(s) {
      rank(abs(s), ties.method="average")})),
    unlist(lapply(cv3.semnig$vb$mpost$beta, function(s) {
      rank(abs(s), ties.method="average")})),
    unlist(lapply(cv4.semnig$vb$mpost$beta, function(s) {
      rank(abs(s), ties.method="average")})),
    unlist(lapply(1:D, function(d) {
      s <- cv5.semnig[[d]]$fit$par[
        names(cv5.semnig[[d]]$fit$par) %in% paste0("beta[", 1:p[d], "]")];
      rank(abs(s), ties.method="average")})),
    unlist(sapply(cv1.lasso, function(s) {
      rank(abs(as.numeric(coef(s, s="lambda.min"))[-1]), 
           ties.method="average")})),
    unlist(sapply(cv1.ridge, function(s) {
      rank(abs(as.numeric(coef(s, s="lambda.min"))[-1]), 
           ties.method="average")})))
  names(brank) <- paste0(
    "brank.", paste0(rep(methods, each=sum(p)), ".",
                     rep(paste0(rep(1:D, times=p), ".",
                                unlist(sapply(p, function(s) {
                                  return(1:s)}))), length(methods))))
  
  list(pmse=pmse, elbo=elbo, elbot=elbot, lpml=lpml, brank=brank)
}
stopCluster(cl=cl)

errorid <- sapply(res, function(s) {is.numeric(s[[1]])})
res2 <- res[errorid]

# prepare and save results table
pmse <- t(sapply(res2, "[[", "pmse"))
elbo <- t(sapply(res2, "[[", "elbo"))
elbot <- t(sapply(res2, "[[", "elbot"))
lpml <- t(sapply(res2, "[[", "lpml"))
brank <- t(sapply(res2, "[[", "brank"))

# Euclidian distance between rank vectors
brankdist <- sapply(methods, function(s) {
  dist(brank[, substr(colnames(brank), 7, nchar(s) + 6)==s])})

# combine tables and save
res <- rbind(pmse, cbind(elbo, NA, NA, NA, NA), cbind(elbot, NA, NA, NA, NA), 
             cbind(lpml, NA, NA, NA, NA), cbind(brankdist, NA))
colnames(res) <- c(methods, "null")
rownames(res) <- c(rep("pmse", nfolds), rep("elbo", nfolds), 
                   rep("elbot", nfolds), rep("lpml", nfolds), 
                   rep("brankdist", nrow(brankdist)))
write.table(res, file="results/analysis_gdsc_res5.txt")


################################## analysis 6 ##################################
### data preparation
# only retain the cell lines that are available in both response and expression
feat.prep <- read.table(file="results/analysis_ccle_res1.txt")
expr.prep <- expr[rownames(expr) %in% rownames(resp), ]
resp.prep <- resp[rownames(resp) %in% rownames(expr), ]

# select features with largest variance separately for in and out pathway genes
psel <- 100
o <- order(-apply(expr.prep, 2, sd))
idsel <- o[1:psel]
expr.sel <- expr.prep[, idsel]
feat.sel <- feat.prep[rownames(feat.prep) %in% colnames(expr.sel), ]

# create data objects used in fitting
D <- ncol(resp.prep)
x <- lapply(1:D, function(d) {scale(expr.sel, scale=FALSE)})
y <- scale(resp.prep, scale=FALSE)
p <- sapply(x, ncol)
n <- nrow(y)

### model fitting
# setting seed
set.seed(2020)

# estimation settings
control <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, epsilon.vb=1e-3,
                maxit.eb=200, maxit.vb=1, maxit.post=100, maxit.block=0)
methods <- c("NIG1", "NIG2", "NIG3", "NIG4", "lasso", "ridge")

# fit external data
C <- lapply(1:D, function(d) {
  s <- matrix(1, nrow=p[d]); colnames(s) <- c("intercept"); return(s)})
Z <- matrix(1, nrow=D)
colnames(Z) <- c("intercept")
C1 <- lapply(1:D, function(d) {
  s <- cbind(1, feat.sel$chifisher); 
  colnames(s) <- c("intercept", "extern"); return(s)})
C2 <- lapply(C1, function(s) {
  s <- cbind(1, feat.sel$pfisher); 
  colnames(s) <- c("intercept", "extern"); return(s)})
C3 <- lapply(1:D, function(d) {
  s <- cbind(1, feat.sel$harmonic); 
  colnames(s) <- c("intercept", "extern"); return(s)})

# model without external covariates
fit1.semnig <- semnig(x=x, y=y, C=C, Z=Z, unpenalized=NULL,
                      standardize=FALSE, intercept=FALSE, fixed.eb="none",
                      full.post=TRUE, init=NULL, control=control)

# models with external covariates
fit2.semnig <- semnig(x=x, y=y, C=C1, Z=Z, unpenalized=NULL,
                      standardize=FALSE, intercept=FALSE, fixed.eb=FALSE,
                      full.post=TRUE, init=NULL, control=control)
fit3.semnig <- semnig(x=x, y=y, C=C2, Z=Z, unpenalized=NULL,
                      standardize=FALSE, intercept=FALSE, fixed.eb=FALSE,
                      full.post=TRUE, init=NULL, control=control)
fit4.semnig <- semnig(x=x, y=y, C=C3, Z=Z, unpenalized=NULL,
                      standardize=FALSE, intercept=FALSE, fixed.eb=FALSE,
                      full.post=TRUE, init=NULL, control=control)

