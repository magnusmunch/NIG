#!/usr/bin/env Rscript

################################################################################
# regression gene expression on SNP data
################################################################################

# number of cores to use
ncores <- 50

### libraries
packages <- c("foreach", "doParallel", "cambridge", "glmnet", "pInc", "rstan")
sapply(packages, library, character.only=TRUE)

### load data
load(file="data/data_eqtl_dat1.Rdata")

################################## analysis 1 ##################################
### data preparation
y <- scale(expr.prep, scale=FALSE)
x <- lapply(snps.prep, scale, scale=FALSE)
G <- length(x)
p <- sapply(x, ncol)
n <- nrow(y)

### model fitting
# setting seed
set.seed(2019)

# estimation settings
control <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, epsilon.vb=1e-3,
                maxit.eb=200, maxit.vb=1, maxit.post=100)
methods <- c("NIG1", "NIG2", "lasso", "ridge", "bSEM")

# model with external covariates
Z <- model.matrix(~ length, data=gene.prep)
colnames(Z) <- c("intercept", "length")
C <- sapply(1:G, function(g) {
  mat <- cbind(model.matrix( ~ factor(feat.prep$exon[[g]], 
                                      levels=c("inter", "intron", "exon"))), 
               feat.prep$distance[[g]],
               feat.prep$minoraf[[g]])
  colnames(mat) <- c("intercept", "intron", "exon", "distance", "minoraf")
  return(mat)})
fit.semnig1 <- semnig(x=x, y=y, C=C, Z=Z, unpenalized=NULL,
                      standardize=FALSE, intercept=FALSE, fixed.eb="none",
                      full.post=TRUE, init=NULL, control=control)

# model without external covariates
Z <- matrix(1, nrow=G)
colnames(Z) <- c("intercept")
C <- lapply(p, function(s) {
  mat <- matrix(1, nrow=s, ncol=1); colnames(mat) <- "intercept"; return(mat)})
fit.semnig2 <- semnig(x=x, y=y, C=C, Z=Z, unpenalized=NULL,
                      standardize=FALSE, intercept=FALSE, fixed.eb="none",
                      full.post=TRUE, init=NULL, control=control)

# bSEM model
fit.bSEM1 <- bSEM(x, split(y, 1:ncol(y)), 
                  lapply(lapply(feat.prep$exon, "!=", "inter"), "+", 1),
                  control=list(maxit=200, trace=TRUE, epsilon=1e-3))

# penalized regression models
fit.lasso1 <- lapply(1:G, function(g) {
  cv.glmnet(x[[g]], y[, g], intercept=FALSE)})
fit.ridge1 <- lapply(1:G, function(g) {
  cv.glmnet(x[[g]], y[, g], alpha=0, intercept=FALSE)})

# saving fitted model objects
save(fit.semnig1, fit.semnig2, fit.lasso1, fit.ridge1, fit.bSEM1, 
     file="results/analysis_eqtl_fit1.Rdata")

# EB estimates
tab <- rbind(c(fit.semnig1$eb$alphaf, fit.semnig1$eb$alphad), 
             c(fit.semnig2$eb$alphaf, rep(NA, 4), 
               fit.semnig2$eb$alphad, NA),
             c(fit.bSEM1$eb$a[nrow(fit.bSEM1$eb$a), 1]/
                 fit.bSEM1$eb$b[nrow(fit.bSEM1$eb$b), 1],
               fit.bSEM1$eb$a[nrow(fit.bSEM1$eb$a), 2]/
                 fit.bSEM1$eb$b[nrow(fit.bSEM1$eb$b), 2] -
                 fit.bSEM1$eb$a[nrow(fit.bSEM1$eb$a), 1]/
                 fit.bSEM1$eb$b[nrow(fit.bSEM1$eb$b), 1], rep(NA, 5)))
rownames(tab) <- methods[c(1:2, 5)]
write.table(tab, file="results/analysis_eqtl_fit1.txt")

### splitting data
set.seed(2019)
ntrain <- floor(n/2)
id.train <- sample(1:n, ntrain)
ytrain <- y[id.train, ]
ytest <- y[-id.train, ]
xtrain <- sapply(x, function(s) {s[id.train, ]})
xtest <- sapply(x, function(s) {s[-id.train, ]})

# model with external covariates
Z <- model.matrix(~ length, data=gene.prep)
colnames(Z) <- c("intercept", "length")
C <- sapply(1:G, function(g) {
  mat <- cbind(model.matrix( ~ factor(feat.prep$exon[[g]], 
                                      levels=c("inter", "intron", "exon"))), 
               feat.prep$distance[[g]],
               feat.prep$minoraf[[g]])
  colnames(mat) <- c("intercept", "intron", "exon", "distance", "minoraf")
  return(mat)})
cv.semnig1 <- semnig(x=xtrain, y=ytrain, C=C, Z=Z, unpenalized=NULL,
                     standardize=FALSE, intercept=FALSE, fixed.eb="none",
                     full.post=TRUE, init=NULL, control=control)

# model without external covariates
Z <- matrix(1, nrow=G)
colnames(Z) <- c("intercept")
C <- lapply(p, function(s) {
  mat <- matrix(1, nrow=s, ncol=1); colnames(mat) <- "intercept"; return(mat)})
cv.semnig2 <- semnig(x=xtrain, y=ytrain, C=C, Z=Z, unpenalized=NULL,
                     standardize=FALSE, intercept=FALSE, fixed.eb="none",
                     full.post=TRUE, init=NULL, control=control)

# bSEM model
cv.bSEM1 <- bSEM(xtrain, split(ytrain, 1:ncol(ytrain)), 
                 lapply(lapply(feat.prep$exon, "!=", "inter"), "+", 1),
                 control=list(maxit=200, trace=TRUE, epsilon=1e-3))

# penalized regression models
cv.lasso1 <- lapply(1:G, function(g) {
  cv.glmnet(xtrain[[g]], ytrain[, g], intercept=FALSE, standardize=FALSE)})
cv.ridge1 <- lapply(1:G, function(g) {
  cv.glmnet(xtrain[[g]], ytrain[, g], alpha=0, intercept=FALSE, 
            standardize=FALSE)})

# calculating average prediction mean squared error
pmse <- c(mean(sapply(1:G, function(d) {
  mean((ytest[, d] - xtest[[d]] %*% cv.semnig1$vb$mpost$beta[[d]])^2)})),
  mean(sapply(1:G, function(d) {
    mean((ytest[, d] - xtest[[d]] %*% cv.semnig2$vb$mpost$beta[[d]])^2)})),
  mean((ytest - sapply(1:G, function(d) {
    predict(cv.lasso1[[d]], xtest[[d]], s="lambda.min")}))^2),
  mean((ytest - sapply(1:G, function(d) {
    predict(cv.ridge1[[d]], xtest[[d]], s="lambda.min")}))^2),
  mean(sapply(1:G, function(d) {
    mean((ytest[, d] - xtest[[d]] %*% cv.bSEM1$vb$beta[[d]][, 1])^2)})),
  mean(apply(ytest, 1, "-", colMeans(ytrain))^2))

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
lpml <- c(sum(logcpo(xtest, ytest, ntrain, cv1.semnig), na.rm=TRUE),
          sum(logcpo(xtest, ytest, ntrain, cv2.semnig), na.rm=TRUE),
          sum(logcpo(xtest, ytest, ntrain, cv3.semnig), na.rm=TRUE),
          sum(logcpo(xtest, ytest, ntrain, cv4.semnig), na.rm=TRUE))





### fitting models
control <- list(conv.post=TRUE, trace=TRUE,
                epsilon.eb=1e-3, epsilon.vb=1e-3, 
                epsilon.opt=sqrt(.Machine$double.eps),
                maxit.eb=100, maxit.vb=2, maxit.opt=200,
                maxit.post=100)

set.seed(2019)
cv.semnig1 <- semnig(x=xtrain, y=ytrain, C=NULL, unpenalized=NULL, 
                      standardize=FALSE, intercept=FALSE, intercept.eb=TRUE,
                      mult.lambda=TRUE, fixed.eb="lambda", full.post=TRUE, 
                      init=list(lambda=rep(1, G)), control=control, 
                      devel=list(cpp=TRUE, model="semnig2", 
                                 fixed.eb.intercept=FALSE))
fit2.semnig <- semnig(x=xtrain, y=ytrain, C=C2, unpenalized=NULL, 
                      standardize=FALSE, intercept=FALSE, intercept.eb=TRUE,
                      mult.lambda=TRUE, fixed.eb="lambda", full.post=TRUE, 
                      init=list(lambda=rep(1, G)), control=control, 
                      devel=list(cpp=TRUE, model="semnig2", 
                                 fixed.eb.intercept=FALSE))
fit3.semnig <- semnig(x=xtrain, y=ytrain, C=C3, unpenalized=NULL, 
                      standardize=FALSE, intercept=FALSE, intercept.eb=TRUE,
                      mult.lambda=TRUE, fixed.eb="lambda", full.post=TRUE, 
                      init=list(lambda=rep(1, G)), control=control, 
                      devel=list(cpp=TRUE, model="semnig2", 
                                 fixed.eb.intercept=FALSE))
fit4.semnig <- semnig(x=xtrain, y=ytrain, C=C3, unpenalized=NULL, 
                      standardize=FALSE, intercept=FALSE, intercept.eb=FALSE,
                      mult.lambda=TRUE, fixed.eb="lambda", full.post=TRUE, 
                      init=list(lambda=rep(1, G)), control=control, 
                      devel=list(cpp=TRUE, model="semnig2", 
                                 fixed.eb.intercept=FALSE))
fit5.semnig <- semnig(x=xtrain, y=ytrain, C=C2, unpenalized=NULL, 
                      standardize=FALSE, intercept=FALSE, intercept.eb=TRUE,
                      mult.lambda=TRUE, fixed.eb="none", full.post=TRUE, 
                      init=list(lambda=rep(1, G), alpha=c(100, 0)), 
                      control=control, 
                      devel=list(cpp=TRUE, model="semnig3", 
                                 fixed.eb.intercept=FALSE, print.eb=TRUE))
fit6.semnig <- semnig(x=xtrain, y=ytrain, C=C2, unpenalized=NULL, 
                      standardize=FALSE, intercept=FALSE, intercept.eb=TRUE,
                      mult.lambda=TRUE, fixed.eb="lambda", full.post=TRUE, 
                      init=list(lambda=rep(0.1, G), alpha=), control=control, 
                      devel=list(cpp=TRUE, model="semnig3", 
                                 fixed.eb.intercept=FALSE, print.eb=TRUE))
lambda <- 0.1
alpha <- unname(c(fit1.pInc$alpha0[1]/fit1.pInc$alpha0[2] - 1/lambda, 
                  fit1.pInc$alpha1[1]/fit1.pInc$alpha1[2] -
                    fit1.pInc$alpha0[1]/fit1.pInc$alpha0[2]))
fit7.semnig <- semnig(x=xtrain, y=ytrain, C=C2, unpenalized=NULL, 
                      standardize=FALSE, intercept=FALSE, intercept.eb=TRUE,
                      mult.lambda=TRUE, fixed.eb="lambda", full.post=TRUE, 
                      init=list(lambda=rep(lambda, G), alpha=alpha), 
                      control=control, 
                      devel=list(cpp=TRUE, model="semnig3", 
                                 fixed.eb.intercept=FALSE, print.eb=TRUE))
fit8.semnig <- semnig(x=xtrain, y=ytrain, C=C2, unpenalized=NULL, 
                      standardize=FALSE, intercept=FALSE, intercept.eb=TRUE,
                      mult.lambda=FALSE, fixed.eb="none", full.post=TRUE, 
                      control=control, 
                      devel=list(cpp=TRUE, model="semnig4", 
                                 fixed.eb.intercept=FALSE, print.eb=TRUE))



### Ginos model
library(graphics)
library(stats)
library(stats4)
library(devtools)
library(grDevices)
library(datasets)
library(Matrix)
library(mvtnorm)
library(rags2ridges)
library(gplots)
library(foreach)
library(glmnet)
library(gsl)
library(corpcor)
library(ROCR)

source("misc/Incorporating_prior_information_and_borrowing_information_in_high-dimensional_sparse_regression_using_the_horseshoe_and_variational_Bayes/code/functions/pInc.R")                   # function to analyse multiple related high-dimensional and complex datasets
source("misc/Incorporating_prior_information_and_borrowing_information_in_high-dimensional_sparse_regression_using_the_horseshoe_and_variational_Bayes/code/functions/spLInit3.R")
source("misc/Incorporating_prior_information_and_borrowing_information_in_high-dimensional_sparse_regression_using_the_horseshoe_and_variational_Bayes/code/functions/spOneIterVBi3.R")
source("misc/Incorporating_prior_information_and_borrowing_information_in_high-dimensional_sparse_regression_using_the_horseshoe_and_variational_Bayes/code/functions/spGetMargLiki3.R")
source("misc/Incorporating_prior_information_and_borrowing_information_in_high-dimensional_sparse_regression_using_the_horseshoe_and_variational_Bayes/code/functions/spGetPostMeani3.R")

fit1.pInc <- pInc(Gtrain, varSel="threshold", maxiter=200)

### ridge model
fit1.ridge <- lapply(1:G, function(g) {
  svd.x <- svd(xtrain[[g]], nu=n)
  opt <- optim(par=c(0, 0), fn=ridge.logmll, method="BFGS",
               uty=as.numeric(t(svd.x$u) %*% ytrain[, g]),
               sv=svd.x$d, n=ntrain, p=p[g], control=list(fnscale=-1))
  fit <- glmnet(xtrain[[g]], ytrain[, g], family="gaussian", alpha=0,
                lambda=1/(ntrain*exp(opt$par[2])), intercept=FALSE, 
                standardize=FALSE)
  fit$opt <- opt
  return(fit)})
fit2.ridge <- lapply(1:G, function(g) {
  fit <- cv.glmnet(x[[g]], y[, g], family="gaussian", alpha=0, 
                   standardize=FALSE, intercept=FALSE)
  return(fit)})

### save fitted models
# save(fit1.semnig, fit2.semnig, fit3.semnig, fit4.semnig, fit5.semnig,
#      fit6.semnig, fit7.semnig, fit8.semnig, fit1.pInc,
#      fit1.ridge, fit2.ridge, file="results/eQTL_fit1.Rdata")
load(file="results/eQTL_fit1.Rdata")

### prediction mean squared error
pmse <- cbind(
  sapply(1:G, function(g) {mean((ytest[, g])^2)}),
  sapply(1:G, function(g) {
    mean((as.numeric(xtest[[g]] %*% fit1.semnig$vb$mu[[g]]) - ytest[, g])^2)}),
  sapply(1:G, function(g) {
    mean((as.numeric(xtest[[g]] %*% fit2.semnig$vb$mu[[g]]) - ytest[, g])^2)}),
  sapply(1:G, function(g) {
    mean((as.numeric(xtest[[g]] %*% fit3.semnig$vb$mu[[g]]) - ytest[, g])^2)}),
  sapply(1:G, function(g) {
    mean((as.numeric(xtest[[g]] %*% fit4.semnig$vb$mu[[g]]) - ytest[, g])^2)}),
  sapply(1:G, function(g) {
    mean((as.numeric(xtest[[g]] %*% fit5.semnig$vb$mu[[g]]) - ytest[, g])^2)}),
  sapply(1:G, function(g) {
    mean((as.numeric(xtest[[g]] %*% fit6.semnig$vb$mu[[g]]) - ytest[, g])^2)}),
  sapply(1:G, function(g) {
    mean((as.numeric(xtest[[g]] %*% fit7.semnig$vb$mu[[g]]) - ytest[, g])^2)}),
  sapply(1:G, function(g) {
    mean(as.numeric(xtest[[g]] %*% fit1.pInc$Betas[[g]][, 1] - 
                      ytest[, g])^2)}),
  sapply(1:G, function(g) {
    mean((predict(fit1.ridge[[g]], xtest[[g]]) - ytest[, g])^2)}),
  sapply(1:G, function(g) {
    mean((predict(fit2.ridge[[g]], xtest[[g]], s="lambda.min") - 
            ytest[, g])^2)}))
colnames(pmse) <- c("null", "intercept", "in gene", "full with EB intercept", 
                    "full without EB intercept", "horseshoe", "horseshoe2",
                    "horseshoe3",
                    "pInc", "MLL ridge", "CV ridge")

pairs(pmse, lower.panel=function(x, y, ...) {
  points(x, y, ...); abline(a=0, b=1, ...)}, 
  upper.panel=function(x, y, ...) {
    points(x, y, ...); abline(a=0, b=1, ...)},
  diag.panel=function(x, ...) {
    lx <- sum(range(x))/2; ly <- min(x) + diff(range(x))/4
    text(lx, ly, paste0("average PMSE \n", round(mean(x), 4), ...))},
  label.pos=2/3)



### convergence
library(sp)
col <- bpy.colors(ncol(C3[[1]]), cutoff.tail=0.3)
opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.1, 1, 1))
layout(matrix(c(rep(c(1:2), each=2, times=2), rep(c(3:4), each=2, times=2),
                rep(c(5:6), each=2, times=2)), 
              byrow=TRUE, ncol=4, nrow=6))
plot(fit1.semnig$seq.eb$alpha[, 1], type="l", xlab="iteration",
     ylab=expression(alpha), main="(a)", ylim=range(fit1.semnig$seq.eb$alpha),
     col=col[1])
legend("topleft", col=col, legend=names(fit3.semnig$eb$alpha), lty=1, seg.len=1)
plot(fit2.semnig$seq.eb$alpha[, 1], type="l", xlab="iteration",
     ylab=expression(alpha), main="(b)", ylim=range(fit2.semnig$seq.eb$alpha),
     col=col[1])
lines(fit2.semnig$seq.eb$alpha[, 2], col=col[2])
plot(fit3.semnig$seq.eb$alpha[, 1], type="l", xlab="iteration",
     ylab=expression(alpha), main="(c)", ylim=range(fit3.semnig$seq.eb$alpha),
     col=col[1])
lines(fit3.semnig$seq.eb$alpha[, 2], col=col[2])
lines(fit3.semnig$seq.eb$alpha[, 3], col=col[3])
lines(fit3.semnig$seq.eb$alpha[, 4], col=col[4])
lines(fit3.semnig$seq.eb$alpha[, 5], col=col[5])
plot(fit4.semnig$seq.eb$alpha[, 1], type="l", xlab="iteration",
     ylab=expression(alpha), main="(d)", ylim=range(fit4.semnig$seq.eb$alpha),
     col=col[2])
lines(fit4.semnig$seq.eb$alpha[, 2], col=col[3])
lines(fit4.semnig$seq.eb$alpha[, 3], col=col[4])
lines(fit4.semnig$seq.eb$alpha[, 4], col=col[5])
plot(fit5.semnig$seq.eb$alpha[, 1], type="l", xlab="iteration",
     ylab=expression(alpha), main="(e)", ylim=range(fit5.semnig$seq.eb$alpha),
     col=col[1])
lines(fit5.semnig$seq.eb$alpha[, 2], col=col[2])
par(opar)

### investigating estimates
library(sp)
g <- 15
col <- bpy.colors(length(unique(inorout[[g]])), cutoff.tail=0.3)
ylim <- range(c(fit1.semnig$vb$mu[[g]], fit2.semnig$vb$mu[[g]],
                fit3.semnig$vb$mu[[g]], fit4.semnig$vb$mu[[g]],
                fit1.pInc$Betas[[g]][, 1], coef(fit1.ridge[[g]])[-1, 1]))
quartz()
opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.1, 1, 1))
layout(matrix(c(rep(c(1:2), each=2, times=2), rep(c(3:4), each=2, times=2),
                rep(c(5:6), each=2, times=2), rep(c(7:8), each=2, times=2),
                rep(c(9:10), each=2, times=2)), 
              byrow=TRUE, ncol=4, nrow=10))
plot(fit1.semnig$vb$mu[[g]], col=col[inorout[[g]] + 1], xlab="SNP",
     ylab=bquote(hat(E)*(beta*"|"*bold(y))), main=colnames(pmse)[-1][1], 
     ylim=ylim)
plot(fit2.semnig$vb$mu[[g]], col=col[inorout[[g]] + 1], xlab="SNP",
     ylab=bquote(hat(E)*(beta*"|"*bold(y))), main=colnames(pmse)[-1][2], 
     ylim=ylim)
plot(fit3.semnig$vb$mu[[g]], col=col[inorout[[g]] + 1], xlab="SNP",
     ylab=bquote(hat(E)*(beta*"|"*bold(y))), main=colnames(pmse)[-1][3], 
     ylim=ylim)
plot(fit4.semnig$vb$mu[[g]], col=col[inorout[[g]] + 1], xlab="SNP",
     ylab=bquote(hat(E)*(beta*"|"*bold(y))), main=colnames(pmse)[-1][4], 
     ylim=ylim)
plot(fit5.semnig$vb$mu[[g]], col=col[inorout[[g]] + 1], xlab="SNP",
     ylab=bquote(hat(E)*(beta*"|"*bold(y))), main=colnames(pmse)[-1][5], 
     ylim=ylim)
plot(fit6.semnig$vb$mu[[g]], col=col[inorout[[g]] + 1], xlab="SNP",
     ylab=bquote(hat(E)*(beta*"|"*bold(y))), main=colnames(pmse)[-1][6], 
     ylim=ylim)
plot(fit7.semnig$vb$mu[[g]], col=col[inorout[[g]] + 1], xlab="SNP",
     ylab=bquote(hat(E)*(beta*"|"*bold(y))), main=colnames(pmse)[-1][7], 
     ylim=ylim)
plot(fit1.pInc$Betas[[g]][, 1], col=col[inorout[[g]] + 1], xlab="SNP",
     ylab=bquote(hat(E)*(beta*"|"*bold(y))), main=colnames(pmse)[-1][8], 
     ylim=ylim)
plot(coef(fit1.ridge[[g]])[-1, 1], col=col[inorout[[g]] + 1], xlab="SNP",
     ylab=bquote(hat(E)*(beta*"|"*bold(y))), main=colnames(pmse)[-1][9], 
     ylim=ylim)
plot(coef(fit2.ridge[[g]], s="lambda.min")[-1, 1], col=col[inorout[[g]] + 1], 
     xlab="SNP", ylab=bquote(hat(E)*(beta*"|"*bold(y))), 
     main=colnames(pmse)[-1][10], ylim=ylim)
par(opar)


var(fit1.semnig$vb$mu[[g]][inorout[[g]]==0])
var(fit1.semnig$vb$mu[[g]][inorout[[g]]==1])

ylim <- range(c(unlist(fit1.semnig$vb$mu), unlist(fit4.semnig$vb$mu)))
quartz()
opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.1, 1, 1))
layout(matrix(c(rep(c(1:2), each=2, times=2)), 
              byrow=TRUE, ncol=4, nrow=2))
plot(unlist(minoraf), unlist(fit1.semnig$vb$mu), ylim=ylim)
plot(unlist(minoraf), unlist(fit4.semnig$vb$mu), ylim=ylim)
par(opar)

genes1 <- c(15, 40, 48, 50, 51, 61)
genes2 <- c(85, 86, 93, 96, 98, 1)
ylim <- c(-0.01, 0.01)
quartz()
opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.1, 1, 1))
layout(matrix(c(rep(c(1:2), each=2, times=2), rep(c(3:4), each=2, times=2),
                rep(c(5:6), each=2, times=2)), 
              byrow=TRUE, ncol=4, nrow=6))
for(g in 1:length(genes1)) {
  plot(minoraf[[g]], fit1.semnig$vb$mu[[g]])
}
par(opar)
