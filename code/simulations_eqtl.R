#!/usr/bin/env Rscript

### installation of packages
if(!("devtools" %in% installed.packages())) {
  install.packages("devtools")
}
library(devtools)
install_github("magnusmunch/cambridge/rpackage", local=FALSE,
               auth_token=Sys.getenv("GITHUB_PAT"))

### libraries
library(cambridge)
library(statmod)
library(glmnet)

### load data
load("data/P38MAPKpathwayKb100K.RData")
rm(G1, G2)

### data preparation
y <- sapply(Gn, "[[", 1)
x <- lapply(sapply(Gn, "[[", 2), scale)
G <- length(Gn)
p <- sapply(x, ncol)
n <- nrow(y)

### creating co-data
inorout <- lapply(Gn, "[[", 4)

### simulation settings
msigma.sq <- 1
vsigma.sq <- 0
mgamma.sq <- c(1, 2)
vgamma.sq <- rep(0, G)

### simulate model parameters
set.seed(2019)
if(vsigma.sq==0) {
  sigma.sq <- rep(msigma.sq, G)
} else {
  sigma.sq <- 1/rgamma(G, msigma.sq^2/vsigma.sq + 2, 
                       msigma.sq^3/vsigma.sq + msigma.sq)  
}
gamma.sq <- sapply(1:G, function(g) {
  if(vgamma.sq[g]==0) {
    mgamma.sq[inorout[[g]] + 1]
  } else {
    rinvgauss(p[g], mgamma.sq[inorout[[g]] + 1], 
              mgamma.sq[inorout[[g]] + 1]^3/vgamma.sq[g])
  }})


beta <- sapply(1:G, function(g) {
  rnorm(p[g], 0, sigma.sq[g]*gamma.sq[[g]])})

### simulate data
set.seed(2019)
ntrain <- floor(n/2)
ntest <- n - ntrain
id.train <- sample(1:n, ntrain)
xtrain <- sapply(x, function(s) {scale(s[id.train, ])})
xtest <- sapply(x, function(s) {scale(s[-id.train, ])})
ytrain <- scale(sapply(1:G, function(g) {
  rnorm(ntrain, as.numeric(xtrain[[g]] %*% beta[[g]]), sqrt(sigma.sq[[g]]))}))
ytest <- scale(sapply(1:G, function(g) {
  rnorm(ntest, as.numeric(xtest[[g]] %*% beta[[g]]), sqrt(sigma.sq[[g]]))}))

### fitting models
control <- list(conv.post=TRUE, trace=TRUE,
                epsilon.eb=1e-3, epsilon.vb=1e-3, 
                epsilon.opt=sqrt(.Machine$double.eps),
                maxit.eb=500, maxit.vb=2, maxit.opt=100,
                maxit.post=100)

C <- sapply(1:G, function(g) {
  as.matrix(inorout[[g]])})

# set.seed(2019)
# fit1.enig <- semnig(xtrain, ytrain, C=NULL, unpenalized=NULL, standardize=FALSE,
#                     intercept=FALSE, intercept.eb=TRUE, mult.lambda=TRUE,
#                     fixed.eb="none", full.post=TRUE, init=NULL, control=control)
# fit2.enig <- semnig(xtrain, ytrain, C=C, unpenalized=NULL, standardize=FALSE,
#                     intercept=FALSE, intercept.eb=TRUE, mult.lambda=TRUE,
#                     fixed.eb="none", full.post=TRUE, init=NULL, control=control)
# fit3.enig <- semnig(xtrain, ytrain, C=NULL, unpenalized=NULL, standardize=FALSE,
#                     intercept=FALSE, intercept.eb=TRUE, mult.lambda=TRUE,
#                     fixed.eb="lambda", full.post=TRUE, init=NULL,
#                     control=control)
# fit4.enig <- semnig(xtrain, ytrain, C=C, unpenalized=NULL, standardize=FALSE,
#                     intercept=FALSE, intercept.eb=TRUE, mult.lambda=TRUE,
#                     fixed.eb="lambda", full.post=TRUE, init=NULL,
#                     control=control)
# fit5.enig <- semnig(xtrain, ytrain, C=NULL, unpenalized=NULL, standardize=FALSE,
#                     intercept=FALSE, intercept.eb=TRUE, mult.lambda=TRUE,
#                     fixed.eb="lambda", full.post=TRUE,
#                     init=list(lambda=rep(0.1, G)),
#                     control=control)
# fit6.enig <- semnig(xtrain, ytrain, C=C, unpenalized=NULL, standardize=FALSE,
#                     intercept=FALSE, intercept.eb=TRUE, mult.lambda=TRUE,
#                     fixed.eb="lambda", full.post=TRUE,
#                     init=list(lambda=rep(0.1, G)),
#                     control=control)
# fit1.glmnet <- lapply(1:G, function(g) {
#   cv.glmnet(xtrain[[g]], ytrain[, g], family="gaussian", alpha=0, 
#             standardize=FALSE, intercept=FALSE)})
# fit2.glmnet <- lapply(1:G, function(g) {
#   cv.glmnet(xtrain[[g]], ytrain[, g], family="gaussian", alpha=1, 
#             standardize=FALSE, intercept=FALSE)})
# save(fit1.enig, fit2.enig, fit3.enig, fit4.enig, fit5.enig, fit6.enig, 
#      fit1.glmnet, fit2.glmnet, file="results/simulations_eqtl_fit1.Rdata")
load(file="results/simulations_eqtl_fit1.Rdata")

### alpha estimates, EMSE, and PMSE
library(knitr)
emse <- cbind(
  sapply(c(1:G), function(g) {mean((beta[[g]] - fit1.enig$vb$mu[[g]])^2)}),
  sapply(c(1:G), function(g) {mean((beta[[g]] - fit2.enig$vb$mu[[g]])^2)}),
  sapply(c(1:G), function(g) {mean((beta[[g]] - fit3.enig$vb$mu[[g]])^2)}),
  sapply(c(1:G), function(g) {mean((beta[[g]] - fit4.enig$vb$mu[[g]])^2)}),
  sapply(c(1:G), function(g) {mean((beta[[g]] - fit5.enig$vb$mu[[g]])^2)}),
  sapply(c(1:G), function(g) {mean((beta[[g]] - fit6.enig$vb$mu[[g]])^2)}),
  sapply(c(1:G), function(g) {
    mean((beta[[g]] - coef(fit1.glmnet[[g]], "lambda.min")[-1, ])^2)}),
  sapply(c(1:G), function(g) {
    mean((beta[[g]] - coef(fit2.glmnet[[g]], "lambda.min")[-1, ])^2)}))
colnames(emse) <- c(paste0("model ", 1:6), "cv ridge", "cv lasso")

pmse <- cbind(
  sapply(1:G, function(g) {
    mean((as.numeric(xtest[[g]] %*% fit1.enig$vb$mu[[g]]) - ytest[, g])^2)}),
  sapply(1:G, function(g) {
    mean((as.numeric(xtest[[g]] %*% fit2.enig$vb$mu[[g]]) - ytest[, g])^2)}),
  sapply(1:G, function(g) {
    mean((as.numeric(xtest[[g]] %*% fit3.enig$vb$mu[[g]]) - ytest[, g])^2)}),
  sapply(1:G, function(g) {
    mean((as.numeric(xtest[[g]] %*% fit4.enig$vb$mu[[g]]) - ytest[, g])^2)}),
  sapply(1:G, function(g) {
    mean((as.numeric(xtest[[g]] %*% fit5.enig$vb$mu[[g]]) - ytest[, g])^2)}),
  sapply(1:G, function(g) {
    mean((as.numeric(xtest[[g]] %*% fit6.enig$vb$mu[[g]]) - ytest[, g])^2)}),
  sapply(c(1:G), function(g) {
    mean((predict(fit1.glmnet[[g]], xtest[[g]],  "lambda.min") - 
            ytest[, g])^2)}),
  sapply(c(1:G), function(g) {
    mean((predict(fit2.glmnet[[g]], xtest[[g]],  "lambda.min") - 
            ytest[, g])^2)}))
colnames(pmse) <- c(paste0("model ", 1:6), "cv ridge", "cv lasso")

pairs(emse, lower.panel=function(x, y, ...) {
  points(x, y, ...); abline(a=0, b=1, ...)}, 
  upper.panel=function(x, y, ...) {
    points(x, y, ...); abline(a=0, b=1, ...)},
  diag.panel=function(x, ...) {
    lx <- sum(range(x))/2; ly <- min(x) + diff(range(x))/4
    text(lx, ly, paste0("average EMSE \n", round(mean(x), 4), ...))},
  label.pos=2/3)

pairs(pmse, lower.panel=function(x, y, ...) {
  points(x, y, ...); abline(a=0, b=1, ...)}, 
  upper.panel=function(x, y, ...) {
    points(x, y, ...); abline(a=0, b=1, ...)},
  diag.panel=function(x, ...) {
    lx <- sum(range(x))/2; ly <- min(x) + diff(range(x))/4
    text(lx, ly, paste0("average PMSE \n", round(mean(x), 4), ...))},
  label.pos=2/3)

tab <- cbind(round(rbind(c(fit1.enig$eb$alpha, NA),
                         fit2.enig$eb$alpha,
                         c(fit3.enig$eb$alpha, NA),
                         fit4.enig$eb$alpha,
                         c(fit5.enig$eb$alpha, NA),
                         fit6.enig$eb$alpha), 2), 
             round(apply(pmse, 2, mean), 3), round(apply(pmse, 2, mean), 3))
colnames(tab) <- c("intercept", "in gene", "EMSE", "PMSE")

tab[which.min(tab[, 3]), 3] <- kableExtra::cell_spec(
  tab[which.min(tab[, 3]), 3], format="latex", bold=TRUE)
tab[which.min(tab[, 4]), 4] <- kableExtra::cell_spec(
  tab[which.min(tab[, 4]), 4], format="latex", bold=TRUE)

rownames(tab) <- c("intercept only", "in gene", 
                   "intercept only, lambda=1", "in gene, lambda=1",
                   "intercept only, lambda=0.1", 
                   "in gene, lambda=0.1")
options(knitr.kable.NA = '')
kableExtra::kable_styling(kableExtra::column_spec(
  knitr::kable(tab, align=rep("r", ncol(tab)),
               col.names=c("intercept", "in gene", "EMSE", "PMSE"),
               caption="Estimated $\\hat{\\bm{\\alpha}}$ hyperparameters, 
               MSE (lowest in bold), and LPML (highest in bold).", 
               "latex", booktabs=TRUE, escape=FALSE, 
               digits=c(rep(2, 6), 4, 4)), 8, border_left=TRUE),
  latex_options=c("HOLD_position"))

### histogram of prior
nbreaks <- 30
h1 <- hist(unlist(fit1.enig$eb$mprior)*rep(fit1.enig$vb$mpost$sigma.sq, 
                                           times=p), plot=FALSE, breaks=nbreaks)
h2 <- hist(unlist(fit2.enig$eb$mprior)*rep(fit2.enig$vb$mpost$sigma.sq, 
                                           times=p), plot=FALSE, breaks=nbreaks)
h3 <- hist(unlist(fit3.enig$eb$mprior)*rep(fit3.enig$vb$mpost$sigma.sq, 
                                           times=p), plot=FALSE, breaks=nbreaks)
h4 <- hist(unlist(fit4.enig$eb$mprior)*rep(fit4.enig$vb$mpost$sigma.sq, 
                                           times=p), plot=FALSE, breaks=nbreaks)
h5 <- hist(unlist(fit5.enig$eb$mprior)*rep(fit5.enig$vb$mpost$sigma.sq, 
                                           times=p), plot=FALSE, breaks=nbreaks)
h6 <- hist(unlist(fit6.enig$eb$mprior)*rep(fit6.enig$vb$mpost$sigma.sq, 
                                           times=p), plot=FALSE, breaks=nbreaks)

breaks <- seq(min(c(h1$breaks, h2$breaks, h3$breaks, h4$breaks, h5$breaks,
                    h6$breaks)),
              max(c(h1$breaks, h2$breaks, h3$breaks, h4$breaks, h5$breaks,
                    h6$breaks)), 
              length.out=nbreaks)
h1 <- hist(unlist(fit1.enig$eb$mprior)*rep(fit1.enig$vb$mpost$sigma.sq, 
                                           times=p), 
           plot=FALSE, breaks=breaks)
h2 <- hist(unlist(fit2.enig$eb$mprior)*rep(fit2.enig$vb$mpost$sigma.sq, 
                                           times=p), 
           plot=FALSE, breaks=breaks)
h3 <- hist(unlist(fit3.enig$eb$mprior)*rep(fit3.enig$vb$mpost$sigma.sq, 
                                           times=p), 
           plot=FALSE, breaks=breaks)
h4 <- hist(unlist(fit4.enig$eb$mprior)*rep(fit4.enig$vb$mpost$sigma.sq, 
                                           times=p), 
           plot=FALSE, breaks=breaks)
h5 <- hist(unlist(fit5.enig$eb$mprior)*rep(fit5.enig$vb$mpost$sigma.sq, 
                                           times=p), 
           plot=FALSE, breaks=breaks)
h6 <- hist(unlist(fit6.enig$eb$mprior)*rep(fit6.enig$vb$mpost$sigma.sq, 
                                           times=p), 
           plot=FALSE, breaks=breaks)

library(sp)
labels <- c("model 1", "model 2", "model 3", "model 4", "model 5", "model 6")
methods <- labels
col <- bpy.colors(length(methods), cutoff.tail=0.3)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.1, 1, 1))
barplot(height=rbind(h1$density, h2$density, h3$density, h4$density, 
                     h5$density, h6$density), 
        beside=TRUE, xlab=expression(V(beta~"|"~hat(sigma)^2)), ylab="Density", 
        names.arg=round(h4$mids, 3), col=col, legend.text=labels, 
        args.legend=list(x="topright", fill=col, border=col), xlim=c(0, 70))
par(opar)

### investigating estimates
library(sp)
g <- 1
col <- bpy.colors(length(unique(inorout[[g]])), cutoff.tail=0.3)
ylim <- range(c(fit5.enig$vb$mu[[g]], fit6.enig$vb$mu[[g]]))
opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.1, 1, 1))
layout(matrix(rep(c(1:2), each=2, times=2), byrow=TRUE, ncol=4, nrow=2))
plot(fit5.enig$vb$mu[[g]], col=col[inorout[[g]] + 1], xlab="SNP",
     ylab=bquote(hat(E)*(beta*"|"*bold(y))), ylim=ylim)
plot(fit6.enig$vb$mu[[g]], col=col[inorout[[g]] + 1], xlab="SNP",
     ylab=bquote(hat(E)*(beta*"|"*bold(y))), ylim=ylim)
par(opar)

### convergence plots
opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.1, 1, 1))
layout(matrix(rep(c(1:2), each=2, times=2), byrow=TRUE, ncol=4, nrow=2))
plot(fit1.enig$seq.eb$alpha, type="l", ylab=expression(hat(alpha)),
     xlab="Iteration", main="(a)")
plot(fit1.enig$seq.eb$lambda[, 1], ylim=range(fit1.enig$seq.eb$lambda), 
     type="l", ylab=expression(hat(lambda)[d]), xlab="Iteration", main="(b)")
for(g in 2:ncol(fit1.enig$seq.eb$lambda)) {
  lines(fit1.enig$seq.eb$lambda[, g], col=g)
}
par(opar)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.1, 1, 1))
layout(matrix(rep(c(1:2), each=2, times=2), byrow=TRUE, ncol=4, nrow=2))
plot(fit2.enig$seq.eb$alpha[, 1], type="l", ylab=expression(hat(alpha)),
     xlab="Iteration", main="(a)", ylim=range(fit2.enig$seq.eb$alpha))
lines(fit2.enig$seq.eb$alpha[, 2], col=2)
plot(fit2.enig$seq.eb$lambda[, 1], ylim=range(fit2.enig$seq.eb$lambda), 
     type="l", ylab=expression(hat(lambda)[d]), xlab="Iteration", main="(b)")
for(g in 2:ncol(fit2.enig$seq.eb$lambda)) {
  lines(fit2.enig$seq.eb$lambda[, g], col=g)
}
par(opar)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.1, 1, 1))
plot(fit3.enig$seq.eb$alpha, type="l", ylab=expression(hat(alpha)),
     xlab="Iteration")
par(opar)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.1, 1, 1))
plot(fit4.enig$seq.eb$alpha[, 1], type="l", ylab=expression(hat(alpha)),
     xlab="Iteration", ylim=range(fit4.enig$seq.eb$alpha))
lines(fit4.enig$seq.eb$alpha[, 2], col=2)
par(opar)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.1, 1, 1))
plot(fit5.enig$seq.eb$alpha, type="l", ylab=expression(hat(alpha)),
     xlab="Iteration")
par(opar)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.1, 1, 1))
plot(fit6.enig$seq.eb$alpha[, 1], type="l", ylab=expression(hat(alpha)),
     xlab="Iteration", ylim=range(fit6.enig$seq.eb$alpha))
lines(fit6.enig$seq.eb$alpha[, 2], col=2)
par(opar)


