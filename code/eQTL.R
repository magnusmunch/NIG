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
library(glmnet)

### load data
load("data/P38MAPKpathwayKb100K.RData")

### data preparation
y <- sapply(Gn, "[[", 1)
x <- lapply(sapply(Gn, "[[", 2), scale)
G <- length(Gn)
p <- sapply(x, ncol)
n <- nrow(y)

### creating co-data
inorout <- lapply(Gn, "[[", 4)
distance <- lapply(Gn, function(g) {  
  sapply(g[[5]], function(s) {
    (sign(prod(s - g[[3]][c(1:2)]))!=1)*
      min(abs(s - g[[3]][c(1:2)]))})})
genelength <- lapply(1:G, function(g) {
  rep(diff(as.numeric(Gn[[g]][[3]][c(1, 2)])), p[g])})
minoraf <- lapply(x, function(g) {f <- apply(g, 2, function(s) {
  f <- mean(s); unname(min(f, 1 - f))})})
# library(biomaRt)
# ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# # View(listFilters(ensembl))
# gb <- getBM(attributes=c('ensembl_exon_id', "exon_chrom_start", 
#                          "exon_chrom_end"), filters=c("hgnc_symbol"), 
#             values=Gn[[1]][[7]]$gene_name, mart=ensembl)
C <- sapply(1:G, function(g) {
  cbind(inorout=inorout[[g]], distance=distance[[g]], 
        genelength=genelength[[g]], minoraf=minoraf[[g]])})

### fitting models
# ENIG model
fit.enig <- enig(x, y, C, mult.lambda=FALSE, fixed.eb="none", init=NULL, 
                 control=list(epsilon.eb=1e-3, epsilon.vb=1e-3, 
                              epsilon.opt=sqrt(.Machine$double.eps), 
                              maxit.eb=30, maxit.vb=2, maxit.opt=100, 
                              trace=TRUE))

# CV penalized likelihood
fit.cv.glmnet <- lapply(1:G, function(g) {
  cv.glmnet(x[[g]], y[, g], intercept=FALSE, alpha=0.5)})

# MML penalized likelihood
opt.res <- lapply(1:G, function(g) {
  svd.x <- svd(x[[g]])
  optim(par=c(0, 0), ridge.mll, sv=svd.x$d, uty=t(svd.x$u) %*% y[, g], 
        yty=sum(y[, g]^2), n=n, p=p[g], 
        method="BFGS", control=list(fnscale=-1, maxit=1000))})
opt.par <- sqrt(exp(t(sapply(opt.res, "[[", 1))))
colnames(opt.par) <- c("sigma", "gamma")
opt.conv <- sapply(opt.res, "[[", 4)
opt.iter <- t(sapply(opt.res, "[[", 3))
lambda <- 1/(n*opt.par[, "gamma"]^2)
fit.mml.glmnet <- lapply(1:G, function(g) {
  glmnet(x[[g]], y[, g], intercept=FALSE, alpha=0, lambda=lambda[g])})

### point estimates
best.enig <- fit.enig$vb$mu
best.cv.glmnet <- sapply(fit.cv.glmnet, function(g) {
  as.numeric(coef(g, s="lambda.min"))[-1]})
best.mml.glmnet <- sapply(fit.mml.glmnet, function(g) {
  as.numeric(coef(g))[-1]})


plot(unlist(best.glmnet), unlist(best.enig))

### prior estimates
scatter.smooth(unlist(distance), unlist(fit.enig$eb$mprior), col=2)
plot(fit.loess, col=3)
fit.enig$vb$mu
sapply(genelength, "[[", 1)
plot(unlist(distance), unlist(fit.enig$eb$mprior))
fit.enig$eb$alpha["genelength"]


     
     