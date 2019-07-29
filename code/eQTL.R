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
library(biomaRt)
library(glmnet)
library(MASS)

### load data
load("data/P38MAPKpathwayKb100K.RData")

### data preparation
y <- sapply(Gn, "[[", 1)
x <- lapply(sapply(Gn, "[[", 2), scale)
G <- length(Gn)
p <- sapply(x, ncol)
n <- nrow(y)

set.seed(2019)
id.train <- sample(1:n, floor(n/3))
ytrain <- y[id.train, ]
ytest <- y[-id.train, ]

xtrain <- sapply(x, function(s) {s[id.train, ]})
xtest <- sapply(x, function(s) {s[-id.train, ]})

### creating co-data
inorout <- lapply(Gn, "[[", 4)
distance <- lapply(Gn, function(g) {  
  sapply(g[[5]], function(s) {
    (sign(prod(s - g[[3]][c(1:2)]))!=1)*
      min(abs(s - g[[3]][c(1:2)]))})})
genelength <- lapply(1:G, function(g) {
  rep(diff(as.numeric(Gn[[g]][[3]][c(1, 2)])), p[g])})
minoraf <- lapply(sapply(Gn, "[[", 2), function(g) {
  f <- apply(g, 2, function(s) {f <- mean(s); unname(min(f, 1 - f))})})

### retrieve exon status of SNPs
# View(listEnsembl())
geneid <- sapply(Gn, function(g) {
  strsplit(g[[7]]$gene_id, "\\.")[[1]][1]})
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",
                      GRCh=37)
# View(listFilters(ensembl))
gb <- getBM(attributes=c("ensembl_gene_id", "ensembl_exon_id", 
                         "exon_chrom_start", "exon_chrom_end"), 
            filters=c("ensembl_gene_id"), 
            values=geneid, 
            mart=ensembl)
exon <- sapply(1:G, function(g) {as.numeric(sapply(Gn[[g]][[5]], function(s) {
  any(apply(gb[gb$ensembl_gene_id %in% geneid[[g]], c(3, 4)], 1, function(e) {
    (s >= e[1]) & (s <= e[2])}))}))})

### transforming co-data
bc.genelength <- boxcox(lm(genelength ~ 1, 
                           data=data.frame(genelength=unlist(genelength))))
tgenelength <- (unlist(genelength)^
                  bc.genelength$x[which.max(bc.genelength$y)] - 1)/
  bc.genelength$x[which.max(bc.genelength$y)]
tgenelength <- unname(split((tgenelength - min(tgenelength))/
                              max(tgenelength - min(tgenelength)), 
                            rep(1:length(p), times=p)))
bc.minoraf <- boxcox(lm(minoraf ~ 1, data=data.frame(minoraf=unlist(minoraf))))
tminoraf <- (unlist(minoraf)^
               bc.minoraf$x[which.max(bc.minoraf$y)] - 1)/
  bc.minoraf$x[which.max(bc.minoraf$y)]
tminoraf <- unname(split((tminoraf - min(tminoraf))/
                           max(tminoraf - min(tminoraf)), 
                         rep(1:length(p), times=p)))

### creating co-data matrix
C <- sapply(1:G, function(g) {
  cbind(inorout=inorout[[g]], 
        exon=exon[[g]],
        genelength=tgenelength[[g]], 
        minoraf=tminoraf[[g]])})

### fitting models
control <- list(conv.post=TRUE, trace=TRUE,
                epsilon.eb=1e-3, epsilon.vb=1e-3, 
                epsilon.opt=sqrt(.Machine$double.eps),
                maxit.eb=50, maxit.vb=2, maxit.opt=100,
                maxit.post=100)

set.seed(2019)
fit.enig <- enig(xtrain, ytrain, C, mult.lambda=FALSE, intercept.eb=TRUE, 
                 fixed.eb="none", full.post=TRUE, init=NULL, control=control)
fit.mult.enig <- enig(xtrain, ytrain, C, mult.lambda=TRUE, intercept.eb=TRUE, 
                      fixed.eb="none", full.post=TRUE, init=NULL, 
                      control=control)
fit.int.enig <- enig(xtrain, ytrain, NULL, mult.lambda=FALSE, intercept.eb=TRUE,
                     fixed.eb="none", full.post=TRUE, init=NULL, 
                     control=control)
fit.mult.int.enig <- enig(xtrain, ytrain, NULL, mult.lambda=TRUE, 
                          intercept.eb=TRUE, fixed.eb="none", full.post=TRUE, 
                          init=NULL, control=control)
fit.glmnet <- lapply(1:G, function(g) {
  cv.glmnet(xtrain[[g]], ytrain[, g], intercept=FALSE, alpha=0)})

save(fit.enig, fit.mult.enig, fit.int.enig, fit.mult.int.enig, fit.glmnet, 
     file="results/eQTL_fit1.Rdata")
# load("results/eQTL_fit1.Rdata")

### point estimates
best.enig <- fit.enig$vb$mu
best.mult.enig <- fit.mult.enig$vb$mu
best.int.enig <- fit.int.enig$vb$mu
best.mult.int.enig <- fit.mult.int.enig$vb$mu
best.glmnet <- sapply(fit.glmnet, function(g) {
  as.numeric(coef(g, s="lambda.min"))[-1]})

pairs(~ enig + mult.enig + int.enig + 
        mult.int.enig + glmnet, 
      data=data.frame(enig=unlist(best.enig), 
                      mult.enig=unlist(best.mult.enig),
                      int.enig=unlist(best.int.enig), 
                      mult.int.enig=unlist(best.mult.int.enig),
                      glmnet=unlist(best.glmnet)))

# heritabilities
hsq.enig <- sapply(1:G, function(g) {
  sum(fit.enig$vb$mpost$gamma.sq[[g]])/
    (sum(fit.enig$vb$mpost$gamma.sq[[g]]) + 1) - 
    sum(fit.enig$vb$vpost$gamma.sq[[g]])/
    (sum(fit.enig$vb$mpost$gamma.sq[[g]]) + 1)^3})
hsq.mult.enig <- sapply(1:G, function(g) {
  sum(fit.mult.enig$vb$mpost$gamma.sq[[g]])/
    (sum(fit.mult.enig$vb$mpost$gamma.sq[[g]]) + 1) - 
    sum(fit.mult.enig$vb$vpost$gamma.sq[[g]])/
    (sum(fit.mult.enig$vb$mpost$gamma.sq[[g]]) + 1)^3})
hsq.int.enig <- sapply(1:G, function(g) {
  sum(fit.int.enig$vb$mpost$gamma.sq[[g]])/
    (sum(fit.int.enig$vb$mpost$gamma.sq[[g]]) + 1) - 
    sum(fit.int.enig$vb$vpost$gamma.sq[[g]])/
    (sum(fit.int.enig$vb$mpost$gamma.sq[[g]]) + 1)^3})
hsq.mult.int.enig <- sapply(1:G, function(g) {
  sum(fit.mult.int.enig$vb$mpost$gamma.sq[[g]])/
    (sum(fit.mult.int.enig$vb$mpost$gamma.sq[[g]]) + 1) - 
    sum(fit.mult.int.enig$vb$vpost$gamma.sq[[g]])/
    (sum(fit.mult.int.enig$vb$mpost$gamma.sq[[g]]) + 1)^3})

hsq.enig.sel <- hsq.mult.enig.sel <- hsq.int.enig.sel <- 
  hsq.mult.int.enig.sel <- lapply(1:G, function(g) {numeric(p[g])})
for(g in 1:G) {
  for(k in 1:p[g]) {
    id <- which(order(fit.enig$vb$mpost$gamma.sq[[g]]) %in% c(1:k))  
    hsq.enig.sel[[g]][k] <- sum(fit.enig$vb$mpost$gamma.sq[[g]][id])/
      (sum(fit.enig$vb$mpost$gamma.sq[[g]][id]) + 1) - 
      sum(fit.enig$vb$vpost$gamma.sq[[g]][id])/
      (sum(fit.enig$vb$mpost$gamma.sq[[g]][id]) + 1)^3
    
    id <- which(order(fit.mult.enig$vb$mpost$gamma.sq[[g]]) %in% c(1:k))  
    hsq.mult.enig.sel[[g]][k] <- sum(fit.mult.enig$vb$mpost$gamma.sq[[g]][id])/
      (sum(fit.mult.enig$vb$mpost$gamma.sq[[g]][id]) + 1) - 
      sum(fit.mult.enig$vb$vpost$gamma.sq[[g]][id])/
      (sum(fit.mult.enig$vb$mpost$gamma.sq[[g]][id]) + 1)^3
    
    id <- which(order(fit.int.enig$vb$mpost$gamma.sq[[g]]) %in% c(1:k))  
    hsq.int.enig.sel[[g]][k] <- sum(fit.int.enig$vb$mpost$gamma.sq[[g]][id])/
      (sum(fit.int.enig$vb$mpost$gamma.sq[[g]][id]) + 1) - 
      sum(fit.int.enig$vb$vpost$gamma.sq[[g]][id])/
      (sum(fit.int.enig$vb$mpost$gamma.sq[[g]][id]) + 1)^3
    
    id <- which(order(fit.mult.int.enig$vb$mpost$gamma.sq[[g]]) %in% c(1:k))  
    hsq.mult.int.enig.sel[[g]][k] <- 
      sum(fit.mult.int.enig$vb$mpost$gamma.sq[[g]][id])/
      (sum(fit.mult.int.enig$vb$mpost$gamma.sq[[g]][id]) + 1) - 
      sum(fit.mult.int.enig$vb$vpost$gamma.sq[[g]][id])/
      (sum(fit.mult.int.enig$vb$mpost$gamma.sq[[g]][id]) + 1)^3
  }
  hsq.enig.sel[[g]] <- c(0, hsq.enig.sel[[g]])
  hsq.mult.enig.sel[[g]] <- c(0, hsq.mult.enig.sel[[g]])
  hsq.int.enig.sel[[g]] <- c(0, hsq.int.enig.sel[[g]])
  hsq.mult.int.enig.sel[[g]] <- c(0, hsq.mult.int.enig.sel[[g]])
}

g <- 6
col <- 1:2
lty <- 1:2
plot(c(0:p[g])/p[g], hsq.enig.sel[[g]]/hsq.enig[[g]], type="l", ylim=c(0, 1), 
     xlab="Fraction of selected SNPs", ylab="Explained heritability", 
     col=col[1], lty=lty[1])
lines(c(0:p[g])/p[g], hsq.mult.enig.sel[[g]]/hsq.mult.enig[[g]], col=col[1],
      lty=lty[2])
lines(c(0:p[g])/p[g], hsq.int.enig.sel[[g]]/hsq.int.enig[[g]], col=col[2], 
      lty=lty[1])
lines(c(0:p[g])/p[g], hsq.mult.int.enig.sel[[g]]/hsq.mult.int.enig[[g]], 
      col=col[2], lty=lty[2])
legend("bottomright", lty=1, col=col, 
       legend=c("Standard", "Intercept only"))

hsq.mult.enig.sel[[1]]
sapply(hsq.mult.enig.sel, function(s) {any(s < 0)})
sapply(hsq.mult.trans.enig.sel, function(s) {any(s < 0)})
sapply(hsq.mult.int.enig.sel, function(s) {any(s < 0)})
hsq.mult.enig.sel[[6]]
hsq.mult.int.enig.sel[[6]]
hsq.mult.trans.enig.sel[[6]]

library(pracma)
# exclude regression where we have negative variance estimates
auc <- t(sapply(c(1:5, 7:G), function(g) {
  c(enig=trapz(c(0:p[g])/p[g], hsq.enig.sel[[g]]/hsq.enig[[g]]),
    mult.enig=trapz(c(0:p[g])/p[g], hsq.mult.enig.sel[[g]]/hsq.mult.enig[[g]]),
    trans.enig=trapz(c(0:p[g])/p[g], hsq.trans.enig.sel[[g]]/
                       hsq.trans.enig[[g]]),
    mult.trans.enig=trapz(c(0:p[g])/p[g], hsq.mult.trans.enig.sel[[g]]/
                            hsq.mult.trans.enig[[g]]),
    int.enig=trapz(c(0:p[g])/p[g], 
                   hsq.int.enig.sel[[g]]/hsq.int.enig[[g]]),
    mult.int.enig=trapz(c(0:p[g])/p[g], hsq.mult.int.enig.sel[[g]]/
                          hsq.mult.int.enig[[g]]))}))



colMeans(auc)
boxplot(auc)

### prior estimates
# standardize alpha estimates
alpha <- cbind(enig=fit.enig$eb$alpha*c(1, apply(Reduce("rbind", C), 2, sd)),
               mult.enig=fit.mult.enig$eb$alpha*
                 c(1, apply(Reduce("rbind", C), 2, sd)),
               trans.enig=fit.trans.enig$eb$alpha*
                 c(1, apply(Reduce("rbind", tC), 2, sd)),
               mult.trans.enig=fit.mult.trans.enig$eb$alpha*
                 c(1, apply(Reduce("rbind", tC), 2, sd)),
               int.enig=c(fit.int.enig$eb$alpha, 
                          setNames(rep(0, 3), colnames(C[[1]]))),
               mult.int.enig=c(fit.mult.int.enig$eb$alpha, 
                          setNames(rep(0, 3), colnames(C[[1]]))))

lambda <- cbind(enig=rep(fit.enig$eb$lambda, G),
                mult.enig=fit.mult.enig$eb$lambda,
                trans.enig=rep(fit.trans.enig$eb$lambda, G),
                mult.trans.enig=fit.mult.trans.enig$eb$lambda,
                int.enig=rep(fit.int.enig$eb$lambda, G),
                mult.int.enig=fit.mult.int.enig$eb$lambda)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.1, 1, 1))
layout(matrix(c(1:3), nrow=1, ncol=3, byrow=TRUE))
hist(lambda[, 2])
abline(v=lambda[1, 1], col=2)

hist(lambda[, 4])
abline(v=lambda[1, 3], col=2)

hist(lambda[, 6])
abline(v=lambda[1, 5], col=2)
par(opar)

### error variance
plot(unlist(fit.mult.int.enig$vb$mpost$sigma.sq), 
     unlist(fit.mult.trans.enig$vb$mpost$sigma.sq))
abline(a=0, b=1, col=2)

plot(unlist(fit.mult.int.enig$eb$mprior), unlist(fit.mult.trans.enig$eb$mprior))
plot(unlist(fit.mult.int.enig$vb$mpost$beta), 
     unlist(fit.mult.trans.enig$vb$mpost$beta))
abline(a=0, b=1, col=2)

# prior gamma mean versus co-data
plot(unlist(tgenelength), unlist(fit.mult.trans.enig$eb$mprior))
boxplot(unlist(fit.mult.trans.enig$eb$mprior) ~ unlist(inorout))
plot(unlist(tminoraf), unlist(fit.mult.trans.enig$eb$mprior))

plot(unlist(tgenelength), unlist(fit.mult.trans.enig$vb$mpost$beta))
abline(h=0, col=2, lty=2)
boxplot(unlist(fit.mult.trans.enig$vb$mpost$beta) ~ unlist(inorout))
abline(h=0, col=2, lty=2)
plot(unlist(tminoraf), unlist(fit.mult.trans.enig$vb$mpost$beta))
abline(h=0, col=2, lty=2)

str(fit.mult.trans.enig$vb$mpost)
length(unlist(fit.mult.trans.enig$vb$mpost$beta))
length(unlist(fit.mult.trans.enig$eb$mprior))

### convergence
col <- 1:ncol(fit.enig$seq.eb$alpha)
plot(fit.enig$seq.eb$alpha[, 1], type="l", ylim=range(fit.enig$seq.eb$alpha),
     col=col[1])
lines(fit.enig$seq.eb$alpha[, 2], col=col[2])
lines(fit.enig$seq.eb$alpha[, 3], col=col[3])
lines(fit.enig$seq.eb$alpha[, 4], col=col[4])
legend("topleft", col=col, legend=names(fit.enig$eb$alpha), lty=1)

col <- 1:ncol(fit.enig$seq.eb$lambda)
plot(fit.enig$seq.eb$lambda[, 1], type="l", ylim=range(fit.enig$seq.eb$lambda),
     col=col[1])
for(i in 2:ncol(fit.enig$seq.eb$lambda)) {
  lines(fit.enig$seq.eb$lambda[, i], col=col[i])  
}

col <- 1:ncol(fit.mult.enig$seq.eb$alpha)
plot(fit.mult.enig$seq.eb$alpha[, 1], type="l", 
     ylim=range(fit.mult.enig$seq.eb$alpha), col=col[1])
lines(fit.mult.enig$seq.eb$alpha[, 2], col=col[2])
lines(fit.mult.enig$seq.eb$alpha[, 3], col=col[3])
lines(fit.mult.enig$seq.eb$alpha[, 4], col=col[4])
legend("topleft", col=col, legend=names(fit.mult.enig$eb$alpha), lty=1)

col <- 1:ncol(fit.enig$seq.eb$lambda)
plot(fit.mult.enig$seq.eb$lambda[, 1], type="l", 
     ylim=range(fit.mult.enig$seq.eb$lambda), col=col[1])
for(i in 2:ncol(fit.mult.enig$seq.eb$lambda)) {
  lines(fit.mult.enig$seq.eb$lambda[, i], col=col[i])  
}


### predictions
pred.enig <- sapply(1:G, function(g) {
  as.numeric(xtest[[g]] %*% fit.enig$vb$mpost$beta[[g]])})
pred.mult.enig <- sapply(1:G, function(g) {
  as.numeric(xtest[[g]] %*% fit.mult.enig$vb$mpost$beta[[g]])})
pred.int.enig <- sapply(1:G, function(g) {
  as.numeric(xtest[[g]] %*% fit.int.enig$vb$mpost$beta[[g]])})
pred.mult.int.enig <- sapply(1:G, function(g) {
  as.numeric(xtest[[g]] %*% fit.mult.int.enig$vb$mpost$beta[[g]])})
pred.glmnet <- sapply(1:G, function(g) {
  predict(fit.glmnet[[g]], xtest[[g]], s="lambda.min")})
pred <- list(enig=pred.enig, mult.enig=pred.mult.enig, 
             int.enig=pred.int.enig, mult.int.enig=pred.mult.int.enig,
             glmnet=pred.glmnet)

corpred <- sapply(pred, function(s) {cor(as.numeric(s), as.numeric(ytest))})
brier <- sapply(pred, function(s) {mean((s - ytest)^2)})

plot(ytest, pred$enig)
abline(a=0, b=1)
plot(ytest, pred$glmnet)
length(pred$mult.trans.enig)

### CPO
f.int.cpo <- function(val, y, zeta, n, p, xtmu, xtSigmax) {
  (val^2/(2*zeta) + 1)^(-(n - length(xtmu) + p + 2)/2)*
    exp(-(val - xtmu + y)^2/(2*xtSigmax))
}

cpoi <- function(model, x, y, n, p, D) {
  out <- lapply(1:D, function(d) {
    zeta <- model$vb$zeta[d]
    xtmu <- as.numeric(x[[d]] %*% model$vb$mu[[d]])
    xtSigmax <- rowSums((x[[d]] %*% model$vb$Sigma[[d]])*x[[d]])
    int <- sapply(1:nrow(matrix(y, ncol=D)), function(i) {
      unlist(integrate(f.int.cpo, -Inf, Inf, y=matrix(y, ncol=D)[i, d], 
                       zeta=zeta, n=n, p=p[d], xtmu=xtmu[i], 
                       xtSigmax=xtSigmax[i])[c(1:2)])})
    
    
    list((2*pi)^(-1)*exp(lgamma((n- nrow(matrix(y, ncol=D)) + p[d] + 2)/2) - 
                            lgamma((n - nrow(matrix(y, ncol=D)) + p[d] + 1)/2))*
            (zeta*xtSigmax)^(-1/2)*int[1, ], int[2, ])})

  return(list(value=unname(sapply(out, "[[", 1)), 
              abs.error=unname(sapply(out, "[[", 2))))
}

lpml <- function(x, y, C, mult.lambda=FALSE, intercept.eb=TRUE,
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
    fit <- enig(x=xtrain, y=ytrain, C=C, mult.lambda=mult.lambda, 
                intercept.eb=intercept.eb, fixed.eb=fixed.eb, 
                full.post=full.post, init=init, control=control)
    
    # estimate CPOi on the test data
    est <- cpoi(fit, xtest, ytest, n, p, D)
    val[foldid==k, ] <- est$value
    err[foldid==k, ] <- est$abs.error
  
  }
  
  out <- list(lpml=colSums(log(val))/n, cpoi=val, abs.error=err)
  return(out)
  
}

set.seed(2019)
nfolds <- 5
foldid <- sample(rep(1:nfolds, times=round(c(rep(
  n %/% nfolds + as.numeric((n %% nfolds)!=0), times=n %% nfolds),
  rep(n %/% nfolds, times=nfolds - n %% nfolds)))))
lpml.enig <- lpml(x, y, C, mult.lambda=FALSE, 
                  intercept.eb=TRUE, fixed.eb="none", full.post=TRUE, 
                  init=NULL, control=control, foldid=foldid)
lmpl.mult.enig <- lpml(x, y, C, mult.lambda=TRUE, 
                       intercept.eb=TRUE, fixed.eb="none", full.post=TRUE, 
                       init=NULL, control=control, foldid=foldid)
lpml.int.enig <- lpml(x, y, NULL, mult.lambda=FALSE, 
                      intercept.eb=TRUE, fixed.eb="none", full.post=TRUE, 
                      init=NULL, control=control, foldid=foldid)
lpml.mult.int.enig <- lpml(x, y, NULL, mult.lambda=TRUE, 
                           intercept.eb=TRUE, fixed.eb="none", full.post=TRUE, 
                           init=NULL, control=control, foldid=foldid)
save(lpml.enig, lmpl.mult.enig, lpml.int.enig, lpml.mult.int.enig, 
     file="eQTL_lpml1.Rdata")

addline <- function(x, y, ...){
  points(x, y, ...)
  abline(a=0, b=1, ...)
}
pairs(cbind(enig=lpml.enig$lpml, mult=lmpl.mult.enig$lpml, 
            trans=lpml.trans.enig$lpml, mult.trans=lpml.mult.trans.enig$lpml, 
            int=lpml.int.enig$lpml, mult.int=lpml.mult.int.enig$lpml),
      lower.panel=addline, upper.panel=addline)



plot(lpml.int.enig, lpml.trans.enig)
plot(lpml.mult.int.enig, lpml.mult.trans.enig)
plot(lpml.int.enig, lpml.mult.int.enig)
abline(a=0, b=1, col=2)
