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
library(rstan)
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
ntrain <- floor(n/5)
id.train <- sample(1:n, ntrain)
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
                         "exon_chrom_start", "exon_chrom_end",
                         "chromosome_name"), 
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

### fitting models
control <- list(conv.post=TRUE, trace=TRUE,
                epsilon.eb=1e-3, epsilon.vb=1e-3, 
                epsilon.opt=sqrt(.Machine$double.eps),
                maxit.eb=100, maxit.vb=2, maxit.opt=100,
                maxit.post=100)

C1 <- sapply(1:G, function(g) {
  cbind(inorout=inorout[[g]],
        exon=exon[[g]],
        genelength=tgenelength[[g]], 
        minoraf=tminoraf[[g]])})
C2 <- sapply(1:G, function(g) {
  cbind(inorout=inorout[[g]],
        genelength=tgenelength[[g]], 
        minoraf=tminoraf[[g]])})
C3 <- sapply(1:G, function(g) {
  cbind(exon=exon[[g]],
        genelength=tgenelength[[g]], 
        minoraf=tminoraf[[g]])})
C4 <- sapply(1:G, function(g) {
  cbind(exon=exon[[g]],
        genelength=tgenelength[[g]], 
        minoraf=tminoraf[[g]],
        minorafsq=tminoraf[[g]]^2)})

set.seed(2019)
fit1.enig <- enig(xtrain, ytrain, C1, mult.lambda=TRUE, intercept.eb=TRUE, 
                  fixed.eb="none", full.post=TRUE, init=NULL, 
                  control=control)
fit2.enig <- enig(xtrain, ytrain, C2, mult.lambda=TRUE, intercept.eb=TRUE, 
                  fixed.eb="none", full.post=TRUE, init=NULL, 
                  control=control)
fit3.enig <- enig(xtrain, ytrain, C3, mult.lambda=TRUE, intercept.eb=TRUE, 
                  fixed.eb="none", full.post=TRUE, init=NULL, 
                  control=control)
fit4.enig <- enig(xtrain, ytrain, C4, mult.lambda=TRUE, intercept.eb=TRUE, 
                  fixed.eb="none", full.post=TRUE, init=NULL, 
                  control=control)
fit.ienig <- enig(xtrain, ytrain, NULL, mult.lambda=TRUE, intercept.eb=TRUE,
                     fixed.eb="none", full.post=TRUE, init=NULL, 
                     control=control)

### estimate LPML
nfolds <- 5
set.seed(2019)
foldid <- sample(rep(1:nfolds, times=round(c(rep(
  n %/% nfolds + as.numeric((n %% nfolds)!=0), times=n %% nfolds),
  rep(n %/% nfolds, times=nfolds - n %% nfolds)))))
lpml1.enig <- lpml(x, y, C1, mult.lambda=TRUE, intercept.eb=TRUE, 
                   fixed.eb="none", full.post=TRUE, init=NULL, control=control, 
                   foldid=foldid)
lpml2.enig <- lpml(x, y, C2, mult.lambda=TRUE, intercept.eb=TRUE, 
                   fixed.eb="none", full.post=TRUE, init=NULL, control=control, 
                   foldid=foldid)
lpml3.enig <- lpml(x, y, C3, mult.lambda=TRUE, intercept.eb=TRUE, 
                   fixed.eb="none", full.post=TRUE, init=NULL, control=control, 
                   foldid=foldid)
lpml4.enig <- lpml(x, y, C4, mult.lambda=TRUE, intercept.eb=TRUE, 
                   fixed.eb="none", full.post=TRUE, init=NULL, control=control, 
                   foldid=foldid)
lpml.ienig <- lpml(x, y, NULL, mult.lambda=TRUE, intercept.eb=TRUE, 
                   fixed.eb="none", full.post=TRUE, init=NULL, control=control, 
                   foldid=foldid)
save(fit1.enig, fit2.enig, fit3.enig, fit4.enig, fit.ienig, lpml1.enig,
     lpml2.enig, lpml3.enig, lpml4.enig, lpml.ienig,
     file="results/eQTL_fit1.Rdata")

### sampling from posteriors
stan.enig <- stan_model("code/enig.stan")
nsamples <- 1000
nwarmup <- 1000
mcmc1.enig <- sapply(1:G, function(g) {
  sampling(stan.enig, chains=1, warmup=nwarmup, iter=nsamples + nwarmup, 
           cores=1, refresh=0, control=list(adapt_delta=0.8, max_treedepth=12),
           data=list(n=ntrain, p=p[g], y=ytrain[, g], x=xtrain[[g]], 
                     ctalphainv=fit1.enig$eb$mprior[[g]], 
                     lambda=fit1.enig$eb$lambda[g]))})
mcmc2.enig <- sapply(1:G, function(g) {
  sampling(stan.enig, chains=1, warmup=nwarmup, iter=nsamples + nwarmup, 
           cores=1, refresh=0, control=list(adapt_delta=0.8, max_treedepth=12),
           data=list(n=ntrain, p=p[g], y=ytrain[, g], x=xtrain[[g]], 
                     ctalphainv=fit2.enig$eb$mprior[[g]], 
                     lambda=fit2.enig$eb$lambda[g]))})
mcmc3.enig <- sapply(1:G, function(g) {
  sampling(stan.enig, chains=1, warmup=nwarmup, iter=nsamples + nwarmup, 
           cores=1, refresh=0, control=list(adapt_delta=0.8, max_treedepth=12),
           data=list(n=ntrain, p=p[g], y=ytrain[, g], x=xtrain[[g]], 
                     ctalphainv=fit3.enig$eb$mprior[[g]], 
                     lambda=fit3.enig$eb$lambda[g]))})
mcmc4.enig <- sapply(1:G, function(g) {
  sampling(stan.enig, chains=1, warmup=nwarmup, iter=nsamples + nwarmup, 
           cores=1, refresh=0, control=list(adapt_delta=0.8, max_treedepth=12),
           data=list(n=ntrain, p=p[g], y=ytrain[, g], x=xtrain[[g]], 
                     ctalphainv=fit4.enig$eb$mprior[[g]], 
                     lambda=fit4.enig$eb$lambda[g]))})
mcmc.ienig <- sapply(1:G, function(g) {
  sampling(stan.enig, chains=1, warmup=nwarmup, iter=nsamples + nwarmup, 
           cores=1, refresh=0, control=list(adapt_delta=0.8, max_treedepth=12),
           data=list(n=ntrain, p=p[g], y=ytrain[, g], x=xtrain[[g]], 
                     ctalphainv=fit.ienig$eb$mprior[[g]], 
                     lambda=fit.ienig$eb$lambda[g]))})
# save(mcmc1.enig, mcmc2.enig, mcmc3.enig, mcmc4.enig, mcmc.ienig,
#      file="results/eQTL_mcmc1.Rdata")
load("results/eQTL_mcmc1.Rdata")

psel <- c(1, 5, 10, 20, 50)
hsq.ienig <- sapply(1:length(psel), function(csel) {
  sapply(1:G, function(g) {
    s <- extract(mcmc.ienig[[g]])
    id <- which(order(abs(colMeans(s$beta)), decreasing=TRUE) %in% c(1:csel))
    gsq <- s$gammasq[, id, drop=FALSE]
    mean(rowSums(gsq)/(rowSums(gsq) + 1))})})
hsq1.enig <- sapply(1:length(psel), function(csel) {
  sapply(1:G, function(g) {
    s <- extract(mcmc1.enig[[g]])
    id <- which(order(abs(colMeans(s$beta)), decreasing=TRUE) %in% c(1:csel))
    gsq <- s$gammasq[, id, drop=FALSE]
    mean(rowSums(gsq)/(rowSums(gsq) + 1))})})
hsq2.enig <- sapply(1:length(psel), function(csel) {
  sapply(1:G, function(g) {
    s <- extract(mcmc2.enig[[g]])
    id <- which(order(abs(colMeans(s$beta)), decreasing=TRUE) %in% c(1:csel))
    gsq <- s$gammasq[, id, drop=FALSE]
    mean(rowSums(gsq)/(rowSums(gsq) + 1))})})
hsq3.enig <- sapply(1:length(psel), function(csel) {
  sapply(1:G, function(g) {
    s <- extract(mcmc3.enig[[g]])
    id <- which(order(abs(colMeans(s$beta)), decreasing=TRUE) %in% c(1:csel))
    gsq <- s$gammasq[, id, drop=FALSE]
    mean(rowSums(gsq)/(rowSums(gsq) + 1))})})
hsq4.enig <- sapply(1:length(psel), function(csel) {
  sapply(1:G, function(g) {
    s <- extract(mcmc4.enig[[g]])
    id <- which(order(abs(colMeans(s$beta)), decreasing=TRUE) %in% c(1:csel))
    gsq <- s$gammasq[, id, drop=FALSE]
    mean(rowSums(gsq)/(rowSums(gsq) + 1))})})
hsq <- data.frame(intercept=hsq.ienig, model1=hsq1.enig, model2=hsq2.enig,
                  model3=hsq3.enig, model4=hsq4.enig)
write.table(hsq, file="results/eQTL_hsq.csv")

best1.enig <- sapply(mcmc1.enig, function(m) {
  get_posterior_mean(m, pars=c("beta0", "beta"))[, 1]})
best2.enig <- sapply(mcmc2.enig, function(m) {
  get_posterior_mean(m, pars=c("beta0", "beta"))[, 1]})
best3.enig <- sapply(mcmc3.enig, function(m) {
  get_posterior_mean(m, pars=c("beta0", "beta"))[, 1]})
best4.enig <- sapply(mcmc4.enig, function(m) {
  get_posterior_mean(m, pars=c("beta0", "beta"))[, 1]})
best.ienig <- sapply(mcmc.ienig, function(m) {
  get_posterior_mean(m, pars=c("beta0", "beta"))[, 1]})

bvar1.enig <- sapply(mcmc1.enig, function(m) {
  s <- summary(m)$summary; s[substr(rownames(s), 1, 4)=="beta", 3]^2})
bvar2.enig <- sapply(mcmc2.enig, function(m) {
  s <- summary(m)$summary; s[substr(rownames(s), 1, 4)=="beta", 3]^2})
bvar3.enig <- sapply(mcmc3.enig, function(m) {
  s <- summary(m)$summary; s[substr(rownames(s), 1, 4)=="beta", 3]^2})
bvar4.enig <- sapply(mcmc4.enig, function(m) {
  s <- summary(m)$summary; s[substr(rownames(s), 1, 4)=="beta", 3]^2})
bvar.ienig <- sapply(mcmc.ienig, function(m) {
  s <- summary(m)$summary; s[substr(rownames(s), 1, 4)=="beta", 3]^2})
# save(best1.enig, best2.enig, best3.enig, best4.enig, best.ienig, bvar1.enig, 
#      bvar2.enig, bvar3.enig, bvar4.enig, bvar.ienig,
#      file="results/eQTL_best1.Rdata")
load("results/eQTL_best1.Rdata")

mse1.enig <- sapply(1:G, function(g) {
  mean((ytest[, g] - as.numeric(cbind(1, xtest[[g]]) %*% best1.enig[[g]]))^2)})
mse2.enig <- sapply(1:G, function(g) {
  mean((ytest[, g] - as.numeric(cbind(1, xtest[[g]]) %*% best2.enig[[g]]))^2)})
mse3.enig <- sapply(1:G, function(g) {
  mean((ytest[, g] - as.numeric(cbind(1, xtest[[g]]) %*% best3.enig[[g]]))^2)})
mse4.enig <- sapply(1:G, function(g) {
  mean((ytest[, g] - as.numeric(cbind(1, xtest[[g]]) %*% best4.enig[[g]]))^2)})
mse.ienig <- sapply(1:G, function(g) {
  mean((ytest[, g] - as.numeric(cbind(1, xtest[[g]]) %*% best.ienig[[g]]))^2)})

mbvar1.enig <- sapply(bvar1.enig, mean)
mbvar2.enig <- sapply(bvar2.enig, mean)
mbvar3.enig <- sapply(bvar3.enig, mean)
mbvar4.enig <- sapply(bvar4.enig, mean)
mbvar.ienig <- sapply(bvar.ienig, mean)

### LPML plot
labels <- c("in gene + \n exon + \n genelength + \n min. allele freq.", 
            "in gene + \n genelength + \n min. allele freq.", 
            "in exon + \n genelength + \n min. allele freq.", 
            "in exon + \n genelength + \n min. allele freq. squared", 
            "intercept only")
pairs(cbind(enig1=lpml1.enig$lpml, enig2=lpml2.enig$lpml, enig3=lpml3.enig$lpml, 
            enig4=lpml4.enig$lpml, ienig=lpml.ienig$lpml),
      lower.panel=function(x, y, ...) {
        points(x, y, ...); abline(a=0, b=1, ...)}, 
      upper.panel=function(x, y, ...) {
        points(x, y, ...); abline(a=0, b=1, ...)}, 
      diag.panel=function(x, ...) {
        lx <- sum(range(x))/2; ly <- min(x) + diff(range(x))/4
        text(lx, ly, paste0("total LPML \n", round(sum(x), 4), ...))},
      label.pos=2/3, labels=labels, main="LPML's for the different equations")

pairs(cbind(enig1=mse1.enig, enig2=mse2.enig, enig3=mse3.enig, enig4=mse4.enig, 
            ienig=mse.ienig),
      lower.panel=function(x, y, ...) {
        points(x, y, ...); abline(a=0, b=1, ...)}, 
      upper.panel=function(x, y, ...) {
        points(x, y, ...); abline(a=0, b=1, ...)}, 
      diag.panel=function(x, ...) {
        lx <- sum(range(x))/2; ly <- min(x) + diff(range(x))/4
        text(lx, ly, paste0("average MSE \n", round(mean(x), 4), ...))},
      label.pos=2/3, labels=labels, main="MSE's for the different equations")

pairs(cbind(enig1=mbvar1.enig, enig2=mbvar2.enig, enig3=mbvar3.enig, 
            enig4=mbvar4.enig, ienig=mbvar.ienig),
      lower.panel=function(x, y, ...) {
        points(x, y, ...); abline(a=0, b=1, ...)},
      upper.panel=function(x, y, ...) {
        points(x, y, ...); abline(a=0, b=1, ...)},
      diag.panel=function(x, ...) {
        lx <- sum(range(x))/2; ly <- min(x) + diff(range(x))/4
        text(lx, ly, paste0("average MSE \n", round(mean(x), 4), ...))},
      label.pos=2/3, labels=labels, 
      main="average posterior variance for the different equations")



### point estimates


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
alpha <- cbind(intercept=fit.ienig$eb$alpha,
               model1=fit1.mult.enig$eb$alpha*
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
# lines(fit.enig$seq.eb$alpha[, 5], col=col[5])
legend("topleft", col=col, legend=names(fit.enig$eb$alpha), lty=1)

col <- 1:ncol(fit.mult.enig$seq.eb$alpha)
plot(fit.mult.enig$seq.eb$alpha[, 1], type="l", 
     ylim=range(fit.mult.enig$seq.eb$alpha), col=col[1])
lines(fit.mult.enig$seq.eb$alpha[, 2], col=col[2])
lines(fit.mult.enig$seq.eb$alpha[, 3], col=col[3])
lines(fit.mult.enig$seq.eb$alpha[, 4], col=col[4])
lines(fit.mult.enig$seq.eb$alpha[, 5], col=col[5])
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





