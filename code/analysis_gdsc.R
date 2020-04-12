#!/usr/bin/env Rscript

################################################################################
# regression GDSC log IC50 values on gene expression/mutations of cell lines
################################################################################

### libraries
packages <- c("foreach", "doParallel", "rstan", "glmnet", "pInc", "cambridge", 
              "xtune")
sapply(packages, library, character.only=TRUE)

# number of CV folds
nfolds <- 10

# number of cores to use
# ncores <- min(nfolds, 100)
ncores <- min(100, nfolds)

### load data
load(file="data/data_gdsc_dat1.Rdata")

### functions
coef.dss <- function(object, psel) {
  best <- matrix(nrow=object$glmnet.fit$dim[1] + 1, ncol=length(psel))
  for(p in 1:length(psel)) {
    ind <- which.min(abs(colSums(coef(object, s=object$lambda)[-1, ]!=0) - 
                           psel[p]))
    best[, p] <- coef(object, s=object$lambda[ind])[, 1]
  }
  return(best)
}

################################################################################
################################## analysis 1 ##################################
################################################################################
# expression data with pvalues from mutationas as external data

### external data analysis
# calculate pvalues for mutation data
mut.pvalues <- lapply(resp, function(s) {
  resp.prep <- s[names(s) %in% rownames(mut)]
  mut.prep <- mut[rownames(mut) %in% names(s), ]
  pvalues <- apply(mut.prep, 2, function(m) {
    if(sum(m) > 1) {
      t.test(resp.prep[m==1], resp.prep[m==0])$p.value
    } else {
      NA
    }})
  pvalues[!is.na(pvalues)]})

### data preparation
# select data
# resp.sel <- resp[names(resp) %in% drug$name]
resp.sel <- sapply(resp, function(s) {
  s <- s[names(s) %in% rownames(expr$expr)]
  s <- s[order(names(s))]})
idcellsel <- lapply(resp.sel, function(s) {
  match(names(s), rownames(expr$expr))})
match.mut <- lapply(1:length(mut.pvalues), function(d) {
  mutnames <- strsplit(names(mut.pvalues[[d]]), ",")
  idsel <- which(colnames(expr$expr) %in% unlist(mutnames))
  pvalues <- sapply(colnames(expr$expr)[idsel], function(s) {
    ids <- which(sapply(mutnames, function(l) {any(s==l)}))
    length(ids)/sum(1/mut.pvalues[[d]][ids])})
  return(list(pvalues=pvalues, idsel=idsel))
  })
idfeatsel <- lapply(match.mut, "[[", "idsel")
pvalues <- lapply(match.mut, "[[", "pvalues")

D <- length(resp.sel)
n <- sapply(resp.sel, length)

x <- lapply(1:D, function(d) {scale(expr$expr[idcellsel[[d]], idfeatsel[[d]]])})
y <- lapply(1:D, function(d) {scale(resp.sel[[d]])[, 1]})

# create co-data
C <- lapply(1:D, function(d) {
  s <- cbind(1, pvalues[[d]]); colnames(s) <- c("intercept", "pvalue"); s})

### model fitting
# settings
control.semnig <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, 
                       epsilon.vb=1e-3, maxit.eb=500, maxit.vb=1, 
                       maxit.post=100, maxit.block=0)
control.ebridge <-list(epsilon=sqrt(.Machine$double.eps), maxit=500, 
                       trace=TRUE, glmnet.fit2=FALSE, beta2=TRUE)

# creating folds
set.seed(2020)
foldid <- lapply(1:D, function(d) {
  sample(c(rep(1:nfolds, each=n[[d]] %/% nfolds),
           rep(1:(n[[d]] %% nfolds), (n[[d]] %% nfolds)!=0)))})

# semnig models
fit.semnig1 <- semnig(x=x, y=y, 
                      C=sapply(x, function(s) {matrix(rep(1, ncol(s)))}), 
                      Z=matrix(1, nrow=D), full.post=TRUE, 
                      control=control.semnig)
fit.semnig2 <- semnig(x=x, y=y, C=C, Z=matrix(1, nrow=D),
                      full.post=TRUE, control=control.semnig)

# penalized regression models
fit.ridge1 <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[[d]], alpha=0, foldid=foldid[[d]], intercept=FALSE)})
fit.lasso1 <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[[d]], foldid=foldid[[d]], intercept=FALSE)})

# fit xtune
fit.xtune1 <- lapply(1:D, function(d) {
  xtune(x[[d]], y[[d]], C[[d]][, -1], family="linear", method="ridge", 
        message=FALSE, control=list(intercept=FALSE))})

# fit ebridge model
fit.ebridge1 <- ebridge(x, y, lapply(C, function(s) {matrix(s[, -1])}), NULL,
                        foldid=foldid, hyper=list(lambda=1/sqrt(n*sapply(
                          fit.ridge1, "[[", "lambda.min")), zeta=0, nu=0),
                        control=control.ebridge)

save(fit.semnig1, fit.semnig2, fit.ridge1, fit.lasso1, fit.xtune1, fit.ebridge1, 
     file="results/analysis_gdsc_fit1.Rdata")

# eb estimates
est <- cbind(c(fit.semnig1$eb$alphaf, NA, fit.semnig1$eb$alphad, 
               fit.semnig1$eb$lambdaf, fit.semnig1$eb$lambdad),
             c(fit.semnig2$eb$alphaf, fit.semnig2$eb$alphad, 
               fit.semnig2$eb$lambdaf, fit.semnig2$eb$lambdad),
             rep(NA, 5), rep(NA, 5),
             c(colMeans(-t(sapply(fit.xtune1, "[[", "alpha.est"))), 
               rep(NA, 3)),
             c(NA, fit.ebridge1$alphaf, NA, fit.ebridge1$alphad, rep(NA, 2)))
rownames(est) <- c("interceptf", "pvalue", "interceptd", "lambdaf", "lambdad")
colnames(est) <- c("semnig1", "semnig2", "ridge1", "lasso1", "xtune1", 
                   "ebridge1")
write.table(est, file="results/analysis_gdsc_fit1.txt")

# log cpo values
lcpo.semnig1 <- logcpo(x, y, n, fit.semnig1)
lcpo.semnig2 <- logcpo(x, y, n, fit.semnig2)
save(lcpo.semnig1, lcpo.semnig2, file="results/analysis_gdsc_cpo1.Rdata")

### predictions and mse
# settings
psel.dss <- 10
control.semnig <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-3, 
                       epsilon.vb=1e-3, maxit.eb=500, maxit.vb=1, 
                       maxit.post=100, maxit.block=0)
control.ebridge <-list(epsilon=sqrt(.Machine$double.eps), maxit=500, 
                       trace=FALSE, glmnet.fit2=FALSE, beta2=FALSE)

set.seed(2020)
foldid <- lapply(1:D, function(d) {
  sample(c(rep(1:nfolds, each=n[[d]] %/% nfolds),
           rep(1:(n[[d]] %% nfolds), (n[[d]] %% nfolds)!=0)))})

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
res <- foreach(r=1:nfolds, .packages=packages, .errorhandling="pass") %dopar% {
  cat("\r", "fold", r)
  set.seed(2020 + r)
  
  # splitting data
  xtrain <- lapply(1:D, function(d) {
    scale(expr$expr[idcellsel[[d]], idfeatsel[[d]]][foldid[[d]]!=r, ])})
  xtest <- lapply(1:D, function(d) {
    scale(expr$expr[idcellsel[[d]], idfeatsel[[d]]][foldid[[d]]==r, ])})
  
  ytrain <- lapply(1:D, function(d) {scale(resp.sel[[d]][foldid[[d]]!=r])[, 1]})
  ytest <- lapply(1:D, function(d) {scale(resp.sel[[d]][foldid[[d]]==r])[, 1]})
  
  ntrain <- sapply(foldid, function(s) {sum(s!=r)})
  foldid2 <- lapply(1:D, function(d) {
    sample(c(rep(1:nfolds, each=ntrain[d] %/% nfolds),
             rep(1:(ntrain[d] %% nfolds), (ntrain[d] %% nfolds)!=0)))})
  
  # create co-data
  C <- lapply(1:D, function(d) {
    s <- cbind(1, pvalues[[d]]); colnames(s) <- c("intercept", "pvalue"); s})
  
  # semnig models
  cv.semnig1 <- semnig(x=xtrain, y=ytrain, 
                       C=sapply(xtrain, function(s) {matrix(rep(1, ncol(s)))}), 
                       Z=matrix(1, nrow=D), full.post=TRUE, 
                       control=control.semnig)
  cv.semnig2 <- semnig(x=xtrain, y=ytrain, C=C, Z=matrix(1, nrow=D),
                       full.post=TRUE, control=control.semnig)
  
  # penalized regression models
  cv.ridge1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[[d]], alpha=0, foldid=foldid2[[d]], 
              intercept=FALSE)})
  cv.lasso1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[[d]], foldid=foldid2[[d]], intercept=FALSE)})
  
  # fit xtune
  cv.xtune1 <- lapply(1:D, function(d) {
    xtune(xtrain[[d]], ytrain[[d]], C[[d]][, -1], family="linear", 
          method="ridge", message=FALSE, control=list(intercept=FALSE))})
  
  # fit ebridge model
  cv.ebridge1 <- ebridge(xtrain, ytrain,
                         lapply(C, function(s) {matrix(s[, -1])}), NULL,
                         foldid=foldid2,
                         hyper=list(lambda=1/sqrt(n*sapply(
                           cv.ridge1, "[[", "lambda.min")), zeta=0, nu=0),
                         control=control.ebridge)
  
  # post selection
  cv.dss.semnig1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], as.numeric(xtrain[[d]] %*% cv.semnig1$vb$mu[[d]]),
              penalty.factor=0.5/abs(cv.semnig1$vb$mu[[d]]), intercept=FALSE)})
  cv.dss.ridge1 <- lapply(1:D, function(d) {
    b <- coef(cv.ridge1[[d]], s="lambda.min")[-1, ]
    cv.glmnet(xtrain[[d]], as.numeric(xtrain[[d]] %*% b),
              penalty.factor=0.5/abs(b), intercept=FALSE)})
  
  # estimates
  best <- list(semnig1=cv.semnig1$vb$mpost$beta, 
               semnig2=cv.semnig2$vb$mpost$beta,
               ridge1=lapply(cv.ridge1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}),
               lasso1=lapply(cv.lasso1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}),
               xtune1=lapply(cv.xtune1, function(s) {
                 unname(s$beta.est[-1, ])}),
               ebridge1=cv.ebridge1$beta1)
  best$dss.semnig1 <- lapply(cv.dss.semnig1, function(s) {
    coef(s, s="lambda.min")[-1, 1]})
  best$dss.semnig2 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.semnig1[[d]], psel=sum(best$lasso1[[d]]!=0))[-1, ]})
  best$dss.semnig3 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.semnig1[[d]], psel=psel.dss)[-1, ]})
  best$dss.lasso1 <- lapply(1:D, function(d) {
    coef.dss(cv.lasso1[[d]], psel=psel.dss)[-1, ]})
  best$dss.ridge1 <- lapply(cv.dss.ridge1, function(s) {
    coef(s, s="lambda.min")[-1, 1]})
  best$dss.ridge2 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.ridge1[[d]], psel=sum(best$lasso1[[d]]!=0))[-1, ]})
  best$dss.ridge3 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.ridge1[[d]], psel=psel.dss)[-1, ]})
  
  pred <- lapply(best, function(b) {
    sapply(1:D, function(d) {as.numeric(xtest[[d]] %*% b[[d]])})})
  pmse <- sapply(pred, function(s) {
    sapply(1:D, function(d) {mean((s[[d]] - ytest[[d]])^2)})})
  psel <- sapply(best, function(b) {
    sapply(b, function(s) {sum(s!=0)})})
  
  rownames(pmse) <- paste0("drug", 1:D)
  rownames(psel) <- paste0("drug", 1:D)
  list(pmse=pmse, psel=psel)
}
stopCluster(cl=cl)

pmse <- Reduce("rbind", lapply(res, "[[", "pmse"))
rownames(pmse) <- paste0("pmse.", rownames(pmse))
psel <- Reduce("rbind", lapply(res, "[[", "psel"))
rownames(psel) <- paste0("psel.", rownames(psel))

res <- rbind(pmse, psel)
write.table(res, file="results/analysis_gdsc_res1.txt")


################################################################################
################################## analysis 2 ##################################
################################################################################
# mutation data as extra data with feature type as co-data

### data preparation
expr.psel <- 300
order.sd <- order(-apply(expr$expr, 2, sd))
expr.sel <- expr$expr[rownames(expr$expr) %in% rownames(mut), 
                      order.sd[1:expr.psel]]
expr.sel <- expr.sel[order(rownames(expr.sel)), ]
mut.sel <- mut[rownames(mut) %in% rownames(expr.sel), ]
mut.sel <- mut.sel[order(rownames(mut.sel)), apply(mut.sel, 2, sd)!=0]
mut.sel <- unique(mut.sel, MARGIN=2)
feat <- cbind(expr.sel, mut.sel)
feat <- feat[order(rownames(feat)), ]
resp.sel <- sapply(resp, function(s) {
  s <- s[names(s) %in% rownames(feat)]
  s <- s[order(names(s))]})
idsel <- lapply(resp.sel, function(s) {
  which(rownames(feat) %in% names(s))})
mutation <- rep(c(0, 1), times=c(ncol(expr.sel), ncol(mut.sel)))

# split into training and test data
set.seed(2020)
D <- length(resp.sel)
n <- sapply(resp.sel, length)

x <- lapply(1:D, function(d) {
  cbind(scale(feat[idsel[[d]], mutation==0]), feat[idsel[[d]], mutation==1])})
y <- lapply(1:D, function(d) {scale(resp.sel[[d]])[, 1]})

# create co-data
C <- cbind(1, mutation); colnames(C) <- c("intercept", "mutation")
C <- replicate(D, list(C))

### model fitting
# settings
control.semnig <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, 
                       epsilon.vb=1e-3, maxit.eb=500, maxit.vb=1, 
                       maxit.post=100, maxit.block=0)
control.ebridge <-list(epsilon=sqrt(.Machine$double.eps), maxit=500, 
                       trace=TRUE, glmnet.fit2=FALSE, beta2=TRUE)

# creating folds
set.seed(2020)
foldid <- lapply(1:D, function(d) {
  sample(c(rep(1:nfolds, each=n[[d]] %/% nfolds),
           rep(1:(n[[d]] %% nfolds), (n[[d]] %% nfolds)!=0)))})

# semnig models
fit.semnig1 <- semnig(x=x, y=y, 
                      C=lapply(x, function(s) {matrix(rep(1, ncol(s)))}), 
                      Z=matrix(1, nrow=D), full.post=TRUE, 
                      control=control.semnig)
fit.semnig2 <- semnig(x=x, y=y, C=C, Z=matrix(1, nrow=D),
                      full.post=TRUE, control=control.semnig)

# penalized regression models
fit.ridge1 <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[[d]], alpha=0, foldid=foldid[[d]], intercept=FALSE)})
fit.lasso1 <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[[d]], foldid=foldid[[d]], intercept=FALSE)})

# fit xtune
fit.xtune1 <- lapply(1:D, function(d) {
  xtune(x[[d]], y[[d]], C[[d]][, -1], family="linear", method="ridge", 
        control=list(intercept=FALSE))})

# fit ebridge model
fit.ebridge1 <- ebridge(x, y, lapply(C, function(s) {matrix(s[, -1])}), NULL,
                        foldid=foldid,
                        hyper=list(lambda=1/sqrt(n*sapply(
                          fit.ridge1, "[[", "lambda.min")), zeta=0, nu=0),
                        control=control.ebridge)

# bSEM model
fit.bSEM1 <- bSEM(x, y, lapply(C, function(s) {s[, -1] + 1}),
                  control=list(maxit=500, trace=TRUE, epsilon=1e-3))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         

# saving fitted model objects
save(fit.semnig1, fit.semnig2, fit.lasso1, fit.ridge1, fit.xtune1, fit.ebridge1, 
     fit.bSEM1, file="results/analysis_gdsc_fit2.Rdata")

est <- cbind(c(fit.semnig1$eb$alphaf, NA, fit.semnig1$eb$alphad, 
               fit.semnig1$eb$lambdaf, fit.semnig1$eb$lambdad),
             c(fit.semnig2$eb$alphaf, fit.semnig2$eb$alphad, 
               fit.semnig2$eb$lambdaf, fit.semnig2$eb$lambdad),
             rep(NA, 5), rep(NA, 5),
             c(colMeans(-t(sapply(fit.xtune1, "[[", "alpha.est"))), 
               rep(NA, 3)),
             c(NA, fit.ebridge1$alphaf, NA, fit.ebridge1$alphad, rep(NA, 2)),
             c(fit.bSEM1$eb$a[nrow(fit.bSEM1$eb$a), 1]/
                 fit.bSEM1$eb$b[nrow(fit.bSEM1$eb$b), 1],
               fit.bSEM1$eb$a[nrow(fit.bSEM1$eb$a), 2]/
                 fit.bSEM1$eb$b[nrow(fit.bSEM1$eb$b), 2] -
                 fit.bSEM1$eb$a[nrow(fit.bSEM1$eb$a), 1]/
                 fit.bSEM1$eb$b[nrow(fit.bSEM1$eb$b), 1], rep(NA, 3)))
rownames(est) <- c("interceptf", "mutation", "interceptd", "lambdaf", "lambdad")
colnames(est) <- c("semnig1", "semnig2", "ridge1", "lasso1", "xtune1", 
                   "ebridge1", "bSEM1")
write.table(est, file="results/analysis_gdsc_fit2.txt")

# log cpo values
lcpo.semnig1 <- logcpo(x, y, n, fit.semnig1)
lcpo.semnig2 <- logcpo(x, y, n, fit.semnig2)
save(lcpo.semnig1, lcpo.semnig2, file="results/analysis_gdsc_cpo2.Rdata")

# predictions and mse
### predictions and mse
# settings
psel.dss <- 10
control.semnig <- list(conv.post=FALSE, trace=FALSE, epsilon.eb=1e-3, 
                       epsilon.vb=1e-3, maxit.eb=500, maxit.vb=1, 
                       maxit.post=100, maxit.block=0)
control.ebridge <-list(epsilon=sqrt(.Machine$double.eps), maxit=500, 
                       trace=FALSE, glmnet.fit2=FALSE, beta2=TRUE)

set.seed(2020)
foldid <- lapply(1:D, function(d) {
  sample(c(rep(1:nfolds, each=n[[d]] %/% nfolds),
           rep(1:(n[[d]] %% nfolds), (n[[d]] %% nfolds)!=0)))})

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
res <- foreach(r=1:nfolds, .packages=packages, .errorhandling="pass") %dopar% {
  cat("\r", "fold", r)
  set.seed(2020 + r)
  
  # splitting data
  xtrain <- lapply(1:D, function(d) {
    cbind(scale(feat[idsel[[d]], mutation==0][foldid[[d]]!=r, ]),
          feat[idsel[[d]], mutation==1][foldid[[d]]!=r, ])})
  xtest <- lapply(1:D, function(d) {
    cbind(scale(feat[idsel[[d]], mutation==0][foldid[[d]]==r, ]),
          feat[idsel[[d]], mutation==1][foldid[[d]]==r, ])})
  duptrain <- lapply(xtrain, function(s) {which(duplicated(s, MARGIN=2))})
  xtrain <- lapply(1:D, function(d) {
    if(length(duptrain[[d]])==0) {
      xtrain[[d]]
    } else {
      xtrain[[d]][, -duptrain[[d]]]
    }})
  
  ytrain <- lapply(1:D, function(d) {scale(resp.sel[[d]][foldid[[d]]!=r])[, 1]})
  ytest <- lapply(1:D, function(d) {scale(resp.sel[[d]][foldid[[d]]==r])[, 1]})
  
  ntrain <- sapply(foldid, function(s) {sum(s!=r)})
  foldid2 <- lapply(1:D, function(d) {
    sample(c(rep(1:nfolds, each=ntrain[d] %/% nfolds),
             rep(1:(ntrain[d] %% nfolds), (ntrain[d] %% nfolds)!=0)))})
  
  # create co-data
  C <- lapply(duptrain, function(dup) {
    if(length(dup)==0) {
      s <- cbind(1, mutation)
    } else {
      s <- cbind(1, mutation[-dup])
    }
    colnames(s) <- c("intercept", "mutation"); s})
  
  # semnig models
  cv.semnig1 <- semnig(x=xtrain, y=ytrain, 
                       C=sapply(xtrain, function(s) {matrix(rep(1, ncol(s)))}), 
                       Z=matrix(1, nrow=D), full.post=TRUE, 
                       control=control.semnig)
  cv.semnig2 <- semnig(x=xtrain, y=ytrain, C=C, Z=matrix(1, nrow=D),
                       full.post=TRUE, control=control.semnig)
  
  # penalized regression models
  cv.ridge1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[[d]], alpha=0, foldid=foldid2[[d]], 
              intercept=FALSE)})
  cv.lasso1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[[d]], foldid=foldid2[[d]], intercept=FALSE)})
  
  # fit xtune
  cv.xtune1 <- lapply(1:D, function(d) {
    xtune(xtrain[[d]], ytrain[[d]], C[[d]][, -1], family="linear", 
          method="ridge", message=FALSE, control=list(intercept=FALSE))})
  
  # fit ebridge model
  cv.ebridge1 <- ebridge(xtrain, ytrain,
                         lapply(C, function(s) {matrix(s[, -1])}), NULL,
                         foldid=foldid2,
                         hyper=list(lambda=1/sqrt(n*sapply(
                           cv.ridge1, "[[", "lambda.min")), zeta=0, nu=0),
                         control=control.ebridge)
  
  # bSEM model
  cv.bSEM1 <- bSEM(xtrain, ytrain, lapply(C, function(s) {s[, -1] + 1}),
                    control=list(maxit=500, trace=FALSE, epsilon=1e-3))
  
  # post selection
  cv.dss.semnig1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], as.numeric(xtrain[[d]] %*% cv.semnig1$vb$mu[[d]]),
              penalty.factor=0.5/abs(cv.semnig1$vb$mu[[d]]), intercept=FALSE)})
  cv.dss.ridge1 <- lapply(1:D, function(d) {
    b <- coef(cv.ridge1[[d]], s="lambda.min")[-1, ]
    cv.glmnet(xtrain[[d]], as.numeric(xtrain[[d]] %*% b),
      penalty.factor=0.5/abs(b), intercept=FALSE)})
  cv.dss.bSEM1 <- lapply(1:D, function(d) {
    b <- cv.bSEM1$vb$beta[[d]][, "mu"]
    cv.glmnet(xtrain[[d]], as.numeric(xtrain[[d]] %*% b),
              penalty.factor=0.5/abs(b), intercept=FALSE)})
  
  # estimates
  best <- list(semnig1=cv.semnig1$vb$mpost$beta, 
               semnig2=cv.semnig2$vb$mpost$beta,
               ridge1=lapply(cv.ridge1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}),
               lasso1=lapply(cv.lasso1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}),
               xtune1=lapply(cv.xtune1, function(s) {
                 unname(s$beta.est[-1, ])}),
               ebridge1=cv.ebridge1$beta1,
               bSEM1=lapply(cv.bSEM1$vb$beta, function(s) {s[, "mu"]}))
  best$dss.semnig1 <- lapply(cv.dss.semnig1, function(s) {
    coef(s, s="lambda.min")[-1, 1]})
  best$dss.semnig2 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.semnig1[[d]], psel=sum(best$lasso1[[d]]!=0))[-1, ]})
  best$dss.semnig3 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.semnig1[[d]], psel=psel.dss)[-1, ]})
  best$dss.lasso1 <- lapply(1:D, function(d) {
    coef.dss(cv.lasso1[[d]], psel=psel.dss)[-1, ]})
  best$dss.ridge1 <- lapply(cv.dss.ridge1, function(s) {
    coef(s, s="lambda.min")[-1, 1]})
  best$dss.ridge2 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.ridge1[[d]], psel=sum(best$lasso1[[d]]!=0))[-1, ]})
  best$dss.ridge3 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.ridge1[[d]], psel=psel.dss)[-1, ]})
  best$dss.bSEM1 <- lapply(cv.dss.bSEM1, function(s) {
    coef(s, s="lambda.min")[-1, 1]})
  best$dss.bSEM2 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.bSEM1[[d]], psel=sum(best$lasso1[[d]]!=0))[-1, ]})
  best$dss.bSEM3 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.bSEM1[[d]], psel=psel.dss)[-1, ]})
  
  # add duplicated variables
  best <- lapply(best, function(b) {
    lapply(1:D, function(d) {
      if(length(duptrain[[d]])==0) {
        b[[d]]
      } else {
        replace(rep(0, ncol(feat)), -duptrain[[d]], b[[d]])  
      }})})
  
  pred <- lapply(best, function(b) {
    sapply(1:D, function(d) {as.numeric(xtest[[d]] %*% b[[d]])})})
  pmse <- sapply(pred, function(s) {
    sapply(1:D, function(d) {mean((s[[d]] - ytest[[d]])^2)})})
  psel <- sapply(best, function(b) {
    sapply(b, function(s) {sum(s!=0)})})
  
  rownames(pmse) <- paste0("drug", 1:D)
  rownames(psel) <- paste0("drug", 1:D)
  list(pmse=pmse, psel=psel)
}
stopCluster(cl=cl)

pmse <- Reduce("rbind", lapply(res, "[[", "pmse"))
rownames(pmse) <- paste0("pmse.", rownames(pmse))
psel <- Reduce("rbind", lapply(res, "[[", "psel"))
rownames(psel) <- paste0("psel.", rownames(psel))

res <- rbind(pmse, psel)
write.table(res, file="results/analysis_gdsc_res2.txt")

################################################################################
################################## analysis 3 ##################################
################################################################################
# pvalues from ccle data
ccle <- read.table(file="results/analysis_ccle_res1.txt", header=TRUE)

### data preparation
psel <- 500
expr.sel <- expr$expr[, colnames(expr$expr) %in% rownames(ccle)]
order.sd <- order(-apply(expr.sel, 2, sd))
expr.sel <- expr.sel[, order.sd[1:psel]]
ccle.sel <- ccle[na.omit(match(colnames(expr.sel), rownames(ccle))), ]
pvalues <- ccle.sel$harmonic^(1/13)

resp.sel <- sapply(resp, function(s) {
  s <- s[names(s) %in% rownames(expr.sel)]
  s <- s[order(names(s))]})
idsel <- lapply(resp.sel, function(s) {
  which(rownames(expr.sel) %in% names(s))})

D <- length(resp.sel)
n <- sapply(resp.sel, length)

x <- lapply(1:D, function(d) {scale(expr.sel[idsel[[d]], ])})
y <- lapply(1:D, function(d) {scale(resp.sel[[d]])[, 1]})

# create co-data
C <- replicate(D, list(cbind(intercept=1, pvalue=pvalues)))

### model fitting
# settings
control.semnig <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, 
                       epsilon.vb=1e-3, maxit.eb=500, maxit.vb=1, 
                       maxit.post=100, maxit.block=0)
control.ebridge <-list(epsilon=sqrt(.Machine$double.eps), maxit=500, 
                       trace=TRUE, glmnet.fit2=FALSE, beta2=FALSE)

# creating folds
set.seed(2020)
foldid <- lapply(1:D, function(d) {
  sample(c(rep(1:nfolds, each=n[[d]] %/% nfolds),
           rep(1:(n[[d]] %% nfolds), (n[[d]] %% nfolds)!=0)))})

# semnig models
fit.semnig1 <- semnig(x=x, y=y, 
                      C=lapply(x, function(s) {matrix(rep(1, ncol(s)))}), 
                      Z=matrix(1, nrow=D), full.post=TRUE, 
                      control=control.semnig)
fit.semnig2 <- semnig(x=x, y=y, C=C, Z=matrix(1, nrow=D),
                      full.post=TRUE, control=control.semnig)

# penalized regression models
fit.ridge1 <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[[d]], alpha=0, foldid=foldid[[d]], intercept=FALSE)})
fit.lasso1 <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[[d]], foldid=foldid[[d]], intercept=FALSE)})

# fit xtune
fit.xtune1 <- lapply(1:D, function(d) {
  xtune(x[[d]], y[[d]], C[[d]][, -1], family="linear", method="ridge", 
        control=list(intercept=FALSE))})

# fit ebridge model
fit.ebridge1 <- ebridge(x, y, lapply(C, function(s) {matrix(s[, -1])}), NULL,
                        foldid=foldid,
                        hyper=list(lambda=1/sqrt(n*sapply(
                          fit.ridge1, "[[", "lambda.min")), zeta=0, nu=0),
                        control=control.ebridge)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                

# saving fitted model objects
save(fit.semnig1, fit.semnig2, fit.lasso1, fit.ridge1, fit.xtune1, fit.ebridge1, 
     file="results/analysis_gdsc_fit3.Rdata")

est <- cbind(c(fit.semnig1$eb$alphaf, NA, fit.semnig1$eb$alphad, 
               fit.semnig1$eb$lambdaf, fit.semnig1$eb$lambdad),
             c(fit.semnig2$eb$alphaf, fit.semnig2$eb$alphad, 
               fit.semnig2$eb$lambdaf, fit.semnig2$eb$lambdad),
             rep(NA, 5), rep(NA, 5),
             c(colMeans(-t(sapply(fit.xtune1, "[[", "alpha.est"))), 
               rep(NA, 3)),
             c(NA, fit.ebridge1$alphaf, NA, fit.ebridge1$alphad, rep(NA, 2)))
rownames(est) <- c("interceptf", "pvalue", "interceptd", "lambdaf", "lambdad")
colnames(est) <- c("semnig1", "semnig2", "ridge1", "lasso1", "xtune1", 
                   "ebridge1")
write.table(est, file="results/analysis_gdsc_fit3.txt")

# log cpo values
lcpo.semnig1 <- logcpo(x, y, n, fit.semnig1)
lcpo.semnig2 <- logcpo(x, y, n, fit.semnig2)
save(lcpo.semnig1, lcpo.semnig2, file="results/analysis_gdsc_cpo3.Rdata")

### predictions and mse
# settings
psel.dss <- 10
control.semnig <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-3, 
                       epsilon.vb=1e-3, maxit.eb=500, maxit.vb=1, 
                       maxit.post=100, maxit.block=0)
control.ebridge <-list(epsilon=sqrt(.Machine$double.eps), maxit=500, 
                       trace=FALSE, glmnet.fit2=FALSE, beta2=FALSE)

set.seed(2020)
foldid <- lapply(1:D, function(d) {
  sample(c(rep(1:nfolds, each=n[[d]] %/% nfolds),
           rep(1:(n[[d]] %% nfolds), (n[[d]] %% nfolds)!=0)))})

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
res <- foreach(r=1:nfolds, .packages=packages, .errorhandling="pass") %dopar% {
  cat("\r", "fold", r)
  set.seed(2020 + r)
  
  ### pmse
  xtrain <- lapply(1:D, function(d) {
    scale(expr.sel[idsel[[d]], ][foldid[[d]]!=r, ])})
  xtest <- lapply(1:D, function(d) {
    scale(expr.sel[idsel[[d]], ][foldid[[d]]==r, ])})
  
  ytrain <- lapply(1:D, function(d) {scale(resp.sel[[d]][foldid[[d]]!=r])[, 1]})
  ytest <- lapply(1:D, function(d) {scale(resp.sel[[d]][foldid[[d]]==r])[, 1]})
  
  
  ntrain <- sapply(foldid, function(s) {sum(s!=r)})
  foldid2 <- lapply(1:D, function(d) {
    sample(c(rep(1:nfolds, each=ntrain[d] %/% nfolds),
             rep(1:(ntrain[d] %% nfolds), (ntrain[d] %% nfolds)!=0)))})
  
  # create co-data
  C <- replicate(D, list(cbind(intercept=1, pvalue=pvalues)))
  
  # semnig models
  cv.semnig1 <- semnig(x=xtrain, y=ytrain, 
                       C=lapply(xtrain, function(s) {matrix(rep(1, ncol(s)))}), 
                       Z=matrix(1, nrow=D), full.post=TRUE, 
                       control=control.semnig)
  cv.semnig2 <- semnig(x=xtrain, y=ytrain, C=C, Z=matrix(1, nrow=D),
                       full.post=TRUE, control=control.semnig)
  
  # penalized regression models
  cv.ridge1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[[d]], alpha=0, foldid=foldid2[[d]], 
              intercept=FALSE)})
  cv.lasso1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[[d]], foldid=foldid2[[d]], intercept=FALSE)})
  
  # fit xtune
  cv.xtune1 <- lapply(1:D, function(d) {
    xtune(xtrain[[d]], ytrain[[d]], C[[d]][, -1], family="linear", 
          method="ridge", message=FALSE, control=list(intercept=FALSE))})
  
  # fit ebridge model
  cv.ebridge1 <- ebridge(xtrain, ytrain,
                         lapply(C, function(s) {matrix(s[, -1])}), NULL,
                         foldid=foldid2,
                         hyper=list(lambda=1/sqrt(n*sapply(
                           cv.ridge1, "[[", "lambda.min")), zeta=0, nu=0),
                         control=control.ebridge)
  
  # post selection
  cv.dss.semnig1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], as.numeric(xtrain[[d]] %*% cv.semnig1$vb$mu[[d]]),
              penalty.factor=0.5/abs(cv.semnig1$vb$mu[[d]]), intercept=FALSE)})
  cv.dss.ridge1 <- lapply(1:D, function(d) {
    b <- coef(cv.ridge1[[d]], s="lambda.min")[-1, ]
    cv.glmnet(xtrain[[d]], as.numeric(xtrain[[d]] %*% b),
              penalty.factor=0.5/abs(b), intercept=FALSE)})
  
  # estimates
  best <- list(semnig1=cv.semnig1$vb$mpost$beta, 
               semnig2=cv.semnig2$vb$mpost$beta,
               ridge1=lapply(cv.ridge1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}),
               lasso1=lapply(cv.lasso1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}),
               xtune1=lapply(cv.xtune1, function(s) {
                 unname(s$beta.est[-1, ])}),
               ebridge1=cv.ebridge1$beta1)
  best$dss.semnig1 <- lapply(cv.dss.semnig1, function(s) {
    coef(s, s="lambda.min")[-1, 1]})
  best$dss.semnig2 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.semnig1[[d]], psel=sum(best$lasso1[[d]]!=0))[-1, ]})
  best$dss.semnig3 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.semnig1[[d]], psel=psel.dss)[-1, ]})
  best$dss.lasso1 <- lapply(1:D, function(d) {
    coef.dss(cv.lasso1[[d]], psel=psel.dss)[-1, ]})
  best$dss.ridge1 <- lapply(cv.dss.ridge1, function(s) {
    coef(s, s="lambda.min")[-1, 1]})
  best$dss.ridge2 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.ridge1[[d]], psel=sum(best$lasso1[[d]]!=0))[-1, ]})
  best$dss.ridge3 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.ridge1[[d]], psel=psel.dss)[-1, ]})
  
  pred <- lapply(best, function(b) {
    sapply(1:D, function(d) {as.numeric(xtest[[d]] %*% b[[d]])})})
  pmse <- sapply(pred, function(s) {
    sapply(1:D, function(d) {mean((s[[d]] - ytest[[d]])^2)})})
  psel <- sapply(best, function(b) {
    sapply(b, function(s) {sum(s!=0)})})
  
  rownames(pmse) <- paste0("drug", 1:D)
  rownames(psel) <- paste0("drug", 1:D)
  list(pmse=pmse, psel=psel)
}
stopCluster(cl=cl)

pmse <- Reduce("rbind", lapply(res, "[[", "pmse"))
rownames(pmse) <- paste0("pmse.", rownames(pmse))
psel <- Reduce("rbind", lapply(res, "[[", "psel"))
rownames(psel) <- paste0("psel.", rownames(psel))

res <- rbind(pmse, psel)
write.table(res, file="results/analysis_gdsc_res3.txt")

################################################################################
################################## analysis 4 ##################################
################################################################################
# expression with regular external data

### data preparation
psel <- 500
order.sd <- order(-apply(expr$expr, 2, sd))
expr.sel <- expr$expr[, order.sd[1:psel]]
resp.sel <- sapply(resp, function(s) {
  s <- s[names(s) %in% rownames(expr.sel)]
  s <- s[order(names(s))]})
resp.sel <- resp.sel[names(resp.sel) %in% drug$name]
drug.sel <- drug[na.omit(match(names(resp), drug$name)), ]
in_pathway.sel <- lapply(expr$in_pathway, function(s) {
  s[na.omit(match(colnames(expr.sel), names(s)))]})
in_pathway.sel <- in_pathway.sel[na.omit(match(names(resp.sel), 
                                               names(in_pathway.sel)))]

idsel <- lapply(resp.sel, function(s) {
  which(rownames(expr.sel) %in% names(s))})

# split into training and test data
D <- length(resp.sel)
n <- sapply(resp.sel, length)

x <- lapply(1:D, function(d) {scale(expr.sel[idsel[[d]], ])})
y <- lapply(1:D, function(d) {scale(resp.sel[[d]])[, 1]})

# create co-data
Z <- model.matrix(~ replace(drug.sel$stage, is.na(drug.sel$stage),
                            "experimental") + replace(drug.sel$action, is.na(
                              drug.sel$action), "unknown"))
colnames(Z) <- c("intercept", "experimental", "in clinical development", "targeted", 
                 "unknown")
C <- lapply(in_pathway.sel, function(s) {
  s <- cbind(1, s); colnames(s) <- c("intercept", "in pathway"); s})

### model fitting
# settings
control.semnig <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, 
                       epsilon.vb=1e-3, maxit.eb=500, maxit.vb=1, 
                       maxit.post=100, maxit.block=0)
control.ebridge <-list(epsilon=sqrt(.Machine$double.eps), maxit=500, 
                       trace=TRUE, glmnet.fit2=FALSE, beta2=TRUE)

# creating folds
set.seed(2020)
foldid <- lapply(1:D, function(d) {
  sample(c(rep(1:nfolds, each=n[[d]] %/% nfolds),
           rep(1:(n[[d]] %% nfolds), (n[[d]] %% nfolds)!=0)))})

# semnig models
fit.semnig1 <- semnig(x=x, y=y, 
                      C=lapply(x, function(s) {matrix(rep(1, ncol(s)))}), 
                      Z=matrix(1, nrow=D), full.post=TRUE, 
                      control=control.semnig)
fit.semnig2 <- semnig(x=x, y=y, C=C, Z=Z,
                      full.post=TRUE, control=control.semnig)

# penalized regression models
fit.ridge1 <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[[d]], alpha=0, foldid=foldid[[d]], intercept=FALSE)})
fit.lasso1 <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[[d]], foldid=foldid[[d]], intercept=FALSE)})

# fit xtune
fit.xtune1 <- lapply(1:D, function(d) {
  xtune(x[[d]], y[[d]], C[[d]][, -1], family="linear", method="ridge", 
        control=list(intercept=FALSE))})

# fit ebridge model
fit.ebridge1 <- ebridge(x, y, lapply(C, function(s) {matrix(s[, -1])}), Z[, -1],
                        foldid=foldid,
                        hyper=list(lambda=1/sqrt(n*sapply(
                          fit.ridge1, "[[", "lambda.min")), zeta=0, nu=0),
                        control=control.ebridge)

# bSEM model
fit.bSEM1 <- bSEM(x, y, lapply(C, function(s) {s[, -1] + 1}),
                  control=list(maxit=500, trace=TRUE, epsilon=1e-3))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         

# saving fitted model objects
save(fit.semnig1, fit.semnig2, fit.lasso1, fit.ridge1, fit.xtune1, fit.ebridge1, 
     fit.bSEM1, file="results/analysis_gdsc_fit4.Rdata")

est <- cbind(c(fit.semnig1$eb$alphaf, NA, fit.semnig1$eb$alphad, rep(NA, 4),
               fit.semnig1$eb$lambdaf, fit.semnig1$eb$lambdad),
             c(fit.semnig2$eb$alphaf, fit.semnig2$eb$alphad, 
               fit.semnig2$eb$lambdaf, fit.semnig2$eb$lambdad),
             rep(NA, 9), rep(NA, 9),
             c(colMeans(-t(sapply(fit.xtune1, "[[", "alpha.est"))), 
               rep(NA, 7)),
             c(NA, fit.ebridge1$alphaf, NA, fit.ebridge1$alphad, rep(NA, 2)),
             c(fit.bSEM1$eb$a[nrow(fit.bSEM1$eb$a), 1]/
                 fit.bSEM1$eb$b[nrow(fit.bSEM1$eb$b), 1],
               fit.bSEM1$eb$a[nrow(fit.bSEM1$eb$a), 2]/
                 fit.bSEM1$eb$b[nrow(fit.bSEM1$eb$b), 2] -
                 fit.bSEM1$eb$a[nrow(fit.bSEM1$eb$a), 1]/
                 fit.bSEM1$eb$b[nrow(fit.bSEM1$eb$b), 1], rep(NA, 7)))
rownames(est) <- c("interceptf", "inpathway", "interceptd", "experimental", 
                   "in clinical development", "targeted", "unknown", "lambdaf", 
                   "lambdad")
colnames(est) <- c("semnig1", "semnig2", "ridge1", "lasso1", "xtune1", 
                   "ebridge1", "bSEM1")
write.table(est, file="results/analysis_gdsc_fit4.txt")

# log cpo values
lcpo.semnig1 <- logcpo(x, y, n, fit.semnig1)
lcpo.semnig2 <- logcpo(x, y, n, fit.semnig2)
save(lcpo.semnig1, lcpo.semnig2, file="results/analysis_gdsc_cpo4.Rdata")

# predictions and mse
### predictions and mse
# settings
psel.dss <- 10
control.semnig <- list(conv.post=FALSE, trace=FALSE, epsilon.eb=1e-3, 
                       epsilon.vb=1e-3, maxit.eb=500, maxit.vb=1, 
                       maxit.post=100, maxit.block=0)
control.ebridge <-list(epsilon=sqrt(.Machine$double.eps), maxit=500, 
                       trace=FALSE, glmnet.fit2=FALSE, beta2=TRUE)

set.seed(2020)
foldid <- lapply(1:D, function(d) {
  sample(c(rep(1:nfolds, each=n[[d]] %/% nfolds),
           rep(1:(n[[d]] %% nfolds), (n[[d]] %% nfolds)!=0)))})

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
res <- foreach(r=1:nfolds, .packages=packages, .errorhandling="pass") %dopar% {
  cat("\r", "fold", r)
  set.seed(2020 + r)
  
  ### pmse
  xtrain <- lapply(1:D, function(d) {
    scale(expr.sel[idsel[[d]], ][foldid[[d]]!=r, ])})
  xtest <- lapply(1:D, function(d) {
    scale(expr.sel[idsel[[d]], ][foldid[[d]]==r, ])})
  
  ytrain <- lapply(1:D, function(d) {scale(resp.sel[[d]][foldid[[d]]!=r])[, 1]})
  ytest <- lapply(1:D, function(d) {scale(resp.sel[[d]][foldid[[d]]==r])[, 1]})
  
  ntrain <- sapply(foldid, function(s) {sum(s!=r)})
  foldid2 <- lapply(1:D, function(d) {
    sample(c(rep(1:nfolds, each=ntrain[d] %/% nfolds),
             rep(1:(ntrain[d] %% nfolds), (ntrain[d] %% nfolds)!=0)))})
  
  # semnig models
  cv.semnig1 <- semnig(x=xtrain, y=ytrain,
                       C=lapply(xtrain, function(s) {matrix(rep(1, ncol(s)))}),
                       Z=matrix(1, nrow=D), full.post=TRUE,
                       control=control.semnig)
  cv.semnig2 <- semnig(x=xtrain, y=ytrain, C=C, Z=Z,
                       full.post=TRUE, control=control.semnig)

  # penalized regression models
  cv.ridge1 <- lapply(c(1:D), function(d) {
    cv.glmnet(xtrain[[d]], ytrain[[d]], alpha=0, foldid=foldid2[[d]],
              intercept=FALSE)})
  cv.lasso1 <- lapply(c(1:D), function(d) {
    cv.glmnet(xtrain[[d]], ytrain[[d]], foldid=foldid2[[d]], intercept=FALSE)})

  # fit xtune
  cv.xtune1 <- lapply(c(1:D), function(d) {
    xtune(xtrain[[d]], ytrain[[d]], C[[d]][, -1], family="linear",
          method="ridge", message=FALSE, control=list(intercept=FALSE))})

  # fit ebridge model
  cv.ebridge1 <- ebridge(xtrain, ytrain,
                         lapply(C, function(s) {matrix(s[, -1])}), Z[, -1],
                         foldid=foldid2,
                         hyper=list(lambda=1/sqrt(n*sapply(
                           cv.ridge1, "[[", "lambda.min")), zeta=0, nu=0),
                         control=control.ebridge)
  
  # bSEM model
  cv.bSEM1 <- bSEM(xtrain, ytrain, lapply(C, function(s) {s[, -1] + 1}),
                   control=list(maxit=500, trace=FALSE, epsilon=1e-3))
  
  # post selection
  cv.dss.semnig1 <- lapply(c(1:D), function(d) {
    cv.glmnet(xtrain[[d]], as.numeric(xtrain[[d]] %*% cv.semnig1$vb$mu[[d]]),
              penalty.factor=0.5/abs(cv.semnig1$vb$mu[[d]]), intercept=FALSE)})
  cv.dss.ridge1 <- lapply(c(1:D), function(d) {
    b <- coef(cv.ridge1[[d]], s="lambda.min")[-1, ]
    cv.glmnet(xtrain[[d]], as.numeric(xtrain[[d]] %*% b),
              penalty.factor=0.5/abs(b), intercept=FALSE)})
  cv.dss.bSEM1 <- lapply(c(1:D), function(d) {
    b <- cv.bSEM1$vb$beta[[d]][, "mu"]
    cv.glmnet(xtrain[[d]], as.numeric(xtrain[[d]] %*% b),
              penalty.factor=0.5/abs(b), intercept=FALSE)})
  
  # estimates
  best <- list(semnig1=cv.semnig1$vb$mpost$beta,
               semnig2=cv.semnig2$vb$mpost$beta,
               ridge1=lapply(cv.ridge1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}),
               lasso1=lapply(cv.lasso1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}),
               xtune1=lapply(cv.xtune1, function(s) {
                 unname(s$beta.est[-1, ])}),
               ebridge1=cv.ebridge1$beta1,
               bSEM1=lapply(cv.bSEM1$vb$beta, function(s) {s[, "mu"]}))
  best$dss.semnig1 <- lapply(cv.dss.semnig1, function(s) {
    coef(s, s="lambda.min")[-1, 1]})
  best$dss.semnig2 <- lapply(c(1:D), function(d) {
    coef.dss(cv.dss.semnig1[[d]], psel=sum(best$lasso1[[d]]!=0))[-1, ]})
  best$dss.semnig3 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.semnig1[[d]], psel=psel.dss)[-1, ]})
  best$dss.lasso1 <- lapply(1:D, function(d) {
    coef.dss(cv.lasso1[[d]], psel=psel.dss)[-1, ]})
  best$dss.ridge1 <- lapply(cv.dss.ridge1, function(s) {
    coef(s, s="lambda.min")[-1, 1]})
  best$dss.ridge2 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.ridge1[[d]], psel=sum(best$lasso1[[d]]!=0))[-1, ]})
  best$dss.ridge3 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.ridge1[[d]], psel=psel.dss)[-1, ]})
  best$dss.bSEM1 <- lapply(cv.dss.bSEM1, function(s) {
    coef(s, s="lambda.min")[-1, 1]})
  best$dss.bSEM2 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.bSEM1[[d]], psel=sum(best$lasso1[[d]]!=0))[-1, ]})
  best$dss.bSEM3 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.bSEM1[[d]], psel=psel.dss)[-1, ]})
  
  pred <- lapply(best, function(b) {
    sapply(1:D, function(d) {as.numeric(xtest[[d]] %*% b[[d]])})})
  pmse <- sapply(pred, function(s) {
    sapply(1:D, function(d) {mean((s[[d]] - ytest[[d]])^2)})})
  psel <- sapply(best, function(b) {
    sapply(b, function(s) {sum(s!=0)})})
  
  rownames(pmse) <- paste0("drug", 1:D)
  rownames(psel) <- paste0("drug", 1:D)
  list(pmse=pmse, psel=psel)
}
stopCluster(cl=cl)

pmse <- Reduce("rbind", lapply(res, "[[", "pmse"))
rownames(pmse) <- paste0("pmse.", rownames(pmse))
psel <- Reduce("rbind", lapply(res, "[[", "psel"))
rownames(psel) <- paste0("psel.", rownames(psel))

res <- rbind(pmse, psel)
write.table(res, file="results/analysis_gdsc_res4.txt")







res4.1 <- read.table("results/analysis_gdsc_res4.1.txt", row.names=NULL)
temp <- res4.1[, 1]
res4.1 <- as.matrix(res4.1[, -1])
rownames(res4.1) <- temp

res4.2 <- read.table("results/analysis_gdsc_res4.2.txt", row.names=NULL)
temp <- res4.2[, 1]
res4.2 <- as.matrix(res4.2[, -1])
rownames(res4.2) <- temp

res4.3 <- read.table("results/analysis_gdsc_res4.3.txt", row.names=NULL)
temp <- res4.3[, 1]
res4.3 <- as.matrix(res4.3[, -1])
rownames(res4.3) <- temp

nrow(res4.1)
nrow(res4.2)
nrow(res4.3)

methods <- c("semnig1", "semnig2", "ridge1", "lasso1", "xtune1", "ebridge1", 
             "bSEM1", "dss.semnig1", "dss.semnig2", "dss.semnig3", "dss.lasso1", 
             "dss.ridge1", "dss.ridge2", "dss.ridge3", "dss.bSEM1", "dss.bSEM2", 
             "dss.bSEM3")


res4 <- cbind(res4.1, res4.2, res4.3)
colnames(res4)[17] <- "bSEM1"
res <- res4[, methods]





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





# correlate mutation pvalues with regression coefficients from expression
expr.sel <- expr$expr[, colnames(expr$expr) %in%
                        unlist(strsplit(colnames(mut), ","))]
fit1.ridge <- lapply(1:length(resp), function(d) {
  resp.prep <- scale(resp[[d]][names(resp[[d]]) %in% rownames(expr.sel)])
  expr.prep <- scale(expr.sel[rownames(expr.sel) %in% names(resp[[d]]), ])
  cv.glmnet(expr.prep, resp.prep, alpha=0, intercept=FALSE)})
best1.ridge <- sapply(fit1.ridge, function(s) {coef(s, s="lambda.min")[-1, ]})
best1.ridge <- best1.ridge[rownames(best1.ridge) %in% rownames(mut.pvalues), ]
best1.ridge <- best1.ridge[order(rownames(best1.ridge)), ]
mut.pvalues <- mut.pvalues[rownames(mut.pvalues) %in% rownames(best1.ridge), ]
mut.pvalues <- mut.pvalues[order(rownames(mut.pvalues)), ]
cor.test(mut.pvalues, best1.ridge)