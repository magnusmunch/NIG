#!/usr/bin/env Rscript

################################################################################
# regression GDSC log IC50 values on gene expression/mutations of cell lines
################################################################################

### libraries
packages <- c("foreach", "doParallel", "glmnet", "NIG")
sapply(packages, library, character.only=TRUE)

# number of CV folds
nfolds <- 10

# number of cores to use
ncores <- 1
# ncores <- min(100, nfolds)

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
  s <- cbind(1, log(pvalues[[d]])); colnames(s) <- c("intercept", "pvalue"); s})

### model fitting
# settings
control <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, 
                epsilon.vb=1e-3, maxit.eb=500, maxit.vb=1, 
                maxit.post=100, maxit.block=0)

# creating folds
set.seed(2020)
foldid <- lapply(1:D, function(d) {
  sample(c(rep(1:nfolds, each=n[[d]] %/% nfolds),
           rep(1:(n[[d]] %% nfolds), (n[[d]] %% nfolds)!=0)))})

# nig models
fit.nig1 <- nig(x=x, y=y, C=sapply(x, function(s) {matrix(rep(1, ncol(s)))}), 
                Z=matrix(1, nrow=D), full.post=TRUE, control=control)
fit.nig2 <- nig(x=x, y=y, C=C, Z=matrix(1, nrow=D), full.post=TRUE, 
                control=control)

# penalized regression models
fit.ridge1 <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[[d]], alpha=0, foldid=foldid[[d]], intercept=FALSE)})
fit.lasso1 <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[[d]], foldid=foldid[[d]], intercept=FALSE)})

save(fit.nig1, fit.nig2, fit.ridge1, fit.lasso1,
     file="results/analysis_gdsc_fit1.Rdata")

# eb estimates
est <- cbind(c(fit.nig1$eb$alphaf, NA, fit.nig1$eb$alphad, 
               fit.nig1$eb$lambdaf, fit.nig1$eb$lambdad),
             c(fit.nig2$eb$alphaf, fit.nig2$eb$alphad, 
               fit.nig2$eb$lambdaf, fit.nig2$eb$lambdad),
             rep(NA, 5), rep(NA, 5))
rownames(est) <- c("interceptf", "pvalue", "interceptd", "lambdaf", "lambdad")
colnames(est) <- c("nig1", "nig2", "ridge1", "lasso1")
write.table(est, file="results/analysis_gdsc_fit1.txt")

# log cpo values
lcpo.nig1 <- logcpo(x, y, n, fit.nig1)
lcpo.nig2 <- logcpo(x, y, n, fit.nig2)
save(lcpo.nig1, lcpo.nig2, y, file="results/analysis_gdsc_cpo1.Rdata")
load(file="results/analysis_gdsc_cpo1.Rdata")

### predictions and mse
# settings
psel.dss <- c(25, 100)
control <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-3, epsilon.vb=1e-3, 
                maxit.eb=500, maxit.vb=1, maxit.post=100, maxit.block=0)

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
    s <- cbind(1, log(pvalues[[d]])); colnames(s) <- c("intercept", "pvalue"); s})
  
  # nig models
  cv.nig1 <- nig(x=xtrain, y=ytrain, 
                 C=sapply(xtrain, function(s) {matrix(rep(1, ncol(s)))}), 
                 Z=matrix(1, nrow=D), full.post=TRUE, control=control)
  cv.nig2 <- nig(x=xtrain, y=ytrain, C=C, Z=matrix(1, nrow=D), full.post=TRUE, 
                 control=control)
  
  # penalized regression models
  cv.ridge1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[[d]], alpha=0, foldid=foldid2[[d]], 
              intercept=FALSE)})
  cv.lasso1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[[d]], foldid=foldid2[[d]], intercept=FALSE)})
  
  # post selection
  cv.dss.nig1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], as.numeric(xtrain[[d]] %*% cv.nig1$vb$mu[[d]]),
              penalty.factor=0.5/abs(cv.nig1$vb$mu[[d]]), intercept=FALSE)})
  cv.dss.nig2 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], as.numeric(xtrain[[d]] %*% cv.nig2$vb$mu[[d]]),
              penalty.factor=0.5/abs(cv.nig2$vb$mu[[d]]), intercept=FALSE)})
  
  # estimates
  best <- list(nig1=cv.nig1$vb$mpost$beta, 
               nig2=cv.nig2$vb$mpost$beta,
               ridge1=lapply(cv.ridge1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}),
               lasso1=lapply(cv.lasso1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}))
  best$dss.nig1.1 <- lapply(cv.dss.nig1, function(s) {
    coef(s, s="lambda.min")[-1, 1]})
  best$dss.nig1.2 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.nig1[[d]], psel=sum(best$lasso1[[d]]!=0))[-1, ]})
  best$dss.nig1.3 <-  lapply(1:D, function(d) {
    coef.dss(cv.dss.nig1[[d]], psel=psel.dss[1])[-1, ]})
  best$dss.nig1.4 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.nig1[[d]], psel=psel.dss[2])[-1, ]})
  best$dss.nig2.1 <- lapply(cv.dss.nig2, function(s) {
    coef(s, s="lambda.min")[-1, 1]})
  best$dss.nig2.2 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.nig2[[d]], psel=sum(best$lasso1[[d]]!=0))[-1, ]})
  best$dss.nig2.3 <-  lapply(1:D, function(d) {
    coef.dss(cv.dss.nig2[[d]], psel=psel.dss[1])[-1, ]})
  best$dss.nig2.4 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.nig2[[d]], psel=psel.dss[2])[-1, ]})
  best$dss.lasso1.1 <- lapply(1:D, function(d) {
    coef.dss(cv.lasso1[[d]], psel=psel.dss[1])[-1, ]})
  best$dss.lasso1.2 <- lapply(1:D, function(d) {
    coef.dss(cv.lasso1[[d]], psel=psel.dss[2])[-1, ]})
  
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
control <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, 
                epsilon.vb=1e-3, maxit.eb=500, maxit.vb=1, 
                maxit.post=100, maxit.block=0)

# creating folds
set.seed(2020)
foldid <- lapply(1:D, function(d) {
  sample(c(rep(1:nfolds, each=n[[d]] %/% nfolds),
           rep(1:(n[[d]] %% nfolds), (n[[d]] %% nfolds)!=0)))})

# nig models
fit.nig1 <- nig(x=x, y=y, C=lapply(x, function(s) {matrix(rep(1, ncol(s)))}), 
                Z=matrix(1, nrow=D), full.post=TRUE, control=control)
fit.nig2 <- nig(x=x, y=y, C=C, Z=matrix(1, nrow=D), full.post=TRUE, 
                control=control)

# penalized regression models
fit.ridge1 <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[[d]], alpha=0, foldid=foldid[[d]], intercept=FALSE)})
fit.lasso1 <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[[d]], foldid=foldid[[d]], intercept=FALSE)})

# saving fitted model objects
save(fit.nig1, fit.nig2, fit.lasso1, fit.ridge1, 
     file="results/analysis_gdsc_fit2.Rdata")

est <- cbind(c(fit.nig1$eb$alphaf, NA, fit.nig1$eb$alphad, 
               fit.nig1$eb$lambdaf, fit.nig1$eb$lambdad),
             c(fit.nig2$eb$alphaf, fit.nig2$eb$alphad, 
               fit.nig2$eb$lambdaf, fit.nig2$eb$lambdad),
             rep(NA, 5), rep(NA, 5))
rownames(est) <- c("interceptf", "mutation", "interceptd", "lambdaf", "lambdad")
colnames(est) <- c("nig1", "nig2", "ridge1", "lasso1")
write.table(est, file="results/analysis_gdsc_fit2.txt")

# log cpo values
lcpo.nig1 <- logcpo(x, y, n, fit.nig1)
lcpo.nig2 <- logcpo(x, y, n, fit.nig2)
save(lcpo.nig1, lcpo.nig2, file="results/analysis_gdsc_cpo2.Rdata")

# predictions and mse
### predictions and mse
# settings
psel.dss <- c(25, 100)
control <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-3, epsilon.vb=1e-3, 
                maxit.eb=500, maxit.vb=1, maxit.post=100, maxit.block=0)

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
  
  # nig models
  cv.nig1 <- nig(x=xtrain, y=ytrain, 
                       C=sapply(xtrain, function(s) {matrix(rep(1, ncol(s)))}), 
                       Z=matrix(1, nrow=D), full.post=TRUE, 
                       control=control)
  cv.nig2 <- nig(x=xtrain, y=ytrain, C=C, Z=matrix(1, nrow=D),
                       full.post=TRUE, control=control)
  
  # penalized regression models
  cv.ridge1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[[d]], alpha=0, foldid=foldid2[[d]], 
              intercept=FALSE)})
  cv.lasso1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[[d]], foldid=foldid2[[d]], intercept=FALSE)})
  
  # post selection
  cv.dss.nig1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], as.numeric(xtrain[[d]] %*% cv.nig1$vb$mu[[d]]),
              penalty.factor=0.5/abs(cv.nig1$vb$mu[[d]]), intercept=FALSE)})
  cv.dss.nig2 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], as.numeric(xtrain[[d]] %*% cv.nig2$vb$mu[[d]]),
              penalty.factor=0.5/abs(cv.nig2$vb$mu[[d]]), intercept=FALSE)})
  
  # estimates
  best <- list(nig1=cv.nig1$vb$mpost$beta, 
               nig2=cv.nig2$vb$mpost$beta,
               ridge1=lapply(cv.ridge1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}),
               lasso1=lapply(cv.lasso1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}))
  best$dss.nig1.1 <- lapply(cv.dss.nig1, function(s) {
    coef(s, s="lambda.min")[-1, 1]})
  best$dss.nig1.2 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.nig1[[d]], psel=sum(best$lasso1[[d]]!=0))[-1, ]})
  best$dss.nig1.3 <-  lapply(1:D, function(d) {
    coef.dss(cv.dss.nig1[[d]], psel=psel.dss[1])[-1, ]})
  best$dss.nig1.4 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.nig1[[d]], psel=psel.dss[2])[-1, ]})
  best$dss.nig2.1 <- lapply(cv.dss.nig2, function(s) {
    coef(s, s="lambda.min")[-1, 1]})
  best$dss.nig2.2 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.nig2[[d]], psel=sum(best$lasso1[[d]]!=0))[-1, ]})
  best$dss.nig2.3 <-  lapply(1:D, function(d) {
    coef.dss(cv.dss.nig2[[d]], psel=psel.dss[1])[-1, ]})
  best$dss.nig2.4 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.nig2[[d]], psel=psel.dss[2])[-1, ]})
  best$dss.lasso1.1 <- lapply(1:D, function(d) {
    coef.dss(cv.lasso1[[d]], psel=psel.dss[1])[-1, ]})
  best$dss.lasso1.2 <- lapply(1:D, function(d) {
    coef.dss(cv.lasso1[[d]], psel=psel.dss[2])[-1, ]})
  
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
control <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, 
                epsilon.vb=1e-3, maxit.eb=500, maxit.vb=1, 
                maxit.post=100, maxit.block=0)

# creating folds
set.seed(2020)
foldid <- lapply(1:D, function(d) {
  sample(c(rep(1:nfolds, each=n[[d]] %/% nfolds),
           rep(1:(n[[d]] %% nfolds), (n[[d]] %% nfolds)!=0)))})

# nig models
fit.nig1 <- nig(x=x, y=y, C=lapply(x, function(s) {matrix(rep(1, ncol(s)))}), 
                Z=matrix(1, nrow=D), full.post=TRUE, control=control)
fit.nig2 <- nig(x=x, y=y, C=C, Z=matrix(1, nrow=D), full.post=TRUE, 
                control=control)

# penalized regression models
fit.ridge1 <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[[d]], alpha=0, foldid=foldid[[d]], intercept=FALSE)})
fit.lasso1 <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[[d]], foldid=foldid[[d]], intercept=FALSE)})

# saving fitted model objects
save(fit.nig1, fit.nig2, fit.lasso1, fit.ridge1,
     file="results/analysis_gdsc_fit3.Rdata")

est <- cbind(c(fit.nig1$eb$alphaf, NA, fit.nig1$eb$alphad, 
               fit.nig1$eb$lambdaf, fit.nig1$eb$lambdad),
             c(fit.nig2$eb$alphaf, fit.nig2$eb$alphad, 
               fit.nig2$eb$lambdaf, fit.nig2$eb$lambdad),
             rep(NA, 5), rep(NA, 5))
rownames(est) <- c("interceptf", "pvalue", "interceptd", "lambdaf", "lambdad")
colnames(est) <- c("nig1", "nig2", "ridge1", "lasso1")
write.table(est, file="results/analysis_gdsc_fit3.txt")

# log cpo values
lcpo.nig1 <- logcpo(x, y, n, fit.nig1)
lcpo.nig2 <- logcpo(x, y, n, fit.nig2)
save(lcpo.nig1, lcpo.nig2, file="results/analysis_gdsc_cpo3.Rdata")

### predictions and mse
# settings
psel.dss <- c(25, 100)
control <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-3, 
                epsilon.vb=1e-3, maxit.eb=500, maxit.vb=1, 
                maxit.post=100, maxit.block=0)

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
  
  # nig models
  cv.nig1 <- nig(x=xtrain, y=ytrain, 
                 C=lapply(xtrain, function(s) {matrix(rep(1, ncol(s)))}), 
                 Z=matrix(1, nrow=D), full.post=TRUE, control=control)
  cv.nig2 <- nig(x=xtrain, y=ytrain, C=C, Z=matrix(1, nrow=D),
                 full.post=TRUE, control=control)
  
  # penalized regression models
  cv.ridge1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[[d]], alpha=0, foldid=foldid2[[d]], 
              intercept=FALSE)})
  cv.lasso1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[[d]], foldid=foldid2[[d]], intercept=FALSE)})
  
  # post selection
  cv.dss.nig1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], as.numeric(xtrain[[d]] %*% cv.nig1$vb$mu[[d]]),
              penalty.factor=0.5/abs(cv.nig1$vb$mu[[d]]), intercept=FALSE)})
  cv.dss.nig2 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], as.numeric(xtrain[[d]] %*% cv.nig2$vb$mu[[d]]),
              penalty.factor=0.5/abs(cv.nig2$vb$mu[[d]]), intercept=FALSE)})
  
  # estimates
  best <- list(nig1=cv.nig1$vb$mpost$beta, 
               nig2=cv.nig2$vb$mpost$beta,
               ridge1=lapply(cv.ridge1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}),
               lasso1=lapply(cv.lasso1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}))
  best$dss.nig1.1 <- lapply(cv.dss.nig1, function(s) {
    coef(s, s="lambda.min")[-1, 1]})
  best$dss.nig1.2 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.nig1[[d]], psel=sum(best$lasso1[[d]]!=0))[-1, ]})
  best$dss.nig1.3 <-  lapply(1:D, function(d) {
    coef.dss(cv.dss.nig1[[d]], psel=psel.dss[1])[-1, ]})
  best$dss.nig1.4 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.nig1[[d]], psel=psel.dss[2])[-1, ]})
  best$dss.nig2.1 <- lapply(cv.dss.nig2, function(s) {
    coef(s, s="lambda.min")[-1, 1]})
  best$dss.nig2.2 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.nig2[[d]], psel=sum(best$lasso1[[d]]!=0))[-1, ]})
  best$dss.nig2.3 <-  lapply(1:D, function(d) {
    coef.dss(cv.dss.nig2[[d]], psel=psel.dss[1])[-1, ]})
  best$dss.nig2.4 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.nig2[[d]], psel=psel.dss[2])[-1, ]})
  best$dss.lasso1.1 <- lapply(1:D, function(d) {
    coef.dss(cv.lasso1[[d]], psel=psel.dss[1])[-1, ]})
  best$dss.lasso1.2 <- lapply(1:D, function(d) {
    coef.dss(cv.lasso1[[d]], psel=psel.dss[2])[-1, ]})
  
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
control <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, 
                epsilon.vb=1e-3, maxit.eb=500, maxit.vb=1, 
                maxit.post=100, maxit.block=0)

# creating folds
set.seed(2020)
foldid <- lapply(1:D, function(d) {
  sample(c(rep(1:nfolds, each=n[[d]] %/% nfolds),
           rep(1:(n[[d]] %% nfolds), (n[[d]] %% nfolds)!=0)))})

# nig models
fit.nig1 <- nig(x=x, y=y, C=lapply(x, function(s) {matrix(rep(1, ncol(s)))}), 
                Z=matrix(1, nrow=D), full.post=TRUE, control=control)
fit.nig2 <- nig(x=x, y=y, C=C, Z=Z, full.post=TRUE, control=control)

# penalized regression models
fit.ridge1 <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[[d]], alpha=0, foldid=foldid[[d]], intercept=FALSE)})
fit.lasso1 <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[[d]], foldid=foldid[[d]], intercept=FALSE)})

# saving fitted model objects
save(fit.nig1, fit.nig2, fit.lasso1, fit.ridge1, 
     file="results/analysis_gdsc_fit4.Rdata")

est <- cbind(c(fit.nig1$eb$alphaf, NA, fit.nig1$eb$alphad, rep(NA, 4),
               fit.nig1$eb$lambdaf, fit.nig1$eb$lambdad),
             c(fit.nig2$eb$alphaf, fit.nig2$eb$alphad, 
               fit.nig2$eb$lambdaf, fit.nig2$eb$lambdad),
             rep(NA, 9), rep(NA, 9))
rownames(est) <- c("interceptf", "inpathway", "interceptd", "experimental", 
                   "in clinical development", "targeted", "unknown", "lambdaf", 
                   "lambdad")
colnames(est) <- c("nig1", "nig2", "ridge1", "lasso1")
write.table(est, file="results/analysis_gdsc_fit4.txt")

# log cpo values
lcpo.nig1 <- logcpo(x, y, n, fit.nig1)
lcpo.nig2 <- logcpo(x, y, n, fit.nig2)
save(lcpo.nig1, lcpo.nig2, file="results/analysis_gdsc_cpo4.Rdata")

# predictions and mse
### predictions and mse
# settings
psel.dss <- c(25, 100)
control <- list(conv.post=FALSE, trace=FALSE, epsilon.eb=1e-3, 
                epsilon.vb=1e-3, maxit.eb=500, maxit.vb=1, 
                maxit.post=100, maxit.block=0)

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
  
  # nig models
  cv.nig1 <- nig(x=xtrain, y=ytrain,
                 C=lapply(xtrain, function(s) {matrix(rep(1, ncol(s)))}),
                 Z=matrix(1, nrow=D), full.post=TRUE, control=control)
  cv.nig2 <- nig(x=xtrain, y=ytrain, C=C, Z=Z, full.post=TRUE, control=control)

  # penalized regression models
  cv.ridge1 <- lapply(c(1:D), function(d) {
    cv.glmnet(xtrain[[d]], ytrain[[d]], alpha=0, foldid=foldid2[[d]],
              intercept=FALSE)})
  cv.lasso1 <- lapply(c(1:D), function(d) {
    cv.glmnet(xtrain[[d]], ytrain[[d]], foldid=foldid2[[d]], intercept=FALSE)})

  # post selection
  cv.dss.nig1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], as.numeric(xtrain[[d]] %*% cv.nig1$vb$mu[[d]]),
              penalty.factor=0.5/abs(cv.nig1$vb$mu[[d]]), intercept=FALSE)})
  cv.dss.nig2 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], as.numeric(xtrain[[d]] %*% cv.nig2$vb$mu[[d]]),
              penalty.factor=0.5/abs(cv.nig2$vb$mu[[d]]), intercept=FALSE)})
  
  # estimates
  best <- list(nig1=cv.nig1$vb$mpost$beta, 
               nig2=cv.nig2$vb$mpost$beta,
               ridge1=lapply(cv.ridge1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}),
               lasso1=lapply(cv.lasso1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}))
  best$dss.nig1.1 <- lapply(cv.dss.nig1, function(s) {
    coef(s, s="lambda.min")[-1, 1]})
  best$dss.nig1.2 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.nig1[[d]], psel=sum(best$lasso1[[d]]!=0))[-1, ]})
  best$dss.nig1.3 <-  lapply(1:D, function(d) {
    coef.dss(cv.dss.nig1[[d]], psel=psel.dss[1])[-1, ]})
  best$dss.nig1.4 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.nig1[[d]], psel=psel.dss[2])[-1, ]})
  best$dss.nig2.1 <- lapply(cv.dss.nig2, function(s) {
    coef(s, s="lambda.min")[-1, 1]})
  best$dss.nig2.2 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.nig2[[d]], psel=sum(best$lasso1[[d]]!=0))[-1, ]})
  best$dss.nig2.3 <-  lapply(1:D, function(d) {
    coef.dss(cv.dss.nig2[[d]], psel=psel.dss[1])[-1, ]})
  best$dss.nig2.4 <- lapply(1:D, function(d) {
    coef.dss(cv.dss.nig2[[d]], psel=psel.dss[2])[-1, ]})
  best$dss.lasso1.1 <- lapply(1:D, function(d) {
    coef.dss(cv.lasso1[[d]], psel=psel.dss[1])[-1, ]})
  best$dss.lasso1.2 <- lapply(1:D, function(d) {
    coef.dss(cv.lasso1[[d]], psel=psel.dss[2])[-1, ]})
  
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

################################################################################
############################### computation times ##############################
################################################################################
load(file="results/analysis_gdsc_fit1.Rdata")
time1 <- sum(fit.nig2$time)
load(file="results/analysis_gdsc_fit2.Rdata")
time2 <- sum(fit.nig2$time)
load(file="results/analysis_gdsc_fit3.Rdata")
time3 <- sum(fit.nig2$time)
load(file="results/analysis_gdsc_fit4.Rdata")
time4 <- sum(fit.nig2$time)
write.table(c(time1, time2, time3, time4), 
            file="results/analysis_gdsc_time1.txt")
