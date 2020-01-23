#!/usr/bin/env Rscript

################################################################################
# regression GDSC log IC50 values on gene expression of cell lines
################################################################################

# number of cores to use
ncores <- 100

### libraries
packages <- c("foreach", "doParallel", "cambridge", "glmnet", "pInc")
sapply(packages, library, character.only=TRUE)

### load data
load(file="data/data_gdsc_dat1.Rdata")

# reassign data
resp <- resp.prep
expr <- expr.prep
meth <- meth.prep

################################## analysis 1 ##################################
### data preparation
# only retain the cell lines that are available in both response and expression
expr.prep <- expr[rownames(expr) %in% rownames(resp), ]
resp.prep <- resp[rownames(resp) %in% rownames(expr), ]

# select features with largest variance separately for in and out pathway genes
D <- ncol(resp.prep)
psel <- 100
o <- order(-apply(expr.prep, 2, sd))
idsel <- lapply(feat.prep$inpathway, function(s) {
  m1 <- o[o %in% which(s==1)]; id1 <- head(m1, n=min(length(m1), psel/2))
  m0 <- o[o %in% which(s==0)]; id0 <- head(m0, n=psel - length(id1))
  c(id1, id0)})
expr.sel <- lapply(idsel, function(s) {expr.prep[, s]})
inpathway <- lapply(1:D, function(d) {feat.prep$inpathway[[d]][idsel[[d]]]})

# create data objects used in fitting
x <- lapply(expr.sel, function(s) {scale(s, scale=FALSE)})
y <- scale(resp.prep, scale=FALSE)
p <- sapply(x, ncol)
n <- nrow(y)

### model fitting
# setting seed
set.seed(2019)

# estimation settings
control <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, epsilon.vb=1e-3,
                maxit.eb=200, maxit.vb=1, maxit.post=100)
methods <- c("NIG1", "NIG2", "NIG3", "NIG4", "lasso", "ridge", "bSEM")

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
fit1.bSEM <- bSEM(x, split(y, 1:ncol(y)), lapply(inpathway, "+", 1),
                  control=list(maxit=200, trace=TRUE, epsilon=1e-3))

# penalized regression models
fit1.lasso <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[, d], intercept=FALSE, standardize=FALSE)})
fit1.ridge <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[, d], alpha=0, intercept=FALSE, standardize=FALSE)})

# saving fitted model objects
save(fit1.semnig, fit2.semnig, fit3.semnig, fit4.semnig, fit1.lasso, fit1.ridge,
     fit1.bSEM, file="results/analysis_gdsc_fit1.Rdata")

# EB estimates
tab <- rbind(c(fit1.semnig$eb$alphaf, NA, fit1.semnig$eb$alphad, rep(NA, 4), 
               fit1.semnig$eb$lambdaf, fit1.semnig$eb$lambdad),
             c(fit2.semnig$eb$alphaf, fit2.semnig$eb$alphad, 
               fit2.semnig$eb$lambdaf, fit2.semnig$eb$lambdad),
             c(fit3.semnig$eb$alphaf, rep(NA, 5), fit3.semnig$eb$lambdaf, NA),
             c(rep(NA, 2), fit4.semnig$eb$alphad, NA, fit2.semnig$eb$lambdad),
             c(fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 1]/
                 fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 1],
               fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 2]/
                 fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 2] -
                 fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 1]/
                 fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 1], rep(NA, 7)))
colnames(tab) <- c(colnames(C[[1]]), colnames(Z), "lambdaf", "lambdad")
rownames(tab) <- methods[c(1:4, 7)]
write.table(tab, file="results/analysis_gdsc_fit1.txt")

### cross-validation of performance measures
# estimation settings
control <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-3, epsilon.vb=1e-3, 
                maxit.eb=200, maxit.vb=1, maxit.post=100)
methods <- c("NIG1", "NIG2", "NIG3", "NIG4", "lasso", "ridge", "bSEM")
nreps <- 50
ntrain <- floor(nrow(resp.prep)/2)

# setup cluster
cl <- makeCluster(ncores) 
if(ncores > 1) {
  registerDoParallel(cl)
} else {
  registerDoSEQ()
}
res <- foreach(r=1:nreps, .packages=packages) %dopar% {
  cat("\r", "rep", r)
  set.seed(2019 + r)
  
  # splitting data
  idtrain <- sample(1:nrow(y), ntrain)
  xtrain <- lapply(x, function(s) {scale(s[idtrain, ])})
  ytrain <- scale(y[idtrain, ], scale=FALSE)
  xtest <- lapply(x, function(s) {scale(s[-idtrain, ])})
  ytest <- scale(y[-idtrain, ], scale=FALSE)
  
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
  cv1.bSEM <- bSEM(xtrain, split(ytrain, 1:ncol(ytrain)),
                   lapply(inpathway, "+", 1),
                   control=list(maxit=200, trace=FALSE, epsilon=1e-3))

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
    mean((ytest - sapply(1:D, function(d) {
      predict(cv1.lasso[[d]], xtest[[d]], s="lambda.min")}))^2),
    mean((ytest - sapply(1:D, function(d) {
      predict(cv1.ridge[[d]], xtest[[d]], s="lambda.min")}))^2),
    mean(sapply(1:D, function(d) {
      mean((ytest[, d] - xtest[[d]] %*% cv1.bSEM$vb$beta[[d]][, 1])^2)})),
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
  
  # determine the ranks of the model parameter point estimates
  brank <- c(unlist(lapply(cv1.semnig$vb$mpost$beta, function(s) {
    rank(abs(s), ties.method="average")})),
    unlist(lapply(cv2.semnig$vb$mpost$beta, function(s) {
      rank(abs(s), ties.method="average")})),
    unlist(lapply(cv3.semnig$vb$mpost$beta, function(s) {
      rank(abs(s), ties.method="average")})),
    unlist(lapply(cv4.semnig$vb$mpost$beta, function(s) {
      rank(abs(s), ties.method="average")})),
    unlist(sapply(cv1.lasso, function(s) {
      rank(abs(as.numeric(coef(s, s="lambda.min"))[-1]), 
           ties.method="average")})),
    unlist(sapply(cv1.ridge, function(s) {
      rank(abs(as.numeric(coef(s, s="lambda.min"))[-1]), 
           ties.method="average")})),
    unlist(lapply(cv1.bSEM$vb$beta, function(s) {
      rank(abs(s[, 1]), ties.method="average")})))
  names(brank) <- paste0(
    "brank.", paste0(rep(methods, each=sum(p)), ".",
                     rep(paste0(rep(1:D, times=p), ".",
                                unlist(sapply(p, function(s) {
                                  return(1:s)}))), length(methods))))

  list(pmse=pmse, elbo=elbo, elbot=elbot, lpml=lpml, brank=brank)
}
stopCluster(cl=cl)

# prepare and save results table
pmse <- t(sapply(res, "[[", "pmse"))
elbo <- t(sapply(res, "[[", "elbo"))
elbot <- t(sapply(res, "[[", "elbot"))
lpml <- t(sapply(res, "[[", "lpml"))
brank <- t(sapply(res, "[[", "brank"))

# Euclidian distance between rank vectors
brankdist <- sapply(methods, function(s) {
  dist(brank[, substr(colnames(brank), 7, nchar(s) + 6)==s])})

# combine tables and save
res <- rbind(pmse, cbind(elbo, NA, NA, NA, NA), cbind(elbot, NA, NA, NA), 
             cbind(lpml, NA, NA, NA, NA), cbind(brankdist, NA))
colnames(res) <- c(methods, "null")
rownames(res) <- c(rep("pmse", nreps), rep("elbo", nreps), rep("elbot", nreps),
                   rep("lpml", nreps), rep("brankdist", nrow(brankdist)))
write.table(res, file="results/analysis_gdsc_res1.txt")


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
x <- lapply(expr.sel, function(s) {scale(s)})
y <- scale(resp.prep, scale=FALSE)

# number of equations, features, and observations
p <- sapply(x, ncol)
n <- nrow(y)

### model fitting
# setting seed
set.seed(2019)

# estimation settings
control <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, epsilon.vb=1e-3,
                maxit.eb=200, maxit.vb=1, maxit.post=100)
methods <- c("NIG1", "NIG2", "NIG3", "NIG4", "lasso", "ridge", "bSEM")

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
fit1.bSEM <- bSEM(x, split(y, 1:ncol(y)), lapply(inpathway, "+", 1),
                  control=list(maxit=200, trace=TRUE, epsilon=1e-3))

# penalized regression models
fit1.lasso <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[, d], intercept=FALSE, standardize=FALSE)})
fit1.ridge <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[, d], alpha=0, intercept=FALSE, standardize=FALSE)})


save(fit1.semnig, fit2.semnig, fit3.semnig, fit4.semnig, fit1.lasso, fit1.ridge,
     fit1.bSEM, file="results/analysis_gdsc_fit2.Rdata")

# EB estimates
tab <- rbind(c(fit1.semnig$eb$alphaf, NA, fit1.semnig$eb$alphad, rep(NA, 4), 
               fit1.semnig$eb$lambdaf, fit1.semnig$eb$lambdad),
             c(fit2.semnig$eb$alphaf, fit2.semnig$eb$alphad, 
               fit2.semnig$eb$lambdaf, fit2.semnig$eb$lambdad),
             c(fit3.semnig$eb$alphaf, rep(NA, 5), fit3.semnig$eb$lambdaf, NA),
             c(rep(NA, 2), fit4.semnig$eb$alphad, NA, fit2.semnig$eb$lambdad),
             c(fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 1]/
                 fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 1],
               fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 2]/
                 fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 2] -
                 fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 1]/
                 fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 1], rep(NA, 7)))
colnames(tab) <- c(colnames(C[[1]]), colnames(Z), "lambdaf", "lambdad")
rownames(tab) <- methods[c(1:4, 7)]
write.table(tab, file="results/analysis_gdsc_fit2.txt")

### cross-validation
# estimation settings
control <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-3, epsilon.vb=1e-3, 
                maxit.eb=200, maxit.vb=1, maxit.post=100)

nreps <- 50
ntrain <- floor(nrow(resp.prep)/2)
methods <- c("NIG1", "NIG2", "NIG3", "NIG4", "lasso", "ridge", "bSEM")

# setup cluster
cl <- makeCluster(ncores) 
if(ncores > 1) {
  registerDoParallel(cl)
} else {
  registerDoSEQ()
}
res <- foreach(r=1:nreps, .packages=packages) %dopar% {
# for(r in 1:nreps) {
  cat("\r", "rep", r)
  set.seed(2019 + r)
  
  # splitting data
  idtrain <- sample(1:nrow(y), ntrain)
  xtrain <- lapply(x, function(s) {scale(s[idtrain, ])})
  ytrain <- scale(y[idtrain, ], scale=FALSE)
  xtest <- lapply(x, function(s) {scale(s[-idtrain, ])})
  ytest <- scale(y[-idtrain, ], scale=FALSE)
  
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
  cv1.bSEM <- bSEM(xtrain, split(ytrain, 1:ncol(ytrain)),
                   lapply(inpathway, "+", 1),
                   control=list(maxit=200, trace=FALSE, epsilon=1e-3))
  
  # lasso and ridge
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
    mean((ytest - sapply(1:D, function(d) {
      predict(cv1.lasso[[d]], xtest[[d]], s="lambda.min")}))^2),
    mean((ytest - sapply(1:D, function(d) {
      predict(cv1.ridge[[d]], xtest[[d]], s="lambda.min")}))^2),
    mean(sapply(1:D, function(d) {
      mean((ytest[, d] - xtest[[d]] %*% cv1.bSEM$vb$beta[[d]][, 1])^2)})),
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
  
  # determine the ranks of the model parameter point estimates
  brank <- c(unlist(lapply(cv1.semnig$vb$mpost$beta, function(s) {
    rank(abs(s), ties.method="average")})),
    unlist(lapply(cv2.semnig$vb$mpost$beta, function(s) {
      rank(abs(s), ties.method="average")})),
    unlist(lapply(cv3.semnig$vb$mpost$beta, function(s) {
      rank(abs(s), ties.method="average")})),
    unlist(lapply(cv4.semnig$vb$mpost$beta, function(s) {
      rank(abs(s), ties.method="average")})),
    unlist(sapply(cv1.lasso, function(s) {
      rank(abs(as.numeric(coef(s, s="lambda.min"))[-1]), 
           ties.method="average")})),
    unlist(sapply(cv1.ridge, function(s) {
      rank(abs(as.numeric(coef(s, s="lambda.min"))[-1]), 
           ties.method="average")})),
    unlist(lapply(cv1.bSEM$vb$beta, function(s) {
      rank(abs(s[, 1]), ties.method="average")})))
  names(brank) <- paste0(
    "brank.", paste0(rep(methods, each=sum(p)), ".",
                     rep(paste0(rep(1:D, times=p), ".",
                                unlist(sapply(p, function(s) {
                                  return(1:s)}))), length(methods))))
  
  list(pmse=pmse, elbo=elbo, elbot=elbot, lpml=lpml, brank=brank)
}
stopCluster(cl=cl)

# prepare and save results table
pmse <- t(sapply(res, "[[", "pmse"))
elbo <- t(sapply(res, "[[", "elbo"))
elbot <- t(sapply(res, "[[", "elbot"))
lpml <- t(sapply(res, "[[", "lpml"))
brank <- t(sapply(res, "[[", "brank"))

# Euclidian distance between rank vectors
brankdist <- sapply(methods, function(s) {
  dist(brank[, substr(colnames(brank), 7, nchar(s) + 6)==s])})

# combine tables and save
res <- rbind(pmse, cbind(elbo, NA, NA, NA, NA), cbind(elbot, NA, NA, NA), 
             cbind(lpml, NA, NA, NA, NA), cbind(brankdist, NA))
colnames(res) <- c(methods, "null")
rownames(res) <- c(rep("pmse", nreps), rep("elbo", nreps), rep("elbot", nreps),
                   rep("lpml", nreps), rep("brankdist", nrow(brankdist)))
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
x <- lapply(expr.sel, function(s) {scale(s)})
y <- scale(resp.prep, scale=FALSE)

# number of equations, features, and observations
p <- sapply(x, ncol)
n <- nrow(y)

### model fitting
# setting seed
set.seed(2019)

# estimation settings
control <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, epsilon.vb=1e-3,
                maxit.eb=200, maxit.vb=1, maxit.post=100)
methods <- c("NIG1", "NIG2", "NIG3", "NIG4", "lasso", "ridge", "bSEM")

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
fit1.bSEM <- bSEM(x, split(y, 1:ncol(y)), lapply(inpathway, "+", 1),
                  control=list(maxit=200, trace=TRUE, epsilon=1e-3))

# penalized regression methods
fit1.lasso <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[, d], intercept=FALSE, standardize=FALSE)})
fit1.ridge <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[, d], alpha=0, intercept=FALSE, standardize=FALSE)})

save(fit1.semnig, fit2.semnig, fit3.semnig, fit4.semnig, fit1.lasso, fit1.ridge,
     file="results/analysis_gdsc_fit3.Rdata")

# EB estimates
tab <- rbind(c(fit1.semnig$eb$alphaf, NA, fit1.semnig$eb$alphad, rep(NA, 4), 
               fit1.semnig$eb$lambdaf, fit1.semnig$eb$lambdad),
             c(fit2.semnig$eb$alphaf, fit2.semnig$eb$alphad, 
               fit2.semnig$eb$lambdaf, fit2.semnig$eb$lambdad),
             c(fit3.semnig$eb$alphaf, rep(NA, 5), fit3.semnig$eb$lambdaf, NA),
             c(rep(NA, 2), fit4.semnig$eb$alphad, NA, fit2.semnig$eb$lambdad),
             c(fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 1]/
                 fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 1],
               fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 2]/
                 fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 2] -
                 fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 1]/
                 fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 1], rep(NA, 7)))
colnames(tab) <- c(colnames(C[[1]]), colnames(Z), "lambdaf", "lambdad")
rownames(tab) <- methods[c(1:4, 7)]
write.table(tab, file="results/analysis_gdsc_fit3.txt")

### cross-validation
# estimation settings
control <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-3, epsilon.vb=1e-3, 
                maxit.eb=200, maxit.vb=1, maxit.post=100)

nreps <- 50
ntrain <- floor(nrow(resp.prep)/2)
methods <- c("NIG1", "NIG2", "NIG3", "NIG4", "lasso", "ridge", "bSEM")

# setup cluster
cl <- makeCluster(ncores) 
if(ncores > 1) {
  registerDoParallel(cl)
} else {
  registerDoSEQ()
}
res <- foreach(r=1:nreps, .packages=packages) %dopar% {
  cat("\r", "rep", r)
  set.seed(2019 + r)
  
  # splitting data
  idtrain <- sample(1:nrow(y), ntrain)
  xtrain <- lapply(x, function(s) {scale(s[idtrain, ])})
  ytrain <- scale(y[idtrain, ], scale=FALSE)
  xtest <- lapply(x, function(s) {scale(s[-idtrain, ])})
  ytest <- scale(y[-idtrain, ], scale=FALSE)
  
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
  cv1.bSEM <- bSEM(xtrain, split(ytrain, 1:ncol(ytrain)),
                   lapply(inpathway, "+", 1),
                   control=list(maxit=200, trace=FALSE, epsilon=1e-3))
  
  # lasso and ridge
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
    mean((ytest - sapply(1:D, function(d) {
      predict(cv1.lasso[[d]], xtest[[d]], s="lambda.min")}))^2),
    mean((ytest - sapply(1:D, function(d) {
      predict(cv1.ridge[[d]], xtest[[d]], s="lambda.min")}))^2),
    mean(sapply(1:D, function(d) {
      mean((ytest[, d] - xtest[[d]] %*% cv1.bSEM$vb$beta[[d]][, 1])^2)})),
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
 
  # determine the ranks of the model parameter point estimates
  brank <- c(unlist(lapply(cv1.semnig$vb$mpost$beta, function(s) {
    rank(abs(s), ties.method="average")})),
    unlist(lapply(cv2.semnig$vb$mpost$beta, function(s) {
      rank(abs(s), ties.method="average")})),
    unlist(lapply(cv3.semnig$vb$mpost$beta, function(s) {
      rank(abs(s), ties.method="average")})),
    unlist(lapply(cv4.semnig$vb$mpost$beta, function(s) {
      rank(abs(s), ties.method="average")})),
    unlist(sapply(cv1.lasso, function(s) {
      rank(abs(as.numeric(coef(s, s="lambda.min"))[-1]), 
           ties.method="average")})),
    unlist(sapply(cv1.ridge, function(s) {
      rank(abs(as.numeric(coef(s, s="lambda.min"))[-1]), 
           ties.method="average")})),
    unlist(lapply(cv1.bSEM$vb$beta, function(s) {
      rank(abs(s[, 1]), ties.method="average")})))
  names(brank) <- paste0(
    "brank.", paste0(rep(methods, each=sum(p)), ".",
                     rep(paste0(rep(1:D, times=p), ".",
                                unlist(sapply(p, function(s) {
                                  return(1:s)}))), length(methods))))
  
  list(pmse=pmse, elbo=elbo, elbot=elbot, lpml=lpml, brank=brank)
}
stopCluster(cl=cl)

# prepare and save results table
pmse <- t(sapply(res, "[[", "pmse"))
elbo <- t(sapply(res, "[[", "elbo"))
elbot <- t(sapply(res, "[[", "elbot"))
lpml <- t(sapply(res, "[[", "lpml"))
brank <- t(sapply(res, "[[", "brank"))

# Euclidian distance between rank vectors
brankdist <- sapply(methods, function(s) {
  dist(brank[, substr(colnames(brank), 7, nchar(s) + 6)==s])})

# combine tables and save
res <- rbind(pmse, cbind(elbo, NA, NA, NA, NA), cbind(elbot, NA, NA, NA), 
             cbind(lpml, NA, NA, NA, NA), cbind(brankdist, NA))
colnames(res) <- c(methods, "null")
rownames(res) <- c(rep("pmse", nreps), rep("elbo", nreps), rep("elbot", nreps),
                   rep("lpml", nreps), rep("brankdist", nrow(brankdist)))
write.table(res, file="results/analysis_gdsc_res3.txt")


################################## analysis 4 ##################################
### data preparation
# only retain the cell lines that are available in both response and methylation
meth.prep <- meth[rownames(meth) %in% rownames(resp), ]
resp.prep <- resp[rownames(resp) %in% rownames(meth), ]

# selecting methylation probes with largest variance
psel <- 100
o <- order(-apply(meth.prep, 2, sd))
idsel <- o[1:psel]
meth.sel <- meth.prep[, idsel]

# transform methylation values
x <- rep(list(scale(log(meth.sel/(1 - meth.sel), base=2), scale=FALSE)), 
         ncol(y))
y <- scale(resp.prep, scale=FALSE)
D <- ncol(y)
p <- sapply(x, ncol)
n <- nrow(y)

### model fitting
# setting seed
set.seed(2020)

# estimation settings
control <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, epsilon.vb=1e-3,
                maxit.eb=200, maxit.vb=1, maxit.post=100)
methods <- c("NIG1", "NIG2", "lasso", "ridge", "bSEM")

# model without external covariates
Z <- matrix(1, nrow=D)
colnames(Z) <- c("intercept")
C <- lapply(p, function(s) {
  s <- matrix(1, nrow=s); colnames(s) <- c("intercept"); return(s)})
fit1.semnig <- semnig(x=x, y=y, C=C, Z=Z, unpenalized=NULL,
                      standardize=FALSE, intercept=FALSE, fixed.eb="none",
                      full.post=TRUE, init=NULL, control=control)

# models with external covariates
Z <- model.matrix(~ replace(drug.prep$stage, is.na(drug.prep$stage),
                            "experimental") + replace(drug.prep$action, is.na(
                              drug.prep$action), "unknown"))
colnames(Z) <- c("intercept", "experimental", "in clinical development",
                 "targeted", "unknown")
C <- lapply(p, function(s) {
  s <- matrix(1, nrow=s); colnames(s) <- c("intercept"); return(s)})
fit2.semnig <- semnig(x=x, y=y, C=C, Z=Z, unpenalized=NULL,
                      standardize=FALSE, intercept=FALSE, fixed.eb=FALSE,
                      full.post=TRUE, init=NULL, control=control)

# penalized regression methods
fit1.lasso <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[, d], intercept=FALSE, standardize=FALSE)})
fit1.ridge <- lapply(1:D, function(d) {
  cv.glmnet(x[[d]], y[, d], alpha=0, intercept=FALSE, standardize=FALSE)})

save(fit1.semnig, fit2.semnig, fit1.lasso, fit1.ridge,
     file="results/analysis_gdsc_fit4.Rdata")

# EB estimates
tab <- rbind(c(fit1.semnig$eb$alphaf, fit1.semnig$eb$alphad, rep(NA, 4), 
               fit1.semnig$eb$lambdaf, fit1.semnig$eb$lambdad),
             c(fit2.semnig$eb$alphaf, fit2.semnig$eb$alphad, 
               fit2.semnig$eb$lambdaf, fit2.semnig$eb$lambdad))
colnames(tab) <- c(colnames(C[[1]]), colnames(Z), "lambdaf", "lambdad")
rownames(tab) <- methods[c(1:2)]
write.table(tab, file="results/analysis_gdsc_fit4.txt")

### cross-validation
# estimation settings
control <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-3, epsilon.vb=1e-3, 
                maxit.eb=200, maxit.vb=1, maxit.post=100)

nreps <- 50
ntrain <- floor(nrow(resp.prep)/2)
methods <- c("NIG1", "NIG2", "lasso", "ridge")

# setup cluster
cl <- makeCluster(ncores) 
if(ncores > 1) {
  registerDoParallel(cl)
} else {
  registerDoSEQ()
}
res <- foreach(r=1:nreps, .packages=packages) %dopar% {
  cat("\r", "rep", r)
  set.seed(2019 + r)
  
  # splitting data
  idtrain <- sample(1:nrow(y), ntrain)
  xtrain <- lapply(x, function(s) {scale(s[idtrain, ])})
  ytrain <- scale(y[idtrain, ], scale=FALSE)
  xtest <- lapply(x, function(s) {scale(s[-idtrain, ])})
  ytest <- scale(y[-idtrain, ], scale=FALSE)
  
  # model without external covariates
  Z <- matrix(1, nrow=D)
  colnames(Z) <- c("intercept")
  C <- lapply(p, function(s) {
    s <- matrix(1, nrow=s); colnames(s) <- c("intercept"); return(s)})
  cv1.semnig <- semnig(x=xtrain, y=ytrain, C=C, Z=Z, unpenalized=NULL, 
                       standardize=FALSE, intercept=FALSE, fixed.eb="none", 
                       full.post=TRUE, init=NULL, control=control)
  
  # models with external covariates
  Z <- model.matrix(~ replace(drug.prep$stage, is.na(drug.prep$stage),
                              "experimental") + replace(drug.prep$action, is.na(
                                drug.prep$action), "unknown"))
  colnames(Z) <- c("intercept", "experimental", "in clinical development",
                   "targeted", "unknown")
  C <- lapply(p, function(s) {
    s <- matrix(1, nrow=s); colnames(s) <- c("intercept"); return(s)})
  cv2.semnig <- semnig(x=xtrain, y=ytrain, C=C, Z=Z, unpenalized=NULL, 
                       standardize=FALSE, intercept=FALSE, fixed.eb="none", 
                       full.post=TRUE, init=NULL, control=control)
  
  # # bSEM model
  # cv1.bSEM <- bSEM(xtrain, split(ytrain, 1:ncol(ytrain)),
  #                  lapply(inpathway, "+", 1),
  #                  control=list(maxit=200, trace=FALSE, epsilon=1e-3))
  
  # lasso and ridge
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
    # mean(sapply(1:D, function(d) {
    #   mean((ytest[, d] - xtest[[d]] %*% cv3.semnig$vb$mpost$beta[[d]])^2)})),
    # mean(sapply(1:D, function(d) {
    #   mean((ytest[, d] - xtest[[d]] %*% cv4.semnig$vb$mpost$beta[[d]])^2)})),
    mean((ytest - sapply(1:D, function(d) {
      predict(cv1.lasso[[d]], xtest[[d]], s="lambda.min")}))^2),
    mean((ytest - sapply(1:D, function(d) {
      predict(cv1.ridge[[d]], xtest[[d]], s="lambda.min")}))^2),
    # mean(sapply(1:D, function(d) {
    #   mean((ytest[, d] - xtest[[d]] %*% cv1.bSEM$vb$beta[[d]][, 1])^2)})),
    mean(apply(ytest, 1, "-", colMeans(ytrain))^2))
  
  # calculate ELBO for semnig models on training and test data
  elbot <- c(mean(cv1.semnig$seq.elbo[nrow(cv1.semnig$seq.elbo), ]),
             mean(cv2.semnig$seq.elbo[nrow(cv2.semnig$seq.elbo), ])
             # ,
             # mean(cv3.semnig$seq.elbo[nrow(cv3.semnig$seq.elbo), ]),
             # mean(cv4.semnig$seq.elbo[nrow(cv4.semnig$seq.elbo), ]),
             # mean(cv1.bSEM$mll[, ncol(cv1.bSEM$mll)])
             )
  elbo <- c(mean(new.elbo(cv1.semnig, xtest, ytest)),
            mean(new.elbo(cv2.semnig, xtest, ytest))
            # ,
            # mean(new.elbo(cv3.semnig, xtest, ytest)),
            # mean(new.elbo(cv4.semnig, xtest, ytest))
            )
  lpml <- c(sum(logcpo(xtest, ytest, ntrain, cv1.semnig), na.rm=TRUE),
            sum(logcpo(xtest, ytest, ntrain, cv2.semnig), na.rm=TRUE)
            # ,
            # sum(logcpo(xtest, ytest, ntrain, cv3.semnig), na.rm=TRUE),
            # sum(logcpo(xtest, ytest, ntrain, cv4.semnig), na.rm=TRUE)
            )
  
  # determine the ranks of the model parameter point estimates
  brank <- c(unlist(lapply(cv1.semnig$vb$mpost$beta, function(s) {
    rank(abs(s), ties.method="average")})),
    unlist(lapply(cv2.semnig$vb$mpost$beta, function(s) {
      rank(abs(s), ties.method="average")})),
    # unlist(lapply(cv3.semnig$vb$mpost$beta, function(s) {
    #   rank(abs(s), ties.method="average")})),
    # unlist(lapply(cv4.semnig$vb$mpost$beta, function(s) {
    #   rank(abs(s), ties.method="average")})),
    unlist(sapply(cv1.lasso, function(s) {
      rank(abs(as.numeric(coef(s, s="lambda.min"))[-1]), 
           ties.method="average")})),
    unlist(sapply(cv1.ridge, function(s) {
      rank(abs(as.numeric(coef(s, s="lambda.min"))[-1]), 
           ties.method="average")}))
    # ,
    # unlist(lapply(cv1.bSEM$vb$beta, function(s) {
    #   rank(abs(s[, 1]), ties.method="average")}))
    )
  names(brank) <- paste0(
    "brank.", paste0(rep(methods[c()], each=sum(p)), ".",
                     rep(paste0(rep(1:D, times=p), ".",
                                unlist(sapply(p, function(s) {
                                  return(1:s)}))), length(methods))))
  
  list(pmse=pmse, elbo=elbo, elbot=elbot, lpml=lpml, brank=brank)
}
stopCluster(cl=cl)

# prepare and save results table
pmse <- t(sapply(res, "[[", "pmse"))
elbo <- t(sapply(res, "[[", "elbo"))
elbot <- t(sapply(res, "[[", "elbot"))
lpml <- t(sapply(res, "[[", "lpml"))
brank <- t(sapply(res, "[[", "brank"))

# Euclidian distance between rank vectors
brankdist <- sapply(methods, function(s) {
  dist(brank[, substr(colnames(brank), 7, nchar(s) + 6)==s])})

# combine tables and save
res <- rbind(pmse, cbind(elbo, NA, NA, NA), cbind(elbot, NA, NA, NA), 
             cbind(lpml, NA, NA, NA), cbind(brankdist, NA))
colnames(res) <- c(methods, "null")
rownames(res) <- c(rep("pmse", nreps), rep("elbo", nreps), rep("elbot", nreps),
                   rep("lpml", nreps), rep("brankdist", nrow(brankdist)))
write.table(res, file="results/analysis_gdsc_res4.txt")


