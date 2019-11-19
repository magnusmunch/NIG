#!/usr/bin/env Rscript

### installation of package
# if(!("cambridge" %in% installed.packages())) {
#   if(!("devtools" %in% installed.packages())) {
#     install.packages("devtools")
#   }
#   library(devtools)
#   install_github("magnusmunch/cambridge/code", local=FALSE, 
#                  auth_token=Sys.getenv("GITHUB_PAT"))
# }

### libraries
if("cambridge" %in% (.packages())){
  detach("package:cambridge", unload=TRUE)
}
if("pInc" %in% (.packages())){
  detach("package:pInc", unload=TRUE)
}
library(cambridge)
library(pInc)
library(glmnet)

################################### analysis ###################################
# load prepared data
load(file="data/data_gdsc_dat1.Rdata")


psel <- 500
o <- order(-apply(expr.prep, 2, sd))
# idsel <- lapply(inpathway, function(s) {
#   m1 <- o[o %in% which(s==1)]; id1 <- head(m1, n=min(length(m1), psel/2))
#   m0 <- o[o %in% which(s==0)]; id0 <- head(m0, n=psel - length(id1))
#   c(id1, id0)})
# expr.sel <- lapply(idsel, function(s) {expr.prep[, s]})
idsel <- o[1:psel]
expr.sel <- lapply(1:D, function(d) {expr.prep[, idsel]})
inpathway <- lapply(inpathway, function(s) {s[idsel]})

# number of equations, features, and observations
D <- ncol(resp.prep)
p <- sapply(expr.sel, ncol)
n <- nrow(resp.prep)



### model fitting
# estimation settings
control <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, epsilon.vb=1e-3,
                maxit.eb=200, maxit.vb=1, maxit.post=100)

# estimation
x <- lapply(expr.sel, function(s) {scale(s)})
y <- scale(resp.prep)
Z <- matrix(1, nrow=D)
colnames(Z) <- c("intercept")
C <- lapply(inpathway, function(s) {
  s <- matrix(1, nrow=length(s)); colnames(s) <- c("intercept"); return(s)})
fit1.semnig <- semnig(x=x, y=y, C=C, Z=Z, unpenalized=NULL,
                      standardize=FALSE, intercept=FALSE, fixed.eb="none",
                      full.post=TRUE, init=NULL, control=control)
# create co-data matrices
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
# fit1.bSEM <- bSEM(rep(list(x), D), split(y, 1:ncol(y)),
#                   lapply(inpathway, "+", 1),
#                   control=list(maxit=200, trace=TRUE, epsilon=1e-3))

# save(C, Z, fit1.semnig, fit2.semnig, fit1.bSEM,
#      file="results/data_gdsc_fit2.Rdata")
load("results/data_gdsc_fit2.Rdata")

### cross-validation
# estimation settings
control <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-3, epsilon.vb=1e-3, 
                maxit.eb=200, maxit.vb=1, maxit.post=100)

# multiple splits
nreps <- 50
ntrain <- floor(nrow(resp.prep)/2)
# apmse <- matrix(NA, nrow=nreps, ncol=5)
# colnames(apmse) <- c("NIG1", "NIG2", "bSEM", "lasso", "ridge")
# brank <- matrix(NA, nrow=nreps, ncol=sum(p)*5)
# colnames(brank) <- paste0(rep(c("NIG1", "NIG2", "bSEM", "lasso", "ridge"), 
#                               each=sum(p)), ".", rep(paste0(
#                                 rep(1:D, times=p), ".",
#                                 unlist(sapply(p, function(s) {
#                                   return(1:s)}))), 5))
# elbo <- matrix(NA, ncol=3, nrow=nreps)
# colnames(elbo) <- c("NIG1", "NIG2", "bSEM")
apmse <- matrix(NA, nrow=nreps, ncol=4)
colnames(apmse) <- c("NIG1", "NIG2", "NIG3", "NIG4")
brank <- matrix(NA, nrow=nreps, ncol=sum(p)*4)
colnames(brank) <- paste0(rep(c("NIG1", "NIG2", "NIG3", "NIG4"), 
                              each=sum(p)), ".", rep(paste0(
                                rep(1:D, times=p), ".",
                                unlist(sapply(p, function(s) {
                                  return(1:s)}))), 4))
elbo <- matrix(NA, ncol=4, nrow=nreps)
colnames(elbo) <- c("NIG1", "NIG2", "NIG3", "NIG4")

for(r in 1:nreps) {
  cat("\r", "replication", r)
  set.seed(2019 + r)
  
  idtrain <- sample(1:nrow(resp.prep), ntrain)
  xtrain <- lapply(expr.sel, function(s) {scale(s[idtrain, ])})
  ytrain <- scale(resp.prep[idtrain, ])
  xtest <- lapply(expr.sel, function(s) {scale(s[-idtrain, ])})
  ytest <- scale(resp.prep[-idtrain, ])
  
  Z <- matrix(1, nrow=D)
  colnames(Z) <- c("intercept")
  C <- lapply(inpathway, function(s) {
    s <- matrix(1, nrow=length(s)); colnames(s) <- c("intercept"); return(s)})
  cv1.semnig <- semnig(x=xtrain, y=ytrain, C=C, Z=Z, unpenalized=NULL, 
                       standardize=FALSE, intercept=FALSE, fixed.eb="none", 
                       full.post=TRUE, init=NULL, control=control)
  # create co-data matrices
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
  # cv1.bSEM <- bSEM(xtrain, split(ytrain, 1:ncol(ytrain)), 
  #                  lapply(inpathway, "+", 1), 
  #                   control=list(maxit=200, trace=FALSE, epsilon=1e-3))
  
  # cv1.lasso <- lapply(1:D, function(d) {
  #   cv.glmnet(xtrain[[d]], ytrain[, d], intercept=FALSE, standardize=FALSE)})
  # cv1.ridge <- lapply(1:D, function(d) {
  #   cv.glmnet(xtrain[[d]], ytrain[, d], alpha=0, intercept=FALSE, 
  #             standardize=FALSE)})
  
  # apmse[r, ] <- c(mean(sapply(1:D, function(d) {
  #                   mean((ytest[, d] - xtest[[d]] %*%
  #                           cv1.semnig$vb$mpost$beta[[d]])^2)})),
  #                 mean(sapply(1:D, function(d) {
  #                   mean((ytest[, d] - xtest[[d]] %*%
  #                         cv2.semnig$vb$mpost$beta[[d]])^2)})),
  #                 mean(sapply(1:D, function(d) {
  #                   mean((ytest[, d] - xtest[[d]] %*%
  #                         cv1.bSEM$vb$beta[[d]][, 1])^2)})),
  #                 mean((ytest - sapply(1:D, function(d) {
  #                   predict(cv1.lasso[[d]], xtest[[d]], s="lambda.min")}))^2),
  #                 mean((ytest - sapply(1:D, function(d) {
  #                   predict(cv1.ridge[[d]], xtest[[d]], s="lambda.min")}))^2))
  # 
  # brank[r, 1:sum(p)] <- unlist(lapply(cv1.semnig$vb$mpost$beta, function(s) {
  #   rank(abs(s), ties.method="average")}))
  # brank[r, (sum(p) + 1):(2*sum(p))] <-
  #   unlist(lapply(cv2.semnig$vb$mpost$beta, function(s) {
  #     rank(abs(s), ties.method="average")}))
  # brank[r, (2*sum(p) + 1):(3*sum(p))] <-
  #   unlist(lapply(cv1.bSEM$vb$beta, function(s) {
  #     rank(abs(s[, 1]), ties.method="average")}))
  # brank[r, (3*sum(p) + 1):(4*sum(p))] <- unlist(sapply(cv1.lasso, function(s) {
  #   rank(abs(as.numeric(coef(s, s="lambda.min"))[-1]), ties.method="average")}))
  # brank[r, (4*sum(p) + 1):(5*sum(p))] <- unlist(sapply(cv1.ridge, function(s) {
  #   rank(abs(as.numeric(coef(s, s="lambda.min"))[-1]), ties.method="average")}))
  # 
  # elbo[r, ] <- c(mean(cv1.semnig$seq.elbo[nrow(cv1.semnig$seq.elbo), ]),
  #                mean(cv2.semnig$seq.elbo[nrow(cv2.semnig$seq.elbo), ]),
  #                mean(cv1.bSEM$mll[, ncol(cv1.bSEM$mll)]))
  apmse[r, ] <- c(mean(sapply(1:D, function(d) {
    mean((ytest[, d] - xtest[[d]] %*%
            cv1.semnig$vb$mpost$beta[[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((ytest[, d] - xtest[[d]] %*%
              cv2.semnig$vb$mpost$beta[[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((ytest[, d] - xtest[[d]] %*%
              cv3.semnig$vb$mpost$beta[[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((ytest[, d] - xtest[[d]] %*%
              cv4.semnig$vb$mpost$beta[[d]])^2)})))
  
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
  
  elbo[r, ] <- c(mean(cv1.semnig$seq.elbo[nrow(cv1.semnig$seq.elbo), ]),
                 mean(cv2.semnig$seq.elbo[nrow(cv2.semnig$seq.elbo), ]),
                 mean(cv3.semnig$seq.elbo[nrow(cv3.semnig$seq.elbo), ]),
                 mean(cv4.semnig$seq.elbo[nrow(cv4.semnig$seq.elbo), ]))
  
}

# save(apmse, brank, elbo, file="results/data_gdsc_cv1.Rdata")
save(apmse, brank, elbo, file="results/data_gdsc_cv2.Rdata")
# load("results/data_gdsc_cv1.Rdata")

# EB estimates
tab <- rbind(c(fit1.semnig$eb$alphad, rep(NA, 4), fit1.semnig$eb$alphaf, NA),
             c(fit2.semnig$eb$alphad, fit2.semnig$eb$alphaf),
             c(rep(NA, 5), fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 1]/
                 fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 1],
               fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 2]/
                 fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 2] -
                 fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 1]/
                 fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 1]))
colnames(tab) <- c(colnames(Z), colnames(C[[1]]))
rownames(tab) <- c("NIG1", "NIG2", "bSEM")

# Spearman correlations between ranks
brankcor <- cbind(NIG1=cor(t(brank[, substr(colnames(brank), 1, 4)=="NIG1"]), 
                           method="spearman")[
                             upper.tri(matrix(0, nrow(brank), nrow(brank)))],
                  NIG2=cor(t(brank[, substr(colnames(brank), 1, 4)=="NIG2"]), 
                           method="spearman")[
                             upper.tri(matrix(0, nrow(brank), nrow(brank)))],
                  bSEM=cor(t(brank[, substr(colnames(brank), 1, 4)=="bSEM"]), 
                           method="spearman")[
                             upper.tri(matrix(0, nrow(brank), nrow(brank)))])

save(apmse, elbo, tab, brankcor, file="results/data_gdsc_res1.Rdata")




x=xtrain
y=split(ytrain, 1:ncol(ytrain))
P=lapply(inpathway, "+", 1)
EBid=NULL
control=list(maxit=200, trace=FALSE, epsilon=1e-3)
bSEM <- function(x, y, P, EBid=NULL,
                 control=list(maxit=100, trace=TRUE, epsilon=1e-3)) {
  
  if(is.null(EBid)) {
    EBid <- sort(unique(unlist(P)))
  }
  
  niter <- control$maxit
  listX <- x
  listy <- y
  listP <- P
  
  # Intialization
  priors <- sort(unique(unlist(listP)))
  nbpriors <- length(priors)
  aprior <- bprior <- matrix(NA, 1, nbpriors)
  colnames(aprior) <- paste("a", 1:nbpriors, sep="")
  colnames(bprior) <- paste("b", 1:nbpriors, sep="")
  aprior[1,] <- bprior[1,] <- 0.001
  idxPriorList <- lapply(listP, function(x){sort(unique(x))})
  allMLs <- matrix(NA, length(listX), 1)
  allMLs[, 1] <- -Inf
  
  conv <- FALSE
  iter <- 0
  while(!conv & (iter < control$maxit)) {
    iter <- j <- iter + 1
    if(control$trace) {cat("iteration ", j, "\n")}
    
    # Prior as lists
    inputaList <- lapply(idxPriorList, function(x){aprior[j,x]})
    inputbList <- lapply(idxPriorList, function(x){bprior[j,x]})
    
    if(j==1){
      mydstarvec <- rep(0, length(inputaList))
      mybstarlist <- lapply(idxPriorList, function(x){rep(0, length(x))})
    } else {
      mydstarvec <- unlist(lapply(res$postSigList, function(x){x[2]}))
      mybstarlist <- lapply(res$postRandList, function(x){x[,2]})
    }
    
    # Fit BSEM
    res <- BSEMVarOneIter(ylist=listy, Xlist=listX, Plist=listP, alist=inputaList, blist=inputbList, bstarlist=mybstarlist, cSigma=0.001, dSigma=0.001, dstarvec=mydstarvec, lincomblist=list())
    
    # Marginal likelihood
    allMLs <- cbind(allMLs, res$allmargs)
    
    # Empirical Bayes
    aprior <- rbind(aprior, NA)
    bprior <- rbind(bprior, NA)
    for(ii in 1:nbpriors){
      if(ii%in%EBid){
        # Get posterior shape and rate parameters
        allaRandStar <- sapply(1:length(idxPriorList), function(x){res$postRandList[[x]][idxPriorList[[x]]==ii,1]}, simplify=TRUE)
        if(is.list(allaRandStar)){
          allaRandStar <- unlist(allaRandStar)
        }
        allbRandStar <- sapply(1:length(idxPriorList), function(x){res$postRandList[[x]][idxPriorList[[x]]==ii,2]}, simplify=TRUE)
        if(is.list(allbRandStar)){
          allbRandStar <- unlist(allbRandStar)
        }
        
        # Variational Empirical Bayes using fixed-point iteration as in Valpola and Honkela (2006)
        ab <- c(aprior[j, ii], bprior[j, ii])
        ab <- fixedPointIterEB(initab=ab, myallaRandStar=allaRandStar, myallbRandStar=allbRandStar, mymaxiter=20, myeps=1e-6)
        if(control$trace) {ab}
        aprior[j+1, ii] <- ab[1]
        bprior[j+1, ii] <- ab[2]
      }
    }
    conv <- all(abs(allMLs[, j] - allMLs[, j + 1]) <= control$epsilon)
  }
  
  out <- list(eb=list(a=aprior, b=bprior),
              vb=list(gamma.sq=lapply(res$postRandList, function(s) {
                c(a=s[1, 1], b=s[1, 2])}),
                sigma.sq=lapply(res$postSigList, function(s) {
                  c(c=s[1, 1], d=s[2, 1])}),
                beta=lapply(res$postBetaList, function(s) {
                  cbind(mu=s[, 1], dSigma=s[, 2]^2)})),
              mll=allMLs, conv=conv, iter=iter)
  return(out)
}