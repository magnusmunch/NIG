#!/usr/bin/env Rscript

# number of cores to use
ncores <- 50

### libraries
packages <- c("foreach", "doParallel", "cambridge", "statmod", "glmnet")
sapply(packages, library, character.only=TRUE)

### load data
load(file="data/data_gdsc_dat1.Rdata")

################################# simulation 1 #################################
### data preparation
# select features
D <- ncol(resp.prep)
psel <- 100
o <- order(-apply(expr.prep, 2, sd))
idsel <- rep(list(o[c(1:psel)]), D)
expr.sel <- lapply(idsel, function(s) {expr.prep[, s]})
x <- lapply(expr.sel, function(s) {scale(s)})

### data preparation
p <- sapply(x, ncol)
n <- nrow(expr.prep)
ntrain <- floor(n/2)

### simulation settings
alphaf <- c(1, 1, 3, 7)
lambdaf <- 1
shape <- 3
rate <- 2

nreps <- 50
methods <- c("NIG$_{\text{f}}^-$", "NIG$_{\text{f}}$", "lasso", "ridge")

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

  ### simulate parameters
  C <- lapply(1:D, function(d) {
    cbind(1, t(replicate(p[d], sample(c(1, 0, 0, 0), 4)[-1])))})
  gamma <- lapply(1:D, function(d) {
    sqrt(rinvgauss(p[d], 1/as.numeric(C[[d]] %*% alphaf), lambdaf))})
  sigma <- 1/sqrt(rgamma(D, shape, rate))
  beta <- lapply(1:D, function(d) {rnorm(p[d], 0, gamma[[d]]*sigma[d])})

  # simulate data
  idtrain <- sample(1:n, ntrain)
  xtrain <- lapply(x, function(s) {scale(s[idtrain, ])})
  xtest <- lapply(x, function(s) {scale(s[-idtrain, ])})
  ytrain <- scale(sapply(1:D, function(d) {
    rnorm(ntrain, as.numeric(xtrain[[d]] %*% beta[[d]]), sigma[d])}),
    scale=FALSE)
  ytest <- scale(sapply(1:D, function(d) {
    rnorm(n - ntrain, as.numeric(xtest[[d]] %*% beta[[d]]), sigma[d])}),
    scale=FALSE)

  ### fitting models
  control <- list(conv.post=TRUE, trace=FALSE,
                  epsilon.eb=1e-3, epsilon.vb=1e-3,
                  maxit.eb=100, maxit.vb=1, maxit.post=200)
  fit2.semnig <- semnig(x=xtrain, y=ytrain, C=C, Z=NULL, unpenalized=NULL,
                        standardize=FALSE, intercept=FALSE, fixed.eb="none",
                        full.post=TRUE, init=NULL, control=control)
  C <- lapply(1:D, function(d) {matrix(1, nrow=p[d], ncol=1)})
  fit1.semnig <- semnig(x=xtrain, y=ytrain, C=C, Z=NULL, unpenalized=NULL,
                        standardize=FALSE, intercept=FALSE, fixed.eb="none",
                        full.post=TRUE, init=NULL, control=control)

  fit1.lasso <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[, d], intercept=FALSE, standardize=FALSE)})
  fit1.ridge <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[, d], alpha=0, intercept=FALSE,
              standardize=FALSE)})

  best1.semnig <- fit1.semnig$vb$mpost$beta
  best2.semnig <- fit2.semnig$vb$mpost$beta
  best1.lasso <- lapply(fit1.lasso, function(s) {
    coef(s, s="lambda.min")[-1, 1]})
  best1.ridge <- lapply(fit1.ridge, function(s) {
    coef(s, s="lambda.min")[-1, 1]})

  pred1.semnig <- sapply(1:D, function(d) {xtest[[d]] %*% best1.semnig[[d]]})
  pred2.semnig <- sapply(1:D, function(d) {xtest[[d]] %*% best2.semnig[[d]]})
  pred1.lasso <- sapply(1:D, function(d) {xtest[[d]] %*% best1.lasso[[d]]})
  pred1.ridge <- sapply(1:D, function(d) {xtest[[d]] %*% best1.ridge[[d]]})

  predt1.semnig <- sapply(1:D, function(d) {xtrain[[d]] %*% best1.semnig[[d]]})
  predt2.semnig <- sapply(1:D, function(d) {xtrain[[d]] %*% best2.semnig[[d]]})
  predt1.lasso <- sapply(1:D, function(d) {xtrain[[d]] %*% best1.lasso[[d]]})
  predt1.ridge <- sapply(1:D, function(d) {xtrain[[d]] %*% best1.ridge[[d]]})

  emse <- c(mean(sapply(1:D, function(d) {
    mean((best1.semnig[[d]] - beta[[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((best2.semnig[[d]] - beta[[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((best1.lasso[[d]] - beta[[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((best1.ridge[[d]] - beta[[d]])^2)})),
    mean(unlist(beta)^2))
  
  cutoffs <- sapply(1:D, function(d) {quantile(abs(beta[[d]]), probs=0.9)})
  emseh <- c(mean(sapply(1:D, function(d) {
    mean((best1.semnig[[d]][abs(beta[[d]]) > cutoffs[d]] - 
            beta[[d]][abs(beta[[d]]) > cutoffs[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((best2.semnig[[d]][abs(beta[[d]]) > cutoffs[d]] - 
              beta[[d]][abs(beta[[d]]) > cutoffs[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((best1.lasso[[d]][abs(beta[[d]]) > cutoffs[d]] - 
              beta[[d]][abs(beta[[d]]) > cutoffs[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((best1.ridge[[d]][abs(beta[[d]]) > cutoffs[d]] - 
              beta[[d]][abs(beta[[d]]) > cutoffs[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean(beta[[d]][abs(beta[[d]]) > cutoffs[[d]]]^2)})))
  
  emsel <- c(mean(sapply(1:D, function(d) {
    mean((best1.semnig[[d]][abs(beta[[d]]) <= cutoffs[d]] - 
            beta[[d]][abs(beta[[d]]) <= cutoffs[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((best2.semnig[[d]][abs(beta[[d]]) <= cutoffs[d]] - 
              beta[[d]][abs(beta[[d]]) <= cutoffs[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((best1.lasso[[d]][abs(beta[[d]]) <= cutoffs[d]] - 
              beta[[d]][abs(beta[[d]]) <= cutoffs[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((best1.ridge[[d]][abs(beta[[d]]) <= cutoffs[d]] - 
              beta[[d]][abs(beta[[d]]) <= cutoffs[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean(beta[[d]][abs(beta[[d]]) <= cutoffs[[d]]]^2)})))
  
  pmse <- c(mean(colMeans((pred1.semnig - ytest)^2)),
            mean(colMeans((pred2.semnig - ytest)^2)),
            mean(colMeans((pred1.lasso - ytest)^2)),
            mean(colMeans((pred1.ridge - ytest)^2)),
            mean(apply(ytest, 1, "-", colMeans(ytrain))^2))

  pmset <- c(mean(colMeans((predt1.semnig - ytrain)^2)),
             mean(colMeans((predt2.semnig - ytrain)^2)),
             mean(colMeans((predt1.lasso - ytrain)^2)),
             mean(colMeans((predt1.ridge - ytrain)^2)),
             mean(apply(ytrain, 1, "-", colMeans(ytrain))^2))

  est <- cbind(c(fit1.semnig$eb$alphaf, rep(NA, 3), fit1.semnig$eb$lambdaf),
               c(fit2.semnig$eb$alphaf, fit2.semnig$eb$lambdaf))

  # calculate ELBO for semnig models
  elbot <- c(mean(fit1.semnig$seq.elbo[nrow(fit1.semnig$seq.elbo), ]),
             mean(fit2.semnig$seq.elbo[nrow(fit2.semnig$seq.elbo), ]))

  elbo <- c(mean(new.elbo(fit1.semnig, xtest, ytest)),
            mean(new.elbo(fit2.semnig, xtest, ytest)))

  list(emse=emse, emsel=emsel, emseh=emseh, pmse=pmse, pmset=pmset, est=est, 
       elbo=elbo, elbot=elbot)
}
stopCluster(cl=cl)

# prepare and save results table
emse <- t(sapply(res, "[[", "emse"))
emsel <- t(sapply(res, "[[", "emsel"))
emseh <- t(sapply(res, "[[", "emseh"))
pmse <- t(sapply(res, "[[", "pmse"))
pmset <- t(sapply(res, "[[", "pmset"))
est <- Reduce("rbind", lapply(1:nrow(res[[1]]$est), function(i) {
  t(sapply(res, function(s) {s$est[i, ]}))}))
elbo <- t(sapply(res, "[[", "elbo"))
elbot <- t(sapply(res, "[[", "elbot"))

res <- rbind(emse, emsel, emseh, pmse, pmset, cbind(est, NA, NA, NA), 
             cbind(elbo, NA, NA, NA), cbind(elbot, NA, NA, NA))
colnames(res) <- c(methods, "null")
rownames(res) <- c(rep(c("emse", "emsel", "emseh", "pmse", "pmset"), 
                       each=nreps),
                   rep(c(paste0("alphaf", 0:3), "lambdaf"), each=nreps),
                   rep(c("elbo", "elbot"), each=nreps))
write.table(res, file="results/simulations_gdsc_res1.txt")




################################# simulation 2 #################################
### data preparation
# select features
D <- ncol(resp.prep)
psel <- 100
o <- order(-apply(expr.prep, 2, sd))
idsel <- rep(list(o[c(1:psel)]), D)
expr.sel <- lapply(idsel, function(s) {expr.prep[, s]})
x <- lapply(expr.sel, function(s) {scale(s)})

### data preparation
p <- sapply(x, ncol)
n <- nrow(expr.prep)
ntrain <- floor(n/2)

### simulation settings
# alphad <- c(0.1, 1/5 - 0.1, 1/2 - 0.1, 1 - 0.1)
alphad <- c(1, 1, 3, 7)
lambdad <- 1
shape <- 3
rate <- 2

nreps <- 50
methods <- c("NIG$_{\text{d}}^-$", "NIG$_{\text{d}}$", "lasso", "ridge")

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

  ### simulate parameters
  Z <- cbind(1, t(replicate(D, sample(c(1, 0, 0, 0), 4)[-1])))
  tau <- sqrt(rinvgauss(D, 1/as.numeric(Z %*% alphad), lambdad))
  sigma <- 1/sqrt(rgamma(D, shape, rate))
  beta <- lapply(1:D, function(d) {rnorm(p[d], 0, tau[d]*sigma[d])})

  # boxplot(tau*sigma ~ as.factor(Z[, 1] + Z[, 2] + 2*Z[, 3] + 3*Z[, 4]))

  # simulate data
  idtrain <- sample(1:n, ntrain)
  xtrain <- lapply(x, function(s) {scale(s[idtrain, ])})
  xtest <- lapply(x, function(s) {scale(s[-idtrain, ])})
  ytrain <- scale(sapply(1:D, function(d) {
    rnorm(ntrain, as.numeric(xtrain[[d]] %*% beta[[d]]), sigma[d])}),
    scale=FALSE)
  ytest <- scale(sapply(1:D, function(d) {
    rnorm(n - ntrain, as.numeric(xtest[[d]] %*% beta[[d]]), sigma[d])}),
    scale=FALSE)

  ### fitting models
  control <- list(conv.post=TRUE, trace=FALSE,
                  epsilon.eb=1e-3, epsilon.vb=1e-3,
                  maxit.eb=100, maxit.vb=1, maxit.post=200)
  fit2.semnig <- semnig(x=xtrain, y=ytrain, C=NULL, Z=Z, unpenalized=NULL,
                        standardize=FALSE, intercept=FALSE, fixed.eb="none",
                        full.post=TRUE, init=NULL, control=control)
  Z <- matrix(1, ncol=1, nrow=D)
  fit1.semnig <- semnig(x=xtrain, y=ytrain, C=NULL, Z=Z, unpenalized=NULL,
                        standardize=FALSE, intercept=FALSE, fixed.eb="none",
                        full.post=TRUE, init=NULL, control=control)

  fit1.lasso <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[, d], intercept=FALSE, standardize=FALSE)})
  fit1.ridge <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[, d], alpha=0, intercept=FALSE,
              standardize=FALSE)})

  best1.semnig <- fit1.semnig$vb$mpost$beta
  best2.semnig <- fit2.semnig$vb$mpost$beta
  best1.lasso <- lapply(fit1.lasso, function(s) {
    coef(s, s="lambda.min")[-1, 1]})
  best1.ridge <- lapply(fit1.ridge, function(s) {
    coef(s, s="lambda.min")[-1, 1]})

  pred1.semnig <- sapply(1:D, function(d) {xtest[[d]] %*% best1.semnig[[d]]})
  pred2.semnig <- sapply(1:D, function(d) {xtest[[d]] %*% best2.semnig[[d]]})
  pred1.lasso <- sapply(1:D, function(d) {xtest[[d]] %*% best1.lasso[[d]]})
  pred1.ridge <- sapply(1:D, function(d) {xtest[[d]] %*% best1.ridge[[d]]})

  predt1.semnig <- sapply(1:D, function(d) {xtrain[[d]] %*% best1.semnig[[d]]})
  predt2.semnig <- sapply(1:D, function(d) {xtrain[[d]] %*% best2.semnig[[d]]})
  predt1.lasso <- sapply(1:D, function(d) {xtrain[[d]] %*% best1.lasso[[d]]})
  predt1.ridge <- sapply(1:D, function(d) {xtrain[[d]] %*% best1.ridge[[d]]})

  emse <- c(mean(sapply(1:D, function(d) {
    mean((best1.semnig[[d]] - beta[[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((best2.semnig[[d]] - beta[[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((best1.lasso[[d]] - beta[[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((best1.ridge[[d]] - beta[[d]])^2)})),
    mean(unlist(beta)^2))
  
  cutoffs <- sapply(1:D, function(d) {quantile(abs(beta[[d]]), probs=0.9)})
  emseh <- c(mean(sapply(1:D, function(d) {
    mean((best1.semnig[[d]][abs(beta[[d]]) > cutoffs[d]] - 
            beta[[d]][abs(beta[[d]]) > cutoffs[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((best2.semnig[[d]][abs(beta[[d]]) > cutoffs[d]] - 
              beta[[d]][abs(beta[[d]]) > cutoffs[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((best1.lasso[[d]][abs(beta[[d]]) > cutoffs[d]] - 
              beta[[d]][abs(beta[[d]]) > cutoffs[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((best1.ridge[[d]][abs(beta[[d]]) > cutoffs[d]] - 
              beta[[d]][abs(beta[[d]]) > cutoffs[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean(beta[[d]][abs(beta[[d]]) > cutoffs[[d]]]^2)})))
  
  emsel <- c(mean(sapply(1:D, function(d) {
    mean((best1.semnig[[d]][abs(beta[[d]]) <= cutoffs[d]] - 
            beta[[d]][abs(beta[[d]]) <= cutoffs[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((best2.semnig[[d]][abs(beta[[d]]) <= cutoffs[d]] - 
              beta[[d]][abs(beta[[d]]) <= cutoffs[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((best1.lasso[[d]][abs(beta[[d]]) <= cutoffs[d]] - 
              beta[[d]][abs(beta[[d]]) <= cutoffs[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((best1.ridge[[d]][abs(beta[[d]]) <= cutoffs[d]] - 
              beta[[d]][abs(beta[[d]]) <= cutoffs[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean(beta[[d]][abs(beta[[d]]) <= cutoffs[[d]]]^2)})))

  pmse <- c(mean(colMeans((pred1.semnig - ytest)^2)),
            mean(colMeans((pred2.semnig - ytest)^2)),
            mean(colMeans((pred1.lasso - ytest)^2)),
            mean(colMeans((pred1.ridge - ytest)^2)),
            mean(apply(ytest, 1, "-", colMeans(ytrain))^2))

  pmset <- c(mean(colMeans((predt1.semnig - ytrain)^2)),
             mean(colMeans((predt2.semnig - ytrain)^2)),
             mean(colMeans((predt1.lasso - ytrain)^2)),
             mean(colMeans((predt1.ridge - ytrain)^2)),
             mean(apply(ytrain, 1, "-", colMeans(ytrain))^2))

  est <- cbind(c(fit1.semnig$eb$alphad, rep(NA, 3), fit1.semnig$eb$lambdad),
               c(fit2.semnig$eb$alphad, fit2.semnig$eb$lambdad))

  # calculate ELBO for semnig models
  elbot <- c(mean(fit1.semnig$seq.elbo[nrow(fit1.semnig$seq.elbo), ]),
             mean(fit2.semnig$seq.elbo[nrow(fit2.semnig$seq.elbo), ]))

  elbo <- c(mean(new.elbo(fit1.semnig, xtest, ytest)),
            mean(new.elbo(fit2.semnig, xtest, ytest)))

  list(emse=emse, emsel=emsel, emseh=emseh, pmse=pmse, pmset=pmset, est=est, 
       elbo=elbo, elbot=elbot)
}
stopCluster(cl=cl)

# prepare and save results table
emse <- t(sapply(res, "[[", "emse"))
emsel <- t(sapply(res, "[[", "emsel"))
emseh <- t(sapply(res, "[[", "emseh"))
pmse <- t(sapply(res, "[[", "pmse"))
pmset <- t(sapply(res, "[[", "pmset"))
est <- Reduce("rbind", lapply(1:nrow(res[[1]]$est), function(i) {
  t(sapply(res, function(s) {s$est[i, ]}))}))
elbo <- t(sapply(res, "[[", "elbo"))
elbot <- t(sapply(res, "[[", "elbot"))

res <- rbind(emse, emsel, emseh, pmse, pmset, cbind(est, NA, NA, NA), 
             cbind(elbo, NA, NA, NA), cbind(elbot, NA, NA, NA))
colnames(res) <- c(methods, "null")
rownames(res) <- c(rep(c("emse", "emsel", "emseh", "pmse", "pmset"), 
                       each=nreps),
                   rep(c(paste0("alphad", 0:3), "lambdad"), each=nreps),
                   rep(c("elbo", "elbot"), each=nreps))
write.table(res, file="results/simulations_gdsc_res2.txt")




################################# simulation 3 #################################
### data preparation
# select features
D <- ncol(resp.prep)
psel <- 100
o <- order(-apply(expr.prep, 2, sd))
idsel <- rep(list(o[c(1:psel)]), D)
expr.sel <- lapply(idsel, function(s) {expr.prep[, s]})
x <- lapply(expr.sel, function(s) {scale(s)})

### data preparation
p <- sapply(x, ncol)
n <- nrow(expr.prep)
ntrain <- floor(n/2)

### simulation settings
alphaf <- c(1, 1, 3, 7)
alphad <- c(1, 1, 3, 7)
lambdaf <- 1
lambdad <- 1
shape <- 3
rate <- 2

nreps <- 50
methods <- c("NIG$_{\text{f+d}}^-$", "NIG$_{\text{f+d}}$", "lasso", "ridge")

# setup cluster
cl <- makeCluster(ncores) 
if(ncores > 1) {
  registerDoParallel(cl)
} else {
  registerDoSEQ()
}
res <- foreach(r=1:nreps, .packages=packages) %dopar% {
# for(r in 1:nreps) {
  cat("\r", "replication", r)
  set.seed(2019 + r)

  ### simulate parameters
  C <- lapply(1:D, function(d) {
    cbind(1, t(replicate(p[d], sample(c(1, 0, 0, 0), 4)[-1])))})
  Z <- cbind(1, t(replicate(D, sample(c(1, 0, 0, 0), 4)[-1])))
  gamma <- lapply(1:D, function(d) {
    sqrt(rinvgauss(p, 1/as.numeric(C[[d]] %*% alphaf), lambdaf))})
  tau <- sqrt(rinvgauss(D, 1/as.numeric(Z %*% alphad), lambdad))
  sigma <- 1/sqrt(rgamma(D, shape, rate))
  beta <- lapply(1:D, function(d) {rnorm(p[d], 0, gamma[[d]]*tau[d]*sigma[d])})

  # simulate data
  idtrain <- sample(1:n, ntrain)
  xtrain <- lapply(x, function(s) {scale(s[idtrain, ])})
  xtest <- lapply(x, function(s) {scale(s[-idtrain, ])})
  ytrain <- scale(sapply(1:D, function(d) {
    rnorm(ntrain, as.numeric(xtrain[[d]] %*% beta[[d]]), sigma[d])}),
    scale=FALSE)
  ytest <- scale(sapply(1:D, function(d) {
    rnorm(n - ntrain, as.numeric(xtest[[d]] %*% beta[[d]]), sigma[d])}),
    scale=FALSE)

  ### fitting models
  control <- list(conv.post=TRUE, trace=FALSE,
                  epsilon.eb=1e-3, epsilon.vb=1e-3,
                  maxit.eb=100, maxit.vb=1, maxit.post=200)
  fit2.semnig <- semnig(x=xtrain, y=ytrain, C=C, Z=Z, unpenalized=NULL,
                        standardize=FALSE, intercept=FALSE, fixed.eb="none",
                        full.post=TRUE, init=NULL, control=control)

  C <- lapply(1:D, function(d) {matrix(1, nrow=p[d], ncol=1)})
  Z <- matrix(1, ncol=1, nrow=D)
  fit1.semnig <- semnig(x=xtrain, y=ytrain, C=C, Z=Z, unpenalized=NULL,
                        standardize=FALSE, intercept=FALSE, fixed.eb="none",
                        full.post=TRUE, init=NULL, control=control)

  fit1.lasso <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[, d], intercept=FALSE, standardize=FALSE)})
  fit1.ridge <- lapply(1:D, function(d) {
    cv.glmnet(xtrain[[d]], ytrain[, d], alpha=0, intercept=FALSE,
              standardize=FALSE)})

  best1.semnig <- fit1.semnig$vb$mpost$beta
  best2.semnig <- fit2.semnig$vb$mpost$beta
  best1.lasso <- lapply(fit1.lasso, function(s) {
    coef(s, s="lambda.min")[-1, 1]})
  best1.ridge <- lapply(fit1.ridge, function(s) {
    coef(s, s="lambda.min")[-1, 1]})

  pred1.semnig <- sapply(1:D, function(d) {xtest[[d]] %*% best1.semnig[[d]]})
  pred2.semnig <- sapply(1:D, function(d) {xtest[[d]] %*% best2.semnig[[d]]})
  pred1.lasso <- sapply(1:D, function(d) {xtest[[d]] %*% best1.lasso[[d]]})
  pred1.ridge <- sapply(1:D, function(d) {xtest[[d]] %*% best1.ridge[[d]]})

  predt1.semnig <- sapply(1:D, function(d) {xtrain[[d]] %*% best1.semnig[[d]]})
  predt2.semnig <- sapply(1:D, function(d) {xtrain[[d]] %*% best2.semnig[[d]]})
  predt1.lasso <- sapply(1:D, function(d) {xtrain[[d]] %*% best1.lasso[[d]]})
  predt1.ridge <- sapply(1:D, function(d) {xtrain[[d]] %*% best1.ridge[[d]]})

  emse <- c(mean(sapply(1:D, function(d) {
    mean((best1.semnig[[d]] - beta[[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((best2.semnig[[d]] - beta[[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((best1.lasso[[d]] - beta[[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((best1.ridge[[d]] - beta[[d]])^2)})),
    mean(unlist(beta)^2))
  
  cutoffs <- sapply(1:D, function(d) {quantile(abs(beta[[d]]), probs=0.9)})
  emseh <- c(mean(sapply(1:D, function(d) {
    mean((best1.semnig[[d]][abs(beta[[d]]) > cutoffs[d]] - 
            beta[[d]][abs(beta[[d]]) > cutoffs[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((best2.semnig[[d]][abs(beta[[d]]) > cutoffs[d]] - 
              beta[[d]][abs(beta[[d]]) > cutoffs[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((best1.lasso[[d]][abs(beta[[d]]) > cutoffs[d]] - 
              beta[[d]][abs(beta[[d]]) > cutoffs[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((best1.ridge[[d]][abs(beta[[d]]) > cutoffs[d]] - 
              beta[[d]][abs(beta[[d]]) > cutoffs[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean(beta[[d]][abs(beta[[d]]) > cutoffs[[d]]]^2)})))
  
  emsel <- c(mean(sapply(1:D, function(d) {
    mean((best1.semnig[[d]][abs(beta[[d]]) <= cutoffs[d]] - 
            beta[[d]][abs(beta[[d]]) <= cutoffs[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((best2.semnig[[d]][abs(beta[[d]]) <= cutoffs[d]] - 
              beta[[d]][abs(beta[[d]]) <= cutoffs[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((best1.lasso[[d]][abs(beta[[d]]) <= cutoffs[d]] - 
              beta[[d]][abs(beta[[d]]) <= cutoffs[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((best1.ridge[[d]][abs(beta[[d]]) <= cutoffs[d]] - 
              beta[[d]][abs(beta[[d]]) <= cutoffs[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean(beta[[d]][abs(beta[[d]]) <= cutoffs[[d]]]^2)})))
  
  pmse <- c(mean(colMeans((pred1.semnig - ytest)^2)),
            mean(colMeans((pred2.semnig - ytest)^2)),
            mean(colMeans((pred1.lasso - ytest)^2)),
            mean(colMeans((pred1.ridge - ytest)^2)),
            mean(apply(ytest, 1, "-", colMeans(ytrain))^2))
  
  pmset <- c(mean(colMeans((predt1.semnig - ytrain)^2)),
             mean(colMeans((predt2.semnig - ytrain)^2)),
             mean(colMeans((predt1.lasso - ytrain)^2)),
             mean(colMeans((predt1.ridge - ytrain)^2)),
             mean(apply(ytrain, 1, "-", colMeans(ytrain))^2))
  
  est <- cbind(c(fit1.semnig$eb$alphaf, rep(NA, 3), fit1.semnig$eb$lambdaf,
                 fit1.semnig$eb$alphad, rep(NA, 3), fit1.semnig$eb$lambdad),
               c(fit2.semnig$eb$alphaf, fit2.semnig$eb$lambdaf,
                 fit2.semnig$eb$alphad, fit2.semnig$eb$lambdad))
  
  # calculate ELBO for semnig models
  elbot <- c(mean(fit1.semnig$seq.elbo[nrow(fit1.semnig$seq.elbo), ]),
             mean(fit2.semnig$seq.elbo[nrow(fit2.semnig$seq.elbo), ]))
  
  elbo <- c(mean(new.elbo(fit1.semnig, xtest, ytest)),
            mean(new.elbo(fit2.semnig, xtest, ytest)))
  
  list(emse=emse, emsel=emsel, emseh=emseh, pmse=pmse, pmset=pmset, est=est, 
       elbo=elbo, elbot=elbot)
}  
stopCluster(cl=cl)

# prepare and save results table
emse <- t(sapply(res, "[[", "emse"))
emsel <- t(sapply(res, "[[", "emsel"))
emseh <- t(sapply(res, "[[", "emseh"))
pmse <- t(sapply(res, "[[", "pmse"))
pmset <- t(sapply(res, "[[", "pmset"))
est <- Reduce("rbind", lapply(1:nrow(res[[1]]$est), function(i) {
  t(sapply(res, function(s) {s$est[i, ]}))}))
elbo <- t(sapply(res, "[[", "elbo"))
elbot <- t(sapply(res, "[[", "elbot"))

res <- rbind(emse, emsel, emseh, pmse, pmset, cbind(est, NA, NA, NA), 
             cbind(elbo, NA, NA, NA), cbind(elbot, NA, NA, NA))
colnames(res) <- c(methods, "null")
rownames(res) <- c(rep(c("emse", "emsel", "emseh", "pmse", "pmset"), 
                       each=nreps),
                   rep(c(paste0("alphaf", 0:3), "lambdaf",
                         paste0("alphad", 0:3), "lambdad"), each=nreps),
                   rep(c("elbo", "elbot"), each=nreps))
write.table(res, file="results/simulations_gdsc_res3.txt")
