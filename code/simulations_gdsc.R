#!/usr/bin/env Rscript

# number of cores to use
ncores <- 100

### libraries
packages <- c("foreach", "doParallel", "cambridge", "msm", "statmod", "glmnet",
              "xtune")
sapply(packages, library, character.only=TRUE)

# ### load data
# load(file="data/data_gdsc_dat1.Rdata")

# number of reps
nreps <- 100

################################# simulation 1 #################################
### simulation functions
.cmtnorm <- function(mu, sigma) {
  if(length(mu)==1 & length(sigma)!=1) {
    mu <- rep(mu, length(sigma))
  }
  if(length(sigma)==1 & length(mu)!=1) {
    sigma <- rep(sigma, length(mu))
  }
  ratio <- dnorm(mu/sigma)/pnorm(mu/sigma)
  m1 <- mu + sigma*ratio
  m2 <- ifelse(sigma==0, 0,
               sigma^2*(1 - (mu/sigma)*ratio - ratio^2))
  m3 <- ifelse(sigma==0, 0,
               sigma^3*ratio*((ratio + mu/sigma)^2 + 
                                ratio*(ratio + mu/sigma) - 1))
  return(list(m1=m1, m2=m2, m3=m3))
}
.simlambda <- function(m1, M2, m3, alpha, ESNR) {
  val1 <- sum(alpha*m1)
  val2 <- sum(t(M2*alpha)*alpha)
  val3 <- sum(alpha^2*m3)
  ESNR/(val2*val3 + 3*val1*val2^2 + val1^3*val2)
}
.ratesigma <- function(Ex, Vx, Egammasq, Etausq, ESNR) {
  rate <- ESNR/(ESNR - Etausq*sum(Egammasq*(Vx + Ex^2)))
  return(rate)
}
.sigmax <- function(shape, Etausq, Egammasq, ESNR) {
  sqrt((1 - 1/shape)*ESNR/(Etausq*sum(Egammasq)))
}

### simulation settings
D <- 100
p <- 100
n <- 500
ntest <- 1000
G <- 4
H <- 4
alphaf <- rep(1, G)
alphad <- rep(1, H)
shape <- 3
rate <- 2
muc <- rep(1, G)
sigmac <- c(0, rep(1, G - 1))
# sigmac <- rep(1, G)
muz <- rep(1, H)
sigmaz <- c(0, rep(1, H - 1))
ESNRf <- 100
ESNRd <- 100
ESNR <- 4
momentsc <- .cmtnorm(muc, sigmac)
momentsz <- .cmtnorm(muz, sigmaz)
lambdaf <- .simlambda(momentsc$m1, diag(momentsc$m2), momentsc$m3, alphaf, 
                      ESNRf)
lambdad <- .simlambda(momentsz$m1, diag(momentsz$m2), momentsz$m3, alphad, 
                      ESNRd)
Egammasq <- rep(mean(1/colSums(matrix(rtnorm(G*100000, lower=0), nrow=G, 
                                      ncol=100000)*alphaf)), p)
Etausq <- mean(1/colSums(matrix(rtnorm(H*100000, lower=0), nrow=H, 
                                ncol=100000)*alphad))
sigmax <- .sigmax(shape, Etausq, Egammasq, ESNR)
mux <- 0

# estimation settings
nfolds <- 10
methods <- c("NIG-", "NIG", "ridge", "lasso", "xtune", "ebridge")
control.semnig <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-3, 
                       epsilon.vb=1e-3, maxit.eb=200, maxit.vb=2, 
                       maxit.post=100, maxit.block=0)
control.ebridge <-list(epsilon=sqrt(.Machine$double.eps), maxit=500, 
                       trace=FALSE, glmnet.fit2=FALSE, beta2=FALSE)

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

  ### simulate parameters
  C <- replicate(D, matrix(rtnorm(p*G, rep(muc, each=p), rep(sigmac, each=p), 
                                  0), ncol=G), simplify=FALSE)
  Z <- matrix(rtnorm(D*H, rep(muc, each=D), rep(sigmac, each=D), 0), ncol=H)
  gamma <- sapply(C, function(s) {
    sqrt(rinvgauss(p, 1/as.numeric(s %*% alphaf), lambdaf))})
  tau <- sqrt(rinvgauss(D, 1/as.numeric(Z %*% alphad), lambdad))
  sigma <- 1/sqrt(rgamma(D, shape, rate))
  
  beta <- sapply(1:D, function(d) {rnorm(p, 0, tau[d]*gamma[, d]*sigma[d])})

  # simulate data
  xtrain <- scale(matrix(rnorm(n*p, 0, sigmax), nrow=n, ncol=p,
                         dimnames=list(c(1:n), NULL)))
  xtest <- scale(matrix(rnorm(ntest*p, 0, sigmax), nrow=ntest, ncol=p,
                        dimnames=dimnames(list(c(1:ntest), NULL))))
  ytrain <- lapply(1:D, function(d) {
    s <- as.numeric(scale(rnorm(n, xtrain %*% beta[, d], sigma[d])))
    names(s) <- c(1:length(s))
    return(s)})
  ytest <- lapply(1:D, function(d) {
    s <- as.numeric(scale(rnorm(ntest, xtest %*% beta[, d], sigma[d])))
    names(s) <- c(1:length(s))
    return(s)})
  
  # creating folds
  foldid <- sample(c(rep(1:nfolds, each=n %/% nfolds),
                     rep(1:(n %% nfolds), (n %% nfolds)!=0)))

  # fitting models
  fit.semnig1 <- semnig(x=rep(list(xtrain), D), y=ytrain, 
                        C=replicate(D, list(matrix(rep(1, p)))), 
                        Z=matrix(rep(1, D)), full.post=TRUE, 
                        control=control.semnig)
  fit.semnig2 <- semnig(x=rep(list(xtrain), D), y=ytrain, C=C, Z=Z, 
                        full.post=TRUE, control=control.semnig)

  # standard penalized methods
  fit.lasso1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain, ytrain[[d]], intercept=FALSE, standardize=TRUE,
              foldid=foldid)})
  fit.ridge1 <- lapply(1:D, function(d) {
    cv.glmnet(xtrain, ytrain[[d]], alpha=0, intercept=FALSE, standardize=TRUE,
              foldid=foldid)})
  
  # empirical bayes penalized methods
  fit.xtune1 <- lapply(1:D, function(d) {
    xtune(xtrain, ytrain[[d]], C[[d]][, -1], family="linear", method="ridge", 
          control=list(intercept=FALSE))})
  fit.ebridge1 <- ebridge(xtrain, ytrain, lapply(C, function(s) {s[, -1]}), 
                          hyper=list(lambda=1/sqrt(n*sapply(
                            fit.ridge1, "[[", "lambda.min")), zeta=0, nu=0),
                          Z[, -1], foldid=rep(list(foldid), D), 
                          control=control.ebridge)
  
  # estimates
  best <- list(semnig1=Reduce("cbind", fit.semnig1$vb$mpost$beta), 
               semnig2=Reduce("cbind", fit.semnig2$vb$mpost$beta),
               ridge1=sapply(fit.ridge1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}),
               lasso1=sapply(fit.lasso1, function(s) {
                 coef(s, s="lambda.min")[-1, 1]}),
               xtune1=sapply(fit.xtune1, function(s) {
                 unname(s$beta.est[-1, ])}),
               ebridge1=Reduce("cbind", fit.ebridge1$beta1))
  
  # predictions
  pred <- lapply(best, function(b) {xtest %*% b})
  predt <- lapply(best, function(b) {xtrain %*% b})
  
  # performance measures
  emse <- sapply(best, function(b) {colMeans((b - beta)^2)})
  
  idh <- apply(-abs(beta), 2, order)[1:floor(0.1*p), ]
  emseh <- sapply(best, function(b) {
    sapply(1:D, function(d) {mean((b[idh[, d], d] - beta[idh[, d], d])^2)})})
    
  idl <- apply(abs(beta), 2, order)[1:floor(0.1*p), ]
  emsel <- sapply(best, function(b) {
    sapply(1:D, function(d) {mean((b[idl[, d], d] - beta[idl[, d], d])^2)})})
  
  pmse <- sapply(pred, function(s) {colMeans((s - Reduce("cbind", ytest))^2)})
  
  pmset <- sapply(predt, function(s) {
    colMeans((s - Reduce("cbind", ytrain))^2)})
  
  # save alpha estimates
  est <- cbind(c(fit.semnig1$eb$alphaf, rep(NA, 3), fit.semnig1$eb$alphaf, 
                 rep(NA, 3), fit.semnig1$eb$lambdaf, fit.semnig1$eb$lambdad),
               c(fit.semnig2$eb$alphaf, fit.semnig2$eb$alphaf, 
                 fit.semnig2$eb$lambdaf, fit.semnig2$eb$lambdad),
               rep(NA, 10), rep(NA, 10),
               c(colMeans(-t(sapply(fit.xtune1, "[[", "alpha.est"))), 
                 rep(NA, 6)),
               c(NA, fit.ebridge1$alphaf, NA, fit.ebridge1$alphad, rep(NA, 2)))

  # calculate ELBO and conditional predictive ordinate for semnig models
  elbot <- cbind(fit.semnig1$seq.elbo[nrow(fit.semnig1$seq.elbo), ],
                 fit.semnig2$seq.elbo[nrow(fit.semnig2$seq.elbo), ],
                 rep(NA, D), rep(NA, D), rep(NA, D), rep(NA, D))

  elbo <- cbind(new.elbo(fit.semnig1, replicate(D, list(xtest)), 
                         Reduce("cbind", ytest)),
                new.elbo(fit.semnig1, replicate(D, list(xtest)), 
                         Reduce("cbind", ytest)), rep(NA, D), rep(NA, D), 
                rep(NA, D), rep(NA, D))

  lpml <- cbind(colMeans(logcpo(replicate(D, list(xtest)), 
                                Reduce("cbind", ytest), n, fit.semnig1)),
                colMeans(logcpo(replicate(D, list(xtest)), 
                                Reduce("cbind", ytest), n, fit.semnig2)), 
                rep(NA, D), rep(NA, D), rep(NA, D), rep(NA, D))
  
  list(emse=emse, emsel=emsel, emseh=emseh, pmse=pmse, pmset=pmset, est=est, 
       elbo=elbo, elbot=elbot, lpml=lpml)
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
lpml <- t(sapply(res, "[[", "lpml"))

res <- rbind(emse, emsel, emseh, pmse, pmset, cbind(est, NA, NA, NA), 
             cbind(elbo, NA, NA, NA), cbind(elbot, NA, NA, NA),
             cbind(lpml, NA, NA, NA))
colnames(res) <- c(methods, "null")
rownames(res) <- c(rep(c("emse", "emsel", "emseh", "pmse", "pmset"), 
                       each=nreps),
                   rep(c(paste0("alphaf", 0:3), "lambdaf"), each=nreps),
                   rep(c("elbo", "elbot", "lpml"), each=nreps))
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

methods <- c("NIGd-", "NIGd", "lasso", "ridge")

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
                  maxit.eb=100, maxit.vb=1, maxit.post=200, maxit.block=0)
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
  lpml <- c(sum(logcpo(xtest, ytest, ntrain, fit1.semnig), na.rm=TRUE),
            sum(logcpo(xtest, ytest, ntrain, fit2.semnig), na.rm=TRUE))

  list(emse=emse, emsel=emsel, emseh=emseh, pmse=pmse, pmset=pmset, est=est, 
       elbo=elbo, elbot=elbot, lpml=lpml)
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
lpml <- t(sapply(res, "[[", "lpml"))

res <- rbind(emse, emsel, emseh, pmse, pmset, cbind(est, NA, NA, NA), 
             cbind(elbo, NA, NA, NA), cbind(elbot, NA, NA, NA),
             cbind(lpml, NA, NA, NA))
colnames(res) <- c(methods, "null")
rownames(res) <- c(rep(c("emse", "emsel", "emseh", "pmse", "pmset"), 
                       each=nreps),
                   rep(c(paste0("alphad", 0:3), "lambdad"), each=nreps),
                   rep(c("elbo", "elbot", "lpml"), each=nreps))
write.table(res, file="results/simulations_gdsc_res2.txt")




################################# simulation 3 #################################
### data preparation
# select features
D <- 100
psel <- 400
o <- order(-apply(expr$expr, 2, sd))
idsel <- rep(list(o[c(1:psel)]), D)
expr.sel <- lapply(idsel, function(s) {expr$expr[, s]})
x <- lapply(expr.sel, function(s) {scale(s)})

### data preparation
p <- sapply(x, ncol)
n <- nrow(expr.sel[[1]])
ntrain <- floor(n/2)

### simulation settings
alphaf <- c(1, 1, 3, 7)*4
alphad <- c(1, 1, 3, 7)*8
lambdaf <- 1
lambdad <- 1
shape <- 3
rate <- 2

SNR <- sum(rep(p[1]/4, 4)/
             as.numeric(model.matrix(~as.factor(c(1:4))) %*% alphaf))/
  as.numeric(model.matrix(~as.factor(c(1:4))) %*% alphad)

methods <- c("NIGfd-", "NIGfd", "lasso", "ridge")

# setup cluster
cl <- makeCluster(ncores) 
if(ncores > 1) {
  registerDoParallel(cl)
} else {
  registerDoSEQ()
}
res <- foreach(r=1:nreps, .packages=packages, .errorhandling="pass") %dopar% {
  cat("\r", "replication", r)
  set.seed(2020 + r)

  ### simulate parameters
  C <- lapply(1:D, function(d) {
    unname(model.matrix(~ as.factor(rep(c(1:4), each=p[d]/4))))})
  Z <- unname(model.matrix(~ as.factor(rep(c(1:4), each=D/4))))
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
    rnorm(ntrain, as.numeric(xtrain[[d]] %*% beta[[d]]), sigma[d])}))
  ytest <- scale(sapply(1:D, function(d) {
    rnorm(n - ntrain, as.numeric(xtest[[d]] %*% beta[[d]]), sigma[d])}))

  ### fitting models
  control <- list(conv.post=TRUE, trace=FALSE,
                  epsilon.eb=1e-3, epsilon.vb=1e-3,
                  maxit.eb=100, maxit.vb=1, maxit.post=200, maxit.block=0)
  fit2.semnig <- semnig(x=xtrain, y=ytrain, C=C, Z=Z, unpenalized=NULL,
                        standardize=FALSE, intercept=FALSE, fixed.eb="none",
                        full.post=TRUE, init=NULL, control=control)
  
  # no external covariates
  control <- list(conv.post=TRUE, trace=FALSE,
                  epsilon.eb=1e-3, epsilon.vb=1e-3,
                  maxit.eb=100, maxit.vb=1, maxit.post=200, maxit.block=0)
  C <- lapply(1:D, function(d) {matrix(1, nrow=p[d], ncol=1)})
  Z <- matrix(1, ncol=1, nrow=D)
  fit1.semnig <- semnig(x=xtrain, y=ytrain, C=C, Z=Z, unpenalized=NULL,
                        standardize=FALSE, intercept=FALSE, fixed.eb="none",
                        full.post=TRUE, init=NULL, control=control)
  
  control <- list(conv.post=TRUE, trace=FALSE,
                  epsilon.eb=1e-3, epsilon.vb=1e-3,
                  maxit.eb=100, maxit.vb=1, maxit.post=200, maxit.block=10)
  
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
  
  lpml <- c(sum(logcpo(xtest, ytest, ntrain, fit1.semnig), na.rm=TRUE),
            sum(logcpo(xtest, ytest, ntrain, fit2.semnig), na.rm=TRUE))
  
  list(emse=emse, emsel=emsel, emseh=emseh, pmse=pmse, pmset=pmset, est=est, 
       elbo=elbo, elbot=elbot, lpml=lpml)
}  
stopCluster(cl=cl)

errorid <- sapply(res, function(s) {is.numeric(s[[1]])})
res2 <- res[errorid]

# prepare and save results table
emse <- t(sapply(res2, "[[", "emse"))
emsel <- t(sapply(res2, "[[", "emsel"))
emseh <- t(sapply(res2, "[[", "emseh"))
pmse <- t(sapply(res2, "[[", "pmse"))
pmset <- t(sapply(res2, "[[", "pmset"))
est <- Reduce("rbind", lapply(1:nrow(res2[[1]]$est), function(i) {
  t(sapply(res2, function(s) {s$est[i, ]}))}))
elbo <- t(sapply(res2, "[[", "elbo"))
elbot <- t(sapply(res2, "[[", "elbot"))
lpml <- t(sapply(res2, "[[", "lpml"))

res <- rbind(emse, emsel, emseh, pmse, pmset, cbind(est, NA, NA, NA), 
             cbind(elbo, NA, NA, NA), cbind(elbot, NA, NA, NA),
             cbind(lpml, NA, NA, NA))
colnames(res) <- c(methods, "null")
rownames(res) <- c(rep(c("emse", "emsel", "emseh", "pmse", "pmset"), 
                       each=nreps),
                   rep(c(paste0("alphaf", 0:3), "lambdaf",
                         paste0("alphad", 0:3), "lambdad"), each=nreps),
                   rep(c("elbo", "elbot", "lpml"), each=nreps))
write.table(res, file="results/simulations_gdsc_res3.txt")


################################# simulation 4 #################################
### data preparation
# select features
D <- 100
psel <- 100
o <- order(-apply(expr$expr, 2, sd))
idsel <- rep(list(o[c(1:psel)]), D)
expr.sel <- lapply(idsel, function(s) {expr$expr[, s]})
x <- lapply(expr.sel, function(s) {scale(s)})

### data preparation
p <- sapply(x, ncol)
n <- nrow(expr.sel[[1]])
ntrain <- floor(n/2)

### simulation settings
alphaf <- c(1, 1, 3, 7)*10
alphad <- c(1, 1, 3, 7)

methods <- c("NIGfd-", "NIGfd", "lasso", "ridge")

# setup cluster
cl <- makeCluster(ncores) 
if(ncores > 1) {
  registerDoParallel(cl)
} else {
  registerDoSEQ()
}
res <- foreach(r=1:nreps, .packages=packages, .errorhandling="pass") %dopar% {
  cat("\r", "replication", r)
  set.seed(2020 + r)
  
  ### simulate parameters
  C <- lapply(1:D, function(d) {
    unname(model.matrix(~ as.factor(rep(c(1:4), each=p[d]/4))))})
  Z <- unname(model.matrix(~ as.factor(rep(c(1:4), each=D/4))))
  sigma <- rep(1, D)
  beta <- lapply(1:D, function(d) {
    rnorm(p[d], 0, sigma[d]*sqrt(1/as.numeric(C[[d]] %*% alphaf))*
            sqrt(1/as.numeric(Z %*% alphad)[d]))})
  
  # simulate data
  idtrain <- sample(1:n, ntrain)
  xtrain <- lapply(x, function(s) {scale(s[idtrain, ])})
  xtest <- lapply(x, function(s) {scale(s[-idtrain, ])})
  ytrain <- scale(sapply(1:D, function(d) {
    rnorm(ntrain, as.numeric(xtrain[[d]] %*% beta[[d]]), sigma[d])}))
  ytest <- scale(sapply(1:D, function(d) {
    rnorm(n - ntrain, as.numeric(xtest[[d]] %*% beta[[d]]), sigma[d])}))
  
  ### fitting models
  control <- list(conv.post=TRUE, trace=FALSE,
                  epsilon.eb=1e-3, epsilon.vb=1e-3,
                  maxit.eb=100, maxit.vb=1, maxit.post=200, maxit.block=0)
  fit2.semnig <- semnig(x=xtrain, y=ytrain, C=C, Z=Z, unpenalized=NULL,
                        standardize=FALSE, intercept=FALSE, fixed.eb="none",
                        full.post=TRUE, init=NULL, control=control)
  
  # no external covariates
  control <- list(conv.post=TRUE, trace=FALSE,
                  epsilon.eb=1e-3, epsilon.vb=1e-3,
                  maxit.eb=100, maxit.vb=1, maxit.post=200, maxit.block=0)
  C <- lapply(1:D, function(d) {matrix(1, nrow=p[d], ncol=1)})
  Z <- matrix(1, ncol=1, nrow=D)
  fit1.semnig <- semnig(x=xtrain, y=ytrain, C=C, Z=Z, unpenalized=NULL,
                        standardize=FALSE, intercept=FALSE, fixed.eb="none",
                        full.post=TRUE, init=NULL, control=control)
  
  control <- list(conv.post=TRUE, trace=FALSE,
                  epsilon.eb=1e-3, epsilon.vb=1e-3,
                  maxit.eb=100, maxit.vb=1, maxit.post=200, maxit.block=10)
  
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
  
  lpml <- c(sum(logcpo(xtest, ytest, ntrain, fit1.semnig), na.rm=TRUE),
            sum(logcpo(xtest, ytest, ntrain, fit2.semnig), na.rm=TRUE))
  
  list(emse=emse, emsel=emsel, emseh=emseh, pmse=pmse, pmset=pmset, est=est, 
       elbo=elbo, elbot=elbot, lpml=lpml)
}  
stopCluster(cl=cl)

errorid <- sapply(res, function(s) {is.numeric(s[[1]])})
res2 <- res[errorid]

# prepare and save results table
emse <- t(sapply(res2, "[[", "emse"))
emsel <- t(sapply(res2, "[[", "emsel"))
emseh <- t(sapply(res2, "[[", "emseh"))
pmse <- t(sapply(res2, "[[", "pmse"))
pmset <- t(sapply(res2, "[[", "pmset"))
est <- Reduce("rbind", lapply(1:nrow(res2[[1]]$est), function(i) {
  t(sapply(res2, function(s) {s$est[i, ]}))}))
elbo <- t(sapply(res2, "[[", "elbo"))
elbot <- t(sapply(res2, "[[", "elbot"))
lpml <- t(sapply(res2, "[[", "lpml"))

res <- rbind(emse, emsel, emseh, pmse, pmset, cbind(est, NA, NA, NA), 
             cbind(elbo, NA, NA, NA), cbind(elbot, NA, NA, NA),
             cbind(lpml, NA, NA, NA))
colnames(res) <- c(methods, "null")
rownames(res) <- c(rep(c("emse", "emsel", "emseh", "pmse", "pmset"), 
                       each=nreps),
                   rep(c(paste0("alphaf", 0:3), "lambdaf",
                         paste0("alphad", 0:3), "lambdad"), each=nreps),
                   rep(c("elbo", "elbot", "lpml"), each=nreps))
write.table(res, file="results/simulations_gdsc_res4.txt")
