#!/usr/bin/env Rscript

### installation of packages
if(!("devtools" %in% installed.packages())) {
  install.packages("devtools")
}
library(devtools)
install_github("magnusmunch/cambridge/rpackage", local=FALSE,
               auth_token=Sys.getenv("GITHUB_PAT"))

### parallelisation
parallel <- TRUE

### libraries
library(cambridge)
library(statmod)
library(GeneralizedHyperbolic)
library(rstan)
library(glmnet)
library(mvtnorm)
library(cramer)
library(foreach)
library(doParallel)

# ################################# simulation 1 #################################
# # settings
# nreps <- 100
# n <- 100
# p <- 100
# D <- 100
# nclass <- 4
# 
# # create fixed parameters
# alpha <- c(1:nclass)
# C.inv.gauss <- model.matrix(~ -1 + factor(rep(c(1:nclass), each=D/nclass)))
# C.inv.gamma <- as.numeric(C.inv.gauss %*% c(1:nclass))
# theta <- as.numeric(1/(C.inv.gauss %*% alpha))
# lambda <- rep(1, D)
# sigma <- rep(1, D)
# SNR <- 50
# 
# # store settings in an object
# temp1 <- c(nreps=nreps, n=n, p=p, D=D, nclass=nclass, alpha=alpha,
#            C=C.inv.gauss, theta=theta, lambda=lambda, sigma=sigma, SNR=SNR)
# set1 <- t(matrix(temp1))
# colnames(set1) <- names(temp1)
# 
# # objects to store simulation results
# methods <- c("igauss.conj", "igauss.non.conj", "igamma.conj", "igamma.non.conj")
# params <- c("alpha", "eta", "theta", "lambda", "mu", "dSigma", "delta", "zeta")
# drugs <- paste("drug", c(1:D), sep="")
# covs <- paste("cov", c(1:nclass), sep="")
# res1 <- as.data.frame(
#   matrix(NA, nrow=nreps, ncol=length(methods)*length(covs) +
#            length(methods)*length(params[-1])*length(drugs),
#          dimnames=list(paste("rep", c(1:nreps), sep=""),
#                        c(apply(expand.grid(covs, methods, params[1]), 1, paste,
#                                collapse="_"),
#                          apply(expand.grid(drugs, methods, params[-1]), 1,
#                                paste, collapse="_")))))
# 
# # set initial values and control parameters
# control <- list(epsilon.eb=-1, epsilon.vb=-1,
#                 epsilon.opt=sqrt(sqrt(.Machine$double.eps)),
#                 maxit.eb=20, maxit.vb=2, maxit.opt=100, trace=FALSE)
# init.inv.gauss <- list(a=rep(mean(1/sigma^2), D), b=rep(mean(theta), D),
#                        theta=theta, lambda=rep(1, D))
# init.inv.gamma <- list(a=rep(mean(1/sigma^2), D), b=rep(mean(theta), D),
#                        eta=rep(theta, D), lambda=rep(1, D))
# 
# # simulation
# set.seed(2018)
# for(r in 1:nreps) {
# 
#   cat("\r", "replication", r)
#   # create data
#   gamma <- sqrt(rinvgauss(D, theta, shape=lambda))
#   beta <- sapply(1:D, function(d) {rnorm(p, 0, sigma[d]*gamma[d])})
#   x <- matrix(rnorm(n*p, 0, sqrt(SNR/(mean(theta)*p))), nrow=n, ncol=p)
#   y <- sapply(1:D, function(d) {rnorm(n, x %*% beta[, d], sigma[d])})
# 
#   # fit inverse Gaussian models
#   fit1.igauss.conj <- est.model(x, y, C.inv.gauss, "inv. Gaussian", TRUE,
#                                 init=init.inv.gauss, control=control)
#   fit1.igauss.non.conj <- est.model(x, y, C.inv.gauss, "inv. Gaussian", FALSE,
#                                     init=init.inv.gauss, control=control)
# 
#   # fit independent inverse Gamma models
#   fit1.igamma.conj <- est.model(x, y, C.inv.gamma, "inv. Gamma", TRUE,
#                                 init=init.inv.gamma, control=control)
#   fit1.igamma.non.conj <- est.model(x, y, C.inv.gamma, "inv. Gamma", FALSE,
#                                     init=init.inv.gamma, control=control)
# 
#   # store results
#   temp1 <- paste("fit1.", apply(expand.grid(methods, "seq.eb", params[1]), 1,
#                                 paste, collapse="$"), "[fit1.", methods,
#                  "$iter$eb", ", ]", sep="")
#   temp2 <- paste("fit1.", apply(expand.grid(methods, "seq.eb", params[2:4]), 1,
#                                 paste, collapse="$"), "[fit1.", methods,
#                  "$iter$eb", ", ]", sep="")
#   temp3 <- paste("fit1.", apply(expand.grid(methods, "vb.post", params[5]), 1,
#                                 paste, collapse="$"), sep="")
#   temp4 <- paste("fit1.", apply(expand.grid(methods, "vb.post", params[6]), 1,
#                                 paste, collapse="$"), sep="")
#   temp5 <- paste("fit1.", apply(expand.grid(methods, "vb.post", params[7:8]), 1,
#                                 paste, collapse="$"), sep="")
#   res1[r, ] <- c(c(sapply(temp1, function(s) {eval(parse(text=s))})),
#                  c(sapply(temp2, function(s) {eval(parse(text=s))})),
#                  c(sapply(temp3, function(s) {
#                    colMeans((eval(parse(text=s)) - beta)^2)})),
#                  c(sapply(temp4, function(s) {colMeans(eval(parse(text=s)))})),
#                  c(sapply(temp5, function(s) {eval(parse(text=s))})))
# 
#   # save fitted models and results
#   fit1 <- data.frame(igauss.conj=fit1.igauss.conj$seq.eb,
#                      igauss.non.conj=fit1.igauss.non.conj$seq.eb,
#                      igamma.conj=fit1.igamma.conj$seq.eb,
#                      igamma.non.conj=fit1.igamma.non.conj$seq.eb)
#   write.table(set1, file="results/simulations_igaussian_set1.csv")
#   write.table(res1, file="results/simulations_igaussian_res1.csv")
#   write.table(fit1, file="results/simulations_igaussian_fit1.csv")
# 
# }
#
# ###############################   simulation 2   ###############################
# # settings
# nreps <- 100
# n <- 100
# p <- 100
# D <- 100
# nclass <- 4
# 
# # create fixed parameters
# alpha <- c(1:nclass) + 2
# C.inv.gauss <- model.matrix(~ -1 + factor(rep(c(1:nclass), each=D/nclass)))
# C.inv.gamma <- as.numeric(C.inv.gauss %*% c(1:nclass))
# eta <- as.numeric(C.inv.gauss %*% alpha)
# lambda <- rep(1, D)
# sigma <- rep(1, D)
# SNR <- 50
# 
# # store settings in an object
# temp1 <- c(nreps=nreps, n=n, p=p, D=D, nclass=nclass, alpha=alpha,
#            C=C.inv.gauss, eta=eta, lambda=lambda, sigma=sigma, SNR=SNR)
# set2 <- t(matrix(temp1))
# colnames(set2) <- names(temp1)
# 
# # objects to store simulation results
# methods <- c("igauss.conj", "igauss.non.conj", "igamma.conj", "igamma.non.conj")
# params <- c("alpha", "eta", "theta", "lambda", "mu", "dSigma", "delta", "zeta")
# drugs <- paste("drug", c(1:D), sep="")
# covs <- paste("cov", c(1:nclass), sep="")
# res2 <- as.data.frame(
#   matrix(NA, nrow=nreps, ncol=length(methods)*length(covs) +
#            length(methods)*length(params[-1])*length(drugs),
#          dimnames=list(paste("rep", c(1:nreps), sep=""),
#                        c(apply(expand.grid(covs, methods, params[1]), 1, paste,
#                                collapse="_"),
#                          apply(expand.grid(drugs, methods, params[-1]), 1,
#                                paste, collapse="_")))))
# 
# # set initial values and control parameters
# control <- list(epsilon.eb=-1, epsilon.vb=-1,
#                 epsilon.opt=sqrt(sqrt(.Machine$double.eps)),
#                 maxit.eb=30, maxit.vb=4, maxit.opt=100, trace=FALSE)
# init.inv.gauss <- list(a=rep(mean(1/sigma^2), D),
#                        b=rep(mean(lambda/(eta - 2)), D),
#                        theta=rep(mean(lambda/(eta - 2)), D), lambda=rep(1, D))
# init.inv.gamma <- list(a=rep(mean(1/sigma^2), D),
#                        b=rep(mean(lambda/(eta - 2)), D),
#                        eta=eta, lambda=rep(1, D))
# 
# # simulation
# set.seed(2018)
# for(r in 1:nreps) {
# 
#   cat("\r", "replication", r)
# 
#   # create data
#   gamma <- 1/sqrt(rgamma(D, shape=eta/2, scale=lambda/2))
#   beta <- sapply(1:D, function(d) {rnorm(p, 0, sigma[d]*gamma[d])})
#   x <- matrix(rnorm(n*p, 0, sqrt(SNR/(mean(lambda/(eta - 2))*p))), nrow=n,
#               ncol=p)
#   y <- sapply(1:D, function(d) {rnorm(n, x %*% beta[, d], sigma[d])})
# 
#   # fit inverse Gaussian models
#   fit2.igauss.conj <- est.model(x, y, C.inv.gauss, "inv. Gaussian", TRUE,
#                                 init=init.inv.gauss, control=control)
#   fit2.igauss.non.conj <- est.model(x, y, C.inv.gauss, "inv. Gaussian", FALSE,
#                                     init=init.inv.gauss, control=control)
# 
#   # fit independent inverse Gamma models
#   fit2.igamma.conj <- est.model(x, y, C.inv.gamma, "inv. Gamma", TRUE,
#                                 init=init.inv.gamma, control=control)
#   fit2.igamma.non.conj <- est.model(x, y, C.inv.gamma, "inv. Gamma", FALSE,
#                                     init=init.inv.gamma, control=control)
# 
#   # store results
#   temp1 <- paste("fit2.", apply(expand.grid(methods, "seq.eb", params[1]), 1,
#                                 paste, collapse="$"), "[fit2.", methods,
#                  "$iter$eb", ", ]", sep="")
#   temp2 <- paste("fit2.", apply(expand.grid(methods, "seq.eb", params[2:4]), 1,
#                                 paste, collapse="$"), "[fit2.", methods,
#                  "$iter$eb", ", ]", sep="")
#   temp3 <- paste("fit2.", apply(expand.grid(methods, "vb.post", params[5]), 1,
#                                 paste, collapse="$"), sep="")
#   temp4 <- paste("fit2.", apply(expand.grid(methods, "vb.post", params[6]), 1,
#                                 paste, collapse="$"), sep="")
#   temp5 <- paste("fit2.", apply(expand.grid(methods, "vb.post", params[7:8]), 1,
#                                 paste, collapse="$"), sep="")
#   res2[r, ] <- c(c(sapply(temp1, function(s) {eval(parse(text=s))})),
#                  c(sapply(temp2, function(s) {eval(parse(text=s))})),
#                  c(sapply(temp3, function(s) {
#                    colMeans((eval(parse(text=s)) - beta)^2)})),
#                  c(sapply(temp4, function(s) {colMeans(eval(parse(text=s)))})),
#                  c(sapply(temp5, function(s) {eval(parse(text=s))})))
#   # save fitted models and results
#   fit2 <- data.frame(igauss.conj=fit2.igauss.conj$seq.eb,
#                      igauss.non.conj=fit2.igauss.non.conj$seq.eb,
#                      igamma.conj=fit2.igamma.conj$seq.eb,
#                      igamma.non.conj=fit2.igamma.non.conj$seq.eb)
#   write.table(set2, file="results/simulations_igaussian_set2.csv")
#   write.table(res2, file="results/simulations_igaussian_res2.csv")
#   write.table(fit2, file="results/simulations_igaussian_fit2.csv")
# 
# }
# 
# warnings()
# 
# 
# ###############################   simulation 3   ###############################
# # settings
# nreps <- 100
# n <- 100
# p <- 200
# D <- 1
# 
# # create fixed parameters
# rho <- 0.5
# Sigma <- matrix(rho, ncol=p, nrow=p); diag(Sigma) <- 1
# gamma <- rep(1, D)
# sigma <- rep(1, D)
# lambda <- 0.1
# theta <- 2
# 
# # store settings in an object
# set <- data.frame(nreps=nreps, n=n, p=p, D=D, gamma=t(as.matrix(gamma)),
#                   sigma=t(as.matrix(sigma)),
#                   lambda=lambda, theta=theta)
# rownames(set) <- paste0("set", 1:nrow(set))
# 
# # objects to store simulation results
# methods <- c("MCMC", "ADVI", "VB", "MAP", "glmnet")
# params <- c("beta0", "beta", "sigma", "gamma")
# 
# # set initial values and control parameters
# stan.model <- stan_model("code/igauss.stan")
# nsamples <- 1000
# nwarmup <- 1000
# 
# ### analysis splits in parallel
# ncores <- min(detectCores() - 1, nreps)
# cluster <- makeForkCluster(ncores)
# if(parallel) {
#   registerDoParallel(cluster)
# } else {
#   registerDoSEQ()
# }
# 
# res <- foreach(k=c(1:nreps), .errorhandling="pass") %dopar% {
#   print(paste("rep", k))
#   set.seed(2019 + k)
#   beta <- sapply(1:D, function(d) {rnorm(p, 0, sigma[d]*gamma[d])})
#   x <- rmvnorm(n, rep(0, p), sigma=Sigma)
#   y <- sapply(1:D, function(d) {rnorm(n, x %*% beta[, d], sigma[d])})
# 
#   fit.mcmc <- sapply(1:D, function(d) {
#     sampling(stan.model, data=list(n=n, p=p, y=y[, d], x=x, theta=theta,
#                                    lambda=lambda), chains=1, warmup=nwarmup,
#              iter=nsamples + nwarmup, cores=1, refresh=0,
#              control=list(adapt_delta=0.8, max_treedepth=12))})
#   fit.advi <- sapply(1:D, function(d) {
#     vb(stan.model, data=list(n=n, p=p, y=y[, d], x=x, theta=theta,
#                              lambda=lambda), output_samples=nsamples)})
#   fit.vb <- est.vb(x, y, theta, lambda, intercept=TRUE, init=NULL,
#                    posterior=TRUE,
#                    control=list(epsilon.vb=1e-6, maxit.vb=100, trace=TRUE))
#   fit.glmnet <- sapply(1:D, function(d) {glmnet(x, y[, d], alpha=0,
#                                                 lambda=1/(n*theta*sigma[d]^2))},
#                        simplify=FALSE)
#   fit.map <- sapply(1:D, function(d) {
#     optimizing(stan.model, data=list(n=n, p=p, y=y[, d], x=x, theta=theta,
#                                      lambda=lambda))}, simplify=FALSE)
# 
#   best <- cbind(mcmc=sapply(fit.mcmc, function(s) {
#     c(get_posterior_mean(s, pars="beta0"),
#       get_posterior_mean(s, pars="beta"))}),
#     advi=sapply(fit.advi, function(s) {
#       c(get_posterior_mean(s, pars="beta0"),
#         get_posterior_mean(s, pars="beta"))}),
#     vb=fit.vb$vb.post$mu,
#     glmnet=sapply(fit.glmnet, function(s) {as.numeric(coef(s))}),
#     map=sapply(fit.map, function(s) {s$par[1:(p + 1)]}))
#   bvar <- cbind(mcmc=sapply(fit.mcmc, function(s) {
#     summary(s)$summary[1:(p + 1), 3]^2}),
#     advi=sapply(fit.advi, function(s) {
#       summary(s)$summary[1:(p + 1), 2]^2}),
#     vb=sapply(1:D, function(d) {diag(fit.vb$vb.post$Sigma[[d]])}))
# 
#   samp.advi <- matrix(unlist(fit.advi[[1]]@sim$samples[[1]][1:(p + 1)]),
#                       ncol=p + 1, byrow=FALSE)
#   samp.vb <- rmvnorm(nsamples, mean=as.numeric(fit.vb$vb.post$mu),
#                      sigma=fit.vb$vb.post$Sigma[[1]])
#   samp.mcmc <- cbind(as.matrix(fit.mcmc[[1]], "beta0"),
#                      as.matrix(fit.mcmc[[1]], "beta"))
#   samp <- cbind(samp.mcmc, samp.advi, samp.vb)
# 
#   cram <- c(advi=cramer.test(samp.advi, samp.mcmc,
#                              just.statistic=TRUE)$statistic,
#             vb=cramer.test(samp.vb, samp.mcmc, just.statistic=TRUE)$statistic)
# 
#   corbvar <- c(advi=cor(bvar[, 2], bvar[, 1]),
#                vb=cor(bvar[, 3], bvar[, 1]))
#   msebvar <- c(advi=mean((bvar[, 2] - bvar[, 1])^2),
#                vb=mean((bvar[, 3] - bvar[, 1])^2))
# 
#   corbest <- c(advi=cor(best[, 2], best[, 1]),
#                vb=cor(best[, 3], best[, 1]),
#                glmnet=cor(best[, 4], best[, 1]),
#                map=cor(best[, 5], best[, 1]))
# 
#   msebest <- c(advi=mean((best[, 2] - best[, 1])^2),
#                vb=mean((best[, 3] - best[, 1])^2),
#                glmnet=mean((best[, 4] - best[, 1])^2),
#                map=mean((best[, 5] - best[, 1])^2))
#   out <- list(cram, corbvar, msebvar, corbest, msebest, samp)
# }
# if(parallel) {stopCluster(cluster)}
# save(res, file="results/simulations_igaussian_res3.Rdata")
# 
# which.excl <- 80
# fit <- res[[1]][[6]]
# dimnames(fit) <- list(paste0("sample", 1:nsamples),
#                       c(paste0(paste0("beta", 0:p), ".mcmc"),
#                         paste0(paste0("beta", 0:p), ".advi"),
#                         paste0(paste0("beta", 0:p), ".vb")))
# res <- data.frame(cram=t(sapply(res[-which.excl], function(s) {s[[1]]})),
#                   corbvar=t(sapply(res[-which.excl], function(s) {s[[2]]})),
#                   msebvar=t(sapply(res[-which.excl], function(s) {s[[3]]})),
#                   corbest=t(sapply(res[-which.excl], function(s) {s[[4]]})),
#                   msebest=t(sapply(res[-which.excl], function(s) {s[[5]]})))
# rownames(res) <- paste0("rep", (1:nreps)[-which.excl])
# write.table(fit, file="results/simulations_igaussian_fit3.csv")
# write.table(res, file="results/simulations_igaussian_res3.csv")
# write.table(set, file="results/simulations_igaussian_set3.csv")
# 
# 
# ###############################   simulation 4   ###############################
# # settings
# nreps <- 100
# n <- 100
# p <- 200
# D <- 1
# 
# # create fixed parameters (eta chosen such that marginal beta variances match)
# rho <- 0.5
# Sigma <- matrix(rho, ncol=p, nrow=p); diag(Sigma) <- 1
# gamma <- rep(1, D)
# sigma <- rep(1, D)
# 
# V <- rep(c(1/2, 1, 2), times=2)
# K <- rep(c(4, 10), each=3)
# thetac <- V/sigma^2
# lambdac <- c(3/2, 3, 6, 3/14, 3/7, 6/7)
# ratio <- besselK(lambdac/thetac, 0, expon.scaled=TRUE)/
#   besselK(lambdac/thetac, 1, expon.scaled=TRUE)
# thetaz <- 1/(ratio + 2*thetac/lambdac)/sigma^2
# lambdaz <- (3*thetaz^3 + 12*thetaz^2*thetac)/
#   (K*lambdac - 3*thetaz^2 - 12*thetaz*thetac)
# eta <- 6/(K - 3) + 4
# lambda <- V*(eta - 2)/sigma^2
# 
# # check variances
# vars <- cbind(nig=sigma^2*thetac,
#               nigig=sigma^2*thetac*thetaz*(ratio + 2*thetac/lambdac),
#               studentst=sigma^2*lambda/(eta - 2))
# all.equal(vars, cbind(nig=V, nigig=V, studentst=V))
# kurts <- cbind(nig=3*thetac/lambdac + 3,
#                nigig=(3*thetaz/lambdaz + 3)*
#                  (4*thetac*(ratio + 2*thetac/lambdac) + 1)/
#                  (lambdac*(ratio + 2*thetac/lambdac)^2),
#                studentst=6/(eta - 4) + 3)
# all.equal(kurts, cbind(nig=K, nigig=K, studentst=K))
# 
# params <- cbind(lambdac=lambdac, lambdaz=lambdaz, thetac=thetac, thetaz=thetaz,
#                 eta=eta, lambda)
# 
# # store settings in an object
# set <- data.frame(nreps=nreps, n=n, p=p, D=D, gamma=gamma, sigma=sigma,
#                   lambdac=lambdac, thetac=thetac, lambdaz=lambdaz,
#                   thetaz=thetaz, eta=eta, lambda=lambda, vars=V, kurts=K)
# rownames(set) <- paste0("set", 1:nrow(set))
# 
# # objects to store simulation results
# methods <- c("nig", "studentst", "nigig")
# params <- c("beta0", "beta", "sigma", "gamma", "tau")
# 
# # set initial values and control parameters
# stan.nig <- stan_model("code/igauss.stan")
# stan.nigig <- stan_model("code/nigig.stan")
# stan.studentst <- stan_model("code/igamma.stan")
# nsamples <- 1000
# nwarmup <- 1000
# 
# ### analysis splits in parallel
# ncores <- min(detectCores() - 1, nreps)
# cluster <- makeForkCluster(ncores)
# if(parallel) {
#   registerDoParallel(cluster)
# } else {
#   registerDoSEQ()
# }
# 
# 
# res <- foreach(s=c(1:nrow(set)), .errorhandling="pass") %:%
#   foreach(k=c(1:nreps), .errorhandling="pass") %dopar% {
#     print(paste0("set ", s, ", rep", k))
# 
#     set.seed(2019 + k)
#     beta <- sapply(1:D, function(d) {rnorm(p, 0, sigma[d]*gamma[d])})
#     x <- rmvnorm(n, rep(0, p), sigma=Sigma)
#     y <- sapply(1:D, function(d) {rnorm(n, x %*% beta[, d], sigma[d])})
# 
#     thetac <- set$thetac[s]
#     lambdac <- set$lambdac[s]
#     thetaz <- set$thetaz[s]
#     lambdaz <- set$lambdaz[s]
#     eta <- set$eta[s]
#     lambda <- set$lambda[s]
#     fit.nig <- sapply(1:D, function(d) {
#       sampling(stan.nig, data=list(n=n, p=p, y=y[, d], x=x, theta=thetac,
#                                  lambda=lambdac), chains=1, warmup=nwarmup,
#                iter=nsamples + nwarmup, cores=1, refresh=0,
#                control=list(adapt_delta=0.8, max_treedepth=12))})
#     fit.nigig <- sapply(1:D, function(d) {
#       sampling(stan.nigig, data=list(n=n, p=p, y=y[, d], x=x, thetac=thetac,
#                                      lambdac=lambdac, thetaz=thetaz, lambdaz),
#                chains=1, warmup=nwarmup, iter=nsamples + nwarmup, cores=1,
#                refresh=0, control=list(adapt_delta=0.8, max_treedepth=12))})
#     fit.studentst <- sapply(1:D, function(d) {
#       sampling(stan.studentst, data=list(n=n, p=p, y=y[, d], x=x, eta=eta,
#                                          lambda=lambda), chains=1,
#                warmup=nwarmup, iter=nsamples + nwarmup, cores=1, refresh=0,
#                control=list(adapt_delta=0.8, max_treedepth=12))})
# 
#     best <- data.frame(nig=sapply(fit.nig, function(s) {
#       c(get_posterior_mean(s, pars="beta0"),
#         get_posterior_mean(s, pars="beta"))}),
#       nigig=sapply(fit.nigig, function(s) {
#         c(get_posterior_mean(s, pars="beta0"),
#           get_posterior_mean(s, pars="beta"))}),
#       studentst=sapply(fit.studentst, function(s) {
#         c(get_posterior_mean(s, pars="beta0"),
#           get_posterior_mean(s, pars="beta"))}))
# 
#     samp.nig <- cbind(as.matrix(fit.nig[[1]], "beta0"),
#                       as.matrix(fit.nig[[1]], "beta"))
#     samp.nigig <- cbind(as.matrix(fit.nigig[[1]], "beta0"),
#                         as.matrix(fit.nigig[[1]], "beta"))
#     samp.studentst <- cbind(as.matrix(fit.studentst[[1]], "beta0"),
#                             as.matrix(fit.studentst[[1]], "beta"))
#     cov.true <- sigma^2*solve(t(x) %*% x + gamma^(-2)*diag(p))
#     mean.true <- cov.true %*% t(x) %*% y/sigma^2
#     samp.true <- cbind(0, rmvnorm(nsamples, mean.true, cov.true))
#     samp <- data.frame(nig=samp.nig, nigig=samp.nigig, studentst=samp.studentst)
# 
#     cram <- c(nig=cramer.test(samp.nig, samp.true,
#                               just.statistic=TRUE)$statistic,
#               nigig=cramer.test(samp.nigig, samp.true,
#                                 just.statistic=TRUE)$statistic,
#               studentst=cramer.test(samp.studentst, samp.true,
#                                     just.statistic=TRUE)$statistic)
# 
#     corbest <- c(nig=cor(c(0, beta), best[, 1]),
#                  nigig=cor(c(0, beta), best[, 2]),
#                  studentst=cor(c(0, beta), best[, 3]))
# 
#     msebest <- c(nig=mean((c(0, beta) - best[, 1])^2),
#                  nigig=mean((c(0, beta) - best[, 2])^2),
#                  studentst=mean((c(0, beta) - best[, 3])^2))
# 
#     ci.nig <- summary(fit.nig[[1]], pars=c("beta0", "beta"))$summary[, c(4, 8)]
#     ci.nigig <- summary(fit.nigig[[1]],
#                         pars=c("beta0", "beta"))$summary[, c(4, 8)]
#     ci.studentst <- summary(fit.studentst[[1]],
#                             pars=c("beta0", "beta"))$summary[, c(4, 8)]
# 
#     cover <- c(nig=mean(sapply(1:(p + 1), function(j) {
#       (ci.nig[j, 1] < c(0, beta)[j]) & (c(0, beta)[j] < ci.nig[j, 2])})),
#       nigig=mean(sapply(1:(p + 1), function(j) {
#         (ci.nigig[j, 1] < c(0, beta)[j]) & (c(0, beta)[j] < ci.nigig[j, 2])})),
#       studentst=mean(sapply(1:(p + 1), function(j) {
#         (ci.studentst[j, 1] < c(0, beta)[j]) &
#           (c(0, beta)[j] < ci.studentst[j, 2])})))
# 
#     list(cram=cram, corbest=corbest, msebest=msebest, cover=cover,
#          samp=samp)
#   }
# if(parallel) {stopCluster(cluster)}
# save(res, file="results/simulations_igaussian_res4.2.Rdata")
# 
# fit <- Reduce("cbind", sapply(res, function(l) {as.matrix(l[[1]][[5]])},
#                               simplify=FALSE))
# dimnames(fit) <- list(paste0("sample", 1:nsamples),
#                       paste(rep(paste0("set", 1:nrow(set)),
#                                 each=(p + 1)*length(methods)),
#                             c(paste0(paste0("beta", 0:p), ".nig"),
#                               paste0(paste0("beta", 0:p), ".nigig"),
#                               paste0(paste0("beta", 0:p), ".studentst")),
#                             sep="."))
# 
# res <- Reduce("cbind", sapply(res, function(r) {
#   data.frame(cram=t(sapply(r, function(s) {s[[1]]})),
#              corbest=t(sapply(r, function(s) {s[[2]]})),
#              msebest=t(sapply(r, function(s) {s[[3]]})),
#              cover=t(sapply(r, function(s) {s[[4]]})))}))
# dimnames(res) <- list(paste0("rep", (1:nreps)), paste(rep(paste0(
#   "set", 1:nrow(set)), each=4*length(methods)), rep(apply(expand.grid(
#     methods, c("cram", "corbest", "msebest", "cover")), 1, function(s) {
#       paste(rev(s), collapse=".")}), times=nrow(set)), sep="."))
# 
# write.table(fit, file="results/simulations_igaussian_fit4.csv")
# write.table(res, file="results/simulations_igaussian_res4.csv")
# write.table(set, file="results/simulations_igaussian_set4.csv")
# 
# 
# ###############################   simulation 5   ###############################
# # settings
# nreps <- 100
# n <- 100
# p <- 30
# D <- 1
# 
# # create fixed parameters (eta chosen such that marginal beta variances match)
# rho <- 0.5
# Sigma <- matrix(rho, ncol=p, nrow=p); diag(Sigma) <- 1
# gamma <- rep(1, D)
# sigma <- rep(1, D)
# nu <- 1
# kappa <- c(nu, 9*nu)
# r <- 0.01
# 
# V <- rep(c(1/2, 1, 2), times=2)
# K <- rep(c(4, 10), each=3)
# thetac <- V/sigma^2
# lambdac <- c(3/2, 3, 6, 3/14, 3/7, 6/7)
# ratio <- besselK(lambdac/thetac, 0, expon.scaled=TRUE)/
#   besselK(lambdac/thetac, 1, expon.scaled=TRUE)
# thetaz <- 1/(ratio + 2*thetac/lambdac)/sigma^2
# lambdaz <- (3*thetaz^3 + 12*thetaz^2*thetac)/
#   (K*lambdac - 3*thetaz^2 - 12*thetaz*thetac)
# eta <- 6/(K - 3) + 4
# lambda <- V*(eta - 2)/sigma^2
# 
# # check variances
# vars <- cbind(nig=sigma^2*thetac, 
#               nigig=sigma^2*thetac*thetaz*(ratio + 2*thetac/lambdac),
#               studentst=sigma^2*lambda/(eta - 2))
# all.equal(vars, cbind(nig=V, nigig=V, studentst=V))
# kurts <- cbind(nig=3*thetac/lambdac + 3,
#                nigig=(3*thetaz/lambdaz + 3)*
#                  (4*thetac*(ratio + 2*thetac/lambdac) + 1)/
#                  (lambdac*(ratio + 2*thetac/lambdac)^2),
#                studentst=6/(eta - 4) + 3)
# all.equal(kurts, cbind(nig=K, nigig=K, studentst=K))
# 
# params <- cbind(lambdac=lambdac, lambdaz=lambdaz, thetac=thetac, thetaz=thetaz, 
#                 eta=eta, lambda)
# 
# # store settings in an object
# set <- data.frame(nreps=nreps, n=n, p=p, D=D, gamma=gamma, sigma=sigma, r=r,
#                   nu=nu, kappa=rep(kappa, each=6), lambdac=lambdac, 
#                   thetac=thetac, lambdaz=lambdaz, 
#                   thetaz=thetaz, eta=eta, lambda=lambda, vars=V, kurts=K)
# rownames(set) <- paste0("set", 1:nrow(set))
# 
# # objects to store simulation results
# methods <- c("nig", "studentst", "nigig")
# params <- c("beta0", "beta", "sigma", "gamma", "tau")
# 
# # set initial values and control parameters
# stan.nig <- stan_model("code/igauss.stan")
# stan.nigig <- stan_model("code/nigig.stan")
# stan.studentst <- stan_model("code/igamma.stan")
# nsamples <- 1000
# nwarmup <- 1000
# 
# ### analysis splits in parallel
# # ncores <- min(detectCores() - 1, nreps)
# ncores <- 50
# cluster <- makeForkCluster(ncores)
# if(parallel) {
#   registerDoParallel(cluster)
# } else {
#   registerDoSEQ()
# }
# 
# res <- foreach(s=c(1:nrow(set)), .errorhandling="pass") %:% 
#   foreach(k=c(1:nreps), .errorhandling="pass") %dopar% {
#     print(paste0("set", s, ", rep", k))
# 
#     set.seed(2019 + k)
#     beta <- sapply(1:D, function(d) {
#       rspikeandslab(p, set[s, "nu"], set[s, "kappa"], 
#                     set[s, "sigma"], set[s, "gamma"], set[s, "r"])})
#     x <- rmvnorm(n, rep(0, p), sigma=Sigma)
#     y <- sapply(1:D, function(d) {rnorm(n, x %*% beta[, d], sigma[d])})
#     
#     thetac <- set$thetac[s]
#     lambdac <- set$lambdac[s]
#     thetaz <- set$thetaz[s]
#     lambdaz <- set$lambdaz[s]
#     eta <- set$eta[s]
#     lambda <- set$lambda[s]
#     fit.nig <- sapply(1:D, function(d) {
#       sampling(stan.nig, data=list(n=n, p=p, y=y[, d], x=x, theta=thetac, 
#                                    lambda=lambdac), chains=1, warmup=nwarmup,
#                iter=nsamples + nwarmup, cores=1, refresh=0, 
#                control=list(adapt_delta=0.8, max_treedepth=12))})
#     fit.nigig <- sapply(1:D, function(d) {
#       sampling(stan.nigig, data=list(n=n, p=p, y=y[, d], x=x, thetac=thetac, 
#                                      lambdac=lambdac, thetaz=thetaz, lambdaz), 
#                chains=1, warmup=nwarmup, iter=nsamples + nwarmup, cores=1, 
#                refresh=0, control=list(adapt_delta=0.8, max_treedepth=12))})
#     fit.studentst <- sapply(1:D, function(d) {
#       sampling(stan.studentst, data=list(n=n, p=p, y=y[, d], x=x, eta=eta, 
#                                          lambda=lambda), chains=1, 
#                warmup=nwarmup, iter=nsamples + nwarmup, cores=1, refresh=0, 
#                control=list(adapt_delta=0.8, max_treedepth=12))})
#     
#     best <- data.frame(nig=sapply(fit.nig, function(s) {
#       c(get_posterior_mean(s, pars="beta0"),
#         get_posterior_mean(s, pars="beta"))}),
#       nigig=sapply(fit.nigig, function(s) {
#         c(get_posterior_mean(s, pars="beta0"),
#           get_posterior_mean(s, pars="beta"))}),
#       studentst=sapply(fit.studentst, function(s) {
#         c(get_posterior_mean(s, pars="beta0"),
#           get_posterior_mean(s, pars="beta"))}))
#     
#     samp.nig <- cbind(as.matrix(fit.nig[[1]], "beta0"), 
#                       as.matrix(fit.nig[[1]], "beta"))
#     samp.nigig <- cbind(as.matrix(fit.nigig[[1]], "beta0"), 
#                         as.matrix(fit.nigig[[1]], "beta"))
#     samp.studentst <- cbind(as.matrix(fit.studentst[[1]], "beta0"), 
#                             as.matrix(fit.studentst[[1]], "beta"))
#     samp.true <- cbind(0, t(
#       gibbs.spikeandslab(x, y, set[s, "gamma"], set[s, "sigma"], set[s, "r"], 
#                          set[s, "nu"], set[s, "kappa"], init=NULL, 
#                          control=list(samples=nsamples, warmup=nwarmup))$beta))
#     samp <- data.frame(nig=samp.nig, nigig=samp.nigig, studentst=samp.studentst)
#     
#     cram <- c(nig=cramer.test(samp.nig, samp.true, 
#                               just.statistic=TRUE)$statistic,
#               nigig=cramer.test(samp.nigig, samp.true, 
#                                 just.statistic=TRUE)$statistic,
#               studentst=cramer.test(samp.studentst, samp.true, 
#                                     just.statistic=TRUE)$statistic)
#     
#     corbest <- c(nig=cor(c(0, beta), best[, 1]),
#                  nigig=cor(c(0, beta), best[, 2]),
#                  studentst=cor(c(0, beta), best[, 3]))
#     
#     msebest <- c(nig=mean((c(0, beta) - best[, 1])^2),
#                  nigig=mean((c(0, beta) - best[, 2])^2),
#                  studentst=mean((c(0, beta) - best[, 3])^2))
#     
#     ci.nig <- summary(fit.nig[[1]], pars=c("beta0", "beta"))$summary[, c(4, 8)]
#     ci.nigig <- summary(fit.nigig[[1]], 
#                         pars=c("beta0", "beta"))$summary[, c(4, 8)]
#     ci.studentst <- summary(fit.studentst[[1]], 
#                             pars=c("beta0", "beta"))$summary[, c(4, 8)]
#     
#     cover <- c(nig=mean(sapply(1:(p + 1), function(j) {
#       (ci.nig[j, 1] < c(0, beta)[j]) & (c(0, beta)[j] < ci.nig[j, 2])})),
#       nigig=mean(sapply(1:(p + 1), function(j) {
#         (ci.nigig[j, 1] < c(0, beta)[j]) & (c(0, beta)[j] < ci.nigig[j, 2])})),
#       studentst=mean(sapply(1:(p + 1), function(j) {
#         (ci.studentst[j, 1] < c(0, beta)[j]) & 
#           (c(0, beta)[j] < ci.studentst[j, 2])})))
#     
#     list(cram=cram, corbest=corbest, msebest=msebest, cover=cover, 
#          samp=samp)
#   }
# if(parallel) {stopCluster(cluster)}
# save(res, file="results/simulations_igaussian_res5.Rdata")
# load(file="results/simulations_igaussian_res5.Rdata")
# str(res[[1]][[1]])
# 
# fit <- Reduce("cbind", sapply(res, function(l) {as.matrix(l[[1]][[5]])},
#                               simplify=FALSE))
# dimnames(fit) <- list(paste0("sample", 1:nsamples),
#                       paste(rep(paste0("set", 1:nrow(set)),
#                                 each=(p + 1)*length(methods)),
#                             c(paste0(paste0("beta", 0:p), ".nig"),
#                               paste0(paste0("beta", 0:p), ".nigig"),
#                               paste0(paste0("beta", 0:p), ".studentst")),
#                             sep="."))
# 
# res <- Reduce("cbind", sapply(res, function(r) {
#   data.frame(cram=t(sapply(r, function(s) {s[[1]]})),
#              corbest=t(sapply(r, function(s) {s[[2]]})),
#              msebest=t(sapply(r, function(s) {s[[3]]})),
#              cover=t(sapply(r, function(s) {s[[4]]})))}))
# dimnames(res) <- list(paste0("rep", (1:nreps)), paste(rep(paste0(
#   "set", 1:nrow(set)), each=4*length(methods)), rep(apply(expand.grid(
#     methods, c("cram", "corbest", "msebest", "cover")), 1, function(s) {
#       paste(rev(s), collapse=".")}), times=nrow(set)), sep="."))
# 
# write.table(fit, file="results/simulations_igaussian_fit4.csv")
# write.table(res, file="results/simulations_igaussian_res4.csv")
# write.table(set, file="results/simulations_igaussian_set4.csv")


library(parallel)
###############################   simulation 5   ###############################
# settings
nreps <- 100
n <- 100
p <- 30
D <- 1

# create fixed parameters (eta chosen such that marginal beta variances match)
rho <- 0.5
Sigma <- matrix(rho, ncol=p, nrow=p); diag(Sigma) <- 1
gamma <- rep(1, D)
sigma <- rep(1, D)
nu <- 1
kappa <- c(nu, 9*nu)
r <- 0.01

V <- rep(c(1/2, 1, 2), times=2)
K <- rep(c(4, 10), each=3)
thetac <- V/sigma^2
lambdac <- c(3/2, 3, 6, 3/14, 3/7, 6/7)
ratio <- besselK(lambdac/thetac, 0, expon.scaled=TRUE)/
  besselK(lambdac/thetac, 1, expon.scaled=TRUE)
thetaz <- 1/(ratio + 2*thetac/lambdac)/sigma^2
lambdaz <- (3*thetaz^3 + 12*thetaz^2*thetac)/
  (K*lambdac - 3*thetaz^2 - 12*thetaz*thetac)
eta <- 6/(K - 3) + 4
lambda <- V*(eta - 2)/sigma^2

# check variances
vars <- cbind(nig=sigma^2*thetac, 
              nigig=sigma^2*thetac*thetaz*(ratio + 2*thetac/lambdac),
              studentst=sigma^2*lambda/(eta - 2))
all.equal(vars, cbind(nig=V, nigig=V, studentst=V))
kurts <- cbind(nig=3*thetac/lambdac + 3,
               nigig=(3*thetaz/lambdaz + 3)*
                 (4*thetac*(ratio + 2*thetac/lambdac) + 1)/
                 (lambdac*(ratio + 2*thetac/lambdac)^2),
               studentst=6/(eta - 4) + 3)
all.equal(kurts, cbind(nig=K, nigig=K, studentst=K))

params <- cbind(lambdac=lambdac, lambdaz=lambdaz, thetac=thetac, thetaz=thetaz, 
                eta=eta, lambda)

# store settings in an object
set <- data.frame(nreps=nreps, n=n, p=p, D=D, gamma=gamma, sigma=sigma, r=r,
                  nu=nu, kappa=rep(kappa, each=6), lambdac=lambdac, 
                  thetac=thetac, lambdaz=lambdaz, 
                  thetaz=thetaz, eta=eta, lambda=lambda, vars=V, kurts=K)
rownames(set) <- paste0("set", 1:nrow(set))

# objects to store simulation results
methods <- c("nig", "studentst", "nigig")
params <- c("beta0", "beta", "sigma", "gamma", "tau")

# set initial values and control parameters
nsamples <- 1000
nwarmup <- 1000

### analysis splits in parallel
# ncores <- min(detectCores() - 1, nreps)
ncores <- 50

res <- lapply(1:nrow(set), function(s) {
  mclapply(1:nreps, function(k) {
    print(paste0("set", s, ", rep", k))
    
    stan.nig <- stan_model("code/igauss.stan")
    stan.nigig <- stan_model("code/nigig.stan")
    stan.studentst <- stan_model("code/igamma.stan")
    
    set.seed(2019 + k)
    beta <- sapply(1:D, function(d) {
      rspikeandslab(p, set[s, "nu"], set[s, "kappa"], 
                    set[s, "sigma"], set[s, "gamma"], set[s, "r"])})
    x <- rmvnorm(n, rep(0, p), sigma=Sigma)
    y <- sapply(1:D, function(d) {rnorm(n, x %*% beta[, d], sigma[d])})
    
    thetac <- set$thetac[s]
    lambdac <- set$lambdac[s]
    thetaz <- set$thetaz[s]
    lambdaz <- set$lambdaz[s]
    eta <- set$eta[s]
    lambda <- set$lambda[s]
    fit.nig <- sapply(1:D, function(d) {
      sampling(stan.nig, data=list(n=n, p=p, y=y[, d], x=x, theta=thetac, 
                                   lambda=lambdac), chains=1, warmup=nwarmup,
               iter=nsamples + nwarmup, cores=1, refresh=0, 
               control=list(adapt_delta=0.8, max_treedepth=12))})
    fit.nigig <- sapply(1:D, function(d) {
      sampling(stan.nigig, data=list(n=n, p=p, y=y[, d], x=x, thetac=thetac, 
                                     lambdac=lambdac, thetaz=thetaz, lambdaz), 
               chains=1, warmup=nwarmup, iter=nsamples + nwarmup, cores=1, 
               refresh=0, control=list(adapt_delta=0.8, max_treedepth=12))})
    fit.studentst <- sapply(1:D, function(d) {
      sampling(stan.studentst, data=list(n=n, p=p, y=y[, d], x=x, eta=eta, 
                                         lambda=lambda), chains=1, 
               warmup=nwarmup, iter=nsamples + nwarmup, cores=1, refresh=0, 
               control=list(adapt_delta=0.8, max_treedepth=12))})
    
    best <- data.frame(nig=sapply(fit.nig, function(s) {
      c(get_posterior_mean(s, pars="beta0"),
        get_posterior_mean(s, pars="beta"))}),
      nigig=sapply(fit.nigig, function(s) {
        c(get_posterior_mean(s, pars="beta0"),
          get_posterior_mean(s, pars="beta"))}),
      studentst=sapply(fit.studentst, function(s) {
        c(get_posterior_mean(s, pars="beta0"),
          get_posterior_mean(s, pars="beta"))}))
    
    samp.nig <- cbind(as.matrix(fit.nig[[1]], "beta0"), 
                      as.matrix(fit.nig[[1]], "beta"))
    samp.nigig <- cbind(as.matrix(fit.nigig[[1]], "beta0"), 
                        as.matrix(fit.nigig[[1]], "beta"))
    samp.studentst <- cbind(as.matrix(fit.studentst[[1]], "beta0"), 
                            as.matrix(fit.studentst[[1]], "beta"))
    samp.true <- cbind(0, t(
      gibbs.spikeandslab(x, y, set[s, "gamma"], set[s, "sigma"], set[s, "r"], 
                         set[s, "nu"], set[s, "kappa"], init=NULL, 
                         control=list(samples=nsamples, warmup=nwarmup))$beta))
    samp <- data.frame(nig=samp.nig, nigig=samp.nigig, studentst=samp.studentst)
    
    cram <- c(nig=cramer.test(samp.nig, samp.true, 
                              just.statistic=TRUE)$statistic,
              nigig=cramer.test(samp.nigig, samp.true, 
                                just.statistic=TRUE)$statistic,
              studentst=cramer.test(samp.studentst, samp.true, 
                                    just.statistic=TRUE)$statistic)
    
    corbest <- c(nig=cor(c(0, beta), best[, 1]),
                 nigig=cor(c(0, beta), best[, 2]),
                 studentst=cor(c(0, beta), best[, 3]))
    
    msebest <- c(nig=mean((c(0, beta) - best[, 1])^2),
                 nigig=mean((c(0, beta) - best[, 2])^2),
                 studentst=mean((c(0, beta) - best[, 3])^2))
    
    ci.nig <- summary(fit.nig[[1]], pars=c("beta0", "beta"))$summary[, c(4, 8)]
    ci.nigig <- summary(fit.nigig[[1]], 
                        pars=c("beta0", "beta"))$summary[, c(4, 8)]
    ci.studentst <- summary(fit.studentst[[1]], 
                            pars=c("beta0", "beta"))$summary[, c(4, 8)]
    
    cover <- c(nig=mean(sapply(1:(p + 1), function(j) {
      (ci.nig[j, 1] < c(0, beta)[j]) & (c(0, beta)[j] < ci.nig[j, 2])})),
      nigig=mean(sapply(1:(p + 1), function(j) {
        (ci.nigig[j, 1] < c(0, beta)[j]) & (c(0, beta)[j] < ci.nigig[j, 2])})),
      studentst=mean(sapply(1:(p + 1), function(j) {
        (ci.studentst[j, 1] < c(0, beta)[j]) & 
          (c(0, beta)[j] < ci.studentst[j, 2])})))
    
    list(cram=cram, corbest=corbest, msebest=msebest, cover=cover, 
         samp=samp)
  }, mc.cores=ncores)
})
save(res, file="results/simulations_igaussian_res5.Rdata")
