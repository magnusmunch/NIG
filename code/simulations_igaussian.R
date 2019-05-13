#!/usr/bin/env Rscript

### installation of packages
if(!("devtools" %in% installed.packages())) {
  install.packages("devtools")
}
library(devtools)
install_github("magnusmunch/cambridge/rpackage", local=FALSE,
               auth_token=Sys.getenv("GITHUB_PAT"))
~/Library/R/3.4.2/library
### parallelisation
parallel <- FALSE

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

################################# simulation 1 #################################
# settings
nreps <- 100
n <- 100
p <- 100
D <- 100
nclass <- 4

# create fixed parameters
alpha <- c(1:nclass)
C.inv.gauss <- model.matrix(~ -1 + factor(rep(c(1:nclass), each=D/nclass)))
C.inv.gamma <- as.numeric(C.inv.gauss %*% c(1:nclass))
theta <- as.numeric(1/(C.inv.gauss %*% alpha))
lambda <- rep(1, D)
sigma <- rep(1, D)
SNR <- 50

# store settings in an object
temp1 <- c(nreps=nreps, n=n, p=p, D=D, nclass=nclass, alpha=alpha, 
           C=C.inv.gauss, theta=theta, lambda=lambda, sigma=sigma, SNR=SNR)
set1 <- t(matrix(temp1))
colnames(set1) <- names(temp1)

# objects to store simulation results
methods <- c("igauss.conj", "igauss.non.conj", "igamma.conj", "igamma.non.conj")
params <- c("alpha", "eta", "theta", "lambda", "mu", "dSigma", "delta", "zeta")
drugs <- paste("drug", c(1:D), sep="")
covs <- paste("cov", c(1:nclass), sep="")
res1 <- as.data.frame(
  matrix(NA, nrow=nreps, ncol=length(methods)*length(covs) + 
           length(methods)*length(params[-1])*length(drugs),
         dimnames=list(paste("rep", c(1:nreps), sep=""), 
                       c(apply(expand.grid(covs, methods, params[1]), 1, paste, 
                               collapse="_"), 
                         apply(expand.grid(drugs, methods, params[-1]), 1, 
                               paste, collapse="_")))))

# set initial values and control parameters
control <- list(epsilon.eb=-1, epsilon.vb=-1, 
                epsilon.opt=sqrt(sqrt(.Machine$double.eps)), 
                maxit.eb=20, maxit.vb=2, maxit.opt=100, trace=FALSE)
init.inv.gauss <- list(a=rep(mean(1/sigma^2), D), b=rep(mean(theta), D), 
                       theta=theta, lambda=rep(1, D))
init.inv.gamma <- list(a=rep(mean(1/sigma^2), D), b=rep(mean(theta), D),  
                       eta=rep(theta, D), lambda=rep(1, D))

# simulation
set.seed(2018)
for(r in 1:nreps) {

  cat("\r", "replication", r)
  
  # create data
  gamma <- sqrt(rinvgauss(D, theta, shape=lambda))
  beta <- sapply(1:D, function(d) {rnorm(p, 0, sigma[d]*gamma[d])})
  x <- matrix(rnorm(n*p, 0, sqrt(SNR/(mean(theta)*p))), nrow=n, ncol=p)
  y <- sapply(1:D, function(d) {rnorm(n, x %*% beta[, d], sigma[d])})
  
  # fit inverse Gaussian models
  fit1.igauss.conj <- est.model(x, y, C.inv.gauss, "inv. Gaussian", TRUE, 
                                init=init.inv.gauss, control=control)
  fit1.igauss.non.conj <- est.model(x, y, C.inv.gauss, "inv. Gaussian", FALSE, 
                                    init=init.inv.gauss, control=control)
  
  # fit independent inverse Gamma models
  fit1.igamma.conj <- est.model(x, y, C.inv.gamma, "inv. Gamma", TRUE, 
                                init=init.inv.gamma, control=control)
  fit1.igamma.non.conj <- est.model(x, y, C.inv.gamma, "inv. Gamma", FALSE, 
                                    init=init.inv.gamma, control=control)
  
  # store results
  temp1 <- paste("fit1.", apply(expand.grid(methods, "seq.eb", params[1]), 1, 
                                paste, collapse="$"), "[fit1.", methods, 
                 "$iter$eb", ", ]", sep="")
  temp2 <- paste("fit1.", apply(expand.grid(methods, "seq.eb", params[2:4]), 1, 
                                paste, collapse="$"), "[fit1.", methods, 
                 "$iter$eb", ", ]", sep="")
  temp3 <- paste("fit1.", apply(expand.grid(methods, "vb.post", params[5]), 1, 
                                paste, collapse="$"), sep="")
  temp4 <- paste("fit1.", apply(expand.grid(methods, "vb.post", params[6]), 1, 
                                paste, collapse="$"), sep="")
  temp5 <- paste("fit1.", apply(expand.grid(methods, "vb.post", params[7:8]), 1, 
                                paste, collapse="$"), sep="")
  res1[r, ] <- c(c(sapply(temp1, function(s) {eval(parse(text=s))})), 
                 c(sapply(temp2, function(s) {eval(parse(text=s))})),
                 c(sapply(temp3, function(s) {
                   colMeans((eval(parse(text=s)) - beta)^2)})),
                 c(sapply(temp4, function(s) {colMeans(eval(parse(text=s)))})),
                 c(sapply(temp5, function(s) {eval(parse(text=s))})))

  
  # save fitted models and results
  fit1 <- data.frame(igauss.conj=fit1.igauss.conj$seq.eb,
                     igauss.non.conj=fit1.igauss.non.conj$seq.eb,
                     igamma.conj=fit1.igamma.conj$seq.eb,
                     igamma.non.conj=fit1.igamma.non.conj$seq.eb)
  write.table(set1, file="results/simulations_igaussian_set1.csv")
  write.table(res1, file="results/simulations_igaussian_res1.csv")
  write.table(fit1, file="results/simulations_igaussian_fit1.csv")

}

###############################   simulation 2   ###############################
# settings
nreps <- 100
n <- 100
p <- 100
D <- 100
nclass <- 4

# create fixed parameters
alpha <- c(1:nclass) + 2
C.inv.gauss <- model.matrix(~ -1 + factor(rep(c(1:nclass), each=D/nclass)))
C.inv.gamma <- as.numeric(C.inv.gauss %*% c(1:nclass))
eta <- as.numeric(C.inv.gauss %*% alpha)
lambda <- rep(1, D)
sigma <- rep(1, D)
SNR <- 50

# store settings in an object
temp1 <- c(nreps=nreps, n=n, p=p, D=D, nclass=nclass, alpha=alpha, 
           C=C.inv.gauss, eta=eta, lambda=lambda, sigma=sigma, SNR=SNR)
set2 <- t(matrix(temp1))
colnames(set2) <- names(temp1)

# objects to store simulation results
methods <- c("igauss.conj", "igauss.non.conj", "igamma.conj", "igamma.non.conj")
params <- c("alpha", "eta", "theta", "lambda", "mu", "dSigma", "delta", "zeta")
drugs <- paste("drug", c(1:D), sep="")
covs <- paste("cov", c(1:nclass), sep="")
res2 <- as.data.frame(
  matrix(NA, nrow=nreps, ncol=length(methods)*length(covs) + 
           length(methods)*length(params[-1])*length(drugs),
         dimnames=list(paste("rep", c(1:nreps), sep=""), 
                       c(apply(expand.grid(covs, methods, params[1]), 1, paste, 
                               collapse="_"), 
                         apply(expand.grid(drugs, methods, params[-1]), 1, 
                               paste, collapse="_")))))

# set initial values and control parameters
control <- list(epsilon.eb=-1, epsilon.vb=-1, 
                epsilon.opt=sqrt(sqrt(.Machine$double.eps)), 
                maxit.eb=30, maxit.vb=4, maxit.opt=100, trace=FALSE)
init.inv.gauss <- list(a=rep(mean(1/sigma^2), D), 
                       b=rep(mean(lambda/(eta - 2)), D), 
                       theta=rep(mean(lambda/(eta - 2)), D), lambda=rep(1, D))
init.inv.gamma <- list(a=rep(mean(1/sigma^2), D), 
                       b=rep(mean(lambda/(eta - 2)), D),  
                       eta=eta, lambda=rep(1, D))

# simulation
set.seed(2018)
for(r in 1:nreps) {
  
  cat("\r", "replication", r)
  
  # create data
  gamma <- 1/sqrt(rgamma(D, shape=eta/2, scale=lambda/2))
  beta <- sapply(1:D, function(d) {rnorm(p, 0, sigma[d]*gamma[d])})
  x <- matrix(rnorm(n*p, 0, sqrt(SNR/(mean(lambda/(eta - 2))*p))), nrow=n, 
              ncol=p)
  y <- sapply(1:D, function(d) {rnorm(n, x %*% beta[, d], sigma[d])})
  
  # fit inverse Gaussian models
  fit2.igauss.conj <- est.model(x, y, C.inv.gauss, "inv. Gaussian", TRUE, 
                                init=init.inv.gauss, control=control)
  fit2.igauss.non.conj <- est.model(x, y, C.inv.gauss, "inv. Gaussian", FALSE, 
                                    init=init.inv.gauss, control=control)
  
  # fit independent inverse Gamma models
  fit2.igamma.conj <- est.model(x, y, C.inv.gamma, "inv. Gamma", TRUE, 
                                init=init.inv.gamma, control=control)
  fit2.igamma.non.conj <- est.model(x, y, C.inv.gamma, "inv. Gamma", FALSE, 
                                    init=init.inv.gamma, control=control)
  
  # store results
  temp1 <- paste("fit2.", apply(expand.grid(methods, "seq.eb", params[1]), 1, 
                                paste, collapse="$"), "[fit2.", methods, 
                 "$iter$eb", ", ]", sep="")
  temp2 <- paste("fit2.", apply(expand.grid(methods, "seq.eb", params[2:4]), 1, 
                                paste, collapse="$"), "[fit2.", methods, 
                 "$iter$eb", ", ]", sep="")
  temp3 <- paste("fit2.", apply(expand.grid(methods, "vb.post", params[5]), 1, 
                                paste, collapse="$"), sep="")
  temp4 <- paste("fit2.", apply(expand.grid(methods, "vb.post", params[6]), 1, 
                                paste, collapse="$"), sep="")
  temp5 <- paste("fit2.", apply(expand.grid(methods, "vb.post", params[7:8]), 1, 
                                paste, collapse="$"), sep="")
  res2[r, ] <- c(c(sapply(temp1, function(s) {eval(parse(text=s))})), 
                 c(sapply(temp2, function(s) {eval(parse(text=s))})),
                 c(sapply(temp3, function(s) {
                   colMeans((eval(parse(text=s)) - beta)^2)})),
                 c(sapply(temp4, function(s) {colMeans(eval(parse(text=s)))})),
                 c(sapply(temp5, function(s) {eval(parse(text=s))})))
  
  # save fitted models and results
  fit2 <- data.frame(igauss.conj=fit2.igauss.conj$seq.eb,
                     igauss.non.conj=fit2.igauss.non.conj$seq.eb,
                     igamma.conj=fit2.igamma.conj$seq.eb,
                     igamma.non.conj=fit2.igamma.non.conj$seq.eb)
  write.table(set2, file="results/simulations_igaussian_set2.csv")
  write.table(res2, file="results/simulations_igaussian_res2.csv")
  write.table(fit2, file="results/simulations_igaussian_fit2.csv")

}

warnings()


###############################   simulation 3   ###############################
# settings
nreps <- 100
n <- 100
p <- 200
D <- 1

# create fixed parameters
rho <- 0.5
Sigma <- matrix(rho, ncol=p, nrow=p); diag(Sigma) <- 1
gamma <- rep(1, D)
sigma <- rep(1, D)
lambda <- 0.1
theta <- 2

# store settings in an object
set <- data.frame(nreps=nreps, n=n, p=p, D=D, gamma=t(as.matrix(gamma)), 
                  sigma=t(as.matrix(sigma)), 
                  lambda=lambda, theta=theta)
rownames(set) <- paste0("set", 1:nrow(set))

# objects to store simulation results
methods <- c("MCMC", "ADVI", "VB", "MAP", "glmnet")
params <- c("beta0", "beta", "sigma", "gamma")

# set initial values and control parameters
stan.model <- stan_model("code/igauss.stan")
nsamples <- 1000
nwarmup <- 1000

### analysis splits in parallel
ncores <- min(detectCores() - 1, nreps)
cluster <- makeForkCluster(ncores)
if(parallel) {
  registerDoParallel(cluster)
} else {
  registerDoSEQ()
}

res <- foreach(k=c(1:nreps), .errorhandling="pass") %dopar% {
  print(paste("rep", k))
  
  set.seed(2019 + k)
  beta <- sapply(1:D, function(d) {rnorm(p, 0, sigma[d]*gamma[d])})
  x <- rmvnorm(n, rep(0, p), sigma=Sigma)
  y <- sapply(1:D, function(d) {rnorm(n, x %*% beta[, d], sigma[d])})

  fit.mcmc <- sapply(1:D, function(d) {
    sampling(stan.model, data=list(n=n, p=p, y=y[, d], x=x, theta=theta, 
                                   lambda=lambda), chains=1, warmup=nwarmup,
             iter=nsamples + nwarmup, cores=1, refresh=0, 
             control=list(adapt_delta=0.8, max_treedepth=12))})
  fit.advi <- sapply(1:D, function(d) {
    vb(stan.model, data=list(n=n, p=p, y=y[, d], x=x, theta=theta, 
                             lambda=lambda), output_samples=nsamples)})
  fit.vb <- est.vb(x, y, theta, lambda, intercept=TRUE, init=NULL, 
                   posterior=TRUE,
                   control=list(epsilon.vb=1e-6, maxit.vb=100, trace=TRUE))
  fit.glmnet <- sapply(1:D, function(d) {glmnet(x, y[, d], alpha=0,
                                                lambda=1/(n*theta*sigma[d]^2))}, 
                       simplify=FALSE)
  fit.map <- sapply(1:D, function(d) {
    optimizing(stan.model, data=list(n=n, p=p, y=y[, d], x=x, theta=theta, 
                                     lambda=lambda))}, simplify=FALSE)
  
  best <- cbind(mcmc=sapply(fit.mcmc, function(s) {
    c(get_posterior_mean(s, pars="beta0"),
      get_posterior_mean(s, pars="beta"))}),
    advi=sapply(fit.advi, function(s) {
      c(get_posterior_mean(s, pars="beta0"),
        get_posterior_mean(s, pars="beta"))}),
    vb=fit.vb$vb.post$mu,
    glmnet=sapply(fit.glmnet, function(s) {as.numeric(coef(s))}),
    map=sapply(fit.map, function(s) {s$par[1:(p + 1)]}))
  
  bvar <- cbind(mcmc=sapply(fit.mcmc, function(s) {
    summary(s)$summary[1:(p + 1), 3]^2}),
    advi=sapply(fit.advi, function(s) {
      summary(s)$summary[1:(p + 1), 2]^2}),
    vb=sapply(1:D, function(d) {diag(fit.vb$vb.post$Sigma[[d]])}))

  samp.advi <- matrix(unlist(fit.advi[[1]]@sim$samples[[1]][1:(p + 1)]), 
                      ncol=p + 1, byrow=FALSE)
  samp.vb <- rmvnorm(nsamples, mean=as.numeric(fit.vb$vb.post$mu), 
                     sigma=fit.vb$vb.post$Sigma[[1]])
  samp.mcmc <- cbind(as.matrix(fit.mcmc[[1]], "beta0"), 
                     as.matrix(fit.mcmc[[1]], "beta"))
  samp <- cbind(samp.mcmc, samp.advi, samp.vb)

  cram <- c(advi=cramer.test(samp.advi, samp.mcmc, 
                             just.statistic=TRUE)$statistic,
            vb=cramer.test(samp.vb, samp.mcmc, just.statistic=TRUE)$statistic)

  corbvar <- c(advi=cor(bvar[, 2], bvar[, 1]),
               vb=cor(bvar[, 3], bvar[, 1]))
  
  msebvar <- c(advi=mean((bvar[, 2] - bvar[, 1])^2),
               vb=mean((bvar[, 3] - bvar[, 1])^2))

  corbest <- c(advi=cor(best[, 2], best[, 1]),
               vb=cor(best[, 3], best[, 1]),
               glmnet=cor(best[, 4], best[, 1]),
               map=cor(best[, 5], best[, 1]))

  msebest <- c(advi=mean((best[, 2] - best[, 1])^2),
               vb=mean((best[, 3] - best[, 1])^2),
               glmnet=mean((best[, 4] - best[, 1])^2),
               map=mean((best[, 5] - best[, 1])^2))
  
  out <- list(cram, corbvar, msebvar, corbest, msebest, samp)
}
if(parallel) {stopCluster(cluster)}
save(res, file="results/simulations_igaussian_res3.Rdata")

which.excl <- 80
fit <- res[[1]][[6]]
dimnames(fit) <- list(paste0("sample", 1:nsamples), 
                      c(paste0(paste0("beta", 0:p), ".mcmc"),
                        paste0(paste0("beta", 0:p), ".advi"),
                        paste0(paste0("beta", 0:p), ".vb")))
res <- data.frame(cram=t(sapply(res[-which.excl], function(s) {s[[1]]})),
                  corbvar=t(sapply(res[-which.excl], function(s) {s[[2]]})),
                  msebvar=t(sapply(res[-which.excl], function(s) {s[[3]]})),
                  corbest=t(sapply(res[-which.excl], function(s) {s[[4]]})),
                  msebest=t(sapply(res[-which.excl], function(s) {s[[5]]})))
rownames(res) <- paste0("rep", (1:nreps)[-which.excl])
write.table(fit, file="results/simulations_igaussian_fit3.csv")
write.table(res, file="results/simulations_igaussian_res3.csv")
write.table(set, file="results/simulations_igaussian_set3.csv")


###############################   simulation 4   ###############################
# settings
nreps <- 100
n <- 100
p <- 200
D <- 1

# create fixed parameters (eta chosen such that marginal beta variances match)
rho <- 0.5
Sigma <- matrix(rho, ncol=p, nrow=p); diag(Sigma) <- 1
gamma <- rep(1, D)
sigma <- rep(1, D)
thetac <- 1
thetaz <- ratio_besselK(1/3, 1) + 2/3
lambdac <- 1/3
lambdaz <- (36*ratio_besselK(1/3, 1)^2 + 249*ratio_besselK(1/3, 1) + 150)/
  (ratio_besselK(1/3, 1) + 6)^2
eta <- 7/3

# store settings in an object
set <- data.frame(nreps=nreps, n=n, p=p, D=D, gamma=t(as.matrix(gamma)), 
                  sigma=t(as.matrix(sigma)), t(as.matrix(tau)),
                  lambdac=lambdac, thetac=thetac,
                  lambdaz=lambdaz, thetaz=thetaz, eta=eta)
rownames(set) <- paste0("set", 1:nrow(set))

# objects to store simulation results
methods <- c("nig", "studentst", "nigig")
params <- c("beta0", "beta", "sigma", "gamma", "tau")

# set initial values and control parameters
stan.nig <- stan_model("code/igauss.stan")
stan.nigig <- stan_model("code/nigig.stan")
stan.studentst <- stan_model("code/igamma.stan")
nsamples <- 1000
nwarmup <- 1000

### analysis splits in parallel
ncores <- min(detectCores() - 1, nreps)
cluster <- makeForkCluster(ncores)
if(parallel) {
  registerDoParallel(cluster)
} else {
  registerDoSEQ()
}

res <- foreach(k=c(1:nreps), .errorhandling="pass") %dopar% {
  print(paste("rep", k))
  
  set.seed(2019 + k)
  beta <- sapply(1:D, function(d) {rnorm(p, 0, sigma[d]*gamma[d])})
  x <- rmvnorm(n, rep(0, p), sigma=Sigma)
  y <- sapply(1:D, function(d) {rnorm(n, x %*% beta[, d], sigma[d])})

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
                                       lambda=lambdac), chains=1, 
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
  cov.true <- sigma^2*solve(t(x) %*% x + gamma^(-2)*diag(p))
  mean.true <- cov.true %*% t(x) %*% y/sigma^2
  samp.true <- cbind(0, rmvnorm(nsamples, mean.true, cov.true))
  samp <- data.frame(nig=samp.nig, nigig=samp.nigig, studentst=samp.studentst)

  cram <- c(nig=cramer.test(samp.nig, samp.true, just.statistic=TRUE)$statistic,
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

  out <- list(cram=cram, corbest=corbest, msebest=msebest, cover, samp)
}
if(parallel) {stopCluster(cluster)}
save(res, file="results/simulations_igaussian_res4.Rdata")

fit <- res[[1]][[5]]
dimnames(fit) <- list(paste0("sample", 1:nsamples), 
                      c(paste0(paste0("beta", 0:p), ".igauss"),
                        paste0(paste0("beta", 0:p), ".igamma")))
res <- data.frame(cram=t(sapply(res, function(s) {s[[1]]})),
                  corbest=t(sapply(res, function(s) {s[[2]]})),
                  msebest=t(sapply(res, function(s) {s[[3]]})),
                  cover=t(sapply(res, function(s) {s[[4]]})))
rownames(res) <- paste0("rep", (1:nreps))
write.table(fit, file="results/simulations_igaussian_fit4.csv")
write.table(res, file="results/simulations_igaussian_res4.csv")
write.table(set, file="results/simulations_igaussian_set4.csv")


### spike and slab model
rspikeandslab <- function(n, nu, kappa, sigma, gamma) {
  z <- rbinom(n, 1, nu/(nu + kappa))
  beta <- rnorm(n, 0, sigma*gamma*z)
  return(beta)
}


.samp.z <- function(z, w, x, y, theta, sigma, p) {
  sorder <- sample(1:p)
  for(i in sorder) {
    z0 <- replace(z, i, 0)
    z1 <- replace(z, i, 1)
    x1 <- x[, z1==1]
    svd.x1 <- svd(x1)
    sv1 <- svd.x1$d
    uty1 <- as.numeric(t(svd.x1$u) %*% y)
    if(sum(z0)==0) {
      R <- exp(0.5*sum(log(sv1^2 + 1/theta)) - 
                 0.5*sum(sv1^2*uty1/(sv1^2 + 1/theta))/sigma^2 - p/2*log(theta))  
    } else {
      x0 <- x[, z0==1]
      svd.x0 <- svd(x0)
      sv0 <- svd.x0$d
      uty0 <- as.numeric(t(svd.x0$u) %*% y)
      R <- exp(0.5*sum(log(sv1^2 + 1/theta)) - 0.5*sum(log(sv0^2 + 1/theta)) -
                 0.5*sum(sv1^2*uty1/(sv1^2 + 1/theta))/sigma^2 + 
                 0.5*sum(sv0^2*uty0/(sv0^2 + 1/theta))/sigma^2)  
    }
    z[i] <- rbinom(1, 1, 1/(1 + R*(1 - w)/w))
  }
  return(z)
}

.samp.w <- function(z, nu, kappa, p) {
  w <- rbeta(1, nu + sum(z), kappa + p - sum(z))
  return(w)
}

svvt=t(v)*sv
zp=sum(z)
x=x[, z==1]
.samp.beta <- function(theta, sigma, v, sv, svvt, x, y, n, zp, p) {
  
  if(zp > p) {
    uvec <- rnorm(p, 0, sqrt(theta)*sigma)
    dvec <- rnorm(n, 0, 1)
    beta <- uvec + as.numeric(t(t(v)*(sv/(sv^2 + 1/theta))) %*% 
                                (sigma*y - sigma*svvt %*% uvec + dvec))  
  } else {
    mat <- solve(t(x) %*% x + diag(zp)/theta)
    beta <- as.numeric(rmvnorm(1, mat %*% t(x) %*% y, sigma^2*mat))
  }
  
  return(beta)
}

theta=1; sigma=1; nu=1; kappa=1; init=NULL
control=list(nsamp=100)
gibbs.spikeandslab <- function(x, y, theta, sigma, nu, kappa, init=NULL,
                               control=list(nsamp=1000)) {
  
  
  n <- nrow(x)
  p <- ncol(x)
  
  if(is.null(init$beta)) {
    beta <- rnorm(p)
  } else {
    beta <- init$beta
  }
  if(is.null(init$z)) {
    z <- rep(1, p)  
  } else {
    z <- init$z
  }
  if(is.null(init$w)) {
    w <- 0.5
  } else {
    w <- init$w
  }
  
  out <- list(beta=matrix(NA, nrow=p, ncol=control$nsamp), 
              z=matrix(NA, nrow=p, ncol=control$nsamp), 
              w=numeric(control$nsamp))
  for(m in 1:control$nsamp) {
    z <- .samp.z(z, w, x, y, theta, sigma, p)
    w <- .samp.w(z, nu, kappa, p)
    
    # if the effective dimension of z is smaller than p, sampling is direct
    if(sum(z) > p) {
      svd.x <- svd(x[, z==1])
      v <- svd.x$v
      sv <- svd.x$d
    } else {
      v <- sv <- NULL
    }
    beta <- rep(0, p)
    if(sum(z) > 0) {
      beta[z==1] <- .samp.beta(theta, sigma, v, sv, t(v)*sv, x[, z==1], y, n, 
                               sum(z), p)  
    }
    out$beta[, m] <- beta
    out$z[, m] <- z
    out$w[m] <- w
  }
  
  return(out)
}


set.seed(123)
n <- 100
p <- 20
x <- matrix(rnorm(n*p), nrow=n, ncol=p)
beta <- rnorm(p, sd=10)
y <- rnorm(n, as.numeric(x %*% beta), 1)

test <- gibbs.spikeandslab(x, y, theta=10, sigma=1, nu=1, kappa=10, 
                           control=list(nsamp=1000))
plot(abs(beta), rowMeans(test$z))


### nigig model
fchisq <- function(x, lambdac, lambdaz, thetac, thetaz) {
  sqrt(lambdac*lambdaz/((2*pi)^2*x^3))*exp(lambdac/thetac + lambdaz/thetaz)*
    2*besselK(2*sqrt(lambdac/(2*thetac^2) + lambdaz/(2*x))*
                sqrt(lambdac/2 + lambdaz*x/(2*thetaz^2)), 0)
}

fdint <- function(x, lambda, a, b, sigma, beta, expon.scaled=FALSE) {
  if(expon.scaled) {
    exp(-a*beta^2/(2*lambda*sigma^2*x) - sqrt(b^2 - 2*a + a*(x + 1/x)) - 
          2*log(x))*besselK(sqrt(b^2 - 2*a + a*(x+ 1/x)), 0, expon.scaled=TRUE)
  } else {
    (1/x^2)*exp(-a*beta^2/(2*lambda*sigma^2*x))*
      besselK(sqrt(b^2 - 2*a + a*(x + 1/x)), 0)  
  }
  
}

fvint <- function(x, a, b) {
  1/sqrt(x)*besselK(sqrt(b^2 - 2*a + a*(x + 1/x)), 0) 
}


dnigig <- function(x, lambda, a, b, sigma, expon.scaled=FALSE) {
  a*exp(b)/sqrt(2*sigma^2*lambda*pi^3)*
    sapply(x, function(s) {
      integrate(fdint, 0, Inf, lambda=lambda, a=a, b=b, sigma=sigma, 
                beta=s, expon.scaled=expon.scaled)$value})
}

vnigig <- function(lambda, a, b, sigma) {
  exp(b)*sigma^2*lambda/(pi*sqrt(a))*
    integrate(fvint, 0, Inf, a=a, b=b)$value
}

vnigig <- function(thetaz, thetac, lambdac, sigma) {
  sigma^2*thetaz*thetac*(ratio_besselK(lambdac/thetac, 1) + 2*thetac/lambdac)
}

set.seed(123)
n <- 100
p <- 20
x <- matrix(rnorm(n*p), nrow=n, ncol=p)
beta <- rnorm(p, 0, 1)
y <- rbinom(n, 1, 1/(1 + exp(-x %*% beta)))
lambda <- 1
a <- 1
b <- 1
sigma <- 1
curve(dnigig(x, lambda=lambda, a=a, b=b, sigma=sigma, expon.scaled=TRUE), 
      -10, 10, n=1000)
vnigig(lambda, a, b, sigma)


lambdac=10; lambdaz=4; thetac=3; thetaz=2
curve(fchisq(x, lambdac=lambdac, lambdaz=lambdaz, thetac=thetac,
             thetaz=thetaz), 0.001, 10)
integrate(fchisq, 0, Inf, lambdac=lambdac, lambdaz=lambdaz, 
          thetac=thetac, thetaz=thetaz)




