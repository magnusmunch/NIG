#!/usr/bin/env Rscript

### installation of packages
# if(!("devtools" %in% installed.packages())) {
#   install.packages("devtools")
# }
# library(devtools)
# install_github("magnusmunch/cambridge/rpackage", local=FALSE,
#                auth_token=Sys.getenv("GITHUB_PAT"))

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
lambda <- 0.5
theta <- 1
eta <- lambda/theta + 2

# store settings in an object
set <- data.frame(nreps=nreps, n=n, p=p, D=D, gamma=t(as.matrix(gamma)), 
                  sigma=t(as.matrix(sigma)), 
                  lambda=lambda, theta=theta, eta=eta)
rownames(set) <- paste0("set", 1:nrow(set))

# objects to store simulation results
methods <- c("MCMC", "ADVI", "VB", "MAP", "glmnet")
params <- c("beta0", "beta", "sigma", "gamma")

# set initial values and control parameters
stan.igauss <- stan_model("code/igauss.stan")
stan.igamma <- stan_model("code/igamma.stan")
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

  fit.igauss <- sapply(1:D, function(d) {
    sampling(stan.igauss, data=list(n=n, p=p, y=y[, d], x=x, theta=theta, 
                                    lambda=lambda), chains=1, warmup=nwarmup,
             iter=nsamples + nwarmup, cores=1, refresh=0, 
             control=list(adapt_delta=0.8, max_treedepth=12))})
  fit.igamma <- sapply(1:D, function(d) {
    sampling(stan.igamma, data=list(n=n, p=p, y=y[, d], x=x, eta=eta, 
                                    lambda=lambda), chains=1, warmup=nwarmup,
             iter=nsamples + nwarmup, cores=1, refresh=0, 
             control=list(adapt_delta=0.8, max_treedepth=12))})
  
  best <- data.frame(igauss=sapply(fit.igauss, function(s) {
    c(get_posterior_mean(s, pars="beta0"),
      get_posterior_mean(s, pars="beta"))}),
    igamma=sapply(fit.igamma, function(s) {
      c(get_posterior_mean(s, pars="beta0"),
        get_posterior_mean(s, pars="beta"))}))

  summary(fit.igauss[[1]])$summary[1:(p + 1), 3]^2
  summary(fit.igamma[[1]])$summary[1:(p + 1), 3]^2
  
  bvar <- cbind(igauss=sapply(fit.igauss, function(s) {
    summary(s)$summary[1:(p + 1), 3]^2}),
  igamma=sapply(fit.igamma, function(s) {
    summary(s)$summary[1:(p + 1), 3]^2}))
  
  samp.igauss <- cbind(as.matrix(fit.igauss[[1]], "beta0"), 
                       as.matrix(fit.igauss[[1]], "beta"))
  samp.igamma <- cbind(as.matrix(fit.igamma[[1]], "beta0"), 
                       as.matrix(fit.igamma[[1]], "beta"))
  samp <- data.frame(igauss=samp.igauss, igamma=samp.igamma)

  cram <- cramer.test(samp.igauss, samp.igamma, just.statistic=TRUE)$statistic
  corbvar <- cor(bvar[, 1], bvar[, 2])

  corbest <- c(igauss=cor(c(0, beta), best[, 1]),
               igamma=cor(c(0, beta), best[, 2]))

  msebest <- c(igauss=mean((c(0, beta) - best[, 1])^2),
               igamma=mean((c(0, beta) - best[, 2])^2))

  out <- list(cram=cram, corbvar=corbvar, 
              corbest=corbest, msebest=msebest, samp)
}
if(parallel) {stopCluster(cluster)}

fit <- res[[1]][[5]]
dimnames(fit) <- list(paste0("sample", 1:nsamples), 
                      c(paste0(paste0("beta", 0:p), ".igauss"),
                        paste0(paste0("beta", 0:p), ".igamma")))
res <- data.frame(cram=t(sapply(res, function(s) {s[[1]]})),
                  corbvar=t(sapply(res, function(s) {s[[2]]})),
                  corbest=t(sapply(res, function(s) {s[[3]]})),
                  msebest=t(sapply(res, function(s) {s[[4]]})))
rownames(res) <- paste0("rep", (1:nreps))
write.table(fit, file="results/simulations_igaussian_fit4.csv")
write.table(res, file="results/simulations_igaussian_res4.csv")
write.table(set, file="results/simulations_igaussian_set4.csv")



