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
library(statmod)
library(GeneralizedHyperbolic)
library(rstan)
library(glmnet)
library(mvtnorm)
library(cramer)

################################# simulation 1 #################################
# settings
nreps <- 100
n <- 200
p <- 100
D <- 100
nclass <- 4

# create fixed parameters
alpha <- c(1:nclass)
C <- lapply(1:D, function(d) {
  mat <- model.matrix(~ factor(rep(1:nclass, each=p/nclass)))[, -1]
  colnames(mat) <- paste0("class", 2:nclass)
  return(mat)})
prior.mean <- lapply(C, function(d) {as.numeric(1/(cbind(1, d) %*% alpha))})
lambda <- rep(1, D)
sigma <- rep(1, D)
SNR <- 50

methods <- "enig"
params <- c(paste0("alpha", 1:length(alpha)), 
            # paste0("lambda", 1:length(lambda)))
            "lambda")

# store settings in an object
set <- data.frame(as.list(c(nreps=nreps, n=n, p=p, D=D, nclass=nclass, 
                            alpha=alpha, lambda=lambda, sigma=sigma, SNR=SNR)))
rownames(set) <- paste0("set", 1:nrow(set))

# objects to store simulation results
res <- matrix(NA, nrow=nreps, ncol=nrow(set)*length(methods)*length(params), 
               dimnames=list(paste0("rep", (1:nreps)), 
                             paste0(rownames(set), ".", methods, ".", params)))

# set initial values and control parameters
control <- list(conv.post=FALSE, trace=FALSE,
                epsilon.eb=1e-3, epsilon.vb=1e-3, 
                epsilon.opt=sqrt(.Machine$double.eps),
                maxit.eb=200, maxit.vb=2, maxit.opt=100,
                maxit.post=100)

# simulation
set.seed(2019)
for(r in 1:nreps) {

  cat("\r", "replication", r)
  # create data
  gamma <- lapply(1:D, function(d) {
    sqrt(rinvgauss(p, prior.mean[[d]], shape=lambda[d]))})
  beta <- lapply(1:D, function(d) {rnorm(p, 0, sigma[d]*gamma[[d]])})
  x <- lapply(1:D, function(d) {
    matrix(rnorm(n*p, 0, sqrt(SNR/(mean(prior.mean[[d]])*p))), nrow=n, ncol=p)})
  y <- sapply(1:D, function(d) {rnorm(n, x[[d]] %*% beta[[d]], sigma[d])})

  # fit enig model
  fit.enig <- enig(x, y, C, mult.lambda=FALSE, 
                   intercept.eb=TRUE, fixed.eb=c("none"), 
                   full.post=TRUE, init=NULL, 
                   control=control)
  
  # store results
  res[r, ] <- c(fit.enig$eb$alpha, fit.enig$eb$lambda)

}

# save fitted models and results
write.table(set, file="results/simulations_set1.csv")
write.table(res, file="results/simulations_res1.csv")
save(fit.enig, file="results/simulations_fit1.Rdata")

plot(alpha, as.numeric(model.matrix(~ factor(1:4)) %*% fit.enig$eb$alpha))
plot(fit.enig$seq.eb$alpha[, 1], ylim=range(fit.enig$seq.eb$alpha), type="l")
lines(fit.enig$seq.eb$alpha[, 2], col=2)
lines(fit.enig$seq.eb$alpha[, 3], col=3)
lines(fit.enig$seq.eb$alpha[, 4], col=4)

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


################################################################################
###############################   simulation 4   ###############################
################################################################################
## compares the different prior models using MCMC samples in a dense setting ###
################################################################################

# settings
nreps <- 100
n <- 100
p <- 200

# create fixed parameters (eta chosen such that marginal beta variances match)
rho <- 0.5
Sigma <- matrix(rho, ncol=p, nrow=p); diag(Sigma) <- 1
gamma <- 1
sigma <- 1

V <- rep(c(1/2, 1, 2), times=2)
K <- rep(c(4, 10), each=3)
ctalphainv <- V/sigma^2
lambda <- c(3/2, 3, 6, 3/14, 3/7, 6/7)
ratios <- ratio_besselK(lambda/ctalphainv, 1)
ztthetainv <- 1/(ratios + 2*ctalphainv/lambda)/sigma^2
kappa <- (3*ztthetainv^3 + 12*ztthetainv^2*ctalphainv)/
  (K*lambda - 3*ztthetainv^2 - 12*ztthetainv*ctalphainv)
eta <- 6/(K - 3) + 4
xi <- V*(eta - 2)/sigma^2

# check variances
vars <- cbind(nig=sigma^2*ctalphainv,
              dnig=sigma^2*ctalphainv*ztthetainv*(ratios + 2*ctalphainv/lambda),
              studentst=sigma^2*xi/(eta - 2))
all.equal(vars, cbind(nig=V, dnig=V, studentst=V))
kurts <- cbind(nig=3*ctalphainv/lambda + 3,
               dnig=(3*ztthetainv/kappa + 3)*
                 (4*ctalphainv*(ratios + 2*ctalphainv/lambda) + 1)/
                 (lambda*(ratios + 2*ctalphainv/lambda)^2),
               studentst=6/(eta - 4) + 3)
all.equal(kurts, cbind(nig=K, dnig=K, studentst=K))

# store settings in an object
set <- data.frame(nreps=nreps, n=n, p=p, gamma=gamma, sigma=sigma,
                  ctalphainv=ctalphainv, lambda=lambda, ztthetainv=ztthetainv, 
                  kappa=kappa, eta=eta, xi=xi, vars=V, kurts=K)
rownames(set) <- paste0("set", 1:nrow(set))
methods1 <- c("nig", "enig", "dnig", "studentst")
methods2 <- c("nig.advi", "enig.advi", "dnig.advi", "enig.vb")
measures1 <- c("cram", "corbest", "msebest", "cover")
measures2 <- c("cram", "corbest", "msebest")

# set initial values and control parameters
stan.nig <- stan_model("code/nig.stan")
stan.enig <- stan_model("code/enig.stan")
stan.dnig <- stan_model("code/dnig.stan")
stan.studentst <- stan_model("code/studentst.stan")
nsamples <- 1000
nwarmup <- 1000

# fit1 <- matrix(NA, nrow=nsamples, ncol=nrow(set)*(p + 1)*length(methods1), 
#                dimnames=list(paste0("sample", 1:nsamples),
#                              paste(rep(paste0("set", 1:nrow(set)),
#                                        each=(p + 1)*length(methods1)),
#                                    sapply(methods1, function(m) {
#                                      paste0(paste0("beta", 0:p), ".", m)}),
#                                    sep=".")))
# fit2 <- matrix(NA, nrow=nsamples, ncol=nrow(set)*(p + 1)*length(methods2), 
#                dimnames=list(paste0("sample", 1:nsamples),
#                              paste(rep(paste0("set", 1:nrow(set)),
#                                        each=(p + 1)*length(methods2)),
#                                    sapply(methods2, function(m) {
#                                      paste0(paste0("beta", 0:p), ".", m)}),
#                                    sep=".")))
res1 <- matrix(NA, nrow=nreps, ncol=nrow(set)*length(methods1)*
                 length(measures1), 
               dimnames=list(paste0("rep", (1:nreps)), paste(rep(paste0(
                 "set", 1:nrow(set)), each=length(measures1)*length(methods1)), 
                 rep(apply(expand.grid(methods1, measures1), 1, function(s) {
                   paste(rev(s), collapse=".")}), times=nrow(set)), sep=".")))
res2 <- matrix(NA, nrow=nreps, ncol=nrow(set)*length(methods2)*
                 length(measures2),
               dimnames=list(paste0("rep", (1:nreps)), paste(rep(paste0(
                 "set", 1:nrow(set)), each=length(measures2)*length(methods2)), 
                 rep(apply(expand.grid(methods2, measures2), 1, function(s) {
                   paste(rev(s), collapse=".")}), times=nrow(set)), sep=".")))

for(s in 1:nrow(set)) {
  for(k in 1:nreps) {
    print(paste0("set", s, ", rep", k))

    set.seed(2019 + k)
    beta <- rnorm(p, 0, sigma*gamma)
    x <- scale(rmvnorm(n, rep(0, p), sigma=Sigma))
    y <- as.numeric(scale(rnorm(n, x %*% beta, sigma)))

    ctalphainv <- set$ctalphainv[s]
    lambda <- set$lambda[s]
    ztthetainv <- set$ztthetainv[s]
    kappa <- set$kappa[s]
    eta <- set$eta[s]
    xi <- set$xi[s]
    
    fit.nig <- sampling(stan.nig, chains=1, warmup=nwarmup, 
                        iter=nsamples + nwarmup, cores=1, refresh=0, 
                        control=list(adapt_delta=0.8, max_treedepth=12),
                        data=list(n=n, p=p, y=y, x=x, ctalphainv=ctalphainv, 
                                  lambda=lambda))
    fit.enig <- sampling(stan.enig, chains=1, warmup=nwarmup, 
                         iter=nsamples + nwarmup, cores=1, refresh=0, 
                         control=list(adapt_delta=0.8, max_treedepth=12),
                         data=list(n=n, p=p, y=y, x=x, 
                                   ctalphainv=ctalphainv, lambda=lambda))
    fit.dnig <- sampling(stan.dnig, chains=1, warmup=nwarmup, 
                         iter=nsamples + nwarmup, cores=1, refresh=0, 
                         control=list(adapt_delta=0.8, max_treedepth=12),
                         data=list(n=n, p=p, y=y, x=x, ctalphainv=ctalphainv, 
                                   lambda=lambda, ztthetainv=ztthetainv, 
                                   kappa=kappa))
    fit.studentst <- sampling(stan.studentst, chains=1, warmup=nwarmup, 
                              iter=nsamples + nwarmup, cores=1, refresh=0, 
                              control=list(adapt_delta=0.8, max_treedepth=12),
                              data=list(n=n, p=p, y=y, x=x, eta=eta, xi=xi))
    
    # Error in sampler$call_sampler(c(args, dotlist)) : 
    #   stan::variational::normal_meanfield::calc_grad: The number of dropped 
    #   evaluations has reached its maximum amount (10). Your model may be 
    #   either severely ill-conditioned or misspecified.
    # fit.nig.advi <- vb(stan.nig, output_samples=nsamples, 
    #                    data=list(n=n, p=p, y=y, x=x, ctalphainv=ctalphainv, 
    #                              lambda=lambda))
    # fit.enig.advi <- vb(stan.enig, output_samples=nsamples, 
    #                     data=list(n=n, p=p, y=y, x=x, ctalphainv=ctalphainv, 
    #                               lambda=lambda))
    # fit.dnig.advi <- vb(stan.enig, output_samples=nsamples, 
    #                     data=list(n=n, p=p, y=y, x=x, ctalphainv=ctalphainv, 
    #                               lambda=lambda, ztthetainv=ztthetainv, 
    #                               kappa=kappa))
    
    init <- list(aold=rep(1, 2), bold=list(p, p), lambda=lambda, 
                 alpha=1, Calpha=rep(1/ctalphainv, 2))
    fit.enig.vb <- enig(list(x, x), cbind(y, y), NULL, mult.lambda=FALSE, 
                        intercept.eb=TRUE, fixed.eb=c("both"), 
                        full.post=TRUE, init=init,
                        control=list(conv.post=TRUE, trace=TRUE,
                                     epsilon.eb=1e-3, epsilon.vb=1e-3, 
                                     epsilon.opt=sqrt(.Machine$double.eps),
                                     maxit.eb=10, maxit.vb=2, maxit.opt=100,
                                     maxit.post=100))
              
    # fit.nig.map <- optimizing(stan.nig, 
    #                           data=list(n=n, p=p, y=y, x=x, 
    #                                     ctalphainv=ctalphainv, lambda=lambda))
    # fit.enig.map <- optimizing(stan.enig, 
    #                            data=list(n=n, p=p, y=y, x=x, 
    #                                      ctalphainv=ctalphainv, lambda=lambda))
    # fit.dnig.map <- optimizing(stan.enig,
    #                            data=list(n=n, p=p, y=y, x=x, 
    #                                      ctalphainv=ctalphainv, lambda=lambda, 
    #                                      ztthetainv=ztthetainv, kappa=kappa))   
    # 
    # fit.glmnet <- glmnet(x, y, alpha=0, lambda=1/(n*ctalphainv*sigma^2))
    # fit.cv.glmnet <- cv.glmnet(x, y, alpha=0)
    
    best <- data.frame(nig=c(get_posterior_mean(fit.nig, pars="beta0"),
                             get_posterior_mean(fit.nig, pars="beta")),
                       enig=c(get_posterior_mean(fit.enig, pars="beta0"),
                              get_posterior_mean(fit.enig, pars="beta")),
                       dnig=c(get_posterior_mean(fit.dnig, pars="beta0"),
                              get_posterior_mean(fit.dnig, pars="beta")),
                       studentst=c(get_posterior_mean(fit.studentst, 
                                                      pars="beta0"),
                                   get_posterior_mean(fit.studentst, 
                                                      pars="beta")),
                       # nig.advi=c(get_posterior_mean(fit.nig.advi, 
                       #                               pars="beta0"),
                       #            get_posterior_mean(fit.nig.advi, 
                       #                               pars="beta")), 
                       # enig.advi=c(get_posterior_mean(fit.enig.advi, 
                       #                                pars="beta0"),
                       #             get_posterior_mean(fit.enig.advi, 
                       #                                pars="beta")), 
                       # dnig.advi=c(get_posterior_mean(fit.dnig.advi, 
                       #                                pars="beta0"),
                       #             get_posterior_mean(fit.dnig.advi, 
                       #                                pars="beta")), 
                       nig.advi=NA, enig.advi=NA, dnig.advi=NA,
                       enig.vb=c(0, fit.enig.vb$vb$mpost$beta[, 1]))
                       # nig.map=fit.nig.map$par[1:(p + 1)],
                       # enig.map=fit.enig.map$par[1:(p + 1)], 
                       # dnig.map=fit.dnig.map$par[1:(p + 1)], 
                       # glmnet=as.numeric(coef(fit.glmnet)),
                       # cv.glmnet=as.numeric(coef(fit.cv.glmnet, 
                       #                           s="lambda.min")))

    samp.nig <- cbind(as.matrix(fit.nig, "beta0"), as.matrix(fit.nig, "beta"))
    samp.enig <- cbind(as.matrix(fit.enig, "beta0"), as.matrix(fit.enig, "beta"))
    samp.dnig <- cbind(as.matrix(fit.dnig, "beta0"), as.matrix(fit.dnig, "beta"))
    samp.studentst <- cbind(as.matrix(fit.studentst, "beta0"),
                            as.matrix(fit.studentst, "beta"))
    # samp.nig.advi <- matrix(unlist(fit.nig.advi@sim$samples[[1]][1:(p + 1)]),
    #                         ncol=p + 1, byrow=FALSE)
    # samp.dnig.advi <- matrix(unlist(fit.dnig.advi@sim$samples[[1]][1:(p + 1)]),
    #                          ncol=p + 1, byrow=FALSE)
    # samp.enig.advi <- matrix(unlist(fit.enig.advi@sim$samples[[1]][1:(p + 1)]),
    #                          ncol=p + 1, byrow=FALSE)
    
    samp.enig.vb <- cbind(0, rmvnorm(nsamples, 
                                     mean=as.numeric(fit.enig.vb$vb$mu[, 1]),
                                     sigma=fit.enig.vb$vb$Sigma[[1]]))
    cov.true <- sigma^2*solve(t(x) %*% x + gamma^(-2)*diag(p))
    mean.true <- cov.true %*% t(x) %*% y/sigma^2
    samp.true <- cbind(0, rmvnorm(nsamples, mean.true, cov.true))
    samp <- data.frame(nig=samp.nig, enig=samp.enig, dnig=samp.dnig, 
                       studentst=samp.studentst)

    cram1 <- c(nig=cramer.test(samp.nig, samp.true, 
                               just.statistic=TRUE)$statistic,
               enig=cramer.test(samp.enig, samp.true,
                                just.statistic=TRUE)$statistic,
               dnig=cramer.test(samp.dnig, samp.true,
                                just.statistic=TRUE)$statistic,
               studentst=cramer.test(samp.studentst, samp.true,
                                     just.statistic=TRUE)$statistic)
    # cram2 <- c(nig.advi=cramer.test(samp.nig, samp.nig.advi,
    #                                 just.statistic=TRUE)$statistic,
    #            enig.advi=cramer.test(samp.enig, samp.enig.advi,
    #                                  just.statistic=TRUE)$statistic,
    #            dnig.advi=cramer.test(samp.dnig, samp.dnig.advi,
    #                                  just.statistic=TRUE)$statistic,
    cram2 <- c(nig.advi=NA, enig.advi=NA, dnig.advi=NA,
               enig.vb=cramer.test(samp.enig, samp.enig.vb,
                                   just.statistic=TRUE)$statistic)

    corbest1 <- apply(best[, c(1:4)], 2, cor, y=c(0, beta))
    # corbest2 <- c(nig.advi=cor(best[, 5], best[, 1]),
    #               enig.advi=cor(best[, 6], best[, 2]),
    #               dnig.advi=cor(best[, 7], best[, 3]),
    corbest2 <- c(nig.advi=NA, enig.advi=NA, dnig.advi=NA,
                  enig.vb=cor(best[, 8], best[, 2]))
  
    msebest1 <- apply(best[, c(1:4)], 2, function(b) {
      mean((b - c(0, beta))^2)})
    # msebest2 <- c(nig.advi=mean((best[, 5] - best[, 1])^2),
    #               enig.advi=mean((best[, 6] - best[, 2])^2),
    #               dnig.advi=mean((best[, 7] - best[, 3])^2),
    msebest2 <- c(nig.advi=NA, enig.advi=NA, dnig.advi=NA,
                  enig.vb=mean((best[, 8] - best[, 2])^2))

    ci.nig <- summary(fit.nig, pars=c("beta0", "beta"))$summary[, c(4, 8)]
    ci.enig <- summary(fit.enig, pars=c("beta0", "beta"))$summary[, c(4, 8)]
    ci.dnig <- summary(fit.dnig, pars=c("beta0", "beta"))$summary[, c(4, 8)]
    ci.studentst <- summary(fit.studentst, 
                            pars=c("beta0", "beta"))$summary[, c(4, 8)]

    cover1 <- c(nig=mean((ci.nig[, 1] < c(0, beta)) & 
                           (c(0, beta) < ci.nig[, 2])),
                enig=mean((ci.enig[, 1] < c(0, beta)) & 
                            (c(0, beta) < ci.enig[, 2])),
                dnig=mean((ci.dnig[, 1] < c(0, beta)) & 
                            (c(0, beta) < ci.dnig[, 2])),
                studentst=mean((ci.studentst[, 1] < c(0, beta)) & 
                                 (c(0, beta) < ci.studentst[, 2])))
    
    res1[k, substr(colnames(res1), 1, 4)==paste0("set", s)] <- 
      c(cram1, corbest1, msebest1, cover1)
    res2[k, substr(colnames(res2), 1, 4)==paste0("set", s)] <- 
      c(cram2, corbest2, msebest2)
    
    # store the first replications samples
    # if(k==1) {
    #   fit1[, substr(colnames(fit1), 1, 4)==paste0("set", s)] <- 
    #     samp
    # }
  
    write.table(res1, file="results/simulations_igaussian_res4.1.csv")
    write.table(res2, file="results/simulations_igaussian_res4.2.csv")
    
  }
}

# write.table(fit, file="results/simulations_igaussian_fit4.csv")

write.table(set, file="results/simulations_igaussian_set4.csv")


################################################################################
###############################   simulation 5   ###############################
################################################################################
## compares the different prior models using MCMC samples in a sparse setting ##
################################################################################
# settings
nreps <- 100
n <- 100
p <- 200

# create fixed parameters (eta chosen such that marginal beta variances match)
rho <- 0.5
Sigma <- matrix(rho, ncol=p, nrow=p); diag(Sigma) <- 1
gamma <- 1
sigma <- 1
phi <- 1
chi <- c(phi, 9*phi)

V <- rep(c(1/2, 1, 2), times=2)
K <- rep(c(4, 10), each=3)
ctalphainv <- V/sigma^2
lambda <- c(3/2, 3, 6, 3/14, 3/7, 6/7)
ratios <- ratio_besselK(lambda/ctalphainv, 1)
ztthetainv <- 1/(ratios + 2*ctalphainv/lambda)/sigma^2
kappa <- (3*ztthetainv^3 + 12*ztthetainv^2*ctalphainv)/
  (K*lambda - 3*ztthetainv^2 - 12*ztthetainv*ctalphainv)
eta <- 6/(K - 3) + 4
xi <- V*(eta - 2)/sigma^2

# check variances
vars <- cbind(nig=sigma^2*ctalphainv,
              dnig=sigma^2*ctalphainv*ztthetainv*(ratios + 2*ctalphainv/lambda),
              studentst=sigma^2*xi/(eta - 2))
all.equal(vars, cbind(nig=V, dnig=V, studentst=V))
kurts <- cbind(nig=3*ctalphainv/lambda + 3,
               dnig=(3*ztthetainv/kappa + 3)*
                 (4*ctalphainv*(ratios + 2*ctalphainv/lambda) + 1)/
                 (lambda*(ratios + 2*ctalphainv/lambda)^2),
               studentst=6/(eta - 4) + 3)
all.equal(kurts, cbind(nig=K, dnig=K, studentst=K))

# store settings in an object
set <- data.frame(nreps=nreps, n=n, p=p, gamma=gamma, sigma=sigma, phi=phi,
                  chi=rep(chi, each=length(V)),
                  ctalphainv=ctalphainv, lambda=lambda, ztthetainv=ztthetainv, 
                  kappa=kappa, eta=eta, xi=xi, vars=V, kurts=K)
rownames(set) <- paste0("set", 1:nrow(set))
methods <- c("nig", "enig", "dnig", "studentst")
measures <- c("cram", "corbest", "msebest", "cover")

# set initial values and control parameters
stan.nig <- stan_model("code/nig.stan")
stan.enig <- stan_model("code/enig.stan")
stan.dnig <- stan_model("code/dnig.stan")
stan.studentst <- stan_model("code/studentst.stan")
nsamples <- 1000
nwarmup <- 1000

fit <- matrix(NA, nrow=nsamples, ncol=nrow(set)*(p + 1)*length(methods), 
              dimnames=list(paste0("sample", 1:nsamples),
                            paste(rep(paste0("set", 1:nrow(set)),
                                      each=(p + 1)*length(methods)),
                                  c(paste0(paste0("beta", 0:p), ".nig"),
                                    paste0(paste0("beta", 0:p), ".nigig"),
                                    paste0(paste0("beta", 0:p), ".studentst")),
                                  sep=".")))
res <- matrix(NA, nrow=nreps, ncol=nrow(set)*length(methods)*length(measures),
              dimnames=list(paste0("rep", (1:nreps)), paste(rep(paste0(
                "set", 1:nrow(set)), each=length(measures)*length(methods)), 
                rep(apply(expand.grid(methods, measures), 1, function(s) {
                  paste(rev(s), collapse=".")}), times=nrow(set)), sep=".")))

for(s in 1:nrow(set)) {
  for(k in 1:nreps) {
    print(paste0("set", s, ", rep", k))

    set.seed(2019 + k)
    beta <- rspikeandslab(p, set[s, "phi"], set[s, "chi"],
                          set[s, "sigma"], set[s, "gamma"], 0.01)
    x <- rmvnorm(n, rep(0, p), sigma=Sigma)
    y <- rnorm(n, x %*% beta, sigma)

    ctalphainv <- set$ctalphainv[s]
    lambda <- set$lambda[s]
    ztthetainv <- set$ztthetainv[s]
    kappa <- set$kappa[s]
    eta <- set$eta[s]
    xi <- set$xi[s]
    
    fit.nig <- sampling(stan.nig, chains=1, warmup=nwarmup, 
                        iter=nsamples + nwarmup, cores=1, refresh=0, 
                        control=list(adapt_delta=0.8, max_treedepth=12),
                        data=list(n=n, p=p, y=y, x=x, ctalphainv=ctalphainv, 
                                  lambda=lambda))
    fit.enig <- sampling(stan.enig, chains=1, warmup=nwarmup, 
                         iter=nsamples + nwarmup, cores=1, refresh=0, 
                         control=list(adapt_delta=0.8, max_treedepth=12),
                         data=list(n=n, p=p, y=y, x=x, 
                                   ctalphainv=ctalphainv, lambda=lambda))
    fit.dnig <- sampling(stan.dnig, chains=1, warmup=nwarmup, 
                         iter=nsamples + nwarmup, cores=1, refresh=0, 
                         control=list(adapt_delta=0.8, max_treedepth=12),
                         data=list(n=n, p=p, y=y, x=x, ctalphainv=ctalphainv, 
                                   lambda=lambda, ztthetainv=ztthetainv, 
                                   kappa=kappa))
    fit.studentst <- sampling(stan.studentst, chains=1, warmup=nwarmup, 
                              iter=nsamples + nwarmup, cores=1, refresh=0, 
                              control=list(adapt_delta=0.8, max_treedepth=12),
                              data=list(n=n, p=p, y=y, x=x, eta=eta, xi=xi))
    
    best <- data.frame(nig=c(get_posterior_mean(fit.nig, pars="beta0"),
                             get_posterior_mean(fit.nig, pars="beta")),
                       enig=c(get_posterior_mean(fit.enig, pars="beta0"),
                              get_posterior_mean(fit.enig, pars="beta")),
                       dnig=c(get_posterior_mean(fit.dnig, pars="beta0"),
                              get_posterior_mean(fit.dnig, pars="beta")),
                       studentst=c(get_posterior_mean(fit.studentst, 
                                                      pars="beta0"),
                                   get_posterior_mean(fit.studentst, 
                                                      pars="beta")))
    
    samp.nig <- cbind(as.matrix(fit.nig, "beta0"),as.matrix(fit.nig, "beta"))
    samp.enig <- cbind(as.matrix(fit.enig, "beta0"),as.matrix(fit.enig, "beta"))
    samp.dnig <- cbind(as.matrix(fit.dnig, "beta0"),as.matrix(fit.dnig, "beta"))
    samp.studentst <- cbind(as.matrix(fit.studentst, "beta0"),
                            as.matrix(fit.studentst, "beta"))
    samp.true <- cbind(0, t(
      gibbs.spikeandslab(x, y, set[s, "gamma"], set[s, "sigma"], 0.01,
                         set[s, "phi"], set[s, "chi"], init=NULL,
                         control=list(samples=nsamples, warmup=nwarmup))$beta))
    samp <- data.frame(nig=samp.nig, enig=samp.enig, dnig=samp.dnig, 
                       studentst=samp.studentst)
    
    cram <- c(nig=cramer.test(samp.nig, samp.true, 
                              just.statistic=TRUE)$statistic,
              enig=cramer.test(samp.enig, samp.true,
                               just.statistic=TRUE)$statistic,
              dnig=cramer.test(samp.dnig, samp.true,
                               just.statistic=TRUE)$statistic,
              studentst=cramer.test(samp.studentst, samp.true,
                                    just.statistic=TRUE)$statistic)
    
    corbest <- c(nig=cor(c(0, beta), best[, 1]),
                 enig=cor(c(0, beta), best[, 2]),
                 dnig=cor(c(0, beta), best[, 3]),
                 studentst=cor(c(0, beta), best[, 4]))
    
    msebest <- c(nig=mean((c(0, beta) - best[, 1])^2),
                 enig=mean((c(0, beta) - best[, 2])^2),
                 dnig=mean((c(0, beta) - best[, 3])^2),
                 studentst=mean((c(0, beta) - best[, 4])^2))
    
    ci.nig <- summary(fit.nig, pars=c("beta0", "beta"))$summary[, c(4, 8)]
    ci.enig <- summary(fit.enig, pars=c("beta0", "beta"))$summary[, c(4, 8)]
    ci.dnig <- summary(fit.dnig, pars=c("beta0", "beta"))$summary[, c(4, 8)]
    ci.studentst <- summary(fit.studentst, 
                            pars=c("beta0", "beta"))$summary[, c(4, 8)]
    
    cover <- c(nig=mean(sapply(1:(p + 1), function(j) {
      (ci.nig[j, 1] < c(0, beta)[j]) & (c(0, beta)[j] < ci.nig[j, 2])})),
      enig=mean(sapply(1:(p + 1), function(j) {
        (ci.enig[j, 1] < c(0, beta)[j]) & (c(0, beta)[j] < ci.enig[j, 2])})),
      dnig=mean(sapply(1:(p + 1), function(j) {
        (ci.dnig[j, 1] < c(0, beta)[j]) & (c(0, beta)[j] < ci.dnig[j, 2])})),
      studentst=mean(sapply(1:(p + 1), function(j) {
        (ci.studentst[j, 1] < c(0, beta)[j]) &
          (c(0, beta)[j] < ci.studentst[j, 2])})))
    
    res[k, substr(colnames(res), 1, 4)==paste0("set", s)] <- 
      c(cram, corbest, msebest, cover)
    # store the first replications samples
    if(k==1) {
      fit[, substr(colnames(fit), 1, 4)==paste0("set", s)] <- 
        samp
    }
  }
}

write.table(fit, file="results/simulations_igaussian_fit5.csv")
write.table(res, file="results/simulations_igaussian_res5.csv")
write.table(set, file="results/simulations_igaussian_set5.csv")







###### ENIG model
set.seed(2019)
n <- 200
D <- 100
p <- sample(100:200, D, replace=TRUE)
x <- lapply(1:D, function(d) {matrix(rnorm(p[d]*n), nrow=n, ncol=p[d])})
C <- lapply(1:D, function(d) {cbind(1, rbinom(p[d], 1, 0.5))})
gamma <- sqrt(c(0.5, 2))
beta <- lapply(1:D, function(d) {
  rnorm(p[d], 0, as.numeric(C[[1]] %*% c(gamma[1], diff(gamma))))})
y <- sapply(1:D, function(d) {as.numeric(x[[d]] %*% beta[[d]]) + rnorm(n)})
x, y, C, mult.lambda=FALSE, intercept.eb=TRUE, 
fixed.eb="none", init=NULL
control=list(conv.post=TRUE, trace=TRUE,
             epsilon.eb=1e-3, epsilon.vb=1e-3, 
             epsilon.opt=sqrt(.Machine$double.eps),
             maxit.eb=30, maxit.vb=2, maxit.opt=100,
             maxit.post=100)


fit.test <- enig(x, y, C, mult.lambda=mult.lambda, fixed.eb=fixed.eb, init=init,
                 control=control)



opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.1, 1, 1))
layout(matrix(rep(rep(1:2, each=2), 2), nrow=2, ncol=4, byrow=TRUE))
plot(0, type="n", main="(a)", ylab=expression(hat(E)(gamma[d]^2)), 
     xlab="Iteration", xlim=c(1, fit.test$iter$eb), 
     ylim=range(fit.test$seq.eb$mprior))
for(d in 1:D) {
  for(j in 1:p[d]) {
    lines(fit.test$seq.eb$mprior[[d]][, j], col=j)
  }
}

fit.test$seq.eb$vprior

plot(0, type="n", main="(b)", ylab=expression(hat(V)(gamma[d]^2)), 
     xlab="Iteration", xlim=c(1, fit.test$iter$eb), 
     ylim=range(sapply(fit.test$seq.eb$vprior, function(i) {i[-c(1:5), ]})))
for(d in 1:D) {
  for(j in 1:p[d]) {
    lines(fit.test$seq.eb$vprior[[d]][-c(1:5), j], col=j)
  }
}
par(opar)

unique(unlist(fit.test$eb$mprior))
gamma^2
unique(unlist(fit.test$eb$vprior))


