#!/usr/bin/env Rscript

### installation of packages
# if(!("devtools" %in% installed.packages())) {
#   install.packages("devtools")
# }
# library(devtools)
# install_github("magnusmunch/cambridge/rpackage", local=FALSE,
#                auth_token=Sys.getenv("GITHUB_PAT"))
if(!("statmod" %in% installed.packages())) {
  install.packages("statmod")
}

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
p <- 30
D <- 1

# create fixed parameters
gamma <- rep(1, D)
sigma <- rep(1, D)
lambda <- 0.1
theta <- 2

# store settings in an object
temp1 <- c(nreps=nreps, n=n, p=p, D=D, gamma=gamma, sigma=sigma, lambda=lambda,
           theta=theta)
set1 <- t(matrix(temp1))
colnames(set1) <- names(temp1)

# objects to store simulation results
methods <- c("MCMC", "ADVI", "VB", "MAP", "glmnet")
params <- c("beta0", "beta", "sigma", "gamma")

# set initial values and control parameters
stan.model <- stan_model("code/igauss.stan")



# simulation
set.seed(2019)
beta <- sapply(1:D, function(d) {rnorm(p, 0, sigma[d]*gamma[d])})
x <- matrix(rnorm(n*p, 0, 1), nrow=n, ncol=p)
y <- sapply(1:D, function(d) {rnorm(n, x %*% beta[, d], sigma[d])})




fit.mcmc <- sapply(1:D, function(d) {
  sampling(stan.model, data=list(n=n, p=p, y=y[, d], x=x, theta=theta, 
                                 lambda=lambda), chains=1, warmup=1000,
       iter=2000, cores=1, refresh=0, control=list(adapt_delta=0.8))})
fit.advi <- sapply(1:D, function(d) {
  vb(stan.model, data=list(n=n, p=p, y=y[, d], x=x, theta=theta, 
                                          lambda=lambda))})
fit.pen <- sapply(1:D, function(d) {
  optimizing(stan.model, data=list(n=n, p=p, y=y[, d], x=x, theta=theta, 
                                   lambda=lambda))}, simplify=FALSE)
fit.glmnet <- sapply(1:D, function(d) {glmnet(x, y[, d], alpha=0,
                                                 lambda=1/(n*theta))}, 
                     simplify=FALSE)
fit.vb <- est.vb(x, y, theta, lambda, intercept=TRUE, init=NULL, posterior=TRUE,
                 control=list(epsilon.vb=1e-6, maxit.vb=100, trace=TRUE))

best.mcmc <- sapply(fit.mcmc, function(s) {
  c(get_posterior_mean(s, pars="beta0"),
    get_posterior_mean(s, pars="beta"))})
best.advi <- sapply(fit.advi, function(s) {
  c(get_posterior_mean(s, pars="beta0"),
    get_posterior_mean(s, pars="beta"))})
best.map <- sapply(fit.pen, function(s) {s$par[1:(p + 1)]})
best.glmnet <- sapply(fit.glmnet, function(s) {as.numeric(coef(s))})
best.vb <- fit.vb$vb.post$mu

bvar.mcmc <- sapply(fit.mcmc, function(s) {
  summary(s)$summary[1:(p + 1), 3]^2})
bvar.vb <- sapply(1:D, function(d) {
  diag(fit.vb$vb.post$Sigma[[d]])})
bvar.advi <- sapply(fit.advi, function(s) {
  summary(s)$summary[1:(p + 1), 2]^2})

samp.advi <- matrix(unlist(fit.advi[[1]]@sim$samples[[1]][1:(p + 1)]), 
                    ncol=p + 1, byrow=FALSE)
samp.vb <- rmvnorm(1000, mean=as.numeric(fit.vb$vb.post$mu), 
                   sigma=fit.vb$vb.post$Sigma[[1]])
samp.mcmc <- cbind(as.matrix(fit.mcmc[[1]], "beta0"), 
                   as.matrix(fit.mcmc[[1]], "beta"))

cram.advi <- cramer.test(samp.advi, samp.mcmc, just.statistic=TRUE)$statistic
cram.vb <- cramer.test(samp.vb, samp.mcmc, just.statistic=TRUE)$statistic

corbvar.advi <- cor(bvar.advi, best.mcmc)
corbvar.vb <- cor(bvar.vb, best.mcmc)

msebvar.advi <- mean((bvar.advi - best.mcmc)^2)
msebvar.vb <- mean((bvar.vb - best.mcmc)^2)

corbest.advi <- cor(best.advi, best.mcmc)
corbest.vb <- cor(best.vb, best.mcmc)
corbest.glmnet <- cor(best.glmnet, best.mcmc)
corbest.map <- cor(best.map, best.mcmc)

msebest.advi <- mean((best.advi - best.mcmc)^2)
msebest.vb <- mean((best.vb - best.mcmc)^2)
msebest.glmnet <- mean((best.glmnet - best.mcmc)^2)
msebest.map <- mean((best.map - best.mcmc)^2)

# ---- hist_igaussian_res3_beta_post ---- 
methods <- c("MCMC", "ADVI", "VB")
labels <- c(methods, expression(hat(E)(beta[0]~"|"~bold(y))))
col <- 2:(length(methods) + 1)
which.beta <- c(1, sample(1:(p + 1), 5))

opar <- par(no.readonly=TRUE)
layout(matrix(c(rep(rep(c(1:3), each=2), 2), rep(rep(c(4:6), each=2), 2)), 
              byrow=TRUE, nrow=4, ncol=6), widths=1, heights=1, respect=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
for(j in 1:length(which.beta)) {
  h1 <- hist(samp.mcmc[, which.beta[j]], plot=FALSE, breaks=30)
  h2 <- hist(samp.advi[, which.beta[j]], plot=FALSE, breaks=30)
  ylim <- c(0, max(dnorm(best.vb[which.beta[j], ], 
                         mean=best.vb[which.beta[j], ], 
                         sd=sqrt(bvar.vb[which.beta[j], ])), 
                   h1$density, h2$density))
  plot(h1, freq=FALSE, 
       main=substitute(beta[i], list(i=which.beta[j] - 1)), 
       xlab=substitute(beta[i], list(i=which.beta[j] - 1)),
       ylim=ylim, col=rgb(1, 0, 0, 0.5))
  plot(h2, freq=FALSE, add=TRUE, col=rgb(0, 1, 0, 0.5))
  curve(dnorm(x, mean=best.vb[which.beta[j], ], 
              sd=sqrt(bvar.vb[which.beta[j], ])), add=TRUE, col=col[3])
  abline(v=best.mcmc[which.beta[j], ], col=col[1], lty=2)
  abline(v=best.advi[which.beta[j], ], col=col[2], lty=2)
  abline(v=best.vb[which.beta[j], ], col=col[3], lty=2)
}
legend("topright", legend=labels, 
       lty=c(rep(NA, length(methods)), 2), col=c(col, 1), 
       seg.len=1, border=c(rep(1, length(methods)), NA), 
       fill=c(col, NA), merge=TRUE)
par(opar)



opar <- par(no.readonly=TRUE)
layout(matrix(c(rep(rep(c(1:3), each=2), 2), rep(c(0, 4, 4, 5, 5, 0), 2)), 
              byrow=TRUE, nrow=4, ncol=6), widths=1, heights=1, respect=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
plot(rbind(0, beta), best.glmnet, ylab=expression(hat(beta)), 
     xlab=expression(beta),
     main=paste0("glmnet, MSE=", 
                 round(mean((best.glmnet - rbind(0, beta))^2), 4)))
plot(rbind(0, beta), best.mcmc, ylab=expression(hat(beta)), 
     xlab=expression(beta),
     main=paste0("MCMC, MSE=", 
                 round(mean((best.mcmc - rbind(0, beta))^2), 4)))
plot(rbind(0, beta), best.vb, ylab=expression(hat(beta)), 
     xlab=expression(beta),
     main=paste0("VB, MSE=", 
                 round(mean((best.vb - rbind(0, beta))^2), 4)))
plot(rbind(0, beta), best.advi, ylab=expression(hat(beta)), 
     xlab=expression(beta),
     main=paste0("ADVI, MSE=", 
                 round(mean((best.advi - rbind(0, beta))^2), 4)))
plot(rbind(0, beta), best.map, ylab=expression(hat(beta)), 
     xlab=expression(beta),
     main=paste0("MAP, MSE=", 
                 round(mean((best.map - rbind(0, beta))^2), 4)))
par(opar)

opar <- par(no.readonly=TRUE)
layout(matrix(c(rep(rep(c(1:2), each=2), 2), rep(rep(c(3:4), each=2), 2)),
              byrow=TRUE, nrow=4, ncol=4), widths=1, heights=1, respect=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
plot(best.mcmc, best.vb, ylab=expression(hat(E)(beta["VB"]~"|"~bold(y))),
     xlab=expression(hat(E)(beta["MCMC"]~"|"~bold(y))))
plot(best.mcmc, best.advi, ylab=expression(hat(beta)["ADVI"]),
     xlab=expression(hat(beta)["MCMC"]))
plot(bvar.mcmc, bvar.vb, ylab=expression(hat(V)(beta["ADVI"]~"|"~bold(y))),
     xlab=expression(hat(V)(beta["MCMC"]~"|"~bold(y))))
plot(bvar.mcmc, bvar.advi, ylab=expression(hat(V)(beta["ADVI"]~"|"~bold(y))),
     xlab=expression(hat(V)(beta["MCMC"]~"|"~bold(y))))
par(opar)






