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
                maxit.eb=20, maxit.vb=2, maxit.opt=100, trace=FALSE)
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
