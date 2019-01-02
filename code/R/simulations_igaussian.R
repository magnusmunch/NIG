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
library(devtools)
install_github("magnusmunch/cambridge/code", local=FALSE, 
               auth_token=Sys.getenv("GITHUB_PAT"))

### libraries
library(cambridge)

### simulations
# settings
nreps <- 100
n <- 100
p <- 200
D <- 100
nclass <- 5

# create fixed parameters
alpha <- c(1:nclass)
C <- model.matrix(~ 1 + factor(rep(c(1:nclass), each=D/nclass)))
theta <- as.numeric(1/(C %*% alpha))
gamma <- sqrt(theta)
sigma <- rep(1, D)
SNR <- 10

# store settings in an object
temp <- c(nreps=nreps, n=n, p=p, D=D, nclass=nclass, alpha=alpha, C=C, 
          theta=theta, gamma=gamma, sigma=sigma, SNR=SNR)
set1 <- t(matrix(temp)); colnames(set1) <- names(temp)

# objects to store simulation results
list1 <- sapply(c("inv Gauss", "ind inv Gauss"), function(m) {
  matrix(NA, ncol=D, nrow=nreps, 
         dimnames=list(NULL, paste("drug", c(1:D), sep="")))}, simplify=FALSE)
list2 <- sapply(c("inv Gamma"), function(m) {
  matrix(NA, ncol=D, nrow=nreps, 
         dimnames=list(NULL, paste("drug", c(1:D), sep="")))}, simplify=FALSE)
list3 <- sapply(c("inv Gauss", "ind inv Gauss"), function(m) {
  matrix(NA, ncol=nclass, nrow=nreps, 
         dimnames=list(NULL, paste("class", c(1:nclass), sep="")))}, 
  simplify=FALSE)
list4 <- sapply(c("inv Gamma"), function(m) {
  matrix(NA, ncol=nclass, nrow=nreps, 
         dimnames=list(NULL, paste("class", c(1:nclass), sep="")))}, 
  simplify=FALSE)
list5 <- sapply(c("inv Gauss", "ind inv Gauss"), function(m) {
  rep(NA, nreps)}, simplify=FALSE)
list6 <- sapply(c("inv Gauss", "ind inv Gauss", "inv Gamma"), function(m) {
  rep(NA, nreps)}, simplify=FALSE)
res1 <- data.frame(theta=list1, delta=list1, zeta=list1, apost=list2, 
                   bpost=list2, cpost=list2, dpost=list2, alpha=list3, 
                   aprior=list4, bprior=list4, lambda=list5, mse.mu=list6, 
                   cor.mu=list6, mSNR=numeric(nreps), 
                   row.names=paste("rep", c(1:nreps), sep=""))

# set initial values and control parameters
control <- list(epsilon.eb=1e-3, epsilon.vb=1e-3, maxit.eb=20, maxit.vb=2, 
                trace=FALSE)
init <- list(alpha=c(mean(gamma^2), rep(0, nclass - 1)), 
             lambda=mean(gamma^2)^2, a=rep(mean(1/gamma^2), D), 
             zeta=rep(0.5*mean(sigma^2)*(n + p - 1), D))

# simulation (SNR not correct)
set.seed(2018)
for(r in 1:nreps) {
  
  cat("\r", "replication", r)
  
  # create data
  beta <- sapply(1:D, function(d) {rnorm(p, 0, sigma[d]*gamma[d])})
  x <- matrix(rnorm(n*p, 0, sqrt(SNR/(mean(gamma^2)*p))), nrow=n, ncol=p)
  y <- sapply(1:D, function(d) {rnorm(n, x %*% beta[, d], sigma[d])})
  
  # fit inverse Gaussian models
  fit1.igauss1 <- est.igauss(x, y, C, control=control, init=init,
                             test=list(ind=FALSE))
  fit1.igauss2 <- est.igauss(x, y, C, control=control, init=init,
                             test=list(ind=TRUE))
  
  # fit independent inverse Gamma priors model (vague, but one-point prior)
  fit1.gwen <- est.gwen(x, y, eqid=rep(c(1:nclass), each=D/nclass),
                        control=list(epsilon.eb=1e-3, epsilon.vb=1e-3,
                                     epsilon.opt=1e-6, maxit.eb=20, maxit.vb=2,
                                     maxit.opt=20, conv.vb="elbo",
                                     trace=FALSE),
                        init=list(aprior=0.001/mean(gamma^2), bprior=0.001))
  
  # store results (not correct)
  res1[r, c(1:D)] <- fit1.igauss1$seq.eb$theta[fit1.igauss1$iter$eb, ]
  res1[r, c(101:(2*D))] <- fit1.igauss2$seq.eb$theta[fit1.igauss2$iter$eb, ]
  res1[r, c(201:(3*D))]  <- fit1.igauss1$vb.post$delta
  res1[r, c(301:(4*D))]  <- fit1.igauss2$vb.post$delta
  res1[r, c(401:(5*D))] <- fit1.igauss1$vb.post$zeta
  res1[r, c(501:(6*D))] <- fit1.igauss2$vb.post$zeta
  res1[r, c(601:(7*D))] <- fit1.gwen$vb.post$apost
  res1[r, c(701:(8*D))] <- fit1.gwen$vb.post$bpost
  res1[r, c(801:(9*D))] <- fit1.gwen$vb.post$cpost
  res1[r, c(901:(10*D))] <- fit1.gwen$vb.post$dpost
  res1[r, c(1001:(10*D + 1*nclass))] <- 
    fit1.igauss1$seq.eb$alpha[fit1.igauss1$iter$eb, ]
  res1[r, c((10*D + 1*nclass + 1):(10*D + 2*nclass))] <- 
    fit1.igauss2$seq.eb$alpha[fit1.igauss2$iter$eb, ]
  res1[r, c((10*D + 2*nclass + 1):(10*D + 3*nclass))] <- 
    fit1.gwen$seq.eb$aprior[fit1.gwen$iter$eb, ]
  res1[r, c((10*D + 3*nclass + 1):(10*D + 4*nclass))] <- 
    fit1.gwen$seq.eb$bprior[fit1.gwen$iter$eb, ]
  res1[r, 10*D + 4*nclass + 1] <- 
    fit1.igauss1$seq.eb$lambda[fit1.igauss1$iter$eb]
  res1[r, 10*D + 4*nclass + 2] <- 
    fit1.igauss2$seq.eb$lambda[fit1.igauss2$iter$eb]
  res1[r, 10*D + 4*nclass + 3] <- mean((beta - fit1.igauss1$vb.post$mu)^2)
  res1[r, 10*D + 4*nclass + 4] <- mean((beta - fit1.igauss2$vb.post$mu)^2)
  res1[r, 10*D + 4*nclass + 5] <- mean((beta - fit1.gwen$vb.post$mu)^2)
  res1[r, 10*D + 4*nclass + 6] <- 
    cor(as.numeric(beta), as.numeric(fit1.igauss1$vb.post$mu))
  res1[r, 10*D + 4*nclass + 7] <- 
    cor(as.numeric(beta), as.numeric(fit1.igauss2$vb.post$mu))
  res1[r, 10*D + 4*nclass + 8] <- 
    cor(as.numeric(beta), as.numeric(fit1.gwen$vb.post$mu))
  
  # store actual mean SNR
  res1[r, 10*D + 4*nclass + 9] <- 
    mean(apply(x %*% beta, 2, var)/apply(y, 2, var))
  
  # save the settings, results, and last model fits for convergence examples
  temp1 <- data.frame(inv.Gauss=fit1.igauss1$seq.eb)
  temp1 <- temp1[rep(seq_len(nrow(temp1)), times=fit1.igauss1$iter$vb + 1), ]
  temp2 <- data.frame(ind.inv.Gauss=fit1.igauss2$seq.eb)
  temp2 <- temp2[rep(seq_len(nrow(temp2)), times=fit1.igauss2$iter$vb + 1), ]
  fit1 <- data.frame(inv.Gauss.elbo=fit1.igauss1$seq.elbo, temp1, 
                     ind.inv.Gauss.elbo=fit1.igauss2$seq.elbo, temp2, 
                     inv.Gamma.elbo=fit1.gwen$seq.elbo,
                     inv.Gamma=fit1.gwen$seq.eb)
 
  write.table(set1, file="../../results/simulations_igaussian_set1.csv")
  write.table(res1, file="../../results/simulations_igaussian_res1.csv")
  write.table(fit1, file="../../results/simulations_igaussian_fit1.csv")
  
}

warnings()

