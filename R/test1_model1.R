library(R2OpenBUGS)

path.mod <- "/Users/magnusmunch/Documents/OneDrive/PhD/project_cambridge/code/"        
path.res <- "/Users/magnusmunch/Documents/OneDrive/PhD/project_cambridge/results/"
path.pgm <- "/Users/magnusmunch/.wine/drive_c/Program Files (x86)/OpenBUGS/OpenBUGS323/OpenBUGS.exe"
path.graph <- "/Users/magnusmunch/Documents/OneDrive/PhD/project_cambridge/graphs/"

### create model
sink(paste(path.mod, "test1_model1.txt", sep=""))        
cat("model {
  for(d in 1:Ndrugs) {
    for (i in 1:Nobs) {
      mu[i, d] <- beta0[d] + inprod(x[i, ], beta[, d]) + inprod(z[i, ], u[, d])
      y[i, d] ~ dnorm(mu[i, d], inv.sigma.squared[d])
    }
  
    temp1[d] <- inv.gamma.squared[d]*inv.sigma.squared[d]
    for(j in 1:Nomics) {
      beta[j, d] ~ dnorm(0, temp1[d])
    }

    temp2[d] <- inv.tau.squared[d]*inv.sigma.squared[d]
    for(t in 1:Ntissue) {
      u[t, d] ~ dnorm(0, temp2[d])
    }

    inv.sigma.squared[d] ~ dgamma(0.001, 0.001)
    inv.tau.squared[d] ~ dgamma(0.001, 0.001)
    inv.gamma.squared[d] ~ dgamma(0.001, 0.001)
    beta0[d] ~ dnorm(0, 0.001)

    sigma[d] <- 1/sqrt(inv.sigma.squared[d])
    tau[d] <- 1/sqrt(inv.tau.squared[d])
    gamma[d] <- 1/sqrt(inv.gamma.squared[d])
  }
}", fill=TRUE)
sink()

### simulate data
Nomics <- 10
Nobs <- 100
Ndrugs <- 5
Ntissue <- 5

sigma <- rep(1, Ndrugs)
tau <- exp(seq(-2, 2, length.out=Ndrugs))
gamma <- exp(seq(-2, 2, length.out=Ndrugs))

beta0 <- rep(0, Ndrugs)

x <- matrix(rnorm(Nobs*Nomics, 0, 1), nrow=Nobs, ncol=Nomics)
z <- matrix(rep(c(rep(1, Nobs/Ntissue), rep(0, Nobs)), Ntissue), nrow=Nobs, ncol=Ntissue)
y <- matrix(0, nrow=Nobs, ncol=Ndrugs)
beta <- matrix(0, nrow=Nomics, ncol=Ndrugs)
u <- matrix(0, nrow=Ntissue, ncol=Ndrugs)
for(d in 1:Ndrugs) {
  beta[, d] <- rnorm(Nomics, 0, sigma[d]*gamma[d])
  u[, d] <- rnorm(Ntissue, 0, sigma[d]*tau[d])
  y[, d] <- rnorm(as.numeric(x %*% beta[, d] + z %*% u[, d]), sigma[d])
}


### MCMC sampling
data <- list("Nobs", "Ntissue", "Nomics", "Ndrugs", "x", "z", "y")

inits <- function() {
  list(beta0=rnorm(Ndrugs, 0, 1), 
       beta=matrix(rnorm(Ndrugs*Nomics, 0, 1), nrow=Nomics, ncol=Ndrugs), 
       u=matrix(rnorm(Ndrugs*Ntissue, 0, 1), nrow=Ntissue, ncol=Ndrugs), 
       inv.sigma.squared=rgamma(Ndrugs, 0.1, 0.1),
       inv.tau.squared=rgamma(Ndrugs, 0.1, 0.1),
       inv.gamma.squared=rgamma(Ndrugs, 0.1, 0.1))
}

model1.sim1 <- bugs(data, inits, useWINE=TRUE,
                    model.file=paste(path.mod, "test1_model1.txt", sep=""),
                    parameters=c("beta0", "beta", "u", "gamma", "tau", "sigma"),
                    n.chains=3, n.iter=1000, n.burnin=500,
                    working.directory=path.res, OpenBUGS.pgm=path.pgm)

colnames(model1.sim1$sims.matrix)
### diagnostics, plots, and tables
# trace plots
a <- 91
par(mfrow=c(2, 1))
plot(model1.sim1$sims.array[, 1, a], type="l", xlab="iteration",
     ylab="posterior mean",
     main=dimnames(model1.sim1$sims.array)[[3]][a],
     ylim=range(model1.sim1$sims.array[, , a]))
lines(model1.sim1$sims.array[, 2, a], col=2)
lines(model1.sim1$sims.array[, 3, a], col=3)
acf(model1.sim1$sims.array[, 1, a], main="")
par(mfrow=c(1, 1))

# correlations between parameters
cor.mat <- cor(model1.sim1$sims.matrix)
heatmap(cor.mat, Rowv=NA, Colv=NA, revC=TRUE, labRow=NA, labCol=NA)

# gelman-rubin statistic (close to one is good, higher is worse)
gelman.diag(model1.sim1)

# posterior means against true
png(paste(path.graph, "test1_model1_posterior_means.png", sep=""),
    units="in", width=7, height=7, res=120)
par(mfrow=c(2, 3), mar=c(4, 4, 4, 2) + 0.1)
plot(beta, model1.sim1$summary[c(6:55), 1], xlab=expression(beta),
     ylab=expression(hat(beta)))
plot(beta0, model1.sim1$summary[c(1:5), 1], xlab=expression(beta[0]),
     ylab=expression(hat(beta[0])))
plot(u, model1.sim1$summary[c(56:80), 1], xlab=expression(u),
     ylab=expression(hat(u)))
par(mar=c(5, 4, 3, 2) + 0.1)
plot(gamma, model1.sim1$summary[c(81:85), 1], xlab=expression(gamma),
     ylab=expression(hat(gamma)))
plot(tau, model1.sim1$summary[c(86:90), 1], xlab=expression(tau),
     ylab=expression(hat(tau)))
plot(sigma, model1.sim1$summary[c(91:95), 1], xlab=expression(sigma),
     ylab=expression(hat(sigma)))
par(mfrow=c(1, 1))
dev.off()



