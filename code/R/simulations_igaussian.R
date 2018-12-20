### change these to user-specific values
path.res <- "/Users/magnusmunch/Documents/OneDrive/PhD/cambridge/results/"
auth_token <- "da10f2b37c513e3383c5b2e0aa1300288329c636"

### installation of package
if(!("cambridge" %in% installed.packages())) {
  if(!("devtools" %in% installed.packages())) {
    install.packages("devtools")
  }
  library(devtools)
  install_github("magnusmunch/cambridge/code", local=FALSE, 
                 auth_token=auth_token)
}

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
set1 <- list(nreps=nreps, n=n, p=p, D=D, nclass=nclass, alpha=alpha, C=C, 
             theta=theta, gamma=gamma, sigma=sigma, SNR=SNR)

# objects to store simulation results
seq.mSNR <- numeric(nreps)
seq.theta <- replicate(3, matrix(NA, ncol=D, nrow=nreps), simplify=FALSE)
seq.alpha <- replicate(3, matrix(NA, ncol=nclass, nrow=nreps), simplify=FALSE)
seq.lambda <- replicate(3, rep(NA, nreps), simplify=FALSE)
seq.delta <- replicate(3, matrix(NA, ncol=D, nrow=nreps), simplify=FALSE)
seq.zeta <- replicate(3, matrix(NA, ncol=D, nrow=nreps), simplify=FALSE)
seq.aprior <- replicate(3, matrix(NA, ncol=nclass, nrow=nreps), simplify=FALSE)
seq.bprior <- replicate(3, matrix(NA, ncol=nclass, nrow=nreps), simplify=FALSE)
seq.apost <- replicate(3, matrix(NA, ncol=D, nrow=nreps), simplify=FALSE)
seq.bpost <- replicate(3, matrix(NA, ncol=D, nrow=nreps), simplify=FALSE)
seq.cpost <- replicate(3, matrix(NA, ncol=D, nrow=nreps), simplify=FALSE)
seq.dpost <- replicate(3, matrix(NA, ncol=D, nrow=nreps), simplify=FALSE)
mse.mu <- replicate(3, rep(NA, nreps), simplify=FALSE)
cor.mu <- replicate(3, rep(NA, nreps), simplify=FALSE)
names(seq.theta) <- names(seq.alpha) <- names(seq.lambda) <- names(seq.delta) <-
  names(seq.zeta) <- names(mse.mu) <- names(cor.mu) <- names(seq.aprior)  <- 
  names(seq.bprior) <- names(seq.apost) <- names(seq.bpost) <- 
  names(seq.cpost) <- names(seq.dpost) <- 
  c("inv Gauss", "ind inv Gauss", "inv Gamma")

# simulation (SNR not correct)
set.seed(2018)
for(r in 1:nreps) {
  
  cat("\r", "replication", r)
  
  # create data
  beta <- sapply(1:D, function(d) {rnorm(p, 0, sigma[d]*gamma[d])})
  x <- matrix(rnorm(n*p, 0, sqrt(SNR/(mean(gamma^2)*p))), nrow=n, ncol=p)
  y <- sapply(1:D, function(d) {rnorm(n, x %*% beta[, d], sigma[d])})

  # store actual mean SNR
  seq.mSNR[r] <- mean(apply(x %*% beta, 2, var)/apply(y, 2, var))

  # set initial values and control parameters
  control <- list(epsilon.eb=1e-3, epsilon.vb=1e-3, maxit.eb=20, maxit.vb=2, 
                  trace=FALSE)
  init <- list(alpha=c(mean(gamma^2), rep(0, nclass - 1)), 
               lambda=mean(gamma^2)^2, a=rep(mean(1/gamma^2), D), 
               zeta=rep(0.5*mean(sigma^2)*(n + p - 1), D))
  
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
  
  # store results
  seq.theta[[1]][r, ] <- fit1.igauss1$seq.eb$theta[fit1.igauss1$iter$eb, ]
  seq.theta[[2]][r, ] <- fit1.igauss2$seq.eb$theta[fit1.igauss2$iter$eb, ]
  seq.alpha[[1]][r, ] <- fit1.igauss1$seq.eb$alpha[fit1.igauss1$iter$eb, ]
  seq.alpha[[2]][r, ] <- fit1.igauss2$seq.eb$alpha[fit1.igauss2$iter$eb, ]
  seq.lambda[[1]][r] <- fit1.igauss1$seq.eb$lambda[fit1.igauss1$iter$eb]
  seq.lambda[[2]][r] <- fit1.igauss2$seq.eb$lambda[fit1.igauss2$iter$eb]
  seq.delta[[1]][r, ] <- fit1.igauss1$vb.post$delta
  seq.delta[[2]][r, ] <- fit1.igauss2$vb.post$delta
  seq.zeta[[1]][r, ] <- fit1.igauss1$vb.post$zeta
  seq.zeta[[2]][r, ] <- fit1.igauss2$vb.post$zeta
  mse.mu[[1]][r] <- mean((beta - fit1.igauss1$vb.post$mu)^2)
  mse.mu[[2]][r] <- mean((beta - fit1.igauss2$vb.post$mu)^2)
  mse.mu[[3]][r] <- mean((beta - fit1.gwen$vb.post$mu)^2)
  cor.mu[[1]][r] <- cor(as.numeric(beta), as.numeric(fit1.igauss1$vb.post$mu))
  cor.mu[[2]][r] <- cor(as.numeric(beta), as.numeric(fit1.igauss2$vb.post$mu))
  cor.mu[[3]][r] <- cor(as.numeric(beta), as.numeric(fit1.gwen$vb.post$mu))
  seq.aprior[[3]][r, ] <- fit1.gwen$seq.eb$aprior[fit1.gwen$iter$eb, ]
  seq.bprior[[3]][r, ] <- fit1.gwen$seq.eb$bprior[fit1.gwen$iter$eb, ]
  seq.apost[[3]][r, ] <- fit1.gwen$vb.post$apost
  seq.bpost[[3]][r, ] <- fit1.gwen$vb.post$bpost
  seq.cpost[[3]][r, ] <- fit1.gwen$vb.post$cpost
  seq.dpost[[3]][r, ] <- fit1.gwen$vb.post$dpost
  
  # save the settings, results, and last model fits for convergence examples
  res1 <- list(theta=seq.theta, alpha=seq.alpha, lambda=seq.lambda, 
               delta=seq.delta, zeta=seq.zeta, aprior=seq.aprior, 
               bprior=seq.bprior, apost=seq.apost, bpost=seq.bpost, 
               cpost=seq.cpost, dpost=seq.dpost, mse.mu=mse.mu, cor.mu=cor.mu,
               seq.mSNR)
  fit1 <- list(fit1.igauss1=fit1.igauss1, fit1.igauss2=fit1.igauss2,
               fit1.gwen=fit1.gwen)
  save(set1, res1, fit1, 
       file=paste(path.res, "simulations_igaussian_res1.RData", sep=""))
  
}


write.table(res1, file=paste(path.res, "test.csv", sep=""))
test <- read.table(paste(path.res, "test.csv", sep=""))

colnames(test)
str(test[, c(1:100)])

cbind(res1$theta$`inv Gauss`[1, ],
      as.numeric(test[1, c(1:100)]))


load(paste(path.res, "simulations_igaussian_res1.RData", sep=""))
### convergence
# elbo convergence
omar <- par()$mar
par(mfrow=c(1, 2), mar=omar + c(0, 1, 0, 0))
plot(fit1.igauss1$seq.elbo[, 1], type="l", xlab="iteration", 
     ylab="ELBO", main="a)", ylim=range(fit1.igauss1$seq.elbo))
for(d in 1:ncol(fit1.igauss1$seq.elbo)) {
  lines(fit1.igauss1$seq.elbo[, d], col=d)
}
plot(fit1.igauss2$seq.elbo[, 1], type="l", xlab="iteration", 
     ylab="ELBO", main="b)", ylim=range(fit1.igauss2$seq.elbo))
for(d in 1:ncol(fit1.igauss2$seq.elbo)) {
  lines(fit1.igauss2$seq.elbo[, d], col=d)
}
par(mfrow=c(1, 1), mar=omar)

# igauss alpha convergences
omar <- par()$mar
par(mfrow=c(1, 2), mar=omar + c(0, 1, 0, 0))
plot(fit1.igauss1$seq.eb$alpha[, 1], type="l", xlab="iteration", 
     ylab=expression(alpha), main="a)", ylim=range(fit1.igauss1$seq.eb$alpha))
lines(fit1.igauss1$seq.eb$alpha[, 2], col=2)
lines(fit1.igauss1$seq.eb$alpha[, 3], col=3)
lines(fit1.igauss1$seq.eb$alpha[, 4], col=4)
lines(fit1.igauss1$seq.eb$alpha[, 5], col=5)
plot(fit1.igauss2$seq.eb$alpha[, 1], type="l", xlab="iteration", 
     ylab=expression(alpha), main="b)", ylim=range(fit1.igauss2$seq.eb$alpha))
lines(fit1.igauss2$seq.eb$alpha[, 2], col=2)
lines(fit1.igauss2$seq.eb$alpha[, 3], col=3)
lines(fit1.igauss2$seq.eb$alpha[, 4], col=4)
lines(fit1.igauss2$seq.eb$alpha[, 5], col=5)
legend("topleft", legend=bquote(.(parse(text=paste(
  "alpha", "[", c(1:nclass), "]", sep="")))),
  lty=1, col=c(1:nclass))
par(mfrow=c(1, 1), mar=omar)

# theta convergence
omar <- par()$mar
par(mfrow=c(1, 2), mar=omar + c(0, 1, 0, 0))
plot(fit1.igauss1$seq.eb$theta[, 1], type="l", xlab="iteration", 
     ylab=expression(theta), main="a)", ylim=range(fit1.igauss1$seq.eb$theta))
lines(fit1.igauss1$seq.eb$theta[, D/nclass + 1], col=2)
lines(fit1.igauss1$seq.eb$theta[, 2*D/nclass + 1], col=3)
lines(fit1.igauss1$seq.eb$theta[, 3*D/nclass + 1], col=4)
lines(fit1.igauss1$seq.eb$theta[, 4*D/nclass + 1], col=5)
plot(fit1.igauss2$seq.eb$theta[, 1], type="l", xlab="iteration", 
     ylab=expression(theta), main="b)", ylim=range(fit1.igauss2$seq.eb$theta))
lines(fit1.igauss2$seq.eb$theta[, D/nclass + 1], col=2)
lines(fit1.igauss2$seq.eb$theta[, 2*D/nclass + 1], col=3)
lines(fit1.igauss2$seq.eb$theta[, 3*D/nclass + 1], col=4)
lines(fit1.igauss2$seq.eb$theta[, 4*D/nclass + 1], col=5)
legend("topleft", legend=bquote(.(parse(text=paste(
  "theta", "[", c(1:nclass), "]", sep="")))),
  lty=1, col=c(1:nclass))
par(mfrow=c(1, 1), mar=omar)

# igauss lambda convergence
omar <- par()$mar
par(mfrow=c(1, 2), mar=omar + c(0, 1, 0, 0))
plot(fit1.igauss1$seq.eb$lambda, type="l", xlab="iteration", 
     ylab=expression(lambda), main="a)")
plot(fit1.igauss2$seq.eb$lambda, type="l", xlab="iteration", 
     ylab=expression(lambda), main="b)")
par(mfrow=c(1, 1), mar=omar)

### posterior mean estimates against true values
# sigma^2
omar <- par()$mar
par(mfrow=c(1, 3), mar=omar + c(0, 1, 0, 0))
plot(sigma^2, 2*fit1.igauss1$vb.post$zeta/(n + p - 1), xlab=expression(sigma^2),
     ylab=expression(hat(E(sigma^2 ~ "|" ~ y))), main="a)")
legend("bottomright", paste("MSE=", round(mean((
  sigma^2 - 2*fit1.igauss1$vb.post$zeta/(n + p - 1))^2), 2), sep=""), bty="n")
abline(a=0, b=1, lty=2)
plot(sigma^2, 2*fit1.igauss2$vb.post$zeta/(n - 2.5), xlab=expression(sigma^2),
     ylab=expression(hat(E(sigma^2 ~ "|" ~ y))), main="b)")
legend("bottomright", paste("MSE=", round(mean((
  sigma^2 - 2*fit1.igauss2$vb.post$zeta/(n + p - 1))^2), 2), sep=""), bty="n")
abline(a=0, b=1, lty=2)
plot(sigma^2, fit1.gwen$vb.post$dpost/(fit1.gwen$vb.post$cpost - 1), 
     xlab=expression(sigma^2), ylab=expression(hat(E(sigma^2 ~ "|" ~ y))), 
     main="c)")
legend("bottomright", paste("MSE=", round(mean((
  sigma^2 - fit1.gwen$vb.post$dpost/(fit1.gwen$vb.post$cpost - 1))^2), 2), 
  sep=""), bty="n")
abline(a=0, b=1, lty=2)
par(mfrow=c(1, 1), mar=omar)

# gamma^2
est1 <- sqrt(fit1.igauss1$vb.post$delta*
               fit1.igauss1$seq.eb$theta[fit1.igauss1$iter$eb, ]^2/
               fit1.igauss1$seq.eb$lambda[fit1.igauss1$iter$eb, ])*
  ratio_besselK_cpp(sqrt(fit1.igauss1$vb.post$delta*fit1.igauss1$seq.eb$lambda[
    fit1.igauss1$iter$eb, ]/
      fit1.igauss1$seq.eb$theta[fit1.igauss1$iter$eb, ]^2), p)
est2 <- sqrt(fit1.igauss2$vb.post$delta*
               fit1.igauss2$seq.eb$theta[fit1.igauss2$iter$eb, ]^2/
               fit1.igauss2$seq.eb$lambda[fit1.igauss2$iter$eb, ])*
  ratio_besselK_cpp(sqrt(fit1.igauss2$vb.post$delta*fit1.igauss2$seq.eb$lambda[
    fit1.igauss2$iter$eb, ]/
      fit1.igauss2$seq.eb$theta[fit1.igauss2$iter$eb, ]^2), p)
omar <- par()$mar
par(mfrow=c(1, 3), mar=omar + c(0, 1, 0, 0))
plot(gamma^2, est1, xlab=expression(gamma^2),
     ylab=expression(hat(E(gamma^2 ~ "|" ~ y))), main="a)")
legend("bottomright", paste("MSE=", round(mean((gamma^2 - est1)^2), 2), sep=""), 
       bty="n") 
abline(a=0, b=1, lty=2)
plot(gamma^2, est2, xlab=expression(gamma^2),
     ylab=expression(hat(E(gamma^2 ~ "|" ~ y))), main="b)")
legend("bottomright", paste("MSE=", round(mean((gamma^2 - est2)^2), 2), sep=""), 
       bty="n")
abline(a=0, b=1, lty=2)
plot(gamma^2, fit1.gwen$vb.post$bpost/(fit1.gwen$vb.post$apost - 1),
     xlab=expression(gamma^2), ylab=expression(hat(E(gamma^2 ~ "|" ~ y))), 
     main="c)")
legend("bottomright", paste("MSE=", round(mean((
  gamma^2 - fit1.gwen$vb.post$bpost/(fit1.gwen$vb.post$apost - 1))^2), 2), 
  sep=""), bty="n")
abline(a=0, b=1, lty=2)
par(mfrow=c(1, 1), mar=omar)

# beta
omar <- par()$mar
par(mfrow=c(1, 3), mar=omar + c(0, 1, 0, 0))
plot(beta, fit1.igauss1$vb.post$mu, xlab=expression(beta),
     ylab=expression(hat(E(beta ~ "|" ~ y))), main="a)")
legend("bottomright", paste("MSE=", round(mean((
  beta - fit1.igauss1$vb.post$mu)^2), 2), sep=""), bty="n") 
abline(a=0, b=1, lty=2)
plot(beta, fit1.igauss2$vb.post$mu, xlab=expression(beta),
     ylab=expression(hat(E(beta ~ "|" ~ y))), main="b)")
legend("bottomright", paste("MSE=", round(mean((
  beta - fit1.igauss2$vb.post$mu)^2), 2), sep=""), bty="n") 
abline(a=0, b=1, lty=2)
plot(beta, fit1.gwen$vb.post$mu, xlab=expression(beta),
     ylab=expression(hat(E(beta ~ "|" ~ y))), main="c)")
legend("bottomright", paste("MSE=", round(mean((
  beta - fit1.gwen$vb.post$mu)^2), 2), sep=""), bty="n") 
abline(a=0, b=1, lty=2)
par(mfrow=c(1, 1), mar=omar)

### EB estimation
# prior mean estimates of gamma^2 against true values
omar <- par()$mar
par(mfrow=c(1, 3), mar=omar + c(0, 1, 0, 0))
plot(theta, 1/(C %*% fit1.igauss1$seq.eb$alpha[fit1.igauss1$iter$eb, ]),
     xlab=expression(E(gamma^2)), ylab=expression(hat(E(gamma^2))), main="a)")
legend("bottomright", paste("MSE=", round(mean((
  theta - 1/(C %*% fit1.igauss1$seq.eb$alpha[fit1.igauss1$iter$eb, ]))^2), 2), 
  sep=""), bty="n") 
abline(a=0, b=1, lty=2)
plot(theta, 1/(C %*% fit1.igauss2$seq.eb$alpha[fit1.igauss2$iter$eb, ]),
     xlab=expression(E(gamma^2)), ylab=expression(hat(E(gamma^2))), main="b)")
legend("bottomright", paste("MSE=", round(mean((
  theta - 1/(C %*% fit1.igauss2$seq.eb$alpha[fit1.igauss2$iter$eb, ]))^2), 2), 
  sep=""), bty="n")
abline(a=0, b=1, lty=2)
plot(theta, rep(fit1.gwen$seq.eb$bprior[fit1.gwen$iter$eb, ]/
                  (fit1.gwen$seq.eb$aprior[fit1.gwen$iter$eb, ] - 1), 
                each=D/nclass),
     xlab=expression(E(gamma^2)), ylab=expression(hat(E(gamma^2))), main="c)")
legend("bottomright", paste("MSE=", round(mean((
  theta - rep(fit1.gwen$seq.eb$bprior[fit1.gwen$iter$eb, ]/
                (fit1.gwen$seq.eb$aprior[fit1.gwen$iter$eb, ] - 1), 
              each=D/nclass))^2), 2), sep=""), bty="n")
abline(a=0, b=1, lty=2)
par(mfrow=c(1, 1), mar=omar)

