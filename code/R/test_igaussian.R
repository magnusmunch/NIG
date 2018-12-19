### paths
path.cppcode <- "/Users/magnusmunch/Documents/OneDrive/PhD/cambridge/code/src/"
path.rcode <- "/Users/magnusmunch/Documents/OneDrive/PhD/cambridge/code/R/"

### installation of package
library(devtools)
install_github("magnusmunch/cambridge/code", local=FALSE,
               auth_token="da10f2b37c513e3383c5b2e0aa1300288329c636")

### libraries
library(Rcpp)
library(statmod)
library(gsl)
library(cambridge)

### Compile and source functions
sourceCpp(paste(path.cppcode, "myfunctions.cpp", sep=""))
sourceCpp(paste(path.cppcode, "functions.cpp", sep=""))
source(paste(path.rcode, "functions.R", sep=""))

use_travis("cambridge")

### testing estimation functions
set.seed(567)
n <- 100
p <- 200
D <- 100
nclass <- 5

alpha <- c(1:nclass)
C <- model.matrix(~ 1 + factor(rep(c(1:nclass), each=D/nclass)))
theta <- as.numeric(1/(C %*% alpha))
gamma <- sqrt(theta)
sigma <- rep(1, D)
SNR <- 5

beta <- sapply(1:D, function(d) {rnorm(p, 0, sigma[d]*gamma[d])})
x <- matrix(rnorm(n*p, 0, sqrt(SNR/(mean(gamma^2)*p))), nrow=n, ncol=p)
y <- sapply(1:D, function(d) {rnorm(n, x %*% beta[, d], sigma[d])})

### simulations
set.seed(123)
n <- 100
p <- 200
D <- 100
nclass <- 5

alpha <- c(1:nclass)
C <- model.matrix(~ 1 + factor(rep(c(1:nclass), each=D/nclass)))
theta <- as.numeric(1/(C %*% alpha))
gamma <- sqrt(theta)
sigma <- rep(1, D)
SNR <- 20

beta <- sapply(1:D, function(d) {rnorm(p, 0, sigma[d]*gamma[d])})
x <- matrix(rnorm(n*p, 0, sqrt(SNR/(mean(gamma^2)*p))), nrow=n, ncol=p)
y <- sapply(1:D, function(d) {rnorm(n, x %*% beta[, d], sigma[d])})

fit1.igauss1 <- est.igauss(x, y, C, 
                           control=list(epsilon.eb=1e-3, epsilon.vb=1e-3, 
                                        maxit.eb=300, maxit.vb=2, trace=TRUE), 
                           init=list(alpha=c(1, rep(0, nclass - 1)), 
                                     lambda=0.001, a=rep(0.001, D), 
                                     zeta=rep(1000, D)),
                           test=list(ind=FALSE))


fit1.igauss2 <- est.igauss(x, y, C, 
                           control=list(epsilon.eb=1e-3, epsilon.vb=1e-3, 
                                        maxit.eb=20, maxit.vb=2, trace=TRUE), 
                           init=list(alpha=c(1, rep(0, nclass - 1)), 
                                     lambda=0.001, a=rep(0.001, D), 
                                     zeta=rep(1000, D)),
                           test=list(ind=TRUE))

fit1.gwen <- est.gwen(x, y, eqid=rep(c(1:nclass), each=D/nclass),
                      control=list(epsilon.eb=1e-3, epsilon.vb=1e-3, 
                                   epsilon.opt=1e-6, maxit.eb=20, maxit.vb=2, 
                                   maxit.opt=20, conv.vb="elbo", 
                                   trace=TRUE), 
                      init=list(aprior=0.001, bprior=0.001))

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
    fit1.igauss1$iter$eb, ]/fit1.igauss1$vb.post$delta*
      fit1.igauss1$seq.eb$theta[fit1.igauss1$iter$eb, ]^2), p)
est2 <- sqrt(fit1.igauss2$vb.post$delta*
               fit1.igauss2$seq.eb$theta[fit1.igauss2$iter$eb, ]^2/
               fit1.igauss2$seq.eb$lambda[fit1.igauss2$iter$eb, ])*
  ratio_besselK_cpp(sqrt(fit1.igauss2$vb.post$delta*fit1.igauss2$seq.eb$lambda[
    fit1.igauss2$iter$eb, ]/fit1.igauss2$vb.post$delta*
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
plot(beta, t(fit1.gwen$vb.post$mu), xlab=expression(beta),
     ylab=expression(hat(E(beta ~ "|" ~ y))), main="c)")
legend("bottomright", paste("MSE=", round(mean((
  beta - t(fit1.gwen$vb.post$mu))^2), 2), sep=""), bty="n") 
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
plot(theta, rep(fit1.gwen$seq.eb$bprior[fit1.gwen$iter$eb + 1, ]/
                  (fit1.gwen$seq.eb$aprior[fit1.gwen$iter$eb + 1, ] - 1), 
                each=D/nclass),
     xlab=expression(E(gamma^2)), ylab=expression(hat(E(gamma^2))), main="c)")
legend("bottomright", paste("MSE=", round(mean((
  theta - rep(fit1.gwen$seq.eb$bprior[fit1.gwen$iter$eb + 1, ]/
                (fit1.gwen$seq.eb$aprior[fit1.gwen$iter$eb + 1, ] - 1), 
              each=D/nclass))^2), 2), sep=""), bty="n")
abline(a=0, b=1, lty=2)
par(mfrow=c(1, 1), mar=omar)

