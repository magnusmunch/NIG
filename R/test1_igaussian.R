### paths
path.cppcode <- "/Users/magnusmunch/Documents/OneDrive/PhD/project_cambridge/code/src/"
path.rcode <- "/Users/magnusmunch/Documents/OneDrive/PhD/project_cambridge/code/R/"

### libraries
library(Rcpp)
library(statmod)
library(gsl)

### Compile cpp functions
sourceCpp(paste(path.cppcode, "myfunctions.cpp", sep=""))

source(paste(path.rcode, "functions.R", sep=""))

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
SNR <- 10

beta <- sapply(1:D, function(d) {rnorm(p, 0, sigma[d]*gamma[d])})
x <- matrix(rnorm(n*p, 0, sqrt(SNR/(mean(gamma^2)*p))), nrow=n, ncol=p)
y <- sapply(1:D, function(d) {rnorm(n, x %*% beta[, d], sigma[d])})

fit1.igauss1 <- est.igauss(x, y, C, 
                           control=list(epsilon.eb=1e-3, epsilon.vb=1e-3, 
                                        maxit.eb=20, maxit.vb=2, trace=TRUE), 
                           init <- list(lambda=1, alpha=rep(1, ncol(C)), 
                                        a=rep(0.1, D), zeta=rep(0.1, D)),
                           test=list(ind.var=FALSE))

fit1.igauss2 <- est.igauss(x, y, C, 
                           control=list(epsilon.eb=1e-3, epsilon.vb=1e-3, 
                                        maxit.eb=20, maxit.vb=2, trace=TRUE), 
                           init <- list(lambda=1, alpha=rep(1, ncol(C)), 
                                        a=rep(0.1, D), zeta=rep(0.1, D)),
                           test=list(ind.var=TRUE))

fit1.gwen <- est.gwen(x, y, eqid=rep(c(1:nclass), each=D/nclass),
                      control=list(epsilon.eb=1e-3, epsilon.vb=1e-3, 
                                   epsilon.opt=1e-6, maxit.eb=20, maxit.vb=2, 
                                   maxit.opt=20, conv.vb="elbo", 
                                   trace=TRUE), 
                      init=list(aprior=0.001, bprior=0.001))


# ### convergence
# # ELBOs
# plot(rowSums(fit1.igauss1$elbo.seq), type="l")
# plot(rowSums(fit1.gwen$elbo.seq), type="l")
# 
# # igauss alpha convergences
# plot(fit1.igauss1$eb.seq$alpha[, 1], type="l", ylim=range(fit1.igauss1$eb.seq$alpha))
# lines(fit1.igauss1$eb.seq$alpha[, 2], col=2)
# lines(fit1.igauss1$eb.seq$alpha[, 3], col=3)
# lines(fit1.igauss1$eb.seq$alpha[, 4], col=4)
# lines(fit1.igauss1$eb.seq$alpha[, 5], col=5)
# plot(fit1.igauss2$eb.seq$alpha[, 1], type="l", ylim=range(fit1.igauss2$eb.seq$alpha))
# lines(fit1.igauss2$eb.seq$alpha[, 2], col=2)
# lines(fit1.igauss2$eb.seq$alpha[, 3], col=3)
# lines(fit1.igauss2$eb.seq$alpha[, 4], col=4)
# lines(fit1.igauss2$eb.seq$alpha[, 5], col=5)
# 
# # theta convergence
# plot(fit1.igauss1$eb.seq$theta[, 1], type="l", ylim=range(fit1.igauss1$eb.seq$theta))
# lines(fit1.igauss1$eb.seq$theta[, D/nclass + 1], col=2)
# lines(fit1.igauss1$eb.seq$theta[, 2*D/nclass + 1], col=3)
# lines(fit1.igauss1$eb.seq$theta[, 3*D/nclass + 1], col=4)
# lines(fit1.igauss1$eb.seq$theta[, 4*D/nclass + 1], col=5)
# plot(fit1.igauss2$eb.seq$theta[, 1], type="l", ylim=range(fit1.igauss2$eb.seq$theta))
# lines(fit1.igauss2$eb.seq$theta[, D/nclass + 1], col=2)
# lines(fit1.igauss2$eb.seq$theta[, 2*D/nclass + 1], col=3)
# lines(fit1.igauss2$eb.seq$theta[, 3*D/nclass + 1], col=4)
# lines(fit1.igauss2$eb.seq$theta[, 4*D/nclass + 1], col=5)
# 
# # igauss lambda convergence
# plot(fit1.igauss1$eb.seq$lambda, type="l")
# plot(fit1.igauss2$eb.seq$lambda, type="l")

### posterior mean estimates against true values
# sigma^2
omar <- par()$mar
par(mfrow=c(1, 3), mar=omar + c(0, 1, 0, 0))
plot(sigma^2, 2*fit1.igauss1$vb.post$zeta/(n + p - 1), xlab=expression(sigma^2),
     ylab=expression(E(sigma^2 ~ "|" ~ y)), main="a)")
plot(sigma^2, 2*fit1.igauss2$vb.post$zeta/(n - 2.5), xlab=expression(sigma^2),
     ylab=expression(E(sigma^2 ~ "|" ~ y)), main="b)")
plot(sigma^2, fit1.gwen$vb.post$dpost/(fit1.gwen$vb.post$cpost - 1), 
     xlab=expression(sigma^2), ylab=expression(E(sigma^2 ~ "|" ~ y)), main="c)")
par(mfrow=c(1, 1), mar=omar)

# gamma^2
est1 <- sqrt(fit1.igauss1$vb.post$delta*
               fit1.igauss1$eb.seq$theta[fit1.igauss1$iter$eb, ]^2/
               fit1.igauss1$eb.seq$lambda[fit1.igauss1$iter$eb, ])*
  ratio_besselK_cpp(sqrt(fit1.igauss1$vb.post$delta*fit1.igauss1$eb.seq$lambda[
    fit1.igauss1$iter$eb, ]/fit1.igauss1$vb.post$delta*
      fit1.igauss1$eb.seq$theta[fit1.igauss1$iter$eb, ]^2), p)
est2 <- sqrt(fit1.igauss2$vb.post$delta*
               fit1.igauss2$eb.seq$theta[fit1.igauss2$iter$eb, ]^2/
               fit1.igauss2$eb.seq$lambda[fit1.igauss2$iter$eb, ])*
  ratio_besselK_cpp(sqrt(fit1.igauss2$vb.post$delta*fit1.igauss2$eb.seq$lambda[
    fit1.igauss2$iter$eb, ]/fit1.igauss2$vb.post$delta*
      fit1.igauss2$eb.seq$theta[fit1.igauss2$iter$eb, ]^2), p)
omar <- par()$mar
par(mfrow=c(1, 3), mar=omar + c(0, 1, 0, 0))
plot(gamma^2, est1, xlab=expression(gamma^2),
     ylab=expression(E(gamma^2 ~ "|" ~ y)), main="a)")
abline(a=0, b=1, lty=2)
plot(gamma^2, est2, xlab=expression(gamma^2),
     ylab=expression(E(gamma^2 ~ "|" ~ y)), main="b)")
abline(a=0, b=1, lty=2)
plot(gamma^2, fit1.gwen$vb.post$bpost/(fit1.gwen$vb.post$apost - 1),
     ylab=expression(E(gamma^2 ~ "|" ~ y)), main="c)")
abline(a=0, b=1, lty=2)
par(mfrow=c(1, 1), mar=omar)

# beta
omar <- par()$mar
par(mfrow=c(1, 3), mar=omar + c(0, 1, 0, 0))
plot(beta, fit1.igauss1$vb.post$mu, xlab=expression(beta),
     ylab=expression(E(beta ~ "|" ~ y)), main="a)")
abline(a=0, b=1, lty=2)
plot(beta, fit1.igauss2$vb.post$mu, xlab=expression(beta),
     ylab=expression(E(beta ~ "|" ~ y)), main="b)")
abline(a=0, b=1, lty=2)
plot(beta, t(fit1.gwen$vb.post$mu), xlab=expression(beta),
     ylab=expression(E(beta ~ "|" ~ y)), main="c)")
abline(a=0, b=1, lty=2)
par(mfrow=c(1, 1), mar=omar)

### EB estimation
# prior mean estimates of gamma^2 against true values
omar <- par()$mar
par(mfrow=c(1, 3), mar=omar + c(0, 1, 0, 0))
plot(theta, 1/(C %*% fit1.igauss1$eb.seq$alpha[fit1.igauss1$iter$eb, ]),
     xlab=expression(E(gamma^2)), ylab=expression(hat(E(gamma^2))), main="a)")
abline(a=0, b=1, lty=2)
plot(theta, 1/(C %*% fit1.igauss2$eb.seq$alpha[fit1.igauss2$iter$eb, ]),
     xlab=expression(E(gamma^2)), ylab=expression(hat(E(gamma^2))), main="b)")
abline(a=0, b=1, lty=2)
plot(theta, rep(fit1.gwen$eb.seq$bprior[fit1.gwen$iter$eb + 1, ]/
                  (fit1.gwen$eb.seq$aprior[fit1.gwen$iter$eb + 1, ] - 1), 
                each=D/nclass),
     xlab=expression(E(gamma^2)), ylab=expression(hat(E(gamma^2))), main="c)")
abline(a=0, b=1, lty=2)
par(mfrow=c(1, 1), mar=omar)





# f.root <- function(x, vary, q) {
#   gamma(3/2)*(1 - q - pgamma(0.5*x*3/vary, 3/2))
# }
# sigma.squared.hat <- (3/5)*apply(y, 2, function(c) {
#   uniroot(f.root, interval=c(0.0001, var(c)), vary=var(c), q=0.95)$root})
#
# f.root <- function(x, vary, n, p, prob) {
#   pgamma(1/vary, 0.5*(n + p + 1), x) + prob - 1
# }
# 
# curve(f.root(x, var(y[, 1]), n, p, 0.9), 0.0001, 100)
# apply(y, 2, function(c) {uniroot(f.root, c(0.0001, var(c)), vary=var(c), n=n, 
#                                  p=p, prob=0.9)$root})
# 
# zeta.start <- 
# library(invgamma)
# curve(dinvgamma(x, 0.5*(p + n + 1), 1), 0.001, 1)
# 
# 
# theta.start <- mean(sigma.squared.hat)
# lambda.start <- length(sigma.squared.hat)/sum(1/sigma.squared.hat - 1/theta.start)
# curve(dinvgauss(x, theta.start, lambda.start), 0.01, 10)

# check log besselK function
x <- seq(0.1, 1000, 0.1)
nu <- seq(1, 1000, 0.5)
grid <- expand.grid(x, nu)
library(gsl)
bessel_lnKnu












