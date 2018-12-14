### paths
path.cppcode <- "/Users/magnusmunch/Documents/OneDrive/PhD/project_cambridge/code/src/"
path.rcode <- "/Users/magnusmunch/Documents/OneDrive/PhD/project_cambridge/code/R/"

### libraries
library(Rcpp)
library(statmod)
library(gsl)

### Compile and source functions
sourceCpp(paste(path.cppcode, "myfunctions.cpp", sep=""))
sourceCpp(paste(path.cppcode, "functions.cpp", sep=""))
source(paste(path.rcode, "functions.R", sep=""))

### testing estimation functions
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
thetaold=rep(1, D)
lambdaold=0.001
deltaold <- rep(1, D)
aold=rep(0.001, D)
zetaold=rep(1, D)
vold <- sqrt(deltaold*theta^2/lambda)*
  ratio_besselK_cpp(sqrt(lambda*deltaold)/theta, p)
elboold <- rep(-Inf, D)
elbo.constold <- rep(-Inf, D)
svd.x <- svd(x)
sv=svd.x$d
uty=t(svd.x$u) %*% y # size n x D
yty=colSums(y*y)

new.vb <- rbind(delta=deltaold, zeta=zetaold, a=aold, v=vold, elbo=elboold,
                elbo.const=elbo.constold)
new.eb <- list(lambda=lambdaold, theta=thetaold)

elbo.seq <- numeric(0)
alpha.seq <- numeric(0)
lambda.seq <- numeric(0)
theta.seq <- numeric(0)
ebiter <- 20
vbiter <- 5
for(i in c(1:ebiter)) {
  
  for(j in c(1:vbiter)) {
    new.vb <- sapply(c(1:D), function(d) {
      single.vb.update(new.vb[2, d], new.vb[3, d], new.eb$lambda, 
                       new.eb$theta[d], sv, n, p, uty[, d], yty[d])})
  }
  
  new.eb <- eb.update(new.vb[4, ], new.vb[3, ], new.vb[6, ], C, D)
  elbo.seq <- rbind(elbo.seq, new.eb$elbo)
  alpha.seq <- rbind(alpha.seq, new.eb$alpha)
  lambda.seq <- rbind(lambda.seq, new.eb$lambda)
  theta.seq <- rbind(theta.seq, new.eb$theta)
  
}

cbind((2.2*y[, 1])^3, 2.2^3*y[, 1]^3)

plot(rowSums(elbo.seq), type="l")
plot(lambda.seq, type="l")
plot(alpha.seq[, 1], type="l", ylim=range(alpha.seq))
lines(alpha.seq[, 2], type="l", col=2)
lines(alpha.seq[, 3], type="l", col=3)
lines(alpha.seq[, 4], type="l", col=4)
lines(alpha.seq[, 5], type="l", col=5)


plot(theta.seq[, 1], type="l", ylim=range(theta.seq))
lines(theta.seq[, 1*D/nclass + 1], type="l", col=2)
lines(theta.seq[, 2*D/nclass + 1], type="l", col=3)
lines(theta.seq[, 3*D/nclass + 1], type="l", col=4)
lines(theta.seq[, 4*D/nclass + 1], type="l", col=5)

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
                                        maxit.eb=20, maxit.vb=2, trace=TRUE), 
                           init=list(alpha=c(1, rep(0, nclass - 1)), 
                                     lambda=0.001, a=rep(0.001, D), 
                                     zeta=rep(1000, D)),
                           test=list(ind.var=FALSE))

fit1.igauss2 <- est.igauss(x, y, C, 
                           control=list(epsilon.eb=1e-3, epsilon.vb=1e-3, 
                                        maxit.eb=20, maxit.vb=2, trace=TRUE), 
                           init=list(alpha=c(1, rep(0, nclass - 1)), 
                                     lambda=0.001, a=rep(0.001, D), 
                                     zeta=rep(1000, D)),
                           test=list(ind.var=TRUE))

fit1.gwen <- est.gwen(x, y, eqid=rep(c(1:nclass), each=D/nclass),
                      control=list(epsilon.eb=1e-3, epsilon.vb=1e-3, 
                                   epsilon.opt=1e-6, maxit.eb=20, maxit.vb=2, 
                                   maxit.opt=20, conv.vb="elbo", 
                                   trace=TRUE), 
                      init=list(aprior=0.001, bprior=0.001))


### convergence
# igauss alpha convergences
omar <- par()$mar
par(mfrow=c(1, 2), mar=omar + c(0, 1, 0, 0))
plot(fit1.igauss1$eb.seq$alpha[, 1], type="l", xlab="iteration", 
     ylab=expression(alpha), main="a)", ylim=range(fit1.igauss1$eb.seq$alpha))
lines(fit1.igauss1$eb.seq$alpha[, 2], col=2)
lines(fit1.igauss1$eb.seq$alpha[, 3], col=3)
lines(fit1.igauss1$eb.seq$alpha[, 4], col=4)
lines(fit1.igauss1$eb.seq$alpha[, 5], col=5)
plot(fit1.igauss2$eb.seq$alpha[, 1], type="l", xlab="iteration", 
     ylab=expression(alpha), main="b)", ylim=range(fit1.igauss2$eb.seq$alpha))
lines(fit1.igauss2$eb.seq$alpha[, 2], col=2)
lines(fit1.igauss2$eb.seq$alpha[, 3], col=3)
lines(fit1.igauss2$eb.seq$alpha[, 4], col=4)
lines(fit1.igauss2$eb.seq$alpha[, 5], col=5)
legend("topleft", legend=bquote(.(parse(text=paste(
  "alpha", "[", c(1:nclass), "]", sep="")))),
  lty=1, col=c(1:nclass))
par(mfrow=c(1, 1), mar=omar)

# theta convergence
omar <- par()$mar
par(mfrow=c(1, 2), mar=omar + c(0, 1, 0, 0))
plot(fit1.igauss1$eb.seq$theta[, 1], type="l", xlab="iteration", 
     ylab=expression(theta), main="a)", ylim=range(fit1.igauss1$eb.seq$theta))
lines(fit1.igauss1$eb.seq$theta[, D/nclass + 1], col=2)
lines(fit1.igauss1$eb.seq$theta[, 2*D/nclass + 1], col=3)
lines(fit1.igauss1$eb.seq$theta[, 3*D/nclass + 1], col=4)
lines(fit1.igauss1$eb.seq$theta[, 4*D/nclass + 1], col=5)
plot(fit1.igauss2$eb.seq$theta[, 1], type="l", xlab="iteration", 
     ylab=expression(theta), main="b)", ylim=range(fit1.igauss2$eb.seq$theta))
lines(fit1.igauss2$eb.seq$theta[, D/nclass + 1], col=2)
lines(fit1.igauss2$eb.seq$theta[, 2*D/nclass + 1], col=3)
lines(fit1.igauss2$eb.seq$theta[, 3*D/nclass + 1], col=4)
lines(fit1.igauss2$eb.seq$theta[, 4*D/nclass + 1], col=5)
legend("topleft", legend=bquote(.(parse(text=paste(
  "theta", "[", c(1:nclass), "]", sep="")))),
  lty=1, col=c(1:nclass))
par(mfrow=c(1, 1), mar=omar)

# igauss lambda convergence
omar <- par()$mar
par(mfrow=c(1, 2), mar=omar + c(0, 1, 0, 0))
plot(fit1.igauss1$eb.seq$lambda, type="l", xlab="iteration", 
     ylab=expression(lambda), main="a)")
plot(fit1.igauss2$eb.seq$lambda, type="l", xlab="iteration", 
     ylab=expression(lambda), main="b)")
par(mfrow=c(1, 1), mar=omar)

### posterior mean estimates against true values
# sigma^2
omar <- par()$mar
par(mfrow=c(1, 3), mar=omar + c(0, 1, 0, 0))
plot(sigma^2, 2*fit1.igauss1$vb.post$zeta/(n + p - 1), xlab=expression(sigma^2),
     ylab=expression(hat(E(sigma^2 ~ "|" ~ y))), main="a)")
plot(sigma^2, 2*fit1.igauss2$vb.post$zeta/(n - 2.5), xlab=expression(sigma^2),
     ylab=expression(hat(E(sigma^2 ~ "|" ~ y))), main="b)")
plot(sigma^2, fit1.gwen$vb.post$dpost/(fit1.gwen$vb.post$cpost - 1), 
     xlab=expression(sigma^2), ylab=expression(hat(E(sigma^2 ~ "|" ~ y))), 
     main="c)")
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
     ylab=expression(hat(E(gamma^2 ~ "|" ~ y))), main="a)")
abline(a=0, b=1, lty=2)
plot(gamma^2, est2, xlab=expression(gamma^2),
     ylab=expression(hat(E(gamma^2 ~ "|" ~ y))), main="b)")
abline(a=0, b=1, lty=2)
plot(gamma^2, fit1.gwen$vb.post$bpost/(fit1.gwen$vb.post$apost - 1),
     xlab=expression(gamma^2), ylab=expression(hat(E(gamma^2 ~ "|" ~ y))), 
     main="c)")
abline(a=0, b=1, lty=2)
par(mfrow=c(1, 1), mar=omar)

# beta
omar <- par()$mar
par(mfrow=c(1, 3), mar=omar + c(0, 1, 0, 0))
plot(beta, fit1.igauss1$vb.post$mu, xlab=expression(beta),
     ylab=expression(hat(E(beta ~ "|" ~ y))), main="a)")
abline(a=0, b=1, lty=2)
plot(beta, fit1.igauss2$vb.post$mu, xlab=expression(beta),
     ylab=expression(hat(E(beta ~ "|" ~ y))), main="b)")
abline(a=0, b=1, lty=2)
plot(beta, t(fit1.gwen$vb.post$mu), xlab=expression(beta),
     ylab=expression(hat(E(beta ~ "|" ~ y))), main="c)")
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

# check log besselK function
x <- seq(0.1, 1000, 0.1)
nu <- seq(1, 1000, 0.5)
grid <- expand.grid(x, nu)
library(gsl)
bessel_lnKnu


set.seed(123)
n <- 200
p <- 100
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
y <- matrix(rnorm(n), nrow=n, ncol=1)
svd.x <- svd(x)
V <- svd.x$v
U <- svd.x$u
d <- svd.x$d
ad <- 15.3
uty <- as.numeric(t(U) %*% y)

trace1 <- as.numeric(t(y) %*% x %*% solve(t(x) %*% x + ad*diag(p)) %*% t(x) %*% y)
trace2 <- sum(d^2*uty^2/(d^2 + ad))
all.equal(trace1, trace2)


