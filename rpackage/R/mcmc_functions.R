.draw_beta <- function(par, x, y, n, p) {
  gamma <- par$gamma
  sigma <- par$sigma
  tau <- par$tau
  u <- rnorm(p, 0, (tau*sigma*gamma)^2)
  delta <- rnorm(n, 0, sigma^2)
  beta <- u + (t(x)*gamma^2) %*% solve(x %*% (t(x)*gamma^2) + diag(n)/tau^2) %*%
    (y - x %*% u - delta)  
  return(as.numeric(beta))
}

.draw_gamma <- function(par, C, p) {
  alphaf <- par$alphaf
  lambdaf <- par$lambdaf
  beta <- par$beta
  sigma <- par$sigma
  tau <- par$tau
  gamma <- sapply(1:p, function(j) {
    sqrt(rgig(n=1, chi=beta[j]^2/(sigma*tau)^2 + lambdaf, 
              psi=lambdaf*as.numeric(C %*% alphaf)[j]^2, lambda=-1))})
  return(gamma)
}

.draw_tau <- function(par, Z, p) {
  alphad <- par$alphad
  lambdad <- par$lambdad
  beta <- par$beta
  sigma <- par$sigma
  gamma <- par$gamma
  tau <- sqrt(rgig(1, -(p + 1)/2, chi=sum(beta^2/(gamma*sigma)^2) + lambdad, 
                     psi=lambdad*as.numeric(Z %*% alphad)^2))
  return(tau)
}

d <- 1
par <- c(sapply(par[-c(5:8)], "[[", d), par[c(5:8)])
x <- x[[d]]
y <- y[, d]
p <- p[d]
.draw_sigma <- function(par, x, y, n, p) {
  beta <- par$beta
  gamma <- par$gamma
  tau <- par$tau
  sigma <- 1/sqrt(rgamma(1, (n + p + 1)/2, 
                         0.5*(sum(y^2) -2*sum(y*as.numeric(x %*% beta)) + 
                                sum(as.numeric(x %*% beta)^2) + 
                                sum(beta^2/gamma^2)/tau^2)))
  return(sigma)
}

.draw_alphaf <- function(par, Cmat) {
  gamma <- unlist(par$gamma)
  lambdaf <- par$lambdaf
  mat <- solve(t(Cmat) %*% diag(gamma^2) %*% Cmat)
  alphaf <- rmvnorm(1, rowSums(mat %*% t(Cmat)), mat/lambdaf)
  return(as.numeric(alphaf))
}

.draw_alphad <- function(par, Zmat) {
  tau <- par$tau
  lambdad <- par$lambdad
  mat <- solve(t(Zmat) %*% diag(tau^2) %*% Zmat)
  alphad <- rmvnorm(1, rowSums(mat %*% t(Zmat)), mat/lambdad)
  return(as.numeric(alphad))
}

.draw_lambdaf <- function(par, Cmat, p, D) {
  alphaf <- par$alphaf
  gamma <- unlist(par$gamma)
  lambdaf <- rgamma(1, p*D/2, 0.5*(sum(as.numeric(
    Cmat %*% alphaf)^2*gamma^2) - 2*sum(Cmat %*% alphaf) + sum(1/gamma^2)))
  return(lambdaf)
}

.draw_lambdad <- function(par, Zmat, D) {
  alphad <- par$alphad
  tau <- par$tau
  lambdad <- rgamma(1, D/2, 0.5*(sum(as.numeric(Zmat %*% alphad)^2*tau^2) -
                                     2*sum(Zmat %*% alphad) + sum(1/tau^2)))
  return(lambdad)
}

.draw_post <- function(par, x, y, C, Z, n, p, D) {
  Cmat <- Reduce("rbind", C)
  for(d in 1:D) {
    par$beta[[d]] <- .draw_beta(c(sapply(par[-c(5:8)], "[[", d), par[c(5:8)]), 
                                x=x[[d]], y=y[, d], n=n, p=p[d])
    par$gamma[[d]] <- .draw_gamma(c(sapply(par[-c(5:8)], "[[", d), par[c(5:8)]), 
                                  C[[d]], p[d])
    par$tau[[d]] <- .draw_tau(c(sapply(par[-c(5:8)], "[[", d), par[c(5:8)]), 
                              Z[d, ], p[d])
    par$sigma[[d]] <- .draw_sigma(c(sapply(par[-c(5:8)], "[[", d), par[c(5:8)]), 
                                  x[[d]], y[, d], n, p[d])
  }
  par$alphaf <- .draw_alphaf(par, Cmat)
  par$alphad <- .draw_alphad(par, Z)
  par$lambdaf <- .draw_lambdaf(par, Cmat, sum(p), D)
  par$lambdad <- .draw_lambdad(par, Z, D)
  return(par)
}

library(GeneralizedHyperbolic)
set.seed(123)
n <- 50
D <- 10
p <- rep(100, D)
x <- sapply(1:D, function(d) {matrix(rnorm(n*p[d]), nrow=n, ncol=p[d])}, 
            simplify=FALSE)
beta <- rbind(cbind(matrix(rnorm(p[1]*D/2, 0, 0.1), nrow=p[1]/2, ncol=D/2), 
                    matrix(rnorm(p[1]*D/2, 0, 0.5), nrow=p[1]/2, ncol=D/2)),
              cbind(matrix(rnorm(p[1]*D/2, 0, 0.5), nrow=p[1]/2, ncol=D/2), 
                    matrix(rnorm(p[1]*D/2, 0, 1), nrow=p[1]/2, ncol=D/2)))
y <- sapply(1:D, function(d) {x[[d]] %*% beta[, d] + rnorm(n)})
Z <- cbind(1, rep(c(0, 1), each=D/2))
C <- sapply(1:D, function(d) {cbind(1, rep(c(0, 1), each=p[d]/2))}, 
            simplify=FALSE)
Cmat <- Reduce("rbind", C)
par <- list(beta=sapply(1:D, function(d) {rnorm(p[d])}, simplify=FALSE), 
            gamma=sapply(1:D, function(d) {rchisq(p[d], 1) + 1}, 
                         simplify=FALSE), 
            tau=rchisq(D, 1) + 1, sigma=rchisq(D, 1) + 1, 
            alphaf=rchisq(2, 1), alphad=rchisq(2, 1), 
            lambdaf=rchisq(1, 1) + 1, lambdad=rchisq(1, 1) + 1)



.draw_beta(c(sapply(par[-c(5:8)], "[[", d), par[c(5:8)]), x[[d]], y[, d], n, p[d]) 
.draw_gamma(c(sapply(par[-c(5:8)], "[[", d), par[c(5:8)]), C[[d]], p[d]) 
.draw_tau(c(sapply(par[-c(5:8)], "[[", d), par[c(5:8)]), Z[d, ], p[d]) 
.draw_sigma(c(sapply(par[-c(5:8)], "[[", d), par[c(5:8)]), x[[d]], y[, d], n, p[d])   
.draw_alphaf(par, Cmat) 
.draw_alphad(par, Z) 
.draw_lambdaf(par, Cmat, p, D)
.draw_lambdad(par, Z, D) 
par <- .draw_post(par, x, y, C, Z, n, p, D)


