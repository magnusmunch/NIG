# sample one beta in the conjugate case
.sample.beta.conj <- function(gamma, sigma, v, sv, svvt, y) {
  uvec <- rnorm(p, 0, gamma*sigma)
  dvec <- rnorm(n, 0, 1)
  beta <- uvec + as.numeric(t(t(v)*(sv/(sv^2 + 1/gamma^2))) %*% 
                              (sigma*y - sigma*svvt %*% uvec + dvec))
  return(beta)
}

# sample one gamma in the conjugate case
.sample.gamma.conj <- function(beta, sigma, eta, lambda, theta, p) {
  
  gamma <- sqrt(rgig(1, sum(beta^2)/sigma^2 + lambda, 
                     lambda/theta^2, -0.5*(p + eta)))
  return(gamma)
}

# sample one sigma in the conjugate case
.sample.sigma.conj <- function(beta, gamma, yty, ytx, x, p, n) {
  alpha <- 0.5*(p + n + 1)
  beta <- 0.5*(yty + sum(beta^2)/gamma^2 + sum((x %*% beta)^2) - 2*ytx %*% beta)
  gamma <- 1/sqrt(rgamma(1, alpha, scale=beta))
  return(gamma)
}

# full sample from the posterior
.sample.post.conj <- function(eta, lambda, theta, nsamp, init, v, sv, svvt, y, 
                              yty, ytx, x, p, n) {
  
  # create the output and initial values
  out <- list(beta=matrix(NA, nrow=p, ncol=nsamp), sigma=rep(NA, nsamp),
              gamma=rep(NA, nsamp))
  beta <- init$beta
  sigma <- init$sigma
  
  # run the Gibbs sampler
  for(m in 1:nsamp) {
    gamma <- .sample.gamma.conj(beta, sigma, eta, lambda, theta, p)
    beta <- .sample.beta.conj(gamma, sigma, v, sv, svvt, y)
    sigma <- .sample.sigma.conj(beta, gamma, yty, ytx, x, p, n)
    out$beta[, m] <- beta
    out$sigma[m] <- sigma
    out$gamma[m] <- gamma
  }
  
  return(out)
}

# runs Gibbs sampler
gibbs.model <- function(x, y, hyperprior=c("inv. Gaussian", "inv. Gamma"), 
                        conjugate, init=NULL,
                        control=list(nsamp, parallel=FALSE, ncores=1)) {
  
  # set the number of cores
  ncores <- min(control$ncores, detectCores())
  
  # set fixed parameters
  D <- ncol(y)
  svd.x <- svd(x)
  v <- svd.x$v
  sv <- svd.x$d
  svvt <- t(v)*sv
  yty <- colSums(y^2)
  ytx <- t(y) %*% x
  
  if(hyperprior=="inv. Gaussian") {
    eta <- rep(1, D)
  }
  
  if(conjugate) {
    out <- mclapply(c(1:D), function(d) {
      cinit <- list(beta=init$beta[, d], gamma=init$gamma[d], sigma=init$sigma[d])
      .sample.post.conj(eta[d], lambda[d], theta[d], control$nsamp, cinit, v, sv, 
                        svvt, y[, d], yty[d], ytx[d, ], x, p, n)}, 
      mc.cores=ncores)  
  } else {
    
  }
  
  
  beta <- aperm(sapply(out, function(s) {s$beta}, simplify="array"), c(1, 3, 2))
  sigma <- t(sapply(out, function(s) {s$sigma}))
  gamma <- t(sapply(out, function(s) {s$gamma}))
  
  return(list(beta=beta, sigma=sigma, gamma=gamma))
  
}
