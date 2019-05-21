# sample one beta in the conjugate case
.sample.beta.conj <- function(gamma, sigma, v, sv, svvt, y, n, p) {
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
    beta <- .sample.beta.conj(gamma, sigma, v, sv, svvt, y, n, p)
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

# sample from spike and slab prior
rspikeandslab <- function(n, nu, kappa, sigma, gamma, r) {
  s <- 1 - nu/(nu + kappa)
  a <- gamma/sqrt((1 - s) + s*r)
  b <- gamma/sqrt((1 - s)/r + s)
  z <- rbinom(n, 1, 1 - s)
  beta <- rnorm(n, 0, sigma*(z*a + (1 - z)*b))
  return(beta)
}

# full conditional for latent indicators in spike and slab
.samp.z <- function(w, beta, gamma, sigma, r, nu, kappa, p) {
  s <- 1 - nu/(nu + kappa)
  R <- exp(-((1 - s)/r + r*s - 1)*beta^2/(2*gamma^2*sigma^2))/sqrt(r)
  prob <- 1/(1 + (1 - w)*R/w)
  z <- rbinom(p, 1, prob)
  return(z)
}

# full conditional for indicator probability in spike and slab
.samp.w <- function(z, nu, kappa, p) {
  w <- rbeta(1, nu + sum(z), kappa + p - sum(z))
  return(w)
}

# full conditional for model parameters in spike and slab
.samp.beta <- function(z, gamma, sigma, r, nu, kappa, x, y, n, p) {
  
  s <- 1 - nu/(nu + kappa)
  dg <- (z - r*(z - 1))/((r - 1)*s + 1) 
  if(p > n) {
    uvec <- rnorm(p, 0, gamma*sigma*dg)
    dvec <- rnorm(n, 0, sigma^2)
    mat1 <- t(x)*dg
    mat2 <- x %*% (t(x)*dg); diag(mat2) <- diag(mat2) + 1/gamma^2
    beta <- uvec + as.numeric(mat1 %*% solve(mat2) %*% (y - x %*% uvec - dvec))  
  } else {
    mat <- solve(t(x) %*% x + diag(dg)/gamma^2)
    cmat <- chol(mat)
    beta0 <- rnorm(p, 0, 1)
    beta <- as.numeric(mat %*% t(x) %*% y + sigma*cmat %*% beta0)
  }
  
  return(beta)
}

# Gibbs sampler for spike and slab
gibbs.spikeandslab <- function(x, y, gamma, sigma, r, nu, kappa, init=NULL,
                               control=list(samples=1000, warmup=1000)) {
  
  n <- nrow(x)
  p <- ncol(x)
  
  if(is.null(init$beta)) {
    beta <- rnorm(p)
  } else {
    beta <- init$beta
  }
  if(is.null(init$w)) {
    w <- 0.5
  } else {
    w <- init$w
  }
  if(is.null(init$z)) {
    z <- rbinom(p, 1, 0.5)
  } else {
    z <- init$z
  }
  
  out <- list(beta=matrix(NA, nrow=p, ncol=control$samples), 
              z=matrix(NA, nrow=p, ncol=control$samples), 
              w=numeric(control$samples))
  for(m in 1:(control$samples + control$warmup)) {
    z <- .samp.z(w, beta, gamma, sigma, r, nu, kappa, p)
    w <- .samp.w(z, nu, kappa, p)
    beta <- .samp.beta(z, gamma, sigma, r, nu, kappa, x, y, n, p)
    
    if(m > control$samples) {
      out$beta[, m - control$samples] <- beta
      out$z[, m - control$samples] <- z
      out$w[m - control$samples] <- w  
    }
  }
  
  return(out)
}
