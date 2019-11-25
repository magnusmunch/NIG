# sample beta
.sample.beta <- function(tau, gamma, sigma, v, sv, svvt, y, n, p) {
  uvec <- rnorm(p, 0, tau*gamma*sigma)
  dvec <- rnorm(n, 0, 1)
  beta <- uvec + as.numeric(t(t(v)*(sv/(sv^2 + 1/gamma^2))) %*% 
                              (sigma*y - sigma*svvt %*% uvec + dvec))
  return(beta)
}

# sample gamma
.sample.gamma <- function(beta, sigma, tau, lambdaf, phi, p) {
  gamma <- sqrt(rgig(p, beta^2/(sigma^2*tau^2) + lambdaf, lambdaf/phi^2, -1))
  return(gamma)
}

# sample tau
.sample.tau <- function(beta, sigma, gamma, lambdad, chi, p) {
  tau <- sqrt(rgig(1, sum(beta^2/gamma^2)/sigma^2 + lambdad, lambdad/chi^2, 
                   -0.5*(p + 1)))
  return(tau)
}

# sample sigma
.sample.sigma <- function(beta, gamma, tau, yty, ytx, x, p, n) {
  a <- 0.5*(p + n + 1)
  b <- 0.5*(yty - 2*as.numeric(ytx %*% beta) + sum((x %*% beta)^2) + 
              sum(beta^2/gamma^2)/tau^2)
  sigma <- 1/sqrt(rgamma(1, a, rate=b))
  return(sigma)
}

# full sample from the posterior
.sample <- function(init, hyper, nsamp, v, sv, svvt, y, yty, ytx, x, p, n) {
  
  # create the output and initial values
  out <- list(beta=matrix(NA, nrow=p, ncol=nsamp), sigma=rep(NA, nsamp),
              tau=rep(NA, nsamp), gamma=matrix(NA, nrow=p, ncol=nsamp))
  beta <- init$beta
  sigma <- init$sigma
  tau <- init$tau
  
  # run the Gibbs sampler
  for(m in 1:nsamp) {
    gamma <- .sample.gamma(beta, sigma, tau, hyper$lambdaf, hyper$phi, p)
    tau <- .sample.tau(beta, sigma, gamma, hyper$lambdad, hyper$chi, p)
    beta <- .sample.beta.conj(gamma, sigma, v, sv, svvt, y, n, p)
    sigma <- .sample.sigma(beta, gamma, tau, yty, ytx, x, p, n)
    beta <- .sample.beta(tau, gamma, sigma, v, sv, svvt, y, n, p)
    out$beta[, m] <- beta
    out$gamma[, m] <- gamma
    out$sigma[m] <- sigma
    out$tau[m] <- tau
  }
  return(out)
}

# runs Gibbs sampler
gibbs <- function(x, y, hyper=NULL, init=NULL) {
  
  # set fixed parameters
  D <- ncol(y)
  p <- sapply(x, ncol)
  n <- nrow(y)
  svd.x <- lapply(x, svd)
  v <- lapply(svd.x, "[[", "v")
  sv <- lapply(svd.x, "[[", "d")
  svvt <- lapply(1:D, function(d) {t(v[[d]])*sv[[d]]})
  yty <- colSums(y^2)
  ytx <- t(sapply(1:D, function(d) {t(y[, d]) %*% x[[d]]}))
  
  # complete hyperparameters
  if(is.null(hyper)) {
    hyper <- list(phi=lapply(p, function(s) {rep(1, p)}),
                  chi=rep(1, D), lambdaf=1, lambdad=1)
  }
  if(is.null(init)) {
    init <- list(beta=lapply(p, rnorm), sigma=sqrt(rchisq(D, 0.5)) + 0.5,
                 tau=sqrt(rchisq(D, 0.5) + 0.5))
  }
  
  # sample from posterior
  res <- lapply(1:D, function(d) {
    
  })
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