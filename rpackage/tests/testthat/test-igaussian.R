context("inverse Gaussian model estimation")

################################## VB testing ##################################

###################################
##### have to finish this one #####
###################################
test_that("auxiliary variables calculation with intercept", {
  n <- 20
  p <- 100
  D <- 1
  
  gamma <- rep(1, D)
  sigma <- rep(1, D)
  beta <- sapply(1:D, function(d) {rnorm(p, 0, sigma[d]*gamma[d])})
  x <- matrix(rnorm(n*p, 0, 1), nrow=n, ncol=p)
  y <- sapply(1:D, function(d) {rnorm(n, x %*% beta[, d], sigma[d])})
  lambda <- 0.1
  theta <- 2
  
  y <- y[, 1]
  c <- 1
  a <- 1
  svd.x <- svd(x)
  v <- svd.x$v
  u <- svd.x$u
  d <- svd.x$d
  
  h <- colSums(u)
  g <- as.numeric(t(y) %*% u)
  s <- 1/(n - sum(h^2*d^2/(d^2 + c)))
  
  sumy <- sum(y)
  aux1 <- sum(d^2*h^2/(d^2 + c))
  aux2 <- sum(d^2*g^2/(d^2 + c))
  aux3 <- sum(d^2*h*g/(d^2 + c))
  aux4 <- sum(d^4*h^2/(d^2 + c)^2)
  aux5 <- sum(d^4*g^2/(d^2 + c)^2)
  aux6 <- sum(d^4*h*g/(d^2 + c)^2)
  aux7 <- sum(d^2*h^2/(d^2 + c)^2)
  
  # trSigma (without intercept term)
  trSigma1 <- sum(diag(solve(t(cbind(1, x)) %*% cbind(1, x) + 
                               diag(c(0, rep(c, p))))/a)[-1])
  trSigma2 <- (sum(1/(d^2 + c)) + s*aux7 + max(p - n, 0)/c)/a 
  trSigma3 <- .aux.var.int(c, a, d, g, h, sumy, n, p)$trSigma
  
  # trXtXSigma
  trXtXSigma1 <- sum(diag(t(cbind(1, x)) %*% cbind(1, x) %*% 
                            solve(t(cbind(1, x)) %*% cbind(1, x) + 
                                    diag(c(0, rep(c, p))))/a))
  trXtXSigma2 <- (sum(d^2/(d^2 + c)) + s*(n - 2*aux1 + aux4))/a
  trXtXSigma3 <- .aux.var.int(c, a, d, g, h, sumy, n, p)$trXtXSigma
  
  # mutmu (without intercept term)
  mutmu1 <- as.numeric((t(y) %*% cbind(1, x) %*% solve(t(cbind(1, x)) %*% cbind(
    1, x) + diag(c(0, rep(c, p)))))[, -1] %*% 
      (solve(t(cbind(1, x)) %*% cbind(1, x) + diag(c(0, rep(c, p)))) %*% 
         t(cbind(1, x)) %*% y)[-1, ])
  mutmu2 <- s^2*(sumy^2 - 2*sumy*aux3 + aux3^2)*sum(h^2*d^2/(d^2 + c)^2) +
    2*s*(aux3 - sumy)*sum(h*g*d^2/(d^2 + c)^2) + sum(g^2*d^2/(d^2 + c)^2)
  mutmu3 <- .aux.var.int(c, a, d, g, h, sumy, n, p)$mutmu
  
  # mutxtxmu
  mutXtXmu1 <- as.numeric(t(y) %*% cbind(1, x) %*% solve(t(cbind(1, x)) %*% cbind(
    1, x) + diag(c(0, rep(c, p)))) %*% t(cbind(1, x)) %*% cbind(1, x) %*% 
      solve(t(cbind(1, x)) %*% cbind(1, x) + diag(c(0, rep(c, p)))) %*% 
      t(cbind(1, x)) %*% y)
  mutXtXmu2 <- n*sumy^2*s^2 + sumy*s*aux3*(2 - 2*n*s) - 2*sumy^2*s^2*aux1 + 
    s*aux3^2*(n*s - 2) + 4*sumy*s^2*aux1*aux3 - 2*s^2*aux1*aux3^2 +
    sumy^2*s^2*aux4 - 2*sumy*s*aux6 - 2*sumy*s^2*aux3*aux4 + 2*s*aux3*aux6 + 
    s^2*aux3^2*aux4 + aux5
  mutXtXmu3 <- .aux.var.int(c, a, d, g, h, sumy, n, p)$mutXtXmu
  
  # ytXmu
  ytXmu1 <- as.numeric(t(y) %*% cbind(1, x) %*% solve(t(cbind(1, x)) %*% cbind(
    1, x) + diag(c(0, rep(c, p)))) %*% t(cbind(1, x)) %*% y)
  ytXmu2 <- sumy^2*s - 2*sumy*s*aux3 + aux2 + s*aux3^2
  ytXmu3 <- .aux.var.int(c, a, d, g, h, sumy, n, p)$ytXmu
  
  # logdetSigma
  logdetSigma1 <- as.numeric(determinant(solve(t(cbind(1, x)) %*% cbind(1, x) + 
                                                 diag(c(0, rep(c, p))))/a, log=TRUE)$modulus)
  logdetSigma2 <- log(s) - (p + 1)*log(a) - sum(log(d^2 + c)) - max(p - n, 0)*log(c)
  logdetSigma3 <- .aux.var.int(c, a, d, g, h, sumy, n, p)$logdetSigma
  
  # mu
  mu1 <- as.numeric(solve(t(cbind(1, x)) %*% cbind(1, x) + diag(
    c(0, rep(c, p)))) %*% t(cbind(1, x)) %*% y)
  mu2 <- c(s*(sumy - aux3), s*(aux3 - sumy)*rowSums(v %*% (t(u)*d/(d^2 + c))) +
             (v %*% (t(u)*d/(d^2 + c))) %*% y)
  
  # dSigma
  dSigma1 <- diag(solve(t(cbind(1, x)) %*% cbind(1, x) + diag(c(0, rep(c, p)))))/a
  dSigma2 <- c(s, colSums((t(svd(x, nv=max(n, p))$v)/
                             sqrt(c(d, rep(0, max(p - n, 0)))^2 + c))^2) +
                 s*rowSums(v %*% (t(u)*d/(d^2 + c)))^2)/a
  
  # Sigma
  Sigma1 <- solve(t(cbind(1, x)) %*% cbind(1, x) + diag(c(0, rep(c, p))))/a
  Sigma2 <- matrix(NA, nrow=p + 1, ncol=p + 1)
  Sigma2[1, 1] <- s/a
  Sigma2[1, 2:(p + 1)] <- Sigma2[2:(p + 1), 1] <- 
    -s*as.numeric(h %*% (t(v)*(d/(d^2 + c))))/a
  Sigma2[2:(p + 1), 2:(p + 1)] <- # first part is correct
    (svd(x, nu=n, nv=max(p, n))$v %*% 
       (t(svd(x, nu=n, nv=max(p, n))$v)/(c(d^2, rep(0, max(p - n, 0))) + c)) + 
       s*t(h %*% (t(v)*d/(d^2 + c))) %*% 
       (h %*% (t(v)*d/(d^2 + c))))/a
  
  all.equal(trSigma1, trSigma2, trSigma3)
  all.equal(trXtXSigma1, trXtXSigma2, trXtXSigma3)
  all.equal(mutmu1, mutmu2, mutmu3)
  all.equal(mutXtXmu1, mutXtXmu2, mutXtXmu3)
  all.equal(ytXmu1, ytXmu2, ytXmu3)
  all.equal(logdetSigma1, logdetSigma2, logdetSigma3)
  all.equal(mu1, mu2)
  all.equal(dSigma1, dSigma2)
  all.equal(Sigma1, Sigma2)
})


test_that("high dim. auxiliary variables calculation", {
  p <- 20; n <- 10; aold <- rchisq(1, 1); cold <- rchisq(1, 3)
  x <- matrix(rnorm(p*n), nrow=n, ncol=p); y <- matrix(rnorm(n))
  svd.x <- svd(x); sv <- svd.x$d; uty <- t(svd.x$u) %*% y
  
  # implemented function
  aux <- .aux.var(cold, aold, sv, uty, n, p)
  
  # the actual calculation are this (see supplement):
  Sigma <- solve(t(x) %*% x + cold*diag(p))/aold
  mu <- solve(t(x) %*% x + cold*diag(p)) %*% t(x) %*% y
  trSigma <- sum(diag(Sigma))
  trXtXSigma <- sum(diag(t(x) %*% x %*% Sigma))
  mutmu <- as.numeric(t(mu) %*% mu)
  mutXtXmu <- as.numeric(t(mu) %*% t(x) %*% x %*% mu)
  logdetSigma <- as.numeric(determinant(Sigma)$modulus)
  ytXmu <- as.numeric(t(y) %*% x %*% mu)
  out <- list(trSigma=trSigma, trXtXSigma=trXtXSigma, mutmu=mutmu,
              mutXtXmu=mutXtXmu, logdetSigma=logdetSigma, ytXmu=ytXmu)
  
  # checks
  expect_equal(aux, out)
  
})

test_that("low dim. auxiliary variables calculation", {
  p <- 10; n <- 20; aold <- rchisq(1, 1); cold <- rchisq(1, 3)
  x <- matrix(rnorm(p*n), nrow=n, ncol=p); y <- matrix(rnorm(n))
  svd.x <- svd(x); sv <- svd.x$d; uty <- t(svd.x$u) %*% y
  
  # implemented function
  aux <- .aux.var(cold, aold, sv, uty, n, p)
  
  # the actual calculation are this (see supplement):
  Sigma <- solve(t(x) %*% x + cold*diag(p))/aold
  mu <- solve(t(x) %*% x + cold*diag(p)) %*% t(x) %*% y
  trSigma <- sum(diag(Sigma))
  trXtXSigma <- sum(diag(t(x) %*% x %*% Sigma))
  mutmu <- as.numeric(t(mu) %*% mu)
  mutXtXmu <- as.numeric(t(mu) %*% t(x) %*% x %*% mu)
  logdetSigma <- as.numeric(determinant(Sigma)$modulus)
  ytXmu <- as.numeric(t(y) %*% x %*% mu)
  out <- list(trSigma=trSigma, trXtXSigma=trXtXSigma, mutmu=mutmu,
              mutXtXmu=mutXtXmu, logdetSigma=logdetSigma, ytXmu=ytXmu)
  
  # checks
  expect_equal(aux, out)
  
})

test_that("low dim. inverse Gaussian ELBO calculation", {
  
  p <- 10; n <- 20; df <- p + n; aold <- rchisq(1, 1); cold <- rchisq(1, 3)
  x <- matrix(rnorm(p*n), nrow=n, ncol=p); y <- matrix(rnorm(n))
  svd.x <- svd(x); sv <- svd.x$d; uty <- t(svd.x$u) %*% y
  aux <- .aux.var(cold, aold, sv, uty, n, p); zeta <- rchisq(1, 7)
  delta <- rchisq(1, 5); a <- rchisq(1, 2); b <- rchisq(1, 4)
  c  <- rchisq(1, 1); e <- rchisq(1, 9); lambda <- rchisq(1, 2)
  theta <- rchisq(1, 1); yty <- as.numeric(t(y) %*% y)
  
  # implemented calculation
  new.elbo <- .elbo.inv.gauss(aux, zeta, delta, a, b, c, e, lambda, theta, df, 
                              p, yty)
  
  # the actual calculation are this (see supplement):
  elbo.const <- 0.5*aux$logdetSigma -
    0.5*a*(yty - 2*aux$ytXmu + aux$mutXtXmu + aux$trXtXSigma) +
    0.5*b*(delta - (a + 1 - b/c)*(aux$mutmu + aux$trSigma)) - 
    0.5*(df + 1)*log(zeta) + 0.25*(p + 1)*log(lambda) - 
    0.5*(p + 1)*log(theta) - 0.25*(p + 1)*log(delta) + 
    gsl::bessel_lnKnu(0.5*(p + 1), sqrt(lambda*delta)/theta) + 
    0.5*lambda*e/theta^2
  elbo <- elbo.const + lambda/theta + 0.5*log(lambda) - 0.5*lambda*e/theta^2 -
    0.5*lambda*b
  out <- list(elbo=elbo, elbo.const=elbo.const)
  
  # checks
  expect_equal(new.elbo, out)
  
})

test_that("high dim. inverse Gaussian ELBO calculation", {
  
  p <- 20; n <- 10; df <- p + n; aold <- rchisq(1, 1); cold <- rchisq(1, 3)
  x <- matrix(rnorm(p*n), nrow=n, ncol=p); y <- matrix(rnorm(n))
  svd.x <- svd(x); sv <- svd.x$d; uty <- t(svd.x$u) %*% y
  aux <- .aux.var(cold, aold, sv, uty, n, p); zeta <- rchisq(1, 7)
  delta <- rchisq(1, 5); a <- rchisq(1, 2); b <- rchisq(1, 4)
  c  <- rchisq(1, 1); e <- rchisq(1, 9); lambda <- rchisq(1, 2)
  theta <- rchisq(1, 1); yty <- as.numeric(t(y) %*% y)
  
  # implemented calculation
  new.elbo <- .elbo.inv.gauss(aux, zeta, delta, a, b, c, e, lambda, theta, df, 
                              p, yty)
  
  # the actual calculation are this (see supplement):
  elbo.const <- 0.5*aux$logdetSigma -
    0.5*a*(yty - 2*aux$ytXmu + aux$mutXtXmu + aux$trXtXSigma) +
    0.5*b*(delta - (a + 1 - b/c)*(aux$mutmu + aux$trSigma)) - 
    0.5*(df + 1)*log(zeta) + 0.25*(p + 1)*log(lambda) - 
    0.5*(p + 1)*log(theta) - 0.25*(p + 1)*log(delta) + 
    gsl::bessel_lnKnu(0.5*(p + 1), sqrt(lambda*delta)/theta) + 
    0.5*lambda*e/theta^2
  elbo <- elbo.const + lambda/theta + 0.5*log(lambda) - 0.5*lambda*e/theta^2 -
    0.5*lambda*b
  out <- list(elbo=elbo, elbo.const=elbo.const)
  
  # checks
  expect_equal(new.elbo, out)
  
})

test_that("low dim. inverse Gamma ELBO calculation", {
  
  p <- 10; n <- 20; df <- p + n; aold <- rchisq(1, 1); cold <- rchisq(1, 3)
  x <- matrix(rnorm(p*n), nrow=n, ncol=p); y <- matrix(rnorm(n))
  svd.x <- svd(x); sv <- svd.x$d; uty <- t(svd.x$u) %*% y
  aux <- .aux.var(cold, aold, sv, uty, n, p); zeta <- rchisq(1, 7)
  delta <- rchisq(1, 5); a <- rchisq(1, 2); b <- rchisq(1, 4)
  c  <- rchisq(1, 1); e <- rchisq(1, 9); lambda <- rchisq(1, 2)
  theta <- rchisq(1, 1); yty <- as.numeric(t(y) %*% y); eta <- rchisq(1, 2)
  
  # implemented calculation
  new.elbo <- .elbo.inv.gamma(aux, zeta, delta, a, b, c, e, lambda, eta, df, p, 
                              yty)
  
  # the actual calculation are this (see supplement):
  elbo.const <- 0.5*aux$logdetSigma -
    0.5*a*(yty - 2*aux$ytXmu + aux$mutXtXmu + aux$trXtXSigma) +
    0.5*b*(delta - (a + 1 - b/c)*(aux$mutmu + aux$trSigma)) - 
    0.5*(df + 1)*log(zeta) - 0.25*(p + 1)*log(delta) + 0.5*eta*e + 
    0.5*eta*log(2) + lgamma(0.5*(p + eta))
  elbo <- elbo.const - lgamma(0.5*eta) - 0.5*eta*e - 0.5*eta*log(2) + 
    0.5*eta*log(lambda) - 0.5*lambda*b
  
  out <- list(elbo=elbo, elbo.const=elbo.const)
  
  # checks
  expect_equal(new.elbo, out)
  
})

test_that("high dim. inverse Gamma ELBO calculation", {
  
  p <- 20; n <- 10; df <- p + n; aold <- rchisq(1, 1); cold <- rchisq(1, 3)
  x <- matrix(rnorm(p*n), nrow=n, ncol=p); y <- matrix(rnorm(n))
  svd.x <- svd(x); sv <- svd.x$d; uty <- t(svd.x$u) %*% y
  aux <- .aux.var(cold, aold, sv, uty, n, p); zeta <- rchisq(1, 7)
  delta <- rchisq(1, 5); a <- rchisq(1, 2); b <- rchisq(1, 4)
  c  <- rchisq(1, 1); e <- rchisq(1, 9); lambda <- rchisq(1, 2)
  theta <- rchisq(1, 1); yty <- as.numeric(t(y) %*% y); eta <- rchisq(1, 2)
  
  # implemented calculation
  new.elbo <- .elbo.inv.gamma(aux, zeta, delta, a, b, c, e, lambda, eta, df, p, 
                              yty)
  
  # the actual calculation are this (see supplement):
  elbo.const <- 0.5*aux$logdetSigma -
    0.5*a*(yty - 2*aux$ytXmu + aux$mutXtXmu + aux$trXtXSigma) +
    0.5*b*(delta - (a + 1 - b/c)*(aux$mutmu + aux$trSigma)) - 
    0.5*(df + 1)*log(zeta) - 0.25*(p + 1)*log(delta) + 0.5*eta*e + 
    0.5*eta*log(2) + lgamma(0.5*(p + eta))
  elbo <- elbo.const - lgamma(0.5*eta) - 0.5*eta*e - 0.5*eta*log(2) + 
    0.5*eta*log(lambda) - 0.5*lambda*b
  
  out <- list(elbo=elbo, elbo.const=elbo.const)
  
  # checks
  expect_equal(new.elbo, out)
  
})

test_that("low dim. conjugate inverse Gaussian VB update", {

  # set parameters
  p <- 10; n <- 20; aold <- rchisq(1, 1); bold <- rchisq(1, 3);
  cold <- bold; x <- matrix(rnorm(p*n), nrow=n, ncol=p); y <- matrix(rnorm(n))
  svd.x <- svd(x); sv <- svd.x$d; uty <- t(svd.x$u) %*% y
  yty <- as.numeric(t(y) %*% y); eta <- rchisq(1, 2); lambda <- rchisq(1, 2)
  theta <- rchisq(1, 1); conjugate <- TRUE; hyperprior <- "inv. Gaussian"
  df <- n + p

  new.vb <- .single.vb.update(aold, bold, eta, theta, lambda, sv, n, p, uty,
                              yty, conjugate, hyperprior)

  # the actual calculation are this (see manuscript):
  Sigma <- solve(t(x) %*% x + cold*diag(p))/aold
  mu <- solve(t(x) %*% x + cold*diag(p)) %*% t(x) %*% y
  delta <- c(delta=(bold/cold)*as.numeric(t(mu) %*% mu) + lambda)
  b <- sqrt(lambda/(theta^2*delta))*
    ratio_besselK(sqrt(lambda*delta)/theta, 0.5*(p + eta)) + (p + eta)/delta
  zeta <- 0.5*(yty - 2*t(y) %*% x %*% mu + sum(diag(t(x) %*% x %*% Sigma)) +
                 t(mu) %*% t(x) %*% x %*% mu + 
                 b*(sum(diag(Sigma)) + t(mu) %*% mu))
  a <- 0.5*(df + 1)/zeta
  e <- unname((b - (p + eta)/delta)*delta*theta^2/lambda)

  # checks
  expect_equal(new.vb["delta"], delta)
  expect_equal(new.vb["zeta"], c(zeta=zeta))
  expect_equal(new.vb["a"], c(a=a))
  expect_equal(new.vb["e"], c(e=e))
})

test_that("high dim. conjugate inverse Gaussian VB update", {
  
  # set parameters
  p <- 20; n <- 10; aold <- rchisq(1, 1); bold <- rchisq(1, 3);
  cold <- bold; x <- matrix(rnorm(p*n), nrow=n, ncol=p); y <- matrix(rnorm(n))
  svd.x <- svd(x); sv <- svd.x$d; uty <- t(svd.x$u) %*% y
  yty <- as.numeric(t(y) %*% y); eta <- rchisq(1, 2); lambda <- rchisq(1, 2)
  theta <- rchisq(1, 1); conjugate <- TRUE; hyperprior <- "inv. Gaussian"
  df <- n + p
  
  new.vb <- .single.vb.update(aold, bold, eta, theta, lambda, sv, n, p, uty,
                              yty, conjugate, hyperprior)
  
  # the actual calculation are this (see manuscript):
  Sigma <- solve(t(x) %*% x + cold*diag(p))/aold
  mu <- solve(t(x) %*% x + cold*diag(p)) %*% t(x) %*% y
  delta <- c(delta=(bold/cold)*as.numeric(t(mu) %*% mu) + lambda)
  b <- sqrt(lambda/(theta^2*delta))*
    ratio_besselK(sqrt(lambda*delta)/theta, 0.5*(p + eta)) + (p + eta)/delta
  zeta <- 0.5*(yty - 2*t(y) %*% x %*% mu + sum(diag(t(x) %*% x %*% Sigma)) +
                 t(mu) %*% t(x) %*% x %*% mu + 
                 b*(sum(diag(Sigma)) + t(mu) %*% mu))
  a <- 0.5*(df + 1)/zeta
  e <- unname((b - (p + eta)/delta)*delta*theta^2/lambda)
  
  # checks
  expect_equal(new.vb["delta"], delta)
  expect_equal(new.vb["zeta"], c(zeta=zeta))
  expect_equal(new.vb["a"], c(a=a))
  expect_equal(new.vb["e"], c(e=e))
})

test_that("low dim. conjugate inverse Gamma VB update", {
  
  # set parameters
  p <- 10; n <- 20; aold <- rchisq(1, 1); bold <- rchisq(1, 3);
  cold <- bold; x <- matrix(rnorm(p*n), nrow=n, ncol=p); y <- matrix(rnorm(n))
  svd.x <- svd(x); sv <- svd.x$d; uty <- t(svd.x$u) %*% y
  yty <- as.numeric(t(y) %*% y); eta <- rchisq(1, 2); lambda <- rchisq(1, 2)
  theta <- rchisq(1, 1); conjugate <- TRUE; hyperprior <- "inv. Gamma"
  df <- n + p
  
  new.vb <- .single.vb.update(aold, bold, eta, theta, lambda, sv, n, p, uty,
                              yty, conjugate, hyperprior)
  
  # the actual calculation are this (see manuscript):
  Sigma <- solve(t(x) %*% x + cold*diag(p))/aold
  mu <- solve(t(x) %*% x + cold*diag(p)) %*% t(x) %*% y
  delta <- c(delta=(bold/cold)*as.numeric(t(mu) %*% mu) + lambda)
  b <- (p + eta)/delta
  zeta <- 0.5*(yty - 2*t(y) %*% x %*% mu + sum(diag(t(x) %*% x %*% Sigma)) +
                 t(mu) %*% t(x) %*% x %*% mu + 
                 b*(sum(diag(Sigma)) + t(mu) %*% mu))
  a <- 0.5*(df + 1)/zeta
  e <- unname(log(delta) - digamma(0.5*(p + eta)) - log(2))
  
  # checks
  expect_equal(new.vb["delta"], delta)
  expect_equal(new.vb["zeta"], c(zeta=zeta))
  expect_equal(new.vb["a"], c(a=a))
  expect_equal(new.vb["e"], c(e=e))
})

test_that("high dim. conjugate inverse Gamma VB update", {
  
  # set parameters
  p <- 20; n <- 10; aold <- rchisq(1, 1); bold <- rchisq(1, 3);
  cold <- bold; x <- matrix(rnorm(p*n), nrow=n, ncol=p); y <- matrix(rnorm(n))
  svd.x <- svd(x); sv <- svd.x$d; uty <- t(svd.x$u) %*% y
  yty <- as.numeric(t(y) %*% y); eta <- rchisq(1, 2); lambda <- rchisq(1, 2)
  theta <- rchisq(1, 1); conjugate <- TRUE; hyperprior <- "inv. Gamma"
  df <- n + p
  
  new.vb <- .single.vb.update(aold, bold, eta, theta, lambda, sv, n, p, uty,
                              yty, conjugate, hyperprior)
  
  # the actual calculation are this (see manuscript):
  Sigma <- solve(t(x) %*% x + cold*diag(p))/aold
  mu <- solve(t(x) %*% x + cold*diag(p)) %*% t(x) %*% y
  delta <- c(delta=(bold/cold)*as.numeric(t(mu) %*% mu) + lambda)
  b <- (p + eta)/delta
  zeta <- 0.5*(yty - 2*t(y) %*% x %*% mu + sum(diag(t(x) %*% x %*% Sigma)) +
                 t(mu) %*% t(x) %*% x %*% mu + 
                 b*(sum(diag(Sigma)) + t(mu) %*% mu))
  a <- 0.5*(df + 1)/zeta
  e <- unname(log(delta) - digamma(0.5*(p + eta)) - log(2))
  
  # checks
  expect_equal(new.vb["delta"], delta)
  expect_equal(new.vb["zeta"], c(zeta=zeta))
  expect_equal(new.vb["a"], c(a=a))
  expect_equal(new.vb["e"], c(e=e))
})

################################## EB testing ##################################
test_that("one inverse Gamma estimating equation evaluation", {
  
  # setting parameters
  eta <- rchisq(1, 10); lambda <- rchisq(1, 5); esum <- rchisq(1, 7)
  size <- 1 + rpois(1, 100)
  
  # implemented function
  eval.froot <- .froot.inv.gamma(eta, lambda, esum, size)
  
  # actual calculations (see supplement)
  eval.froot.test <- log(lambda) - log(2) - digamma(eta/2) - esum/size
  
  # checks
  expect_equal(eval.froot, eval.froot.test)
  
})

test_that("inverse Gamma EB update", {
  
  # setting parameters
  nclass <- 10; epsilon <- 1e-3; maxit <- 100; elbo.const <- 10; 
  sclass <- rpois(nclass, 100) + 1; esum <- rchisq(nclass, 19); 
  bsum <- exp(esum/sclass)*sclass + rchisq(nclass, 2); # true by Jensen's ineq.
  D <- sum(sclass); p <- 100; b <- rchisq(D, 1); old.eta <- rchisq(nclass, 3)
  delta <- rchisq(D, 8)
  
  # implemented function
  new.eb <- .eb.update.inv.gamma(esum, bsum, b, delta, old.eta, elbo.const, 
                                 nclass, sclass, p, D, epsilon, maxit)
  
  # actual calculations (see supplement)
  old.eta <- eta.star <- 1/(log(bsum) + esum/sclass - log(sclass))
  old.lambda <- old.eta*sclass/bsum
  conv <- FALSE
  iter <- 0
  while(!(conv | iter >= maxit)) {
    iter <- iter + 1
    eta <- sapply(c(1:nclass), function(c) {
      uniroot(.froot.inv.gamma, c(eta.star[c], 2*eta.star[c]), 
              old.lambda[c], esum[c], sclass[c])$root})
    lambda <- eta*sclass/bsum
    conv <- all(abs(c(eta - old.eta, lambda - old.lambda)) < epsilon)
    old.lambda <- lambda
    old.eta <- eta
  }
  elbo <- elbo.const - 0.5*(p + rep(eta, sclass))*log(delta) - 
    rep(lgamma(0.5*eta), sclass) + 
    0.5*rep(eta*digamma(0.5*(p + old.eta)), sclass) + 
    0.5*rep(eta*log(lambda), sclass) - 0.5*rep(lambda, sclass)*b
  
  # checks
  expect_equal(new.eb["eta"], list(eta=rep(eta, times=sclass)))
  expect_equal(new.eb["lambda"], list(lambda=rep(lambda, times=sclass)))
  expect_equal(new.eb["conv"], list(conv=conv))
  expect_equal(new.eb["elbo"], list(elbo=elbo))
  expect_equal(all(new.eb$lambda > 0), TRUE)
  expect_equal(all(new.eb$eta > 0), TRUE)
  
})

test_that("one multiple lambda inverse Gaussian EB update", {
  
  # setting parameters
  D <- 100; dcov <- 5; epsilon <- 1e-3
  C <- matrix(rchisq(dcov*D, 1), nrow=D, ncol=dcov); e <- rchisq(D, 1)
  b <- 1/e + rchisq(D, 1) # b > 1/e through Jensen's inequality
  old.lambda <- rchisq(D, 1); old.alpha <- rchisq(dcov, 1)
  
  # implemented function
  new.eb <- .eb.update.mult.lambda(old.lambda, old.alpha, e, b, C, epsilon)
  
  # actual calculations (see supplement)
  alpha <- as.numeric(solve(t(C) %*% diag(old.lambda*e) %*% C) %*% t(C) %*% 
                        diag(old.lambda) %*% matrix(1, nrow=D, ncol=1))
  lambda <- 1/(b + e*as.numeric(C %*% alpha)^2 - 2*as.numeric(C %*% alpha))
  conv <- all(abs(c(alpha - old.alpha, lambda - old.lambda)) < epsilon)
  
  # checks
  expect_equal(new.eb["alpha"], list(alpha=alpha))
  expect_equal(new.eb["lambda"], list(lambda=lambda))
  expect_equal(new.eb["conv"], list(conv=conv))
  expect_equal(all(new.eb$lambda > 0), TRUE)

})

test_that("single lambda inverse Gaussian EB update", {
  
  # setting parameters
  D <- 100; dcov <- 5; epsilon <- 1e-3; maxit=100; elbo.const <- 10
  C <- matrix(rchisq(dcov*D, 1), nrow=D, ncol=dcov); e <- rchisq(D, 1)
  b <- 1/e + rchisq(D, 1) # b > 1/e through Jensen's inequality
  
  # implemented function
  new.eb <- .eb.update.inv.gauss(e, b, elbo.const, C, D, mult.lambda=FALSE, 
                                 epsilon, maxit)
  
  # actual calculations (see manuscript)
  alpha <- as.numeric(solve(t(C) %*% diag(e) %*% C) %*% t(C) %*% 
                        matrix(1, nrow=D, ncol=1))
  theta <- 1/as.numeric(C %*% alpha)
  lambda <- D/as.numeric(sum(b) - t(alpha) %*% t(C) %*% 
                           matrix(1, nrow=D, ncol=1))
  elbo <- elbo.const + lambda/theta + 0.5*log(lambda) - 0.5*lambda*e/theta^2 -
    0.5*lambda*b
  
  # checks
  expect_equal(new.eb["alpha"], list(alpha=alpha))
  expect_equal(new.eb["theta"], list(theta=theta))
  expect_equal(new.eb["lambda"], list(lambda=lambda))
  expect_equal(new.eb["elbo"], list(elbo=elbo))
  expect_equal(all(new.eb$lambda > 0), TRUE)
  
})

test_that("multiple lambda inverse Gaussian EB update", {
  
  # setting parameters
  D <- 100; dcov <- 5; epsilon <- 1e-3; maxit=100; elbo.const <- 10
  mult.lambda <- TRUE; C <- matrix(rchisq(dcov*D, 1), nrow=D, ncol=dcov)
  e <- rchisq(D, 1)
  b <- 1/e + rchisq(D, 1) # b > 1/e through Jensen's inequality
  
  # implemented function
  new.eb <- .eb.update.inv.gauss(e, b, elbo.const, C, D, mult.lambda, epsilon, 
                                 maxit)
  
  # actual calculations (see supplement)
  old.alpha <- as.numeric(solve(t(C) %*% diag(e) %*% C) %*% t(C) %*% 
                            matrix(1, nrow=D, ncol=1))
  old.lambda <- rep(D/as.numeric(sum(b) - t(old.alpha) %*% t(C) %*% 
                                   matrix(1, nrow=D, ncol=1)), D)
  conv <- FALSE
  iter <- 0
  while(!(conv | iter >= maxit)) {
    iter <- iter + 1
    new.eb.test <- .eb.update.mult.lambda(old.lambda, old.alpha, e, b, C, 
                                          epsilon)
    alpha <- new.eb.test$alpha
    lambda <- new.eb.test$lambda
    conv <- new.eb.test$conv
    old.lambda <- lambda
    old.alpha <- alpha
  }
  theta <- 1/as.numeric(C %*% alpha)
  elbo <- elbo.const + lambda/theta + 0.5*log(lambda) - 0.5*lambda*e/theta^2 -
    0.5*lambda*b
  
  # checks
  expect_equal(new.eb["alpha"], list(alpha=alpha))
  expect_equal(new.eb["theta"], list(theta=theta))
  expect_equal(new.eb["lambda"], list(lambda=lambda))
  expect_equal(new.eb["elbo"], list(elbo=elbo))
  expect_equal(all(new.eb$lambda > 0), TRUE)
  
})