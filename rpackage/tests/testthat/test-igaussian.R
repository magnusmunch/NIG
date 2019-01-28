context("inverse Gaussian model estimation")

################################## VB testing ##################################
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