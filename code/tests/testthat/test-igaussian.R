context("inverse Gaussian model estimation")

test_that("dependent beta and sigma^2 caluclations", {
  
  n <- 100; p <- 200; set.seed(123); x <- matrix(rnorm(n*p), nrow=n, ncol=p)
  y <- rnorm(n); lambda <- 1; theta <- 1; aold=2.2; zetaold <- 1
  svd.x <- svd(x); sv <- svd.x$d; uty <- t(svd.x$u) %*% y; yty <- sum(y^2)
  new.vb <- single.vb.update(zetaold, aold, lambda, theta, sv, n, p, uty, yty)
  
  # the actual calculation are this (see manuscript):
  Sigma <- 2*zetaold*solve(t(x) %*% x + aold*diag(p))/(n + p + 1)
  mu <- solve(t(x) %*% x + aold*diag(p)) %*% t(x) %*% y
  delta <- 0.5*(n + p + 1)*(sum(diag(Sigma)) + sum(mu^2))/zetaold + lambda
  a <- sqrt(lambda/(theta^2*delta))*
    ratio_besselK_cpp(sqrt(lambda*delta)/theta, p) + (p + 1)/delta
  zeta <- as.numeric(0.5*a*(sum(y^2) - 2*t(y) %*% x %*% mu +
                              sum(diag((t(x) %*% x + diag(p)) %*% Sigma)) +
                              t(mu) %*% (t(x) %*% x + diag(p)) %*% mu))
  v <- sqrt(delta*theta^2/lambda)*
    ratio_besselK_cpp(sqrt(lambda*delta)/theta, p)
  
  expect_equal(new.vb["delta"], c(delta=delta))
  expect_equal(new.vb["zeta"], c(zeta=zeta))
  expect_equal(new.vb["a"], c(a=a))
  expect_equal(new.vb["v"], c(v=v))
})

test_that("independent beta and sigma^2 caluclations", {
  
  n <- 100; p <- 200; set.seed(123); x <- matrix(rnorm(n*p), nrow=n, ncol=p)
  y <- rnorm(n); lambda <- 1; theta <- 1; aold=2.2; zetaold <- 1
  svd.x <- svd(x); sv <- svd.x$d; uty <- t(svd.x$u) %*% y; yty <- sum(y^2)
  new.vb <- single.vb.update.ind(zetaold, aold, lambda, theta, sv, n, p, uty, 
                                 yty)
  
  # the actual calculation are this (see manuscript):
  Sigma <- solve(0.5*(n + 1)*t(x) %*% x/zetaold + aold*diag(p))
  mu <- 0.5*(n + 1)*Sigma %*% t(x) %*% y/zetaold
  delta <- sum(diag(Sigma)) + sum(mu^2) + lambda
  a <- sqrt(lambda/(theta^2*delta))*
    ratio_besselK_cpp(sqrt(lambda*delta)/theta, p) + (p + 1)/delta
  zeta <- as.numeric(0.5*(sum(y^2) - 2*t(y) %*% x %*% mu +
                              sum(diag(t(x) %*% x %*% Sigma)) +
                              t(mu) %*% t(x) %*% x %*% mu))
  v <- sqrt(delta*theta^2/lambda)*
    ratio_besselK_cpp(sqrt(lambda*delta)/theta, p)
  
  expect_equal(new.vb["delta"], c(delta=delta))
  expect_equal(new.vb["zeta"], c(zeta=zeta))
  expect_equal(new.vb["a"], c(a=a))
  expect_equal(new.vb["v"], c(v=v))
})

test_that("empirical Bayes update", {
  
  D <- 100; dcov <- 5; set.seed(123)
  C <- matrix(rchisq(dcov*D, 1), nrow=D, ncol=dcov); v <- rchisq(D, 1)
  a <- rchisq(D, 1); elbo.const <- 10
  new.eb <- eb.update(v, a, elbo.const, C, D)
  
  # actual calculations (see manuscript)
  alpha <- as.numeric(solve(t(C) %*% diag(v) %*% C) %*% t(C) %*% 
                        matrix(1, nrow=D, ncol=1))
  theta <- 1/as.numeric(C %*% alpha)
  lambda <- D/as.numeric(sum(a) - t(alpha) %*% t(C) %*% matrix(1, nrow=D, ncol=1))
  
  expect_equal(new.eb["alpha"], list(alpha=alpha))
  expect_equal(new.eb["theta"], list(theta=theta))
  expect_equal(new.eb["lambda"], list(lambda=lambda))
})