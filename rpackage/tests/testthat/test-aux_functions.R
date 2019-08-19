context("miscellaneous functions")
test_that("ratio of modified Bessel functions", {
  
  # setting parameters
  x <- exp(seq(-10, 10, by=0.01))
  nu <- c(1, 10, 100, 1000)
  grid <- expand.grid(x=x, nu=nu)
  
  # implemented function
  test.ratio <- sapply(nu, function(n) {ratio_besselK(x, n)})
  
  # actual calculation (see supplement)
  ratio <- sapply(nu, function(n) {besselK(x, n - 1)/besselK(x, n)})
  
  # check
  expect_equal(test.ratio[!is.nan(ratio) & ratio!=0], 
               ratio[!is.nan(ratio) & ratio!=0])
  
})

context("aux variables without unpenalized covariates")
test_that("low dim. auxiliary variables calculation", {
  
  # setting parameters
  n <- 100
  p <- 10
  x <- matrix(rnorm(n*p), nrow=n, ncol=p)
  y <- rnorm(n)
  ytx <- as.numeric(t(y) %*% x)
  aold <- rchisq(1, 1) + 1
  bold <- rchisq(p, 1) + 1
  
  # implemented function
  test.aux <- .aux.var(aold, bold, y, x, ytx)
  
  # actual calculation (see supplement)
  H <- diag(bold)
  Sigma <- solve(t(x) %*% x + H)/aold
  mu <- as.numeric(solve(t(x) %*% x + H) %*% t(x) %*% y)
  aux <- list(mu=mu, dSigma=diag(Sigma), ytXmu=as.numeric(t(y) %*% x %*% mu),
              trXtXSigma=sum(diag(t(x) %*% x %*% Sigma)), 
              mutXtXmu=as.numeric(t(mu) %*% t(x) %*% x %*% mu))
  
  # check
  expect_equal(aux, test.aux)
  
})

test_that("high dim. auxiliary variables calculation", {
  
  # setting parameters
  n <- 10
  p <- 100
  x <- matrix(rnorm(n*p), nrow=n, ncol=p)
  y <- rnorm(n)
  ytx <- as.numeric(t(y) %*% x)
  aold <- rchisq(1, 1) + 1
  bold <- rchisq(p, 1) + 1
  
  # implemented function
  test.aux <- .aux.var(aold, bold, y, x, ytx)
  
  # actual calculation (see supplement)
  H <- diag(bold)
  Sigma <- solve(t(x) %*% x + H)/aold
  mu <- as.numeric(solve(t(x) %*% x + H) %*% t(x) %*% y)
  aux <- list(mu=mu, dSigma=diag(Sigma), ytXmu=as.numeric(t(y) %*% x %*% mu),
              trXtXSigma=sum(diag(t(x) %*% x %*% Sigma)), 
              mutXtXmu=as.numeric(t(mu) %*% t(x) %*% x %*% mu))
  
  # check
  expect_equal(aux, test.aux)
  
})

context("aux variables with unpenalized covariates")
test_that("low dim. auxiliary variables calculation", {
  
  # setting parameters
  n <- 100
  u <- 5
  r <- 10
  xu <- matrix(rnorm(n*u), nrow=n, ncol=u)
  xr <- matrix(rnorm(n*r), nrow=n, ncol=r)
  y <- rnorm(n)
  aold <- rchisq(1, 1) + 1
  bold <- rchisq(r, 1) + 1
  
  # implemented function
  test.aux <- .aux.var.unp(aold, bold, y, xu, xr, u, r)
  
  # actual calculation (see supplement)
  H <- diag(c(rep(0, u), bold))
  x <- cbind(xu, xr)
  Sigma <- solve(t(x) %*% x + H)/aold
  mu <- as.numeric(solve(t(x) %*% x + H) %*% t(x) %*% y)
  aux <- list(mu=mu, dSigma=diag(Sigma), ytXmu=as.numeric(t(y) %*% x %*% mu),
              trXtXSigma=sum(diag(t(x) %*% x %*% Sigma)), 
              mutXtXmu=as.numeric(t(mu) %*% t(x) %*% x %*% mu))
  
  # check
  expect_equal(aux, test.aux)
  
})

test_that("high dim. auxiliary variables calculation", {
  
  # setting parameters
  n <- 10
  u <- 5
  r <- 100
  xu <- matrix(rnorm(n*u), nrow=n, ncol=u)
  xr <- matrix(rnorm(n*r), nrow=n, ncol=r)
  y <- rnorm(n)
  aold <- rchisq(1, 1) + 1
  bold <- rchisq(r, 1) + 1
  
  # implemented function
  test.aux <- .aux.var.unp(aold, bold, y, xu, xr, u, r)
  
  # actual calculation (see supplement)
  H <- diag(c(rep(0, u), bold))
  x <- cbind(xu, xr)
  Sigma <- solve(t(x) %*% x + H)/aold
  mu <- as.numeric(solve(t(x) %*% x + H) %*% t(x) %*% y)
  aux <- list(mu=mu, dSigma=diag(Sigma), ytXmu=as.numeric(t(y) %*% x %*% mu),
              trXtXSigma=sum(diag(t(x) %*% x %*% Sigma)), 
              mutXtXmu=as.numeric(t(mu) %*% t(x) %*% x %*% mu))
  
  # check
  expect_equal(aux, test.aux)
  
})

context("full posterior covariance without unpenalized covariates")
test_that("low dim. full posterior calculation", {
  
  # setting parameters
  n <- 100
  p <- 10
  x <- matrix(rnorm(n*p), nrow=n, ncol=p)
  aold <- rchisq(1, 1) + 1
  bold <- rchisq(p, 1) + 1
  
  # implemented function
  test.Sigma <- .Sigma(aold, bold, x)
  
  # actual calculation (see supplement)
  Sigma <- solve(t(x) %*% x + diag(bold))/aold
  
  # check
  expect_equal(Sigma, test.Sigma)
  
})

test_that("high dim. full posterior calculation", {
  
  # setting parameters
  n <- 10
  p <- 100
  x <- matrix(rnorm(n*p), nrow=n, ncol=p)
  aold <- rchisq(1, 1) + 1
  bold <- rchisq(p, 1) + 1
  
  # implemented function
  test.Sigma <- .Sigma(aold, bold, x)
  
  # actual calculation (see supplement)
  Sigma <- solve(t(x) %*% x + diag(bold))/aold
  
  # check
  expect_equal(Sigma, test.Sigma)
  
})

context("full posterior covariance with unpenalized covariates")
test_that("low dim. full posterior calculation", {
  
  # setting parameters
  n <- 100
  u <- 5
  r <- 10
  xu <- matrix(rnorm(n*u), nrow=n, ncol=u)
  xr <- matrix(rnorm(n*r), nrow=n, ncol=r)
  aold <- rchisq(1, 1) + 1
  bold <- rchisq(r, 1) + 1
  
  # implemented function
  test.Sigma <- .Sigma.unp(aold, bold, xu, xr, u, r)
  
  # actual calculation (see supplement)
  x <- cbind(xu, xr)
  H <- diag(c(rep(0, u), bold))
  Sigma <- solve(t(x) %*% x + H)/aold
  
  # check
  expect_equal(Sigma, test.Sigma)
  
})

test_that("high dim. full posterior calculation", {
  
  # setting parameters
  n <- 10
  u <- 5
  r <- 100
  xu <- matrix(rnorm(n*u), nrow=n, ncol=u)
  xr <- matrix(rnorm(n*r), nrow=n, ncol=r)
  aold <- rchisq(1, 1) + 1
  bold <- rchisq(r, 1) + 1
  
  # implemented function
  test.Sigma <- .Sigma.unp(aold, bold, xu, xr, u, r)
  
  # actual calculation (see supplement)
  x <- cbind(xu, xr)
  H <- diag(c(rep(0, u), bold))
  Sigma <- solve(t(x) %*% x + H)/aold
  
  # check
  expect_equal(Sigma, test.Sigma)
  
})