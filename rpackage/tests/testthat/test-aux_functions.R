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

test_that("low dim. aux. variables in general model without intercept", {

  n <- 20
  p <- 10
  y <- rnorm(n)
  x <- matrix(rnorm(n*p), ncol=p, nrow=n)
  aold <- rchisq(1, 1)
  bold <- rchisq(1, 1)
  cold <- rchisq(p, 1)
  ytx <- as.numeric(t(y) %*% x)
  
  # fast calculation
  aux <- .aux.var.gen(aold, bold, cold, x, ytx)

  # actual calculation
  Sigma <- solve(t(x) %*% x + bold*diag(cold))/aold
  mu <- solve(t(x) %*% x + bold*diag(cold)) %*% t(x) %*% y
  mutdiagcmu <- as.numeric(t(mu) %*% diag(cold) %*% mu)
  trdiagcSigma <- sum(diag(diag(cold) %*% Sigma))
  mutmu <- as.numeric(t(mu) %*% mu)
  trSigma <- sum(diag(Sigma))
  ytXmu <- as.numeric(t(y) %*% x %*% mu)
  trXtXSigma <- sum(diag(t(x) %*% x %*% Sigma))
  mutXtXmu <- as.numeric(t(mu) %*% t(x) %*% x %*% mu)

  # check
  expect_equal(aux$mutdiagcmu, mutdiagcmu)
  expect_equal(aux$trdiagcSigma, trdiagcSigma)
  expect_equal(aux$mutmu, mutmu)
  expect_equal(aux$trSigma, trSigma)
  expect_equal(aux$ytXmu, ytXmu)
  expect_equal(aux$trXtXSigma, trXtXSigma)
  expect_equal(aux$mutXtXmu, mutXtXmu)

})

test_that("high dim. aux. variables in general model without intercept", {
  
  n <- 10
  p <- 20
  y <- rnorm(n)
  x <- matrix(rnorm(n*p), ncol=p, nrow=n)
  aold <- rchisq(1, 1)
  bold <- rchisq(1, 1)
  cold <- rchisq(p, 1)
  ytx <- as.numeric(t(y) %*% x)
  
  # fast calculation
  aux <- .aux.var.gen(aold, bold, cold, x, ytx)
  
  # actual calculation
  Sigma <- solve(t(x) %*% x + bold*diag(cold))/aold
  mu <- solve(t(x) %*% x + bold*diag(cold)) %*% t(x) %*% y
  mutdiagcmu <- as.numeric(t(mu) %*% diag(cold) %*% mu)
  trdiagcSigma <- sum(diag(diag(cold) %*% Sigma))
  mutmu <- as.numeric(t(mu) %*% mu)
  trSigma <- sum(diag(Sigma))
  ytXmu <- as.numeric(t(y) %*% x %*% mu)
  trXtXSigma <- sum(diag(t(x) %*% x %*% Sigma))
  mutXtXmu <- as.numeric(t(mu) %*% t(x) %*% x %*% mu)
  
  # check
  expect_equal(aux$mutdiagcmu, mutdiagcmu)
  expect_equal(aux$trdiagcSigma, trdiagcSigma)
  expect_equal(aux$mutmu, mutmu)
  expect_equal(aux$trSigma, trSigma)
  expect_equal(aux$ytXmu, ytXmu)
  expect_equal(aux$trXtXSigma, trXtXSigma)
  expect_equal(aux$mutXtXmu, mutXtXmu)
  
})