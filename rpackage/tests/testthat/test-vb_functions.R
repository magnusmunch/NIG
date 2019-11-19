context("vb update without unpenalized covariates")
test_that("low dim. vb update", {
  
  # setting parameters
  n <- 100
  p <- 10
  x <- matrix(rnorm(p*n), nrow=n, ncol=p)
  y <- rnorm(n)
  aold <- rchisq(1, 1) + 1
  bold <- rchisq(p, 1) + 1
  gold <- rchisq(1, 1) + 1
  Calphaf <- rchisq(p, 1) + 1
  Zalphad <- rchisq(1, 1) + 1
  lambdaf <- rchisq(1, 1) + 1
  lambdad <- rchisq(1, 1) + 1
  ytx <- as.numeric(t(y) %*% x)
  yty <- as.numeric(t(y) %*% y)
  
  # implemented function
  test.vb <- .single.vb.update(aold, bold, gold, Calphaf, Zalphad, lambdaf, 
                               lambdad, y, x, ytx, yty, n, p)
  
  # actual calculation (see supplement)
  Sigma <- solve(t(x) %*% x + gold*diag(bold))/aold
  mu <- as.numeric(solve(t(x) %*% x + gold*diag(bold)) %*% t(x) %*% y)
  delta <- aold*gold*(mu^2 + diag(Sigma)) + lambdaf
  b <- sqrt(lambdaf*Calphaf^2/delta)*
    ratio_besselK(sqrt(lambdaf*delta)*Calphaf, 1) + 2/delta
  eta <- aold*sum(bold*(mu^2 + diag(Sigma))) + lambdad 
  g <- sqrt(lambdad*Zalphad^2/eta)*
    ratio_besselK(sqrt(lambdad*eta)*Zalphad, (p + 1)/2) + (p + 1)/eta
  zeta <- 0.5*as.numeric(
    t(y) %*% y  - 2*t(y) %*% x %*% mu + sum(diag(t(x) %*% x %*% Sigma)) +
      t(mu) %*% t(x) %*% x %*% mu + g*sum(diag(diag(b) %*% Sigma)) + 
      g*t(mu) %*% diag(b) %*% mu)
  a <- (n + p + 1)/(2*zeta)
  e <- (b - 2/delta)*delta/(lambdaf*Calphaf^2)
  f <- (g - (p + 1)/eta)*eta/(lambdad*Zalphad^2)
  aux <- .aux.var(aold, gold*bold, y, x, ytx)
  vb <- list(delta=delta, eta=eta, zeta=zeta, a=a, b=b, g=g, e=e, f=f,
             aux=aux)
  
  # check
  expect_equal(vb, test.vb)
  
})

test_that("high dim. vb update", {
  
  # setting parameters
  n <- 10
  p <- 100
  x <- matrix(rnorm(p*n), nrow=n, ncol=p)
  y <- rnorm(n)
  aold <- rchisq(1, 1) + 1
  bold <- rchisq(p, 1) + 1
  gold <- rchisq(1, 1) + 1
  Calphaf <- rchisq(p, 1) + 1
  Zalphad <- rchisq(1, 1) + 1
  lambdaf <- rchisq(1, 1) + 1
  lambdad <- rchisq(1, 1) + 1
  ytx <- as.numeric(t(y) %*% x)
  yty <- as.numeric(t(y) %*% y)
  
  # implemented function
  test.vb <- .single.vb.update(aold, bold, gold, Calphaf, Zalphad, lambdaf, 
                               lambdad, y, x, ytx, yty, n, p)
  
  # actual calculation (see supplement)
  Sigma <- solve(t(x) %*% x + gold*diag(bold))/aold
  mu <- as.numeric(solve(t(x) %*% x + gold*diag(bold)) %*% t(x) %*% y)
  delta <- aold*gold*(mu^2 + diag(Sigma)) + lambdaf
  b <- sqrt(lambdaf*Calphaf^2/delta)*
    ratio_besselK(sqrt(lambdaf*delta)*Calphaf, 1) + 2/delta
  eta <- aold*sum(bold*(mu^2 + diag(Sigma))) + lambdad 
  g <- sqrt(lambdad*Zalphad^2/eta)*
    ratio_besselK(sqrt(lambdad*eta)*Zalphad, (p + 1)/2) + (p + 1)/eta
  zeta <- 0.5*as.numeric(
    t(y) %*% y  - 2*t(y) %*% x %*% mu + sum(diag(t(x) %*% x %*% Sigma)) +
      t(mu) %*% t(x) %*% x %*% mu + g*sum(diag(diag(b) %*% Sigma)) + 
      g*t(mu) %*% diag(b) %*% mu)
  a <- (n + p + 1)/(2*zeta)
  e <- (b - 2/delta)*delta/(lambdaf*Calphaf^2)
  f <- (g - (p + 1)/eta)*eta/(lambdad*Zalphad^2)
  aux <- .aux.var(aold, gold*bold, y, x, ytx)
  vb <- list(delta=delta, eta=eta, zeta=zeta, a=a, b=b, g=g, e=e, f=f,
             aux=aux)
  
  # check
  expect_equal(vb, test.vb)
  
})

context("vb update with unpenalized covariates")
test_that("low dim. vb update", {
  
  # setting parameters
  n <- 100
  u <- 5
  r <- 10
  xu <- matrix(rnorm(n*u), nrow=n, ncol=u)
  xr <- matrix(rnorm(n*r), nrow=n, ncol=r)
  y <- rnorm(n)
  aold <- rchisq(1, 1) + 1
  bold <- rchisq(r, 1) + 1
  gold <- rchisq(1, 1) + 1
  Calphaf <- rchisq(r, 1) + 1
  Zalphad <- rchisq(1, 1) + 1
  lambdaf <- rchisq(1, 1) + 1
  lambdad <- rchisq(1, 1) + 1
  yty <- as.numeric(t(y) %*% y)
  
  # implemented function
  test.vb <- .single.vb.update.unp(aold, bold, gold, Calphaf, Zalphad, lambdaf, 
                                   lambdad, y, xu, xr, yty, n, u, r)
  
  # actual calculation (see supplement)
  hold <- c(rep(0, u), bold)
  x <- cbind(xu, xr)
  Sigma <- solve(t(x) %*% x + gold*diag(hold))/aold
  mu <- as.numeric(solve(t(x) %*% x + gold*diag(hold)) %*% t(x) %*% y)
  delta <- aold*gold*(mu[-c(1:u)]^2 + diag(Sigma)[-c(1:u)]) + lambdaf
  b <- sqrt(lambdaf*Calphaf^2/delta)*
    ratio_besselK(sqrt(lambdaf*delta)*Calphaf, 1) + 2/delta
  h <- c(rep(0, u), b)
  eta <- aold*sum(bold*(mu[-c(1:u)]^2 + diag(Sigma)[-c(1:u)])) + lambdad 
  g <- sqrt(lambdad*Zalphad^2/eta)*
    ratio_besselK(sqrt(lambdad*eta)*Zalphad, (r + 1)/2) + (r + 1)/eta
  zeta <- 0.5*as.numeric(
    t(y) %*% y  - 2*t(y) %*% x %*% mu + sum(diag(t(x) %*% x %*% Sigma)) +
      t(mu) %*% t(x) %*% x %*% mu + g*sum(diag(diag(h) %*% Sigma)) + 
      g*t(mu) %*% diag(h) %*% mu)
  a <- (n + r + u + 1)/(2*zeta)
  e <- (b - 2/delta)*delta/(lambdaf*Calphaf^2)
  f <- (g - (r + 1)/eta)*eta/(lambdad*Zalphad^2)
  aux <- .aux.var.unp(aold, gold*bold, y, xu, xr, u, r)
  vb <- list(delta=delta, eta=eta, zeta=zeta, a=a, b=b, g=g, e=e, f=f,
             aux=aux)
  
  # check
  expect_equal(vb, test.vb)
  
})

test_that("high dim. vb update", {
  
  # setting parameters
  n <- 10
  u <- 5
  r <- 100
  xu <- matrix(rnorm(n*u), nrow=n, ncol=u)
  xr <- matrix(rnorm(n*r), nrow=n, ncol=r)
  y <- rnorm(n)
  aold <- rchisq(1, 1) + 1
  bold <- rchisq(r, 1) + 1
  gold <- rchisq(1, 1) + 1
  Calphaf <- rchisq(r, 1) + 1
  Zalphad <- rchisq(1, 1) + 1
  lambdaf <- rchisq(1, 1) + 1
  lambdad <- rchisq(1, 1) + 1
  yty <- as.numeric(t(y) %*% y)
  
  # implemented function
  test.vb <- .single.vb.update.unp(aold, bold, gold, Calphaf, Zalphad, lambdaf, 
                                   lambdad, y, xu, xr, yty, n, u, r)
  
  # actual calculation (see supplement)
  hold <- c(rep(0, u), bold)
  x <- cbind(xu, xr)
  Sigma <- solve(t(x) %*% x + gold*diag(hold))/aold
  mu <- as.numeric(solve(t(x) %*% x + gold*diag(hold)) %*% t(x) %*% y)
  delta <- aold*gold*(mu[-c(1:u)]^2 + diag(Sigma)[-c(1:u)]) + lambdaf
  b <- sqrt(lambdaf*Calphaf^2/delta)*
    ratio_besselK(sqrt(lambdaf*delta)*Calphaf, 1) + 2/delta
  h <- c(rep(0, u), b)
  eta <- aold*sum(bold*(mu[-c(1:u)]^2 + diag(Sigma)[-c(1:u)])) + lambdad 
  g <- sqrt(lambdad*Zalphad^2/eta)*
    ratio_besselK(sqrt(lambdad*eta)*Zalphad, (r + 1)/2) + (r + 1)/eta
  zeta <- 0.5*as.numeric(
    t(y) %*% y  - 2*t(y) %*% x %*% mu + sum(diag(t(x) %*% x %*% Sigma)) +
      t(mu) %*% t(x) %*% x %*% mu + g*sum(diag(diag(h) %*% Sigma)) + 
      g*t(mu) %*% diag(h) %*% mu)
  a <- (n + r + u + 1)/(2*zeta)
  e <- (b - 2/delta)*delta/(lambdaf*Calphaf^2)
  f <- (g - (r + 1)/eta)*eta/(lambdad*Zalphad^2)
  aux <- .aux.var.unp(aold, gold*bold, y, xu, xr, u, r)
  vb <- list(delta=delta, eta=eta, zeta=zeta, a=a, b=b, g=g, e=e, f=f,
             aux=aux)
  
  # check
  expect_equal(vb, test.vb)
  
})