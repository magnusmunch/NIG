# single VB update in NIG model with unpenalized covariates
.single.vb.update.unp <- function(aold, bold, gold, Calphaf, Zalphad, lambdaf, 
                                  lambdad, y, xu, xr, yty, n, u, r, var.scale) {
  
  # auxiliary variables involving mu and Sigma
  aux <- .aux.var.unp(aold, gold*bold/var.scale, y, xu, xr, u, r)  
  
  # vb parameters
  delta <- aold*gold*(aux$mu[-c(1:u)]^2 + aux$dSigma[-c(1:u)])/var.scale + 
    lambdaf
  if(is.null(Calphaf)) {
    b <- rep(1, r)  
  } else {
    b <- sqrt(lambdaf*Calphaf^2/delta)*
      ratio_besselK(sqrt(lambdaf*delta)*Calphaf, 1) + 2/delta  
  }
  eta <- aold*sum(bold*(aux$mu[-c(1:u)]^2 + aux$dSigma[-c(1:u)]))/var.scale + 
    lambdad 
  if(is.null(Zalphad)) {
    g <- 1
  } else {
    g <- sqrt(lambdad*Zalphad^2/eta)*
      ratio_besselK(sqrt(lambdad*eta)*Zalphad, (r + 1)/2) + (r + 1)/eta
  }
  zeta <- 0.5*(yty - 2*aux$ytXmu + aux$trXtXSigma + aux$mutXtXmu + 
                 g*sum(aux$dSigma[-c(1:u)]*b)/var.scale + 
                 g*sum(aux$mu[-c(1:u)]^2*b)/var.scale)
  a <- (n + u + r + 1)/(2*zeta)
  e <- (b - 2/delta)*delta/(lambdaf*Calphaf^2)
  f <- (g - (r + 1)/eta)*eta/(lambdad*Zalphad^2)
  
  out <- list(delta=delta, eta=eta, zeta=zeta, a=a, b=b, g=g, e=e, f=f, aux=aux)
  return(out)
}

# single VB update in NIG model
.single.vb.update <- function(aold, bold, gold, Calphaf, Zalphad, lambdaf, 
                              lambdad, y, x, ytx, yty, n, p, var.scale) {
  
  # auxiliary variables involving mu and Sigma
  aux <- .aux.var(aold, gold*bold/var.scale, y, x, ytx)
  
  # vb parameters
  delta <- aold*gold*(aux$mu^2 + aux$dSigma)/var.scale + lambdaf
  if(is.null(Calphaf)) {
    b <- rep(1, p)  
  } else {
    b <- sqrt(lambdaf*Calphaf^2/delta)*
      ratio_besselK(sqrt(lambdaf*delta)*Calphaf, 1) + 2/delta  
  }
  eta <- aold*sum(bold*(aux$mu^2 + aux$dSigma))/var.scale + lambdad 
  if(is.null(Zalphad)) {
    g <- 1
  } else {
    g <- sqrt(lambdad*Zalphad^2/eta)*
      ratio_besselK(sqrt(lambdad*eta)*Zalphad, (p + 1)/2) + (p + 1)/eta
  }
  zeta <- (yty - 2*aux$ytXmu + aux$trXtXSigma + aux$mutXtXmu + 
             sum(aux$dSigma*g*b)/var.scale + sum(aux$mu^2*g*b)/var.scale)/2
  a <- (n + p + 1)/(2*zeta)
  e <- (b - 2/delta)*delta/(lambdaf*Calphaf^2)
  f <- (g - (p + 1)/eta)*eta/(lambdad*Zalphad^2)
  
  out <- list(delta=delta, eta=eta, zeta=zeta, a=a, b=b, g=g, e=e, f=f,
              aux=aux)
  return(out)
}
