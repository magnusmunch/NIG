# single VB update in ENIG model with unpenalized covariates (not tested)
.single.vb.update.unp <- function(aold, bold, Calpha, lambda, y, xu, xr, yty, 
                                  n, u, r) {
  
  # auxiliary variables involving mu and Sigma
  aux <- .aux.var.unp(aold, bold, y, xu, xr, u, r)
  
  # vb parameters
  delta <- aold*(aux$mu[-c(1:u)]^2 + aux$dSigma[-c(1:u)]) + lambda
  b <- sqrt(lambda*Calpha^2/delta)*
    ratio_besselK(sqrt(lambda*delta)*Calpha, 1) + 2/delta
  zeta <- (yty - 2*aux$ytXmu + aux$trXtXSigma + aux$mutXtXmu + 
             sum(aux$dSigma*c(rep(0, u), b)) + sum(aux$mu^2*c(rep(0, u), b)))/2
  a <- (n + u + r + 1)/(2*zeta)
  e <- (b - 2/delta)*delta/(lambda*Calpha^2)
  
  out <- list(delta=delta, zeta=zeta, a=a, b=b, e=e)
  return(out)
}

# single VB update in ENIG model (not tested)
.single.vb.update <- function(aold, bold, Calpha, lambda, y, x, ytx, yty, n, 
                              p) {
  
  # auxiliary variables involving mu and Sigma
  aux <- .aux.var(aold, bold, y, x, ytx)
  
  # vb parameters
  delta <- aold*(aux$mu^2 + aux$dSigma) + lambda
  b <- sqrt(lambda*Calpha^2/delta)*
    ratio_besselK(sqrt(lambda*delta)*Calpha, 1) + 2/delta
  zeta <- (yty - 2*aux$ytXmu + aux$trXtXSigma + aux$mutXtXmu + 
             sum(aux$dSigma*b) + sum(aux$mu^2*b))/2
  a <- (n + p + 1)/(2*zeta)
  e <- (b - 2/delta)*delta/(lambda*Calpha^2)
    
  out <- list(delta=delta, zeta=zeta, a=a, b=b, e=e)
  return(out)
}
