# update EB parameters (not tested)
.eb.updated <- function(f, g, Z, D, lambdad=NULL) {
  
  # eb parameters
  mat <- t(Z*f) %*% Z
  if(!is.finite(determinant(mat)$modulus)) {
    alphad <- rowSums(solve(mat + sqrt(.Machine$double.eps)*diag(nrow(mat))) %*% 
                        t(Z))
  } else {
    alphad <- tryCatch({
      rowSums(solve(mat) %*% t(Z))}, error=function(e) {
        rowSums(solve(mat + sqrt(.Machine$double.eps)*
                        diag(nrow(mat))) %*% t(Z))})
  }
  
  # if lambdad is fixed we do not estimate
  if(is.null(lambdad)) {
    lambdad <- D/(sum(g) + sum(colSums(t(Z*sqrt(f))*alphad)^2) - 
                    2*sum(t(Z)*alphad))    
  }
  
  
  # compute inner product, mean and variance
  Zalphad <- ifelse(as.numeric(Z %*% alphad) < 0, abs(as.numeric(Z %*% alphad)),
                    as.numeric(Z %*% alphad))
  mpriord <- 1/Zalphad
  vpriord <- 1/(Zalphad^3*lambdad)
  
  out <- list(alphad=alphad, lambdad=lambdad, mpriord=mpriord, vpriord=vpriord, 
              Zalphad=Zalphad)
  return(out)
}

# update EB parameters (not tested)
.eb.updatef <- function(e, b, Cmat, p, D, lambdaf=NULL) {
  
  # eb parameters
  mat <- t(Cmat*Reduce("c", e)) %*% Cmat
  if(!is.finite(determinant(mat)$modulus)) {
    alphaf <- rowSums(solve(mat + sqrt(.Machine$double.eps)*diag(nrow(mat))) %*% 
                        t(Cmat))
  } else {
    alphaf <- tryCatch({
      rowSums(solve(mat) %*% t(Cmat))}, error=function(e) {
        rowSums(solve(mat + sqrt(.Machine$double.eps)*
                        diag(nrow(mat))) %*% t(Cmat))})
  }
  if(is.null(lambdaf)) {
    lambdaf <- sum(p)/
      (sum(Reduce("c", b)) + sum(colSums(t(Cmat*sqrt(Reduce("c", e)))*
                                           alphaf)^2) - 2*sum(t(Cmat)*alphaf))    
  }
  
  # compute inner product, mean and variance
  Calphaf <- unname(split(as.numeric(Cmat %*% alphaf), rep(1:D, times=p)))
  Calphaf <- lapply(Calphaf, function(ca) {ifelse(ca < 0, abs(ca), ca)})
  mpriorf <- lapply(1:D, function(d) {1/Calphaf[[d]]})
  vpriorf <- lapply(1:D, function(d) {
    1/(Calphaf[[d]]^3*lambdaf)})
  
  out <- list(alphaf=alphaf, lambdaf=lambdaf, mpriorf=mpriorf, vpriorf=vpriorf, 
              Calphaf=Calphaf)
  return(out)
}
