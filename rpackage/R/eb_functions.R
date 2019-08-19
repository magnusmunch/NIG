# update EB parameters (not tested)
.eb.update <- function(e, b, old.lambda, Cmat, p, D, mult.lambda, fixed.eb,
                       epsilon, maxit) {
  
  # eb parameters
  mat <- t(Cmat*Reduce("c", e)) %*% Cmat
  if(!is.finite(determinant(mat)$modulus)) {
    alpha <- rowSums(solve(mat + sqrt(.Machine$double.eps)*diag(nrow(mat))) %*% 
                       t(Cmat))
  } else {
    alpha <- tryCatch({rowSums(solve(mat) %*% t(Cmat))},
                      error=function(e) {
                        rowSums(solve(mat + sqrt(.Machine$double.eps)*
                                        diag(nrow(mat))) %*% t(Cmat))})
  }
  if(fixed.eb=="lambda") {
    lambda <- old.lambda
  } else {
    lambda <- sum(p)/(sum(Reduce("c", b)) + 
                        sum(colSums(t(Cmat*sqrt(Reduce("c", e)))*alpha)^2) - 
                        2*sum(t(Cmat)*alpha))  
  }
  
  # IRLS
  if(mult.lambda & fixed.eb=="none") {
    new.eb <-  list(alpha=alpha, lambda=rep(lambda, D), conv=FALSE)
    iter <- 0
    while(!(new.eb$conv | iter >= maxit)) {
      
      # update iteration number and EB parameters
      iter <- iter + 1
      new.eb <- .eb.update.mult(new.eb$lambda, new.eb$alpha, e, b, Cmat, 
                                p, D, epsilon)
    }
    alpha <- new.eb$alpha
    lambda <- new.eb$lambda
    conv <- new.eb$conv
  } else {
    conv <- TRUE
    iter <- 1
  }
  
  # compute inner product, mean and variance
  Calpha <- unname(split(as.numeric(Cmat %*% alpha), rep(1:D, times=p)))
  mprior <- lapply(1:D, function(d) {1/Calpha[[d]]})
  vprior <- lapply(1:D, function(d) {
    1/(Calpha[[d]]^3*ifelse(mult.lambda, lambda[d], lambda))})
  
  out <- list(alpha=alpha, lambda=lambda, mprior=mprior, vprior=vprior, 
              Calpha=Calpha, conv=conv, iter=iter)
  return(out)
}

# one EB IRLS iteration (not tested)
.eb.update.mult <- function(old.lambda, old.alpha, e, b, Cmat, p, D, epsilon) {
  
  # create new IRLS updates
  mat <- t(Cmat*Reduce("c", e)*rep(old.lambda, p)) %*% Cmat
  if(!is.finite(determinant(mat)$modulus)) {
    alpha <- rowSums(solve(mat + sqrt(.Machine$double.eps)*diag(nrow(mat))) %*% 
                       t(Cmat*rep(old.lambda, p)))
  } else {
    alpha <- tryCatch({rowSums(solve(mat) %*% t(Cmat*rep(old.lambda, p)))},
                      error=function(e) {
                        rowSums(solve(mat + sqrt(.Machine$double.eps)*
                                        diag(nrow(mat))) %*% 
                                  t(Cmat*rep(old.lambda, p)))})
  }
  Calpha <- unname(split(as.numeric(Cmat %*% alpha), rep(1:D, times=p)))
  lambda <- p/(sapply(b, sum) + sapply(1:D, function(d) {
    sum((sqrt(e[[d]])*Calpha[[d]])^2)}) - 2*sapply(Calpha, sum))
  
  # check convergence of IRLS
  conv <- all(abs(c(alpha - old.alpha, lambda - old.lambda)) < epsilon)
  
  out <- list(alpha=unname(alpha), lambda=unname(lambda), conv=unname(conv))
  return(out)
}
