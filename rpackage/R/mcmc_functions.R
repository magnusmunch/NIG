.draw_beta <- function(par, x, y, n, p) {
  gamma <- par$gamma
  sigma <- par$sigma
  tau <- par$tau
  u <- rnorm(p, 0, tau*sigma*gamma)
  delta <- rnorm(n, 0, sigma)
  beta <- u + (t(x)*gamma^2) %*% solve(x %*% (t(x)*gamma^2) + diag(n)/tau^2) %*%
    (y - x %*% u - delta)  
  return(as.numeric(beta))
}

.draw_gamma <- function(par, C, p) {
  alphaf <- par$alphaf
  lambdaf <- par$lambdaf
  beta <- par$beta
  sigma <- par$sigma
  tau <- par$tau
  gamma <- sapply(1:p, function(j) {
    chi <- max(beta[j]^2/(sigma*tau)^2 + lambdaf, sqrt(.Machine$double.eps))
    psi <- max(lambdaf*as.numeric(C %*% alphaf)[j]^2, sqrt(.Machine$double.eps))
    sqrt(rgig(n=1, chi=chi, psi=psi, lambda=-1))})
  return(gamma)
}

.draw_tau <- function(par, Z, p) {
  alphad <- par$alphad
  lambdad <- par$lambdad
  beta <- par$beta
  sigma <- par$sigma
  gamma <- par$gamma
  chi <- max(sum(beta^2/(gamma*sigma)^2) + lambdad, sqrt(.Machine$double.eps))
  psi <- max(lambdad*as.numeric(Z %*% alphad)^2, sqrt(.Machine$double.eps))
  tau <- sqrt(rgig(1, chi=chi, psi=psi, lambda=-(p + 1)/2))
  return(tau)
}

.draw_sigma <- function(par, x, y, n, p) {
  beta <- par$beta
  gamma <- par$gamma
  tau <- par$tau
  sigma <- 1/sqrt(rgamma(1, (n + p + 1)/2, 
                         0.5*(sum(y^2) -2*sum(y*as.numeric(x %*% beta)) + 
                                sum(as.numeric(x %*% beta)^2) + 
                                sum(beta^2/gamma^2)/tau^2)))
  return(sigma)
}

.draw_alphaf <- function(par, Cmat, hyper) {
  nu0 <- hyper$nu0[1]
  gamma <- unlist(par$gamma)
  lambdaf <- par$lambdaf
  mat <- t(Cmat*gamma^2) %*% Cmat
  diag(mat) <- diag(mat) + 1/(nu0^2)
  mat <- solve(mat)
  alphaf <- rmvnorm(1, rowSums(mat %*% t(Cmat)), mat/lambdaf)
  return(as.numeric(alphaf))
}

.draw_alphad <- function(par, Zmat, hyper) {
  nu0 <- hyper$nu0[2]
  tau <- par$tau
  lambdad <- par$lambdad
  mat <- t(Zmat*tau^2) %*% Zmat
  diag(mat) <- diag(mat) + 1/(nu0^2)
  mat <- solve(mat)
  alphad <- rmvnorm(1, rowSums(mat %*% t(Zmat)), mat/lambdad)
  return(as.numeric(alphad))
}

.draw_lambdaf <- function(par, Cmat, hyper, p) {
  alpha0 <- hyper$alpha0[1]
  beta0 <- hyper$beta0[1]
  alphaf <- par$alphaf
  gamma <- unlist(par$gamma)
  lambdaf <- rgamma(1, p/2 + alpha0, 
                    0.5*sum(gamma^2*(Cmat %*% alphaf - 1/gamma^2)^2) + beta0)
  return(lambdaf)
}

.draw_lambdad <- function(par, Zmat, hyper, D) {
  alpha0 <- hyper$alpha0[2]
  beta0 <- hyper$beta0[2]
  alphad <- par$alphad
  tau <- par$tau
  lambdad <- rgamma(1, D/2 + alpha0, 
                    0.5*sum(tau^2*(Zmat %*% alphad - 1/tau^2)^2) + beta0)
  return(lambdad)
}

.draw_post <- function(par, x, y, C, Z, hyper, n, p, D) {
  Cmat <- Reduce("rbind", C)
  for(d in 1:D) {
    par$beta[[d]] <- .draw_beta(c(sapply(par[-c(5:8)], "[[", d), par[c(5:8)]), 
                                x=x[[d]], y=y[, d], n=n, p=p[d])
    par$gamma[[d]] <- .draw_gamma(c(sapply(par[-c(5:8)], "[[", d), par[c(5:8)]), 
                                  C[[d]], p[d])
    par$tau[[d]] <- .draw_tau(c(sapply(par[-c(5:8)], "[[", d), par[c(5:8)]), 
                              Z[d, ], p[d])
    par$sigma[[d]] <- .draw_sigma(c(sapply(par[-c(5:8)], "[[", d), par[c(5:8)]), 
                                  x[[d]], y[, d], n, p[d])
  }
  par$alphaf <- .draw_alphaf(par, Cmat, hyper)
  par$alphad <- .draw_alphad(par, Z, hyper)
  par$lambdaf <- .draw_lambdaf(par, Cmat, hyper, sum(p))
  par$lambdad <- .draw_lambdad(par, Z, hyper, D)
  return(par)
}

gibbs.semnig <- function(x, y, C, Z, 
                         hyper=list(nu0=c(Inf, Inf), alpha0=c(0, 0), 
                                    beta0=c(0, 0)), 
                         control=list(niter=200, warmup=floor(niter/2), 
                                      trace=TRUE)) {
  
  if(is.null(control$warmup)) {
    control$warmup <- floor(control$niter/2)
  }
  D <- ncol(y)
  p <- sapply(x, ncol)
  n <- nrow(y)
  G <- ncol(C[[1]])
  H <- ncol(Z)
  
  draws <- list(beta=lapply(1:D, function(d) {matrix(nrow=p[d], 
                                                     ncol=control$niter)}), 
                gamma=lapply(1:D, function(d) {matrix(nrow=p[d], 
                                                      ncol=control$niter)}),
                tau=matrix(NA, nrow=D, ncol=control$niter),
                sigma=matrix(NA, nrow=D, ncol=control$niter),
                alphaf=matrix(NA, nrow=G, ncol=control$niter),
                alphad=matrix(NA, nrow=H, ncol=control$niter),
                lambdaf=numeric(control$niter),
                lambdad=numeric(control$niter))
  
  par <- list(beta=sapply(1:D, function(d) {rnorm(p[d])}, simplify=FALSE), 
              gamma=sapply(1:D, function(d) {rchisq(p[d], 1) + 1}, 
                           simplify=FALSE), 
              tau=rchisq(D, 1) + 1, sigma=rchisq(D, 1) + 1, 
              alphaf=rchisq(2, 1), alphad=rchisq(2, 1), 
              lambdaf=rchisq(1, 1) + 1, lambdad=rchisq(1, 1) + 1)
  
  srt <- proc.time()[3]
  for(k in 1:control$warmup) {
    if(control$trace) {
      if((k %% 10)==0) {
        cat("\r", "Warmup, iteration ", k)
      }
    }
    
    par <- .draw_post(par, x, y, C, Z, hyper, n, p, D)
  }
  for(k in 1:control$niter) {
    if(control$trace) {
      if((k %% 10)==0) {
        cat("\r", "Sampling from posterior, iteration ", k)
      }
    }
    
    par <- .draw_post(par, x, y, C, Z, hyper, n, p, D)
    for(d in 1:D) {
      draws$beta[[d]][, k] <- par$beta[[d]]
      draws$gamma[[d]][, k] <- par$gamma[[d]]
    }
    draws$tau[, k] <- par$tau
    draws$sigma[, k] <- par$sigma
    draws$alphaf[, k] <- par$alphaf
    draws$alphad[, k] <- par$alphad
    draws$lambdaf[k] <- par$lambdaf
    draws$lambdad[k] <- par$lambdad
  }
  draw.time <- proc.time()[3] - srt
  if(control$trace) {
    cat("\n", "Generated ", control$niter, " draws from the posterior in ", 
        round(draw.time), " seconds", sep="")  
  }
  out <- list(post=draws, time=draw.time)
  return(out)
}

acf.gibbs.semnig <- function(object, pars=NULL, d=NULL, j=NULL) {
  if(is.null(pars)) {
    pars <- names(object$post)
  }
  if(is.null(d)) {
    d <- 1:length(object$post$beta)
  }
  if(is.null(j)) {
    j <- lapply(object$post$beta, function(b) {1:nrow(b)})
  }
  out <- vector("list", length(pars))
  names(out) <- pars
  
  for(par in pars) {
    if(par %in% c("beta", "gamma")) {
      out[[par]] <- lapply(1:length(d), function(cd) {
        t(apply(object$post[[par]][[d[cd]]][j[[cd]], ], 1, function(s) {
          acf(s, plot=FALSE)$acf}))})
    }
    if(par %in% c("tau", "sigma")) {
      out[[par]] <- t(apply(object$post[[par]], 1, function(s) {
        acf(s, plot=FALSE)$acf}))
    }
    if(par %in% c("alphaf", "alphad")) {
      out[[par]] <- t(apply(object$post[[par]], 1, function(s) {
        acf(s, plot=FALSE)$acf}))
    }
    if(par %in% c("lambdaf", "lambdad")) {
      out[[par]] <- as.numeric(acf(object$post[[par]], plot=FALSE)$acf)
    }
  }
  return(out)
}

.plot.acf <- function(y) {
  plot(0, 0, type="n", xlim=c(1, length(y)), ylim=c(0, 1))
  for(l in 1:length(y)) {
    lines(c(l, l), y=c(0, y[l]))
  }
}

plot.acf.gibbs.semnig <- function(object, pars=NULL) {
  if(is.null(pars)) {
    pars <- names(object)
  }
  for(par in pars) {
    if(par %in% c("beta", "gamma")) {
      for(d in 1:length(object[[par]] )) {
        for(j in 1:nrow(object[[par]][[d]])) {
          .plot.acf(object[[par]][[d]][j, ])
        }
      }
    }
    if(par %in% c("tau", "sigma")) {
      for(d in 1:length(object[[par]] )) {
        .plot.acf(object[[par]][d, ])
      }
    }
    if(par %in% c("alphaf", "alphad")) {
      for(j in 1:nrow(object[[par]])) {
        .plot.acf(object[[par]][j, ])
      }
    }
    if(par %in% c("lambdaf", "lambdad")) {
      .plot.acf(object[[par]])
    }
  }
}
