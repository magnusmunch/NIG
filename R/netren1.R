library(plot3Drgl)
mll <- function(tau, rho, T, D, w1, w2, w3) {
  if(tau > max(rho - T*rho, rho, 0)) {
    out <- D*(T^2 + T)*log(tau)/2 + D*T*tau*log(tau)/2 + 
      D*(T + 1)*log(tau - (1 - T)*rho)/2 + D*tau*log(tau - (1 - T)*rho)/2 + 
      D*(T^2 - 1)*log(tau - rho)/2 + D*(T - 1)*tau*log(tau - rho)/2 - 
      D*sum(lgamma(tau/2 + (c(1:T) + 1)/2)) - tau^2*w1 + tau*rho*w2 + tau*w3
    if(is.infinite(out)) {out <- NA}
  } else {
    out <- NA
  }
  return(out)
}

rho <- seq(-1000, 1000, 1)
tau <- seq(0.01, 1000, 1)
z <- matrix(NA, length(tau), length(rho))
for(i in 1:nrow(z)) {
  for(j in 1:ncol(z)) {
    z[i, j] <- mll(tau[i], rho[j], 20, 20, 10, -9, -10)
  }
}


persp3Drgl(tau, rho, z, xlab=expression(tau), ylab=expression(rho), zlab="mll")
par3d(windowRect=c(20, 30, 800, 800))


tau[which(z==max(z, na.rm=TRUE), arr.ind=TRUE)[1]]
rho[which(z==max(z, na.rm=TRUE), arr.ind=TRUE)[2]]

plot(x, digamma(x), type="l")

