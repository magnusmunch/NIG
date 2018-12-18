library(R2OpenBUGS)

path.mod <- "/Users/magnusmunch/Documents/OneDrive/PhD/project_cambridge/code/"        
path.res <- "/Users/magnusmunch/Documents/OneDrive/PhD/project_cambridge/results/"
path.pgm <- "/Users/magnusmunch/.wine/drive_c/Program Files (x86)/OpenBUGS/OpenBUGS323/OpenBUGS.exe"

### school example
sink(paste(path.mod, "test1_openbugs_model1.txt", sep=""))        
cat("
model {
  for (j in 1:J)
  {
    y[j] ~ dnorm (theta[j], tau.y[j])
    theta[j] ~ dnorm (mu.theta, tau.theta)
    tau.y[j] <- pow(sigma.y[j], -2)
  }
  mu.theta ~ dnorm (0.0, 1.0E-6)
  tau.theta <- pow(sigma.theta, -2)
  sigma.theta ~ dunif (0, 1000)
}
", fill=TRUE)
sink()


data(schools)
J <- nrow(schools)
y <- schools$estimate
sigma.y <- schools$sd
data <- list ("J", "y", "sigma.y")

inits <- function(){
  list(theta=rnorm(J, 0, 100), mu.theta=rnorm(1, 0, 100), 
       sigma.theta=runif(1, 0, 100))
}

schools.sim <- bugs(data, inits, useWINE=TRUE,
                    model.file=paste(path.mod, "test1_openbugs_model1.txt", sep=""),
                    parameters=c("theta", "mu.theta", "sigma.theta"),
                    n.chains=3, n.iter=10000, n.burnin=500,
                    working.directory=path.res,
                    OpenBUGS.pgm=path.pgm)


print(schools.sim, digits=3)
schools.DIC <- schools.sim$DIC

plot(schools.sim) # gives a summary plot of coefficients and credible intervals
par(mfrow=c(1, 1))
plot(schools.sim$sims.array[, 1, 1], type="l", 
     ylim=range(schools.sim$sims.array[, , 1]))
lines(schools.sim$sims.array[, 2, 1], type="l", col=2)
lines(schools.sim$sims.array[, 3, 1], type="l", col=3)

hist(schools.sim$sims.array[, , 1], breaks=80)



T <- 10
rho <- -1
tau <- (1 - T)*rho + 0.1
eigen(matrix(rho, T, T) + (tau - rho)*diag(T))$values
tau - rho
(T-1)*rho + tau

dens <- function(gamma2, a, b, c) {
  exp(-a/gamma2 - (log(gamma2) - b)^2/c)
}
x <- seq(0.01, 10, 0.01)
plot(x, dens(x, 0.1, -1, 1), type="l")
integrate(dens, 0.01, 10, 2, 0, 1)

