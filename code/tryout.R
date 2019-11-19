n <- 100
p <- 150
D <- 100
x <- replicate(D, matrix(rnorm(n*p), nrow=n, ncol=p), simplify=FALSE)
beta <- replicate(D, rnorm(p), simplify=FALSE)
y <- sapply(1:D, function(d) {x[[d]] %*% beta[[d]]})
Z <- NULL
# Z <- cbind(1, rbinom(D, 1, 0.5))
C <- lapply(1:D, function(d) {cbind(1, rbinom(p, 1, 0.5))})
unpenalized=NULL
standardize=FALSE
intercept=FALSE
fixed.eb="none"
full.post=FALSE
init=NULL
control=list(conv.post=TRUE, trace=TRUE,
             epsilon.eb=1e-3, epsilon.vb=1e-3,
             maxit.eb=10, maxit.vb=2, maxit.post=100)
library(Rcpp)
sourceCpp("rpackage/src/aux_functions.cpp")