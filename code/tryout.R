### MM algorithm
# libraries
library(cambridge)

# data generation
set.seed(123)
D <- 100
n <- rep(100, D)
p <- rep(200, D)
x <- matrix(rnorm(n[1]*p[1]), nrow=n[1], ncol=p[1])
rownames(x) <- 1:nrow(x)
alpha <- c(1, 1)
C <- lapply(1:D, function(d) {
  cbind(matrix(rnorm(p[d])))})
Z <- cbind(matrix(rnorm(D)))
G <- ncol(C[[1]])
H <- ncol(Z)
lambda <- rep(1, D)
beta <- lapply(1:D, function(d) {
  rnorm(p[d], 0, lambda[d]^2*exp(colSums(t(
    cbind(C[[d]], matrix(rep(Z[d, ], each=p[d]), ncol=H)))*alpha)))})
y <- lapply(1:D, function(d) {
  vec <- rnorm(n[d], x %*% beta[[d]], 1);
  names(vec) <- 1:n[d];
  return(vec)})

# estimation
test.ebridge <- ebridge2(x, y, C, Z)
test.ebridge$alphad

vprior <- lapply(1:D, function(d) {
  test.ebridge$lambda[d]^2*exp(
    colSums(t(C[[d]])*test.ebridge$alphaf) + sum(Z[d, ]*test.ebridge$alphaf))})
plot(Z, sapply(vprior, function(s) {mean(log(s))}))
plot(unlist(C), log(unlist(vprior)))
hist(log(unlist(vprior)))
  
