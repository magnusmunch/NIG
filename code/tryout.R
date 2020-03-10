# 
# .f.optim(rep(0, H + G), lambda=hyper$lambda, nu=hyper$nu, zeta=hyper$zeta, 
#          Cmat=Cmat, Z=Z, n=n, p=p, D=D, G=G, y=y, x=x, yty=yty)
# alpha=rep(0, H + G); lambda=hyper$lambda; nu=hyper$nu; zeta=hyper$zeta; 
# Cmat=Cmat; Z=Z; n=n; p=p; D=D; G=G; y=y; x=x; yty=yty
# 


library(glmnet)
set.seed(2020)
n <- 80
ntest <- 1000
p <- 100
D <- 50
C <- replicate(D, matrix(rnorm(p)), simplify=FALSE)
Z <- matrix(rnorm(D))
G <- ncol(C[[1]])
H <- ncol(Z)
lambda <- 1
alpha <- c(1, 1)
beta <- sapply(1:D, function(d) {
  rnorm(p, 0, lambda^2*exp(C[[d]]*alpha[1]/2)*exp(Z[d, ]*alpha[2]/2))})
x <- scale(matrix(rnorm(n*p), ncol=p, nrow=n))
xtest <- scale(matrix(rnorm(ntest*p), ncol=p, nrow=ntest))
y <- scale(x %*% beta + matrix(rnorm(n*D), ncol=D, nrow=n))
ytest <- scale(xtest %*% beta + matrix(rnorm(ntest*D), ncol=D, nrow=ntest))

x=xtrain[[1]]; y=ytrain; C=C2; Z=Z2; mult.lambda=TRUE
hyper=list(lambda=NULL, zeta=0, nu=0);
control=list(epsilon=sqrt(.Machine$double.eps), 
             maxit=500, trace=TRUE)



test1 <- ebridge(x, y, C, Z, mult.lambda=FALSE,
                 hyper=list(lambda=1, zeta=0, nu=0),
                 control=list(epsilon=sqrt(.Machine$double.eps), 
                              maxit=500, trace=TRUE))
best1 <- test1$beta1
best2 <- test1$beta2
best3 <- sapply(test1$glmnet.fit2, function(s) {coef(s)[-1, 1]})

pairs(cbind(ebridge1=c(best1), ebridge2=c(best2), ridge=c(best3)))

pred1 <- xtest %*% best1
pred2 <- xtest %*% best2
pred3 <- xtest %*% best3

c(ebridge1=mean((pred1 - ytest)^2), ebridge2=mean((pred2 - ytest)^2),
  ridge=mean((pred3 - ytest)^2), null=mean((ytest - colMeans(ytest))^2))

hist(test1$tau^2)
hist(test1$gamma^2)
hist(test1$tau^2)
summary(c(test1$gamma^2))
test1$alphad
test1$alphaf



### permutation distribution
x <- expr.sel[[1]]
test1 <- Reduce("+", lapply(1:nrow(x), function(i) {
  as.matrix(x[i, ]) %*% t(as.matrix(x[i, ]))}))/nrow(x)
str(test1)
diag(test1)
  
cbind(diag(x[1, ] %*% t(x[1, ])  ), x[1, ]^2)

