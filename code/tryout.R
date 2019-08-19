D <- 10
r <- 200
u <- 5
n <- 100
xu <- matrix(rnorm(u*n), nrow=n, ncol=u)
xr <- matrix(rnorm(r*n), nrow=n, ncol=r)
y <- matrix(rnorm(D*n), nrow=n, ncol=D)
C <- lapply(1:D, function(d) {as.matrix(rbinom(r, 1, 0.5))})

library(cambridge)
test <- semnig(x=rep(list(xr), D), y=y, C=C, unpenalized=rep(list(xu), D), 
               standardize=TRUE, intercept=TRUE, intercept.eb=TRUE, 
               mult.lambda=FALSE, fixed.eb="none", full.post=TRUE, init=NULL,
               control=list(conv.post=TRUE, trace=TRUE,
                            epsilon.eb=1e-3, epsilon.vb=1e-3, 
                            epsilon.opt=sqrt(.Machine$double.eps),
                            maxit.eb=100, maxit.vb=2, maxit.opt=100,
                            maxit.post=100))
plot(test$seq.eb$alpha[, 1], type="l", ylim=range(test$seq.eb$alpha))
lines(test$seq.eb$alpha[, 2], col=2)

test$eb$mprior

str(test)
x=rep(list(xr), D); y=y; C=NULL; unpenalized=rep(list(xu), D); 
standardize=TRUE; intercept=TRUE; intercept.eb=TRUE; 
mult.lambda=FALSE; fixed.eb="both"; full.post=TRUE; init=NULL;
control=list(conv.post=TRUE, trace=TRUE,
             epsilon.eb=1e-3, epsilon.vb=1e-3, 
             epsilon.opt=sqrt(.Machine$double.eps),
             maxit.eb=10, maxit.vb=2, maxit.opt=100,
             maxit.post=100)

mat1 <- solve(cbind(rbind(t(xu) %*% xu, t(xr) %*% xu), 
                    rbind(t(xu) %*% xr, t(xr) %*% xr + diag(h))))

Om <- solve(xr %*% diag(1/h) %*% t(xr) + diag(n))
mat2.11 <- solve(t(xu) %*% Om %*% xu)
mat2.12 <- - mat2.11 %*% t(xu) %*% xr %*% diag(1/h) +
  mat2.11 %*% t(xu) %*% xr %*% diag(1/h) %*% t(xr) %*% Om %*% xr %*% diag(1/h)
mat2.21 <- t(mat2.12)
mat2.22 <- diag(1/h) - diag(1/h) %*% t(xr) %*% Om %*% xr %*% diag(1/h) -
  diag(1/h) %*% t(xr) %*% xu %*% mat2.12 + 
  diag(1/h) %*% t(xr) %*% Om %*% xr %*% diag(1/h) %*% t(xr) %*% xu %*% mat2.12
mat2 <- cbind(rbind(mat2.11, mat2.21), rbind(mat2.12, mat2.22))
mat1[1:5, 1:5]
mat2[1:5, 1:5]

diag1 <- diag(diag(1/h) %*% t(xr) %*% Om %*% xr %*% diag(1/h))
diag2 <- rowSums((diag(1/h) %*% t(xr) %*% Om)*(diag(1/h) %*% t(xr)))
all.equal(diag1, diag2)
