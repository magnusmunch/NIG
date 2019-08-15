r <- 200
u <- 5
n <- 100
xu <- matrix(rnorm(u*n), nrow=n, ncol=u)
xr <- matrix(rnorm(r*n), nrow=n, ncol=r)
h <- rep(1, r)

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

