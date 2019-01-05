# ---- figures ----
# ---- boxplot_igaussian_res1_sigma2 ----
res1 <- read.table("../results/simulations_igaussian_res1.csv")
set1 <- read.table("../results/simulations_igaussian_set1.csv")
names1 <- c("inv Gauss", "ind inv Gauss", "inv Gamma")
temp1 <- list(res1[, substr(colnames(res1), 1, 14)=="zeta.inv.Gauss"],
              res1[, substr(colnames(res1), 1, 18)=="zeta.ind.inv.Gauss"],
              res1[, substr(colnames(res1), 1, 5)=="cpost"],
              res1[, substr(colnames(res1), 1, 5)=="dpost"])
temp2 <- cbind(colMeans((2*temp1[[1]]/(set1$n + set1$p - 1) - as.numeric(set1[, substr(colnames(set1), 1, 5)=="sigma"])^2)^2),
               colMeans((2*temp1[[2]]/(set1$n - 1) - 
                           as.numeric(set1[, substr(colnames(set1), 1, 5)=="sigma"])^2)^2),
               colMeans((temp1[[4]]/(temp1[[3]] - 1) - 
                           as.numeric(set1[, substr(colnames(set1), 1, 5)=="sigma"])^2)^2))
par(cex=1.3)
boxplot(temp2, outline=FALSE, names=names1, ylab="MSE")
legend("topright", title=expression(bold("median MSE")), bty="n", 
       paste(names1, round(apply(temp2, 2, median), 2), sep="="))
par(cex=1)

# ---- boxplot_igaussian_res1_gamma2 ----
res1 <- read.table("../results/simulations_igaussian_res1.csv")
set1 <- read.table("../results/simulations_igaussian_set1.csv")
names1 <- c("inv Gauss", "ind inv Gauss", "inv Gamma")
temp1 <- list(res1[, substr(colnames(res1), 1, 15)=="delta.inv.Gauss"],
              res1[, substr(colnames(res1), 1, 15)=="theta.inv.Gauss"],
              res1[, substr(colnames(res1), 1, 16)=="lambda.inv.Gauss"],
              res1[, substr(colnames(res1), 1, 19)=="delta.ind.inv.Gauss"],
              res1[, substr(colnames(res1), 1, 19)=="theta.ind.inv.Gauss"],
              res1[, substr(colnames(res1), 1, 20)=="lambda.ind.inv.Gauss"],
              res1[, substr(colnames(res1), 1, 15)=="apost.inv.Gamma"],
              res1[, substr(colnames(res1), 1, 15)=="bpost.inv.Gamma"])
temp2 <- list(as.matrix(sqrt(temp1[[1]]*temp1[[2]]^2/temp1[[3]]))*
                ratio_besselK_cpp(as.matrix(sqrt(
                  temp1[[1]]*temp1[[3]])/temp1[[2]]), set1$p),
              as.matrix(sqrt(temp1[[4]]*temp1[[5]]^2/temp1[[6]]))*
                ratio_besselK_cpp(as.matrix(sqrt(
                  temp1[[4]]*temp1[[6]])/temp1[[5]]), set1$p),
              as.matrix(temp1[[8]]/(temp1[[7]] - 1)))
temp3 <- sapply(temp2, function(m) {
  colMeans((t(m) - as.numeric(set1[, substr(colnames(set1), 1, 5)=="gamma"])^2)^2)})
temp4 <- sapply(temp2, function(m) {apply(t(m), 2, function(r) {
  cor(r, as.numeric(set1[, substr(colnames(set1), 1, 5)=="gamma"])^2)})})
par(cex=1.3)
layout(matrix(c(rep(1, 4), rep(2, 4)), 2, 4), widths=1, heights=1, respect=TRUE)
boxplot(temp3, names=names1, ylab="MSE", main="a)")
legend("right", title=expression(bold("median MSE")), bty="n", 
       paste(names1, round(apply(temp3, 2, median), 2), sep="="))
boxplot(temp4, names=names1, ylab="Correlation", main="b)")
legend("bottomleft", title=expression(bold("median correlation")), bty="n", 
       paste(names1, round(apply(temp4, 2, median), 2), sep="="))
par(cex=1)
layout(matrix(1, 1, 1), widths=1, heights=1, respect=TRUE)

# ---- boxplot_igaussian_res1_mu ----
res1 <- read.table("../results/simulations_igaussian_res1.csv")
names1 <- c("inv Gauss", "ind inv Gauss", "inv Gamma")
temp1 <- cbind(res1$mse.mu.inv.Gauss, res1$mse.mu.ind.inv.Gauss, 
               res1$mse.mu.inv.Gamma)
temp2 <- cbind(res1$cor.mu.inv.Gauss, res1$cor.mu.ind.inv.Gauss, 
               res1$cor.mu.inv.Gamma)
par(cex=1.3)
layout(matrix(c(rep(1, 4), rep(2, 4)), 2, 4), widths=1, heights=1, respect=TRUE)
boxplot(temp1, outline=FALSE, names=names1, ylab="MSE", main="a)")
legend("topleft", title=expression(bold("median MSE")), bty="n", 
       paste(names1, round(apply(temp1, 2, median), 2), sep="="))
boxplot(temp2, names=names1, ylab="Correlation", main="b)")
legend("bottomleft", title=expression(bold("median correlation")), bty="n", 
       paste(names1, round(apply(temp2, 2, median), 2), 
             sep="="))
par(cex=1)
layout(matrix(1, 1, 1), widths=1, heights=1, respect=TRUE)

# ---- boxplot_igaussian_res1_theta ----
res1 <- read.table("../results/simulations_igaussian_res1.csv")
set1 <- read.table("../results/simulations_igaussian_set1.csv")
names1 <- c("inv Gauss", "ind inv Gauss", "inv Gamma")
temp1 <- list(res1[, substr(colnames(res1), 1, 15)=="theta.inv.Gauss"],
              res1[, substr(colnames(res1), 1, 19)=="theta.ind.inv.Gauss"],
              res1[, substr(colnames(res1), 1, 16)=="aprior.inv.Gamma"],
              res1[, substr(colnames(res1), 1, 16)=="bprior.inv.Gamma"])
temp2 <- list(colMeans((t(temp1[[1]]) - 
                          as.numeric(set1[, substr(colnames(set1), 1, 5)=="theta"]))^2),
              colMeans((t(temp1[[2]]) - 
                          as.numeric(set1[, substr(colnames(set1), 1, 5)=="theta"]))^2),
              colMeans((t(temp1[[4]]/(temp1[[3]] - 1)) - 
                          unique(as.numeric(set1[, substr(colnames(set1), 1, 5)=="theta"])))^2))
temp3 <- list(apply(t(temp1[[1]]), 2, function(r) {
  cor(r, as.numeric(set1[, substr(colnames(set1), 1, 5)=="theta"]))}),
  apply(t(temp1[[2]]), 2, function(r) {
    cor(r, as.numeric(set1[, substr(colnames(set1), 1, 5)=="theta"]))}),
  apply(t(temp1[[4]]/(temp1[[3]] - 1)), 2, function(r) {
    cor(r, unique(as.numeric(set1[, substr(colnames(set1), 1, 5)=="theta"])))}))
par(cex=1.3)
layout(matrix(c(rep(1, 4), rep(2, 4)), 2, 4), widths=1, heights=1, respect=TRUE)
boxplot(temp2, names=names1, ylab="MSE", main="a)")
legend("bottomright", title=expression(bold("median MSE")), bty="n",
       paste(names1, round(sapply(temp2, median), 2), sep="="))
boxplot(temp3, names=names1, ylab="Correlation", main="b)")
legend("bottomleft", title=expression(bold("median correlation")), bty="n",
       paste(names1, round(sapply(temp3, median), 2), sep="="))
par(cex=1)
layout(matrix(1, 1, 1), widths=1, heights=1, respect=TRUE)

# ---- boxplot_igaussian_res1_alpha ----
res1 <- read.table("../results/simulations_igaussian_res1.csv")
set1 <- read.table("../results/simulations_igaussian_set1.csv")

names1 <- c("inv Gauss", "ind inv Gauss")
temp1 <- list(res1[, substr(colnames(res1), 1, 15)=="alpha.inv.Gauss"],
              res1[, substr(colnames(res1), 1, 19)=="alpha.ind.inv.Gauss"])
temp2 <- list(colMeans((t(temp1[[1]]) - as.numeric(set1[, substr(colnames(set1), 1, 5)=="alpha"]))^2),
              colMeans((t(temp1[[2]]) - as.numeric(set1[, substr(colnames(set1), 1, 5)=="alpha"]))^2))
temp3 <- list(apply(t(temp1[[1]]), 2, function(r) {cor(r, as.numeric(set1[, substr(colnames(set1), 1, 5)=="alpha"]))}),
              apply(t(temp1[[2]]), 2, function(r) {cor(r, as.numeric(set1[, substr(colnames(set1), 1, 5)=="alpha"]))}))
par(cex=1.3)
layout(matrix(c(rep(1, 4), rep(2, 4)), 2, 4), widths=1, heights=1, respect=TRUE)
boxplot(temp2, names=names1, ylab="MSE", main="a)")
legend("topright", title=expression(bold("median MSE")), bty="n",
       paste(names1, round(sapply(temp2, median), 2), sep="="))
boxplot(temp3, names=names1, ylab="Correlation", main="b)")
legend("bottomright", title=expression(bold("median correlation")), bty="n",
       paste(names1, round(sapply(temp3, median), 2), sep="="))
par(cex=1)
layout(matrix(1, 1, 1), widths=1, heights=1, respect=TRUE)

# ---- lines_igaussian_res1_elbo ----
fit1 <- read.table("../results/simulations_igaussian_fit1.csv")
par(cex=1.3)
layout(matrix(c(rep(c(1, 1, 2, 2), times=2), rep(c(0, 3, 3, 0), times=2)), 
              4, 4, byrow=TRUE), widths=1, heights=1, respect=TRUE)
plot(fit1[, substr(colnames(fit1), 1, 14)=="inv.Gauss.elbo"][, 1], type="l",
     ylim=range(fit1[, substr(colnames(fit1), 1, 14)=="inv.Gauss.elbo"]),
     main="a)", xlab="Iteration", ylab="ELBO")
for(i in 2:100) {
  lines(fit1[, substr(colnames(fit1), 1, 14)=="inv.Gauss.elbo"][, i], col=i)
}
plot(fit1[, substr(colnames(fit1), 1, 18)=="ind.inv.Gauss.elbo"][, 1], type="l",
     ylim=range(fit1[, substr(colnames(fit1), 1, 18)=="ind.inv.Gauss.elbo"]),
     main="b)", xlab="Iteration", ylab="ELBO")
for(i in 2:100) {
  lines(fit1[, substr(colnames(fit1), 1, 18)=="ind.inv.Gauss.elbo"][, i], col=i)
}
plot(fit1[, substr(colnames(fit1), 1, 14)=="inv.Gamma.elbo"][, 1], type="l",
     ylim=range(fit1[, substr(colnames(fit1), 1, 14)=="inv.Gamma.elbo"]),
     main="c)", xlab="Iteration", ylab="ELBO")
for(i in 2:100) {
  lines(fit1[, substr(colnames(fit1), 1, 14)=="inv.Gamma.elbo"][, i], col=i)
}
par(cex=1)
layout(matrix(1, 1, 1), widths=1, heights=1, respect=TRUE)

# ---- lines_igaussian_res1_theta ----
fit1 <- read.table("../results/simulations_igaussian_fit1.csv")
set1 <- read.table("../results/simulations_igaussian_set1.csv")

par(cex=1.3)
layout(matrix(c(rep(c(1, 1, 2, 2), times=2), rep(c(0, 3, 3, 0), times=2)), 
              4, 4, byrow=TRUE), widths=1, heights=1, respect=TRUE)
plot(fit1[, substr(colnames(fit1), 1, 15)=="inv.Gauss.theta"][, 1], type="l",
     ylim=range(fit1[, substr(colnames(fit1), 1, 15)=="inv.Gauss.theta"]),
     main="a)", xlab="Iteration", ylab=expression(theta[d]))
for(i in 2:5) {
  lines(fit1[, substr(colnames(fit1), 1, 15)=="inv.Gauss.theta"][
    , (i - 1)*(set1$D/set1$nclass) + 1], col=i)
}
plot(fit1[, substr(colnames(fit1), 1, 19)=="ind.inv.Gauss.theta"][, 1], 
     type="l", 
     ylim=range(fit1[, substr(colnames(fit1), 1, 19)=="ind.inv.Gauss.theta"]),
     main="b)", xlab="Iteration", ylab=expression(theta[d]))
for(i in 2:5) {
  lines(fit1[, substr(colnames(fit1), 1, 19)=="ind.inv.Gauss.theta"][
    , (i - 1)*(set1$D/set1$nclass) + 1], col=i)
}
plot(sort(fit1[, substr(colnames(fit1), 1, 16)=="inv.Gamma.bprior"][, 1]/
            (fit1[, substr(colnames(fit1), 1, 16)=="inv.Gamma.aprior"][
              , 1] - 1)), 
     type="l", 
     ylim=range(fit1[, substr(colnames(fit1), 1, 16)=="inv.Gamma.bprior"]/
                  (fit1[, substr(colnames(fit1), 1, 16)=="inv.Gamma.aprior"] - 
                     1)),
     main="c)", xlab="Iteration", ylab=expression(beta/(alpha - 1)))
for(i in 2:5) {
  lines(sort((fit1[, substr(colnames(fit1), 1, 16)=="inv.Gamma.bprior"][, ]/
                (fit1[, substr(colnames(fit1), 1, 16)=="inv.Gamma.aprior"][
                  , 1] - 1))[, i]), col=i)
}
par(cex=1)
layout(matrix(1, 1, 1), widths=1, heights=1, respect=TRUE)

# ---- boxplot_igaussian_res2_theta ----
res2 <- read.table("../results/simulations_igaussian_res2.csv")
set2 <- read.table("../results/simulations_igaussian_set2.csv")
names1 <- c("inv Gauss", "ind inv Gauss", "inv Gamma")
temp1 <- list(res2[, substr(colnames(res2), 1, 15)=="theta.inv.Gauss"],
              res2[, substr(colnames(res2), 1, 19)=="theta.ind.inv.Gauss"],
              res2[, substr(colnames(res2), 1, 16)=="aprior.inv.Gamma"],
              res2[, substr(colnames(res2), 1, 16)=="bprior.inv.Gamma"])
temp2 <- list(colMeans((t(temp1[[1]]) - 
                          as.numeric(set2[, substr(colnames(set2), 1, 5)=="theta"]))^2),
              colMeans((t(temp1[[2]]) - 
                          as.numeric(set2[, substr(colnames(set2), 1, 5)=="theta"]))^2),
              colMeans((t(temp1[[4]]/(temp1[[3]] - 1)) - 
                          unique(as.numeric(set2[, substr(colnames(set2), 1, 5)=="theta"])))^2))
par(cex=1.3)
boxplot(temp2, names=names1, ylab="MSE", main="a)")
legend("bottomright", title=expression(bold("median MSE")), bty="n",
       paste(names1, round(sapply(temp2, median), 2), sep="="))
par(cex=1)

# ---- lines_igaussian_res2_alpha ----
fit2 <- read.table("../results/simulations_igaussian_fit2.csv")
set2 <- read.table("../results/simulations_igaussian_set2.csv")

par(cex=1.3)
layout(matrix(rep(c(1, 1, 2, 2), times=2), 2, 4, byrow=TRUE), widths=1, 
       heights=1, respect=TRUE)
plot(fit2[, substr(colnames(fit2), 1, 15)=="inv.Gauss.alpha"][, 1], type="l",
     ylim=range(fit2[, substr(colnames(fit2), 1, 15)=="inv.Gauss.alpha"]),
     main="a)", xlab="Iteration", ylab=expression(theta[d]))
for(i in 2:5) {
  lines(fit2[, substr(colnames(fit2), 1, 15)=="inv.Gauss.alpha"][, i], col=i)
}
plot(fit2[, substr(colnames(fit2), 1, 19)=="ind.inv.Gauss.alpha"][, 1], 
     type="l", 
     ylim=range(fit2[, substr(colnames(fit2), 1, 19)=="ind.inv.Gauss.alpha"]),
     main="b)", xlab="Iteration", ylab=expression(theta[d]))
for(i in 2:5) {
  lines(fit2[, substr(colnames(fit2), 1, 19)=="ind.inv.Gauss.alpha"][, i], 
        col=i)
}
par(cex=1)
layout(matrix(1, 1, 1), widths=1, heights=1, respect=TRUE)
