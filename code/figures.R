################################# main document ################################
# ---- figures ----
# ---- dens_kappa ----
library(sp)
labels <- c("NIG", "Student's t", "lasso")
col <- bpy.colors(length(labels), cutoff.tail=0.3)
lty <- c(1:length(labels))
dknigig <- function(x, lambda1, theta1, lambda2, theta2) {
  sqrt(lambda1*lambda2)/pi*exp(lambda1/theta1 + lambda2/theta2)*
    x^(-2)*(1/x - 1)^(-3/2)*
    besselK(sqrt((lambda1/theta1^2 + lambda2/(1/x - 1))*
                   (lambda1 + lambda2*(1/x - 1)/theta2^2)), 0)
}
dknig <- function(x, lambda, theta) {
  sqrt(lambda/(2*pi))*(1/x - 1)^(-3/2)*(1/x^2)*
    exp(-0.5*lambda*(1/x - 1 - theta)^2/(theta^2*(1/x - 1)))
}
dkstudentst <- function(x, a, b) {
  (b^a/gamma(a))*(1/x - 1)^(-a-1)*exp(-b/(1/x - 1))*(1/x^2)
}
dklasso <- function(x, lambda) {
  lambda*exp(lambda)*(1/x^2)*exp(-lambda/x)
}
x.seq <- seq(.00001, 0.99999, length.out=10000)
opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.2, 1, 1), cex=2, cex.axis=1.7, cex.lab=1.7,
    cex.main=1.7, lwd=1.7)
layout(matrix(rep(rep(c(1:3), each=2), 2), nrow=2, ncol=6, byrow=TRUE))
plot(x.seq, dknig(x.seq, 0.1, 10), type="l", xlab=expression(kappa),
     ylab=expression(pi(kappa)), main="(a)", col=col[1], lty=lty[1])
lines(x.seq, dknig(x.seq, 1, 10), col=col[1], lty=lty[2])
lines(x.seq, dknig(x.seq, 1, 1), col=col[1], lty=lty[3])

plot(x.seq, dkstudentst(x.seq, a=0.9, b=1.3), type="l", ylim=c(0, 4), 
     xlab=expression(kappa), ylab=expression(pi(kappa)), main="(b)", col=col[2], 
     lty=lty[1])
lines(x.seq, dkstudentst(x.seq, a=1, b=4), col=col[2], lty=lty[2])
lines(x.seq, dkstudentst(x.seq, a=1, b=0.2), col=col[2], lty=lty[3])

plot(x.seq, dklasso(x.seq, 1), type="l", ylim=c(0, 3.5), 
     xlab=expression(kappa), ylab=expression(pi(kappa)), main="(c)", col=col[3], 
     lty=lty[1])
lines(x.seq, dklasso(x.seq, 0.2), col=col[3], lty=lty[2])
lines(x.seq, dklasso(x.seq, 0.1), col=col[3], lty=lty[3])
par(opar)

# ---- simulations_gdsc_est1 ----
library(sp)
res <- read.table("results/simulations_gdsc_res1.txt", row.names=NULL)
temp <- res[, 1]
res <- as.matrix(res[, -1])
rownames(res) <- temp

alphaf <- c(1, 1, 3, 7)
phi <- 1/as.numeric(alphaf %*% t(cbind(1, rbind(0, diag(3)))))
est1 <- sapply(0:3, function(s) {
  res[rownames(res)==paste0("alphaf" ,s), 2]})
est2 <- 1/(est1 %*% t(cbind(1, rbind(0, diag(3)))))

col <- bpy.colors(length(alphaf), cutoff.tail=0.3)
labels1 <- expression(alpha["feat,0"], alpha["feat,1"], alpha["feat,2"], 
                      alpha["feat,3"])
labels2 <- expression(phi["1"], phi["2"], phi["3"], phi["4"])
ylim1 <- range(c(unlist(boxplot(est1, plot=FALSE)[c("stats", "out")]), alphaf))
ylim2 <- range(c(unlist(boxplot(est2, plot=FALSE)[c("stats", "out")]), phi))
opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(rep(c(1:2), each=2), 2), nrow=2, ncol=4, byrow=TRUE))
boxplot(est1, ylim=ylim1, main="(a)", ylab=expression(hat(alpha)), 
        names=labels1, col=col[c(1:4)], cex=1.5, cex.lab=1.5, cex.axis=1.5)
points(c(1:4), alphaf, pch=2, col=col[c(1:4)], cex=1.5, 
       cex.lab=1.5, cex.axis=1.4)
boxplot(est2, ylim=ylim2, main="(b)", ylab=expression(hat(phi)), 
        names=labels2, col=col[c(1:4)], cex=1.5, cex.lab=1.5, cex.axis=1.5)
points(c(1:4), phi, pch=2, col=col[c(1:4)], cex=1.5, 
       cex.lab=1.5, cex.axis=1.5)
par(opar)

# ---- simulations_gdsc_est2 ----
library(sp)
res <- read.table("results/simulations_gdsc_res2.txt", row.names=NULL)
temp <- res[, 1]
res <- as.matrix(res[, -1])
rownames(res) <- temp

alphad <- c(1, 1, 3, 7)
chi <- 1/as.numeric(alphad %*% t(cbind(1, rbind(0, diag(3)))))
est1 <- sapply(0:3, function(s) {
  res[rownames(res)==paste0("alphad" ,s), 2]})
est2 <- 1/(est1 %*% t(cbind(1, rbind(0, diag(3)))))

col <- bpy.colors(length(alphad), cutoff.tail=0.3)
labels1 <- expression(alpha["drug,0"], alpha["drug,1"], alpha["drug,2"], 
                      alpha["drug,3"])
labels2 <- expression(chi["1"], chi["2"], chi["3"], chi["4"])
ylim1 <- range(c(unlist(boxplot(est1, plot=FALSE)[c("stats", "out")]), alphad))
ylim2 <- range(c(unlist(boxplot(est2, plot=FALSE)[c("stats", "out")]), chi))
opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(rep(c(1:2), each=2), 2), nrow=2, ncol=4, byrow=TRUE))
boxplot(est1, ylim=ylim1, main="(a)", ylab=expression(hat(alpha)), 
        names=labels1, col=col[c(1:4)], cex=1.5, cex.lab=1.5, cex.axis=1.5)
points(c(1:4), alphaf, pch=2, col=col[c(1:4)], cex=1.5, 
       cex.lab=1.5, cex.axis=1.4)
boxplot(est2, ylim=ylim2, main="(b)", ylab=expression(hat(chi)), 
        names=labels2, col=col[c(1:4)], cex=1.5, cex.lab=1.5, cex.axis=1.5)
points(c(1:4), phi, pch=2, col=col[c(1:4)], cex=1.5, 
       cex.lab=1.5, cex.axis=1.5)
par(opar)

# ---- simulations_gdsc_est3 ----
library(sp)
res <- read.table("results/simulations_gdsc_res3.txt", row.names=NULL)
temp <- res[, 1]
res <- as.matrix(res[, -1])
rownames(res) <- temp

alphaf <- c(1, 1, 3, 7)
alphad <- c(1, 1, 3, 7)
phi <- 1/as.numeric(alphaf %*% t(cbind(1, rbind(0, diag(3)))))
chi <- 1/as.numeric(alphad %*% t(cbind(1, rbind(0, diag(3)))))
est1 <- sapply(0:3, function(s) {
  res[rownames(res)==paste0("alphaf" ,s), 2]})
est2 <- 1/(est1 %*% t(cbind(1, rbind(0, diag(3)))))
est3 <- sapply(0:3, function(s) {
  res[rownames(res)==paste0("alphad" ,s), 2]})
est4 <- 1/(est3 %*% t(cbind(1, rbind(0, diag(3)))))

col <- bpy.colors(length(alphad), cutoff.tail=0.3)
labels1 <- expression(alpha["feat,0"], alpha["feat,1"], alpha["feat,2"], 
                      alpha["feat,3"])
labels2 <- expression(phi["1"], phi["2"], phi["3"], phi["4"])
labels3 <- expression(alpha["drug,0"], alpha["drug,1"], alpha["drug,2"], 
                      alpha["drug,3"])
labels4 <- expression(chi["1"], chi["2"], chi["3"], chi["4"])

ylim1 <- range(c(unlist(boxplot(est1, plot=FALSE)[c("stats", "out")]), alphaf))
ylim2 <- range(c(unlist(boxplot(est2, plot=FALSE)[c("stats", "out")]), phi))
ylim3 <- range(c(unlist(boxplot(est3, plot=FALSE)[c("stats", "out")]), alphad))
ylim4 <- range(c(unlist(boxplot(est4, plot=FALSE)[c("stats", "out")]), chi))
opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(c(rep(rep(c(1:2), each=2), 2), rep(rep(c(3:4), each=2), 2)), 
              nrow=4, ncol=4, byrow=TRUE))
boxplot(est1, ylim=ylim1, main="(a)", ylab=expression(hat(alpha)), 
        names=labels1, col=col[c(1:4)], cex=1.5, cex.lab=1.5, cex.axis=1.5)
points(c(1:4), alphaf, pch=2, col=col[c(1:4)], cex=1.5, 
       cex.lab=1.5, cex.axis=1.4)
boxplot(est2, ylim=ylim2, main="(b)", ylab=expression(hat(phi)), 
        names=labels2, col=col[c(1:4)], cex=1.5, cex.lab=1.5, cex.axis=1.5)
points(c(1:4), phi, pch=2, col=col[c(1:4)], cex=1.5, 
       cex.lab=1.5, cex.axis=1.5)
boxplot(est3, ylim=ylim3, main="(c)", ylab=expression(hat(alpha)), 
        names=labels3, col=col[c(1:4)], cex=1.5, cex.lab=1.5, cex.axis=1.5)
points(c(1:4), alphad, pch=2, col=col[c(1:4)], cex=1.5, 
       cex.lab=1.5, cex.axis=1.4)
boxplot(est4, ylim=ylim4, main="(d)", ylab=expression(hat(chi)), 
        names=labels4, col=col[c(1:4)], cex=1.5, cex.lab=1.5, cex.axis=1.5)
points(c(1:4), chi, pch=2, col=col[c(1:4)], cex=1.5, 
       cex.lab=1.5, cex.axis=1.5)
par(opar)

################################# presentation #################################
# ---- dens_beta_prior1 ---- 
library(GeneralizedHyperbolic)
dprior <- function(x, lambda, theta, sigma) {
  dnig(x, 0, sigma*sqrt(lambda), sqrt(lambda/(theta*sigma)), 0)
}
colors <- sp::bpy.colors(4)[-c(1, 4)]
pdf(file="figs/dens_beta_prior1.pdf", width=10)
curve(dprior(x, 0.1, 10, 1), -2, 2, col=colors[1], lwd=4, main="", yaxt="n", xaxt="n", 
      bty="n", ann=FALSE, n=1000)
curve(dprior(x, 10, 0.1, 1), -2, 2, col=colors[2], lwd=4, add=TRUE, main="", yaxt="n", 
      xaxt="n", bty="n", ann=FALSE, n=1000)
legend("topright", legend=c(NA, NA), fill=colors, box.col=NA, cex=3, 
       border=colors)
dev.off()


################################################################################
################################################################################
################################################################################

res <- read.table("results/analysis_gdsc_res1.txt", row.names=NULL)
temp <- res[, 1]
res <- as.matrix(res[, -1])
rownames(res) <- temp
test1 <- Reduce("rbind", by(res[substr(rownames(res), 1, 4)=="pmse", ], 
                            list(as.factor(rownames(res)[substr(rownames(res), 1, 
                                                                4)=="pmse"])), 
                            function(s) {apply(s, 2, median)}))
D <- 168
nfolds <- 10
test2 <- Reduce("rbind", by(res[substr(rownames(res), 1, 4)=="pmse", ], 
                            list(as.factor(rep(1:nfolds, each=D))), 
                            colMeans))

pairs(test1[, 1:9])
boxplot(test2)

res <- read.table("results/analysis_gdsc_res5.txt", row.names=NULL)
temp <- res[, 1]
res <- as.matrix(res[, -1])
rownames(res) <- temp
boxplot(res[substr(rownames(res), 1, 4)=="pmse", ])
sort(apply(res[substr(rownames(res), 1, 4)=="pmse", ], 2, median))

sort(apply(res[rownames(res)=="pmse", ], 2, median))
sort(apply(res[rownames(res)=="elbo", ], 2, median), decreasing=TRUE)
sort(apply(res[rownames(res)=="lpml", ], 2, median), decreasing=TRUE)
sort(apply(res[rownames(res)=="brankdist", ], 2, median))
boxplot(res[rownames(res)=="elbo", ])
boxplot(res[rownames(res)=="lpml", ])
boxplot(res[rownames(res)=="brankdist", ])

load("results/analysis_gdsc_fit1.Rdata")
fit1.semnig$eb$alphaf
fit1.semnig$eb$alphad
fit2.semnig$eb$alphaf
fit2.semnig$eb$alphad
fit3.semnig$eb$alphaf
fit4.semnig$eb$alphad
fit5.semnig$eb$alphaf
fit5.semnig$eb$alphad


test1 <- Reduce("cbind", fit2.semnig$vb$mu)
test2 <- sapply(fit1.ridge, function(s) {as.numeric(coef(s, s="lambda.min"))[-1]})
test3 <- sapply(1:D, function(d) {
  fit1.optim[[d]]$par[names(fit1.optim[[d]]$par) %in% paste0("beta[", 1:p[d], "]")]})
plot(test1, test3)
abline(a=0, b=1, col=2)

plot(colMeans((y - x[[1]] %*% test2)^2),
      colMeans((y - x[[1]] %*% test3)^2), xlab="ridge", ylab="NIG")
abline(a=0, b=1, col=2)

mean(colMeans((y - x[[1]] %*% test2)^2))
mean(colMeans((y - x[[1]] %*% test3)^2))

plot(y - x[[1]] %*% test2, y - x[[1]] %*% test3)
plot(y, x[[1]] %*% test2)
plot(y, x[[1]] %*% test3)
abline(a=0, b=1, col=2)

dim(x[[1]])
100*168
length(fit2.semnig$eb$alphaf) + length(fit2.semnig$eb$alphad) +
  length(fit2.semnig$eb$lambdaf) + length(fit2.semnig$eb$lambdad)


library(plotly)
mat <- matrix(colMeans(test$cvmat), ncol=length(chi), nrow=length(phi))
cv.plot <- plot_ly(x=chi, y=phi, z=mat) %>% 
  add_surface() %>%
  layout(title="Cross-validated mean squared error",
         scene=list(xaxis=list(title="chi"),
                    yaxis=list(title="phi"),
                    zaxis=list(title="MSE")))
cv.plot
