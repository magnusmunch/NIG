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
load("results/simulations_gdsc_set1.Rdata")
col <- bpy.colors(length(set$alphaf), cutoff.tail=0.3)
labels <- expression(alpha["feat,0"], alpha["feat,1"], alpha["feat,2"], 
                     alpha["feat,3"])
ylim1 <- range(c(unlist(boxplot(est[, c(1:4)], plot=FALSE)[c("stats", "out")]),
                 set$alphaf))
ylim2 <- range(c(unlist(boxplot(est[, c(6:9)], plot=FALSE)[c("stats", "out")]),
                 set$alphaf))
opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(rep(c(1:2), each=2), 2), nrow=2, ncol=4, byrow=TRUE))
boxplot(est[, c(1:4)], ylim=ylim1, main="(a)", ylab=expression(hat(alpha)), 
        names=labels, col=col[c(1:4)], cex=1.5, cex.lab=1.5, cex.axis=1.5)
points(c(1:4), set$alphaf, pch=2, col=col[c(1:4)], cex=1.5, 
       cex.lab=1.5, cex.axis=1.4)
boxplot(est[, c(6:9)], ylim=ylim2, main="(b)", ylab=expression(hat(alpha)), 
        names=labels, col=col[c(1:4)], cex=1.5, cex.lab=1.5, cex.axis=1.5)
points(c(1:4), set$alphaf, pch=2, col=col[c(1:4)], cex=1.5, 
       cex.lab=1.5, cex.axis=1.5)
par(opar)

# ---- simulations_gdsc_phi1 ----
library(sp)
load("results/simulations_gdsc_set1.Rdata")
phi <- as.numeric(1/(set$alphaf %*% t(unname(model.matrix(~ factor(1:4))))))
phi1 <- cbind(1/(est[, c(1:4)] %*% 
                   t(unname(model.matrix(~ factor(1:4))))),
              1/(est[, c(6:9)] %*% 
                   t(unname(model.matrix(~ factor(1:4))))))
col <- bpy.colors(length(phi), cutoff.tail=0.3)
labels <- expression(phi["1"], phi["2"], phi["3"], phi["4"])
ylim1 <- range(c(unlist(boxplot(phi1[, c(1:4)], plot=FALSE)[c("stats", "out")]),
                 phi))
ylim2 <- range(c(unlist(boxplot(phi1[, c(5:8)], plot=FALSE)[c("stats", "out")]),
                 phi))
opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(rep(c(1:2), each=2), 2), nrow=2, ncol=4, byrow=TRUE))
boxplot(phi1[, c(1:4)], ylim=ylim1, main="(a)", ylab=expression(hat(phi)), 
        names=labels, col=col[c(1:4)], cex=1.5, cex.lab=1.5, cex.axis=1.5)
points(c(1:4), phi, pch=2, col=col[c(1:4)], cex=1.5, 
       cex.lab=1.5, cex.axis=1.4)
boxplot(phi1[, c(5:8)], ylim=ylim2, main="(b)", ylab=expression(hat(phi)), 
        names=labels, col=col[c(1:4)], cex=1.5, cex.lab=1.5, cex.axis=1.5)
points(c(1:4), phi, pch=2, col=col[c(1:4)], cex=1.5, 
       cex.lab=1.5, cex.axis=1.5)
par(opar)

# ---- simulations_gdsc_bias1 ----
library(sp)
load("results/simulations_gdsc_set1.Rdata")
phi <- as.numeric(1/(set$alphaf %*% t(unname(model.matrix(~ factor(1:4))))))
phi1 <- cbind(1/(est[, c(1:4)] %*% 
                   t(unname(model.matrix(~ factor(1:4))))),
              1/(est[, c(6:9)] %*% 
                   t(unname(model.matrix(~ factor(1:4))))))
bias1 <- t(apply(phi1, 1, "/", rep(phi, 2)))
col <- bpy.colors(length(phi), cutoff.tail=0.3)
labels <- expression(hat(phi)["1"]/phi["1"], hat(phi)["2"]/phi["2"], 
                     hat(phi)["3"]/phi["3"], hat(phi)["4"]/phi["4"])
ylim1 <- range(c(unlist(boxplot(bias1[, c(1:4)], plot=FALSE)[c("stats", 
                                                               "out")])))
ylim2 <- range(c(unlist(boxplot(bias1[, c(5:8)], plot=FALSE)[c("stats", 
                                                               "out")])))
opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(rep(c(1:2), each=2), 2), nrow=2, ncol=4, byrow=TRUE))
boxplot(bias1[, c(1:4)], ylim=ylim1, main="(a)", 
        ylab=expression(hat(phi)/phi), 
        names=labels, col=col[c(1:4)], cex=1.5, cex.lab=1.5, cex.axis=1.5)
boxplot(bias1[, c(5:8)], ylim=ylim2, main="(b)", 
        ylab=expression(hat(phi)/phi), 
        names=labels, col=col[c(1:4)], cex=1.5, cex.lab=1.5, cex.axis=1.5)
par(opar)

# ---- simulations_gdsc_est2 ----
library(sp)
load("results/simulations_gdsc_set2.Rdata")
col <- bpy.colors(length(set$alphaf), cutoff.tail=0.3)
labels1 <- expression(alpha["feat,0"], alpha["feat,1"], alpha["feat,2"], 
                      alpha["feat,3"])
labels2 <- expression(alpha["drug,0"], alpha["drug,1"], alpha["drug,2"], 
                      alpha["drug,3"])
ylim1 <- range(c(unlist(boxplot(est[, c(1:4)], plot=FALSE)[c("stats", "out")]),
                 set$alphaf))
ylim2 <- range(c(unlist(boxplot(est[, c(11:14)], plot=FALSE)[c("stats", "out")]),
                 set$alphaf))
ylim3 <- range(c(unlist(boxplot(est[, c(5:8)], plot=FALSE)[c("stats", "out")]),
                 set$alphad))
ylim4 <- range(c(unlist(boxplot(est[, c(15:18)], plot=FALSE)[c("stats", "out")]),
                 set$alphad))
opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(c(rep(rep(c(1:2), each=2), 2), rep(rep(c(3:4), each=2), 2)), 
              nrow=4, ncol=4, byrow=TRUE))
boxplot(est[, c(1:4)], ylim=ylim1, main="(a)", ylab=expression(hat(alpha)), 
        names=labels1, col=col[c(1:4)], cex=1.5, cex.lab=1.5, cex.axis=1.5)
points(c(1:4), set$alphaf, pch=2, col=col[c(1:4)], cex=1.5, 
       cex.lab=1.5, cex.axis=1.4)
boxplot(est[, c(11:14)], ylim=ylim2, main="(b)", ylab=expression(hat(alpha)), 
        names=labels1, col=col[c(1:4)], cex=1.5, cex.lab=1.5, cex.axis=1.5)
points(c(1:4), set$alphaf, pch=2, col=col[c(1:4)], cex=1.5, 
       cex.lab=1.5, cex.axis=1.5)
boxplot(est[, c(5:8)], ylim=ylim3, main="(c)", ylab=expression(hat(alpha)), 
        names=labels2, col=col[c(1:4)], cex=1.5, cex.lab=1.5, cex.axis=1.5)
points(c(1:4), set$alphad, pch=2, col=col[c(1:4)], cex=1.5, 
       cex.lab=1.5, cex.axis=1.4)
boxplot(est[, c(15:18)], ylim=ylim4, main="(d)", ylab=expression(hat(alpha)), 
        names=labels2, col=col[c(1:4)], cex=1.5, cex.lab=1.5, cex.axis=1.5)
points(c(1:4), set$alphad, pch=2, col=col[c(1:4)], cex=1.5, 
       cex.lab=1.5, cex.axis=1.5)
par(opar)

# ---- simulations_gdsc_phi2 ----
library(sp)
load("results/simulations_gdsc_set2.Rdata")
phi <- as.numeric(1/(set$alphaf %*% t(unname(model.matrix(~ factor(1:4))))))
chi <- as.numeric(1/(set$alphad %*% t(unname(model.matrix(~ factor(1:4))))))
phi1 <- cbind(1/(est[, c(1:4)] %*% 
                   t(unname(model.matrix(~ factor(1:4))))),
              1/(est[, c(11:14)] %*% 
                   t(unname(model.matrix(~ factor(1:4))))))
chi1 <- cbind(1/(est[, c(5:8)] %*% 
                   t(unname(model.matrix(~ factor(1:4))))),
              1/(est[, c(15:18)] %*% 
                   t(unname(model.matrix(~ factor(1:4))))))
col <- bpy.colors(length(phi), cutoff.tail=0.3)
labels1 <- expression(phi["1"], phi["2"], phi["3"], phi["4"])
labels2 <- expression(chi["1"], chi["2"], chi["3"], chi["4"])
ylim1 <- range(c(unlist(boxplot(phi1[, c(1:4)], plot=FALSE)[c("stats", "out")]),
                 phi))
ylim2 <- range(c(unlist(boxplot(phi1[, c(5:8)], plot=FALSE)[c("stats", "out")]),
                 phi))
ylim3 <- range(c(unlist(boxplot(chi1[, c(1:4)], plot=FALSE)[c("stats", "out")]),
                 chi))
ylim4 <- range(c(unlist(boxplot(chi1[, c(5:8)], plot=FALSE)[c("stats", "out")]),
                 chi))
opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(c(rep(rep(c(1:2), each=2), 2), rep(rep(c(3:4), each=2), 2)), 
              nrow=4, ncol=4, byrow=TRUE))
boxplot(phi1[, c(1:4)], ylim=ylim1, main="(a)", ylab=expression(hat(phi)), 
        names=labels1, col=col[c(1:4)], cex=1.5, cex.lab=1.5, cex.axis=1.5)
points(c(1:4), phi, pch=2, col=col[c(1:4)], cex=1.5, 
       cex.lab=1.5, cex.axis=1.4)
boxplot(phi1[, c(5:8)], ylim=ylim2, main="(b)", ylab=expression(hat(phi)), 
        names=labels1, col=col[c(1:4)], cex=1.5, cex.lab=1.5, cex.axis=1.5)
points(c(1:4), phi, pch=2, col=col[c(1:4)], cex=1.5, 
       cex.lab=1.5, cex.axis=1.5)
boxplot(chi1[, c(1:4)], ylim=ylim3, main="(c)", ylab=expression(hat(chi)), 
        names=labels2, col=col[c(1:4)], cex=1.5, cex.lab=1.5, cex.axis=1.5)
points(c(1:4), chi, pch=2, col=col[c(1:4)], cex=1.5, 
       cex.lab=1.5, cex.axis=1.4)
boxplot(chi1[, c(5:8)], ylim=ylim4, main="(d)", ylab=expression(hat(chi)), 
        names=labels2, col=col[c(1:4)], cex=1.5, cex.lab=1.5, cex.axis=1.5)
points(c(1:4), phi, pch=2, col=col[c(1:4)], cex=1.5, 
       cex.lab=1.5, cex.axis=1.5)
par(opar)

# ---- data_gdsc_res1_brank_hist ----
library(sp)
load("results/data_gdsc_cv1.Rdata")
col <- bpy.colors(ncol(brankcor), cutoff.tail=0.3)
opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.1, 1, 1))
boxplot(brankcor, ylab=expression("Spearman's"~rho), col=col)
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



        