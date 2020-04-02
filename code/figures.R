################################# main document ################################
# ---- figures ----
# ---- dens_kappa ----
suppressWarnings(suppressMessages(library(sp)))
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
     lty=lty[2])
lines(x.seq, dkstudentst(x.seq, a=1, b=4), col=col[2], lty=lty[1])
lines(x.seq, dkstudentst(x.seq, a=1, b=0.2), col=col[2], lty=lty[3])

plot(x.seq, dklasso(x.seq, 1), type="l", ylim=c(0, 3.5), 
     xlab=expression(kappa), ylab=expression(pi(kappa)), main="(c)", col=col[3], 
     lty=lty[1])
lines(x.seq, dklasso(x.seq, 0.2), col=col[3], lty=lty[2])
lines(x.seq, dklasso(x.seq, 0.1), col=col[3], lty=lty[3])
par(opar)

# ---- simulations_gdsc_est1 ----
suppressWarnings(suppressMessages(library(sp)))
res <- read.table("results/simulations_gdsc_res1.2.txt", row.names=NULL)
temp <- res[, 1]
res <- as.matrix(res[, -1])
rownames(res) <- temp

pmse <- res[substr(rownames(res), 1, 5)=="pmse.", ]
mpmse <- aggregate(pmse, by=list(substr(rownames(pmse), 6, 1000000L)), 
                   FUN="mean")[, -1]
boxplot(mpmse)


est <- sapply(c(paste0("alphaf", 0:3), paste0("alphad", 0:3), 
                  "lambdaf", "lambdad"), function(s) {
                    s <- res[rownames(res)==s, ]; rownames(s) <- NULL
                    s}, simplify=FALSE)

alpha <- sapply(est, function(s) {s[, 2]})
boxplot(alpha[, -c(5, 9, 10)])
abline(h=1)
test <- alpha$alphad0
boxplot(cbind(test[, c(1, 2)], exp(-test[, 5])), outline=FALSE)

apply(alphaf0, 2, median, na.rm=TRUE)
alphaf0

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
suppressWarnings(suppressMessages(library(sp)))
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
suppressWarnings(suppressMessages(library(sp)))
res <- read.table("results/simulations_gdsc_res3.txt", row.names=NULL)
temp <- res[, 1]
res <- as.matrix(res[, -1])
rownames(res) <- temp

alphaf <- c(1, 1, 3, 7)*10
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

# ---- simulations_gdsc_est4 ----
suppressWarnings(suppressMessages(library(sp)))
res <- read.table("results/simulations_gdsc_res4.txt", row.names=NULL)
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

# ---- analysis_gdsc_res1 ----
suppressWarnings(suppressMessages(library(sp)))
res <- read.table("results/analysis_gdsc_res1.txt", row.names=NULL)
temp <- res[, 1]
res <- as.matrix(res[, -1])
rownames(res) <- temp
pmse <- res[substr(rownames(res), 1, 4)=="pmse", ]

ylim1 <- range(1 - c(res[, "semnig1"], res[, "ridge1"]))
ylim2 <- range(1 - c(res[, "ebridge1"], res[, "ridge1"]))

col <- bpy.colors(3, cutoff.tail=0.2)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1:2), each=2), nrow=2, ncol=2, byrow=TRUE))
plot(1 - pmse[, "ridge1"], xlab="Drug", 
     ylab="MSE reduction", main="(a)",
     cex=0.5, pch=16, ylim=ylim1, col=col[1])
points(1 - pmse[, "semnig1"], cex=0.5, pch=16, col=col[2])
abline(v=order(abs(pmse[, "ridge1"] - pmse[, "semnig1"]), 
               decreasing=TRUE)[c(1:5)], lty=2, 
       col=col[as.numeric((pmse[, "ridge1"] - pmse[, "semnig1"] > 0) + 1)[
         order(abs(pmse[, "ridge1"] - pmse[, "semnig1"]), 
               decreasing=TRUE)[c(1:5)]]])
abline(h=0, lty=3)
legend("topleft", c("ridge + CV", "NIG"), pch=16, col=col[c(1, 2)])

plot(1 - pmse[, "ridge1"], xlab="Drug", 
     ylab="MSE reduction", main="(b)",
     cex=0.5, pch=16, ylim=ylim2, col=col[1])
points(1 - pmse[, "ebridge1"], cex=0.5, pch=16, col=col[3])
abline(v=order(abs(pmse[, "ridge1"] - pmse[, "ebridge1"]), 
               decreasing=TRUE)[c(1:5)], lty=2, 
       col=col[as.numeric((pmse[, "ridge1"] - pmse[, "ebridge1"] > 0)*2 + 1)[
         order(abs(pmse[, "ridge1"] - pmse[, "ebridge1"]), 
               decreasing=TRUE)[c(1:5)]]])
abline(h=0, lty=3)
legend("topleft", c("ridge + CV", "ridge + EB"), pch=16, col=col[c(1, 3)])

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


