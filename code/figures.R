################################################################################
################################# main document ################################
################################################################################
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

# ---- simulations_gdsc_var ----
suppressWarnings(suppressMessages(library(sp)))
res <- read.table("results/simulations_gdsc_res3.txt", row.names=NULL)
temp <- res[, 1]
res <- as.matrix(res[, -1])
rownames(res) <- temp

alphaf <- c(1, 1, 3, 7)
alphad <- c(1, 1, 3, 7)
est.phi <- 1/(sapply(0:3, function(s) {
  res[rownames(res)==paste0("alphaf", s), 2]}) %*% 
    t(cbind(1, rbind(0, diag(3)))))
est.chi <- 1/(sapply(0:3, function(s) {
  res[rownames(res)==paste0("alphad", s), 2]}) %*% 
    t(cbind(1, rbind(0, diag(3)))))
x <- as.numeric(outer(1/(alphaf %*% t(cbind(1, rbind(0, diag(3))))),
                      1/(alphad %*% t(cbind(1, rbind(0, diag(3)))))))
y <- apply(sapply(1:nrow(est.phi), function(s) {
  as.numeric(outer(est.phi[s, ], est.chi[s, ]))}), 1, mean)

col <- bpy.colors(1, cutoff.tail=0.5)
lty <- 3
pch <- 1

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
plot(x, y, xlab=expression(V(beta)), 
     ylab=expression(hat(V)(beta)),
     xlim=range(c(x, y)), ylim=range(c(x, y)), col=col[1], pch=pch[1],
     cex=1.5, cex.lab=1.5, cex.axis=1.5, lwd=2)
abline(a=0, b=1, col=col[1], lty=lty[1], cex=1.5, cex.lab=1.5, cex.axis=1.5,
       lwd=2)
par(opar)

# ---- simulations_gdsc_kurt ----
suppressWarnings(suppressMessages(library(sp)))
res <- read.table("results/simulations_gdsc_res3.txt", row.names=NULL)
temp <- res[, 1]
res <- as.matrix(res[, -1])
rownames(res) <- temp

alphaf <- c(1, 1, 3, 7)
alphad <- c(1, 1, 3, 7)
est.phi <- 1/(sapply(0:3, function(s) {
  res[rownames(res)==paste0("alphaf", s), 2]}) %*% 
    t(cbind(1, rbind(0, diag(3)))))
est.chi <- 1/(sapply(0:3, function(s) {
  res[rownames(res)==paste0("alphad", s), 2]}) %*% 
    t(cbind(1, rbind(0, diag(3)))))
lambdaf <- res[rownames(res)=="lambdaf", 2]
lambdad <- res[rownames(res)=="lambdad", 2]
x <- 3*as.numeric(outer(
  as.numeric(1/(alphaf %*% t(cbind(1, rbind(0, diag(3)))))),
  as.numeric(1/(alphad %*% t(cbind(1, rbind(0, diag(3)))))), 
  FUN=function(X, Y) {X*Y + X + Y})) + 3
y <- rowMeans(sapply(1:nrow(est.phi), function(s) {
  3*as.numeric(outer(est.phi[s, ]/lambdaf[s], est.chi[s, ]/lambdad[s], 
                     function(X, Y) {X*Y + X + Y})) + 3}))

col <- bpy.colors(1, cutoff.tail=0.5)
lty <- 3
pch <- 1

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
plot(x, y, xlab=expression(K(beta)), 
     ylab=expression(hat(K)(beta)),
     xlim=range(c(x, y)), ylim=range(c(x, y)), col=col[1], pch=pch[1],
     cex=1.5, cex.lab=1.5, cex.axis=1.5, lwd=2)
abline(a=0, b=1, col=col[1], lty=lty[1], cex=1.5, cex.lab=1.5, cex.axis=1.5,
       lwd=2)
par(opar)

# ---- simulations_gdsc_est4 ----
suppressWarnings(suppressMessages(library(sp)))
res <- read.table("results/simulations_gdsc_res4.txt", row.names=NULL)
temp <- res[, 1]
res <- as.matrix(res[, -1])
rownames(res) <- temp

fracs <- c(1/10, 1/5, 1/3, 1/2, 2/3, 4/5, 1)
phi <- lapply(fracs, function(q) {
  cq <- paste0("frac", q, ".alphaf")
  apply(1/(sapply(0:3, function(s) {
    res[rownames(res)==paste0(cq, s), 2]}) %*% 
      t(cbind(1, rbind(0, diag(3))))), 2, mean)})
chi <- lapply(fracs, function(q) {
  cq <- paste0("frac", q, ".alphad")
  apply(1/(sapply(0:3, function(s) {
    res[rownames(res)==paste0(cq, s), 2]}) %*% 
      t(cbind(1, rbind(0, diag(3))))), 2, mean)})

suppressWarnings(suppressMessages(library(sp)))
res <- read.table("results/simulations_gdsc_res3.txt", row.names=NULL)
temp <- res[, 1]
res <- as.matrix(res[, -1])
rownames(res) <- temp

phi <- c(list(apply(1/(sapply(0:3, function(s) {
  res[rownames(res)==paste0("alphaf", s), 2]}) %*% 
    t(cbind(1, rbind(0, diag(3))))), 2, mean)), phi)
chi <- c(list(apply(1/(sapply(0:3, function(s) {
  res[rownames(res)==paste0("alphad", s), 2]}) %*% 
    t(cbind(1, rbind(0, diag(3))))), 2, mean)), chi)

lty <- c(1:(1 + length(fracs)))
col <- bpy.colors(1, cutoff.tail=0.3, 
                  alpha=seq(1, 0.4, length.out=length(fracs) + 1))
pch <- c(15:(15 + length(fracs)))
labels1 <- expression(phi["1"], phi["2"], phi["3"], phi["4"])
labels2 <- expression(chi["1"], chi["2"], chi["3"], chi["4"])

ylim1 <- range(unlist(phi))
ylim2 <- range(unlist(chi))

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(rep(c(1:2), each=2), 2), nrow=2, ncol=4, byrow=TRUE))
plot(phi[[1]], ylim=ylim1, main="(a)", ylab=expression(hat(phi)), xlab="", 
     type="l", col=col[1], cex=1.5, cex.lab=1.5, cex.axis=1.5, lty=lty[1], 
     xaxt="n")
axis(1, at=c(1:4), cex=1.5, cex.lab=1.5, cex.axis=1.5, labels=labels1)
for(q in 2:(length(fracs) + 1)) {
  lines(phi[[q]], col=col[q], lty=lty[q], cex=1.5, cex.lab=1.5, cex.axis=1.4)
}
plot(chi[[1]], ylim=ylim2, main="(b)", ylab=expression(hat(chi)), xlab="", 
     type="l", col=col[1], cex=1.5, cex.lab=1.5, cex.axis=1.5, lty=lty[1], 
     xaxt="n")
axis(1, at=c(1:4), cex=1.5, cex.lab=1.5, cex.axis=1.5, labels=labels2)
for(q in 2:(length(fracs) + 1)) {
  lines(chi[[q]], col=col[q], lty=lty[q], cex=1.5, cex.lab=1.5, cex.axis=1.4)
}
legend("topright", lty=lty, col=col, title="Permuted rows",
       legend=paste0(round(100*c(0, fracs), 0), "%"))
par(opar)

################################################################################
################################## supplement ##################################
################################################################################
# ---- simulations_gdsc_pmse4 ----
res <- read.table("results/simulations_gdsc_res3.txt", row.names=NULL)
temp <- res[, 1]
res <- as.matrix(res[, -1])
rownames(res) <- temp

D <- 251
nreps <- 100
fracs <- c(1/10, 1/5, 1/3, 1/2, 2/3, 4/5, 1)
measures <- c("emse", "emseh", "emsel", "pmse")

tabm <- aggregate(res, by=list(substr(rownames(res), 1, 5)), FUN="mean")
plotm <- data.frame(measure=c("emse", "emseh", "emsel", "pmse"), frac=0, 
                    tabm[tabm$Group.1 %in% 
                           c("emse.", "emseh", "emsel", "pmse."), -1])

tabsd <- sapply(c("emse.", "emseh", "emsel", "pmse."), function(s) {
  apply(aggregate(res[substr(rownames(res), 1, 5)==s, ], 
                  list(rep(1:nreps, each=D)), mean), 2, sd)})[-1, ]
plotsd <- data.frame(measure=c("emse", "emseh", "emsel", "pmse"),
                     frac=0, t(tabsd))

suppressWarnings(suppressMessages(library(sp)))
res <- read.table("results/simulations_gdsc_res4.txt", row.names=NULL)
temp <- res[, 1]
res <- as.matrix(res[, -1])
rownames(res) <- temp

combns <- apply(expand.grid(paste0("frac", fracs), paste0(measures, ".")), 1, 
                paste0, collapse=".")

plotm <- rbind(plotm, data.frame(
  measure=rep(measures, each=length(fracs)),
  frac=rep(fracs, length(measures)),
  t(sapply(combns, function(s) {
    colMeans(res[substr(rownames(res), 1, nchar(s))==s, ])}))))

plotsd <- rbind(plotsd, data.frame(
  measure=rep(measures, each=length(fracs)),
  frac=rep(fracs, length(measures)),
  t(sapply(combns, function(s) {
    apply(aggregate(res[substr(rownames(res), 1, nchar(s))==s, ], 
                    list(rep(1:nreps, each=D)), mean)[, -1], 2, sd)}))))

col <- bpy.colors(6, cutoff.tail=0.3)
labels1 <- c("NIG$_{\\text{f+d}}^-$", 
             "NIG$_{\\text{f+d}}$", "ridge", "lasso", "mxtune")

ylim1 <- range(plotm[plotm$measure=="emse", -c(1, 2)])
ylim2 <- range(plotm[plotm$measure=="emsel", -c(1, 2)])
ylim3 <- range(plotm[plotm$measure=="emseh", -c(1, 2)])
ylim4 <- range(plotm[plotm$measure=="pmse", -c(1, 2)])

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(1:4, nrow=2, byrow=TRUE))
plot(plotm$frac[plotm$measure=="emse"], plotm$NIG.[plotm$measure=="emse"], 
     ylim=ylim1, ylab="EMSE", xlab="Proportion of permuted rows", type="l",
     cex=1.5, cex.lab=1.5, cex.axis=1.5, col=col[1], lwd=1.5, main="(a)")
for(m in 2:(ncol(plotm) - 2)) {
  lines(plotm$frac[plotm$measure=="emse"], 
        plotm[plotm$measure=="emse", colnames(plotm)[m + 2]], 
        ylim=ylim1, lwd=1.5, cex=1.5, cex.lab=1.5, cex.axis=1.5, col=col[m])
}
plot(plotm$frac[plotm$measure=="emsel"], plotm$NIG.[plotm$measure=="emsel"], 
     ylim=ylim2, ylab=expression(EMSE[bottom]), 
     xlab="Proportion of permuted rows", type="l",
     cex=1.5, cex.lab=1.5, cex.axis=1.5, col=col[1], lwd=1.5, main="(b)")
for(m in 2:(ncol(plotm) - 2)) {
  lines(plotm$frac[plotm$measure=="emsel"], 
        plotm[plotm$measure=="emsel", colnames(plotm)[m + 2]], 
        lwd=1.5, cex=1.5, col=col[m])
}
plot(plotm$frac[plotm$measure=="emseh"], plotm$NIG.[plotm$measure=="emseh"], 
     ylim=ylim3, ylab=expression(EMSE[top]), xlab="Proportion of permuted rows", 
     type="l", cex=1.5, cex.lab=1.5, cex.axis=1.5, col=col[1], lwd=1.5, 
     main="(c)")
for(m in 2:(ncol(plotm) - 2)) {
  lines(plotm$frac[plotm$measure=="emseh"], 
        plotm[plotm$measure=="emseh", colnames(plotm)[m + 2]], 
        lwd=1.5, cex=1.5, col=col[m])
}
plot(plotm$frac[plotm$measure=="pmse"], plotm$NIG.[plotm$measure=="pmse"], 
     ylim=ylim4, ylab="PMSE", xlab="Proportion of permuted rows", type="l",
     cex=1.5, cex.lab=1.5, cex.axis=1.5, col=col[1], lwd=1.5, main="(d)")
for(m in 2:(ncol(plotm) - 2)) {
  lines(plotm$frac[plotm$measure=="pmse"], 
        plotm[plotm$measure=="pmse", colnames(plotm)[m + 2]], 
        lwd=1.5, cex=1.5, col=col[m])
}
legend("topright", legend=colnames(plotm)[-c(1, 2)], col=col, pch=16)
par(opar)

# ---- simulations_gdsc_post5 ----
load(file="results/simulations_gdsc_res5.Rdata")
D <- 100
p <- 100
H <- 4
G <- 4
opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(c(1:(G*H)), nrow=H, ncol=G, byrow=TRUE))
for(h in 1:H) {
  for(g in 1:G) {
    hi <- hist(post.mcmc1[[h]][[g]], breaks=40, plot=FALSE)
    m <- post.nig1$mu[g, h]
    s <- post.nig1$sigma[g, h]
    ylim <- range(c(hi$density, dnorm(m, m, s)))
    plot(hi, freq=FALSE, ylim=ylim, xlab=expression(beta), 
         ylab=expression(paste("p(", beta, "|", bold(y), ")")), 
         main=paste0("(", letters[(h - 1)*H + g], ")"),
         lwd=1.5, cex=1.5)
    curve(dnorm(x, m, s), add=TRUE, lwd=1.5, col=2)
  }
}
par(opar)

# ---- simulations_gdsc_res5 ----
load(file="results/simulations_gdsc_res5.Rdata")
D <- 100
p <- 100
H <- 4
G <- 4
alphaf <- c(1, 1, 3, 7)
alphad <- c(1, 1, 3, 7)
phi <- 1/as.numeric(cbind(1, rbind(0, diag(3))) %*% alphaf)
chi <- 1/as.numeric(cbind(1, rbind(0, diag(3))) %*% alphad)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(c(1:8), nrow=4, ncol=2, byrow=TRUE))
hist(1/(cbind(1, rbind(0, diag(3))) %*% post.mcmc2$alphaf)[1, ], freq=FALSE,
     main="(a)", xlab=expression(phi[1]), 
     ylab=expression(paste("p(", phi, "|", bold(y), ")")), breaks=40,
     lwd=1.5, cex=1.5, cex.lab=1.5, cex.axis=1.4)
abline(v=phi[1], col=2)
hist(1/(cbind(1, rbind(0, diag(3))) %*% post.mcmc2$alphad)[1, ], freq=FALSE,
     main="(b)", xlab=expression(chi[1]), 
     ylab=expression(paste("p(", chi, "|", bold(y), ")")), breaks=40,
     lwd=1.5, cex=1.5, cex.lab=1.5, cex.axis=1.4)
abline(v=chi[1], col=2)
hist(1/(cbind(1, rbind(0, diag(3))) %*% post.mcmc2$alphaf)[2, ], freq=FALSE,
     main="(c)", xlab=expression(phi[2]), 
     ylab=expression(paste("p(", phi, "|", bold(y), ")")), breaks=40,
     lwd=1.5, cex=1.5, cex.lab=1.5, cex.axis=1.4)
abline(v=phi[2], col=2)
hist(1/(cbind(1, rbind(0, diag(3))) %*% post.mcmc2$alphad)[2, ], freq=FALSE,
     main="(d)", xlab=expression(chi[2]), 
     ylab=expression(paste("p(", chi, "|", bold(y), ")")), breaks=40,
     lwd=1.5, cex=1.5, cex.lab=1.5, cex.axis=1.4)
abline(v=chi[2], col=2)
hist(1/(cbind(1, rbind(0, diag(3))) %*% post.mcmc2$alphaf)[3, ], freq=FALSE,
     main="(e)", xlab=expression(phi[3]), 
     ylab=expression(paste("p(", phi, "|", bold(y), ")")), breaks=40,
     lwd=1.5, cex=1.5, cex.lab=1.5, cex.axis=1.4)
abline(v=phi[3], col=2)
hist(1/(cbind(1, rbind(0, diag(3))) %*% post.mcmc2$alphad)[3, ], freq=FALSE,
     main="(f)", xlab=expression(chi[3]), 
     ylab=expression(paste("p(", chi, "|", bold(y), ")")), breaks=40,
     lwd=1.5, cex=1.5, cex.lab=1.5, cex.axis=1.4)
abline(v=chi[3], col=2)
hist(1/(cbind(1, rbind(0, diag(3))) %*% post.mcmc2$alphaf)[4, ], freq=FALSE,
     main="(g)", xlab=expression(phi[4]), 
     ylab=expression(paste("p(", phi, "|", bold(y), ")")), breaks=40,
     lwd=1.5, cex=1.5, cex.lab=1.5, cex.axis=1.4)
abline(v=phi[4], col=2)
hist(1/(cbind(1, rbind(0, diag(3))) %*% post.mcmc2$alphad)[4, ], freq=FALSE,
     main="(h)", xlab=expression(chi[4]), 
     ylab=expression(paste("p(", chi, "|", bold(y), ")")), breaks=40,
     lwd=1.5, cex=1.5, cex.lab=1.5, cex.axis=1.4)
abline(v=chi[4], col=2)
par(opar)

# ---- analysis_gdsc_cpo ----
load(file="data/data_gdsc_dat1.Rdata")
load("results/analysis_gdsc_cpo1.Rdata")
hist(exp(unlist(lcpo.nig2)), freq=FALSE, breaks=100)
abline(v=0.01)

################################################################################
################################# presentation #################################
################################################################################
# ---- dens_beta_prior1 ---- 
suppressWarnings(suppressMessages(library(GeneralizedHyperbolic)))
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