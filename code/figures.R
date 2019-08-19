################################# main document ################################
# ---- figures ----
# ---- dens_igaussian_marginalbeta ----
library(GeneralizedHyperbolic)
library(metRology)
library(sp)
sigma <- 1
mgammasq <- 1/sigma^2
lambda <- c(0.1, 2)
lambda1 <- sqrt(2)
eta <- lambda/mgammasq + 2
labels <- as.expression(c(bquote("NIG, "~lambda==.(lambda[1])*","~
                                   (bold(c)^"T"*bold(alpha))^-1==.(mgammasq)),
                          bquote("NIG, "~lambda==.(lambda[2])*","~
                                   (bold(c)^"T"*bold(alpha))^-1==.(mgammasq)),
                          bquote("Student's t, "~xi==.(lambda[1])*","~
                                   eta==.(eta[1])),
                          bquote("Student's t, "~xi==.(lambda[2])*","~
                                   eta==.(eta[2])),
                          "ridge", "lasso"))
dprior <- function(x, lambda, mgammasq, sigma) {
  dnig(x, 0, sigma*sqrt(lambda), sqrt(lambda/(mgammasq*sigma)), 0)
}
dlasso <- function(x, lambda1) {
  0.5*lambda1*exp(-lambda1*abs(x))
}

col <- bpy.colors(length(labels), cutoff.tail=0.3)
lty <- c(1:length(labels))

ylim <- c(0, dt.scaled(0, eta[1]/2, 0, sqrt(lambda[1]*sigma^2/eta[1])))
opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.1, 1, 1))
curve(dprior(x, lambda[1], mgammasq, sigma), -3, 3, 
      ylab=expression(p(beta)), 
      xlab=expression(beta), n=1000, ylim=ylim,
      col=col[1], lty=lty[1])
curve(dprior(x, lambda[2], mgammasq, sigma), add=TRUE, n=1000, 
      col=col[2], lty=lty[2])
curve(dt.scaled(x, eta[1]/2, 0, sqrt(lambda[1]*sigma^2/eta[1])), add=TRUE, 
      n=1000, col=col[3], lty=lty[3])
curve(dt.scaled(x, eta[2]/2, 0, sqrt(lambda[2]*sigma^2/eta[2])), add=TRUE, 
      n=1000, col=col[4], lty=lty[4])
curve(dnorm(x, 0, 1), add=TRUE, n=1000, col=col[5], lty=lty[5])
curve(dlasso(x, lambda1), add=TRUE, n=1000, col=col[6], lty=lty[6])
legend("topright", legend=labels, col=col, lty=lty, 
       title="Prior", seg.len=1)
par(opar)

# ---- boxplots_igaussian_res4.1_post1 ----  
library(sp)
res <- read.table("results/simulations_igaussian_res4.1.csv")
set <- read.table("results/simulations_igaussian_set4.csv")

plot.data <- sapply(1:nrow(set), function(s) {
    list(res[, substr(colnames(res), 1, 9)==paste0("set", s, ".cram") &
               substr(colnames(res), nchar(colnames(res)) - 3, 
                      nchar(colnames(res))) %in% c("enig", "ntst")],
         res[, substr(colnames(res), 1, 12)==paste0("set", s, ".corbest") &
               substr(colnames(res), nchar(colnames(res)) - 3, 
                      nchar(colnames(res))) %in% c("enig", "ntst")],
         res[, substr(colnames(res), 1, 12)==paste0("set", s, ".msebest") &
               substr(colnames(res), nchar(colnames(res)) - 3, 
                      nchar(colnames(res))) %in% c("enig", "ntst")],
         res[, substr(colnames(res), 1, 10)==paste0("set", s, ".cover") &
               substr(colnames(res), nchar(colnames(res)) - 3, 
                      nchar(colnames(res))) %in% c("enig", "ntst")])},
    simplify=FALSE)

methods <- c("NIG", "Student's t")
labels1 <- c("NIG", "Student's t")
labels2 <- c("Cramér-von Mises criterion", expression("Corr"(hat(beta), beta)), 
             expression("MSD"(hat(beta), beta)), "95% coverage")
col <- bpy.colors(length(methods), cutoff.tail=0.3)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.1, 1, 1))
layout(matrix(c(rep(rep(c(1:4), each=2), 2),
                rep(rep(c(5:8), each=2), 2),
                rep(rep(c(9:12), each=2), 2)), nrow=6, ncol=8, byrow=TRUE))
for(r in 1:3) {
  for(m in 1:4) {
    boxplot(plot.data[[r]][[m]], main=paste0("(", letters[(r - 1)*4 + m], ")"),
            ylab=labels2[m], las=2, names=labels1, col=col, outline=FALSE)
    if(m==4) {abline(h=0.95)}
  }
}
par(opar)

# ---- boxplots_igaussian_res4.1_post2 ----  
library(sp)
res <- read.table("results/simulations_igaussian_res4.1.csv")
set <- read.table("results/simulations_igaussian_set4.csv")

plot.data <- sapply(1:nrow(set), function(s) {
  list(res[, substr(colnames(res), 1, 9)==paste0("set", s, ".cram") &
             substr(colnames(res), nchar(colnames(res)) - 3, 
                    nchar(colnames(res))) %in% c("enig", "ntst")],
       res[, substr(colnames(res), 1, 12)==paste0("set", s, ".corbest") &
             substr(colnames(res), nchar(colnames(res)) - 3, 
                    nchar(colnames(res))) %in% c("enig", "ntst")],
       res[, substr(colnames(res), 1, 12)==paste0("set", s, ".msebest") &
             substr(colnames(res), nchar(colnames(res)) - 3, 
                    nchar(colnames(res))) %in% c("enig", "ntst")],
       res[, substr(colnames(res), 1, 10)==paste0("set", s, ".cover") &
             substr(colnames(res), nchar(colnames(res)) - 3, 
                    nchar(colnames(res))) %in% c("enig", "ntst")])},
  simplify=FALSE)

methods <- c("NIG", "Student's t")
labels1 <- c("NIG", "Student's t")
labels2 <- c("Cramér-von Mises criterion", expression("Corr"(hat(beta), beta)), 
             expression("MSD"(hat(beta), beta)), "95% coverage")
col <- bpy.colors(length(methods), cutoff.tail=0.3)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.1, 1, 1))
layout(matrix(c(rep(rep(c(1:4), each=2), 2),
                rep(rep(c(5:8), each=2), 2),
                rep(rep(c(9:12), each=2), 2)), nrow=6, ncol=8, byrow=TRUE))
for(r in 4:6) {
  for(m in 1:4) {
    boxplot(plot.data[[r]][[m]], main=paste0("(", letters[(r - 4)*4 + m], ")"),
            ylab=labels2[m], las=2, names=labels1, col=col, outline=FALSE)
    if(m==4) {abline(h=0.95)}
  }
}
par(opar)

# ---- boxplots_igaussian_res5_post1 ----  
library(sp)
res <- read.table("results/simulations_igaussian_res5.csv")
set <- read.table("results/simulations_igaussian_set5.csv")

plot.data <- sapply(1:nrow(set), function(s) {
  list(res[, substr(colnames(res), 1, 9)==paste0("set", s, ".cram") &
             substr(colnames(res), nchar(colnames(res)) - 3, 
                    nchar(colnames(res))) %in% c("enig", "ntst")],
       res[, substr(colnames(res), 1, 12)==paste0("set", s, ".corbest") &
             substr(colnames(res), nchar(colnames(res)) - 3, 
                    nchar(colnames(res))) %in% c("enig", "ntst")],
       res[, substr(colnames(res), 1, 12)==paste0("set", s, ".msebest") &
             substr(colnames(res), nchar(colnames(res)) - 3, 
                    nchar(colnames(res))) %in% c("enig", "ntst")],
       res[, substr(colnames(res), 1, 10)==paste0("set", s, ".cover") &
             substr(colnames(res), nchar(colnames(res)) - 3, 
                    nchar(colnames(res))) %in% c("enig", "ntst")])},
  simplify=FALSE)

methods <- c("NIG", "Student's t")
labels1 <- c("NIG", "Student's t")
labels2 <- c("Cramér-von Mises criterion", expression("Corr"(hat(beta), beta)), 
             expression("MSD"(hat(beta), beta)), "95% coverage")
col <- bpy.colors(length(methods), cutoff.tail=0.3)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.1, 1, 1))
layout(matrix(c(rep(rep(c(1:4), each=2), 2),
                rep(rep(c(5:8), each=2), 2),
                rep(rep(c(9:12), each=2), 2)), nrow=6, ncol=8, byrow=TRUE))
for(r in 1:3) {
  for(m in 1:4) {
    boxplot(plot.data[[r]][[m]], main=paste0("(", letters[(r - 1)*4 + m], ")"),
            ylab=labels2[m], las=2, names=labels1, col=col, outline=FALSE)
    if(m==4) {abline(h=0.95)}
  }
}
par(opar)

# ---- boxplots_igaussian_res5_post2 ----  
library(sp)
res <- read.table("results/simulations_igaussian_res5.csv")
set <- read.table("results/simulations_igaussian_set5.csv")

plot.data <- sapply(1:nrow(set), function(s) {
  list(res[, substr(colnames(res), 1, 9)==paste0("set", s, ".cram") &
             substr(colnames(res), nchar(colnames(res)) - 3, 
                    nchar(colnames(res))) %in% c("enig", "ntst")],
       res[, substr(colnames(res), 1, 12)==paste0("set", s, ".corbest") &
             substr(colnames(res), nchar(colnames(res)) - 3, 
                    nchar(colnames(res))) %in% c("enig", "ntst")],
       res[, substr(colnames(res), 1, 12)==paste0("set", s, ".msebest") &
             substr(colnames(res), nchar(colnames(res)) - 3, 
                    nchar(colnames(res))) %in% c("enig", "ntst")],
       res[, substr(colnames(res), 1, 10)==paste0("set", s, ".cover") &
             substr(colnames(res), nchar(colnames(res)) - 3, 
                    nchar(colnames(res))) %in% c("enig", "ntst")])},
  simplify=FALSE)

methods <- c("NIG", "Student's t")
labels1 <- c("NIG", "Student's t")
labels2 <- c("Cramér-von Mises criterion", expression("Corr"(hat(beta), beta)), 
             expression("MSD"(hat(beta), beta)), "95% coverage")
col <- bpy.colors(length(methods), cutoff.tail=0.3)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.1, 1, 1))
layout(matrix(c(rep(rep(c(1:4), each=2), 2),
                rep(rep(c(5:8), each=2), 2),
                rep(rep(c(9:12), each=2), 2)), nrow=6, ncol=8, byrow=TRUE))
for(r in 4:6) {
  for(m in 1:4) {
    boxplot(plot.data[[r]][[m]], main=paste0("(", letters[(r - 4)*4 + m], ")"),
            ylab=labels2[m], las=2, names=labels1, col=col, outline=FALSE)
    if(m==4) {abline(h=0.95)}
  }
}
par(opar)

# ---- boxplots_igaussian_res4.2_post1 ----  
library(sp)
res <- read.table("results/simulations_igaussian_res4.2.csv")
set <- read.table("results/simulations_igaussian_set4.csv")

plot.data <- sapply(1:nrow(set), function(s) {
  list(res[, substr(colnames(res), 1, 9)==paste0("set", s, ".cram") &
             substr(colnames(res), nchar(colnames(res)) - 1, 
                    nchar(colnames(res)))!="et"],
       res[, substr(colnames(res), 1, 12)==paste0("set", s, ".corbest") &
             substr(colnames(res), nchar(colnames(res)) - 1, 
                    nchar(colnames(res)))!="et"],
       res[, substr(colnames(res), 1, 12)==paste0("set", s, ".msebest") &
             substr(colnames(res), nchar(colnames(res)) - 1, 
                    nchar(colnames(res)))!="et"])},
  simplify=FALSE)

methods <- c("ADVI", "VB", "MAP")
labels1.1 <- c("ADVI", "VB")
labels1.2 <- c("ADVI", "VB", "MAP")
labels2 <- c("Cramér-von Mises criterion", expression("Corr"(hat(beta), beta)), 
             expression("MSD"(hat(beta), beta)))
col <- bpy.colors(length(methods), cutoff.tail=0.3)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.1, 1, 1))
layout(matrix(c(rep(rep(c(1:3), each=2), 2),
                rep(rep(c(4:6), each=2), 2),
                rep(rep(c(7:9), each=2), 2)), nrow=6, ncol=6, byrow=TRUE))
for(r in 1:3) {
  for(m in 1:3) {
    if(m==1) {lab1 <- labels1.1} else {lab1 <- labels1.2}
    boxplot(plot.data[[r]][[m]], main=paste0("(", letters[(r - 1)*3 + m], ")"),
            ylab=labels2[m], las=2, names=lab1, col=col, outline=FALSE)
  }
}
par(opar)

# ---- boxplots_igaussian_res4.2_post2 ----  
library(sp)
res <- read.table("results/simulations_igaussian_res4.2.csv")
set <- read.table("results/simulations_igaussian_set4.csv")

plot.data <- sapply(1:nrow(set), function(s) {
  list(res[, substr(colnames(res), 1, 9)==paste0("set", s, ".cram") &
             substr(colnames(res), nchar(colnames(res)) - 1, 
                    nchar(colnames(res)))!="et"],
       res[, substr(colnames(res), 1, 12)==paste0("set", s, ".corbest") &
             substr(colnames(res), nchar(colnames(res)) - 1, 
                    nchar(colnames(res)))!="et"],
       res[, substr(colnames(res), 1, 12)==paste0("set", s, ".msebest") &
             substr(colnames(res), nchar(colnames(res)) - 1, 
                    nchar(colnames(res)))!="et"])},
  simplify=FALSE)

methods <- c("ADVI", "VB", "MAP")
labels1.1 <- c("ADVI", "VB")
labels1.2 <- c("ADVI", "VB", "MAP")
labels2 <- c("Cramér-von Mises criterion", expression("Corr"(hat(beta), beta)), 
             expression("MSD"(hat(beta), beta)))
col <- bpy.colors(length(methods), cutoff.tail=0.3)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.1, 1, 1))
layout(matrix(c(rep(rep(c(1:3), each=2), 2),
                rep(rep(c(4:6), each=2), 2),
                rep(rep(c(7:9), each=2), 2)), nrow=6, ncol=6, byrow=TRUE))
for(r in 4:6) {
  for(m in 1:3) {
    if(m==1) {lab1 <- labels1.1} else {lab1 <- labels1.2}
    boxplot(plot.data[[r]][[m]], main=paste0("(", letters[(r - 4)*3 + m], ")"),
            ylab=labels2[m], las=2, names=lab1, col=col, outline=FALSE)
  }
}
par(opar)

# ---- hist_igaussian_fit4 ----
library(sp)
load("results/simulations_igaussian_fit4.Rdata")
set <- read.table("results/simulations_igaussian_set4.csv")

beta <- cbind(rstan::extract(fit.enig, pars="beta0")$beta0, 
              rstan::extract(fit.enig, pars="beta")$beta)

labels1 <- c("VB", "ADVI", "MAP")
labels2 <- c("density", "point estimate")
lty <- c(1, 2)
methods <- labels1
col <- bpy.colors(length(methods), cutoff.tail=0.3)

set.seed(2019)
varid <- sample(1:(ncol(beta) - 1), 8)
opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.1, 1, 1), lwd=2)
layout(matrix(c(rep(rep(c(1:3), each=2), 2),
                rep(rep(c(4:6), each=2), 2),
                rep(rep(c(7:9), each=2), 2)), nrow=6, ncol=6, byrow=TRUE))
for(id in 1:length(varid)) {
  j <- varid[id]
  h <- hist(beta[, j + 1], plot=FALSE)
  hist(beta[, j + 1], breaks=30, prob=TRUE, 
       ylim=range(c(dnorm(fit.enig.vb$vb$mpost$beta[[1]][j], 
                          fit.enig.vb$vb$mpost$beta[[1]][j],
                          sqrt(fit.enig.vb$vb$vpost$beta[[1]][j, j])), 
                    h$density)),
       xlab=expression(beta), main=bquote(beta[.(j)]))
  curve(dnorm(x, fit.enig.vb$vb$mpost$beta[[1]][j],
              sqrt(fit.enig.vb$vb$vpost$beta[[1]][j, j])), add=TRUE, col=col[1])
  curve(dnorm(x, summary(fit.enig.advi)$summary[substr(rownames(summary(
    fit.enig.advi)$summary), 1, 4)=="beta" & substr(rownames(summary(
      fit.enig.advi)$summary), 5, 6)!="sd", 1][j + 1], 
    summary(fit.enig.advi)$summary[substr(rownames(summary(
      fit.enig.advi)$summary), 1, 4)=="beta" & substr(rownames(summary(
        fit.enig.advi)$summary), 5, 6)!="sd", 2][j + 1]^2),
    add=TRUE, col=col[2])
  abline(v=summary(fit.enig.advi)$summary[substr(rownames(summary(
    fit.enig.advi)$summary), 1, 4)=="beta" & substr(rownames(summary(
      fit.enig.advi)$summary), 5, 6)!="sd", 1][j + 1], lty=2, col=col[2])
  abline(v=fit.enig.map$par[substr(names(fit.enig.map$par), 1, 4)=="beta"][j + 1], 
         lty=2, col=col[3])
}
plot.new()
legend("center", legend=c(labels1, labels2), fill=c(col, NA, NA), 
       border=c(rep(1, 3), NA, NA), lty=c(rep(NA, 3), lty),
       seg.len=1, merge=TRUE, bg="white")
par(opar)



# ---- boxplots_igaussian_res3_vb ----  
library(sp)
res <- read.table("results/simulations_igaussian_res3.csv")

plot.data1 <- res[, substr(colnames(res), 1, 4)=="cram"]
plot.data2 <- res[, substr(colnames(res), 1, 7)=="corbvar"]
plot.data3 <- res[, substr(colnames(res), 1, 7)=="msebvar"]
methods <- c("ADVI", "VB")
labels <- c("ADVI", "VB")
col <- bpy.colors(length(methods), cutoff.tail=0.3)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2, 3, 3), 2), nrow=2, ncol=6, byrow=TRUE))
boxplot(plot.data1, names=labels, col=col, outline=FALSE, main="(a)",
        ylab="Cramér-von Mises criterion")
boxplot(plot.data2, names=labels, col=col, outline=FALSE, main="(b)",
        ylab=expression("Corr"(hat(V)(beta),hat(V)(beta["MCMC"]))))
boxplot(plot.data3, names=labels, col=col, outline=FALSE, main="(c)",
        ylab=expression("MSD"(hat(V)(beta),hat(V)(beta["MCMC"]))))
par(opar)

# ---- boxplots_igaussian_res3_best ----  
library(sp)
res <- read.table("results/simulations_igaussian_res3.csv")

plot.data1 <- res[, substr(colnames(res), 1, 7)=="corbest"]
plot.data2 <- res[, substr(colnames(res), 1, 7)=="msebest"]
methods <- c("ADVI", "VB", "glmnet", "MAP")
col <- bpy.colors(length(methods), cutoff.tail=0.3)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2), 2), nrow=2, ncol=4, byrow=TRUE))
boxplot(plot.data1, names=labels, col=col, outline=FALSE, main="(a)",
        ylab=expression("Corr"(hat(beta),hat(beta["MCMC"]))))
boxplot(plot.data2, names=labels, col=col, outline=FALSE, main="(b)",
        ylab=expression("MSD"(hat(beta),hat(beta["MCMC"]))))
par(opar)

# ---- boxplots_igaussian_res6_eb ----  
library(sp)
res <- read.table("results/simulations_igaussian_res6.csv")
set <- read.table("results/simulations_igaussian_set6.csv")

plot.data <- res[, substr(colnames(res), 1, 5)=="alpha"]
alpha <- set[1, substr(colnames(set), 1, 5)=="alpha"]

labels <- expression(alpha[0], alpha[1], alpha[2], alpha[3], alpha[4])
col <- bpy.colors(length(labels), cutoff.tail=0.3)

boxplot(plot.data, names=labels, col=col, ylab=expression(hat(alpha)))
lines(c(0.5, 1.5), rep(alpha[1], 2), lty=2, lwd=2, col=col[1])
lines(c(1.5, 2.5), rep(alpha[2], 2), lty=2, lwd=2, col=col[2])
lines(c(2.5, 3.5), rep(alpha[3], 2), lty=2, lwd=2, col=col[3])
lines(c(3.5, 4.5), rep(alpha[4], 2), lty=2, lwd=2, col=col[4])
lines(c(4.5, 5.5), rep(alpha[5], 2), lty=2, lwd=2, col=col[5])
legend("topleft", lty=2, col=1, legend=expression("true"~alpha))

# ---- lines_eqtl_hsq ----
library(sp)
res <- read.table("results/eQTL_hsq.csv")

plot.data1 <- list(intercept=colMeans(res[, substr(colnames(res), 1, 
                                                   6)=="interc"]),
                   model1=colMeans(res[, substr(colnames(res), 1, 
                                                6)=="model1"]),
                   model2=colMeans(res[, substr(colnames(res), 1, 
                                                6)=="model2"]),
                   model3=colMeans(res[, substr(colnames(res), 1, 
                                                6)=="model3"]),
                   model4=colMeans(res[, substr(colnames(res), 1, 
                                                6)=="model4"]))
plot.data2 <- list(intercept=apply(res[, substr(colnames(res), 1, 
                                                   6)=="interc"], 2, sd),
                   model1=apply(res[, substr(colnames(res), 1, 
                                             6)=="model1"], 2, sd),
                   model2=apply(res[, substr(colnames(res), 1, 
                                             6)=="model2"], 2, sd),
                   model3=apply(res[, substr(colnames(res), 1, 
                                             6)=="model3"], 2, sd),
                   model4=apply(res[, substr(colnames(res), 1, 
                                             6)=="model4"], 2, sd))
psel <- c(1, 5, 10, 20, 50)

labels <- c("intercept only", "model 1", "model 2", "model 3", "model 4")
col <- bpy.colors(length(plot.data), cutoff.tail=0.3)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.1, 1, 1))
plot(psel, plot.data1[[1]], pch=19, ylab="mean heritability", 
     xlab="number of selected SNPs",
     ylim=range(c(unlist(plot.data1) - unlist(plot.data2),
                  unlist(plot.data1) + unlist(plot.data2))), col=col[1])
arrows(psel, plot.data1[[1]] - plot.data2[[1]], psel, 
       plot.data1[[1]] + plot.data2[[1]], length=0.05, angle=90, code=3, 
       col=col[1])
points(psel, plot.data1[[2]], col=col[2], pch=19)
arrows(psel, plot.data1[[2]] - plot.data2[[2]], psel, 
       plot.data1[[2]] + plot.data2[[2]], length=0.05, angle=90, code=3, 
       col=col[2])
points(psel, plot.data1[[3]], col=col[3], pch=19)
arrows(psel, plot.data1[[3]] - plot.data2[[3]], psel, 
       plot.data1[[3]] + plot.data2[[3]], length=0.05, angle=90, code=3, 
       col=col[3])
points(psel, plot.data1[[4]], col=col[4], pch=19)
arrows(psel, plot.data1[[4]] - plot.data2[[4]], psel, 
       plot.data1[[4]] + plot.data2[[4]], length=0.05, angle=90, code=3, 
       col=col[4])
points(psel, plot.data1[[5]], col=col[5], pch=19)
arrows(psel, plot.data1[[5]] - plot.data2[[5]], psel, 
       plot.data1[[5]] + plot.data2[[5]], length=0.05, angle=90, code=3, 
       col=col[5])
legend("bottomright", lty=1, legend=labels, col=col)
par(opar)

# # ---- boxplot_igaussian_res1_prior_mean ----  
# res1 <- as.matrix(read.table("results/simulations_igaussian_res1.csv"))
# set1 <- read.table("results/simulations_igaussian_set1.csv")
# labels1 <- c("conj. inv. \n Gaussian", "non-conj. inv. \n Gaussian",
#               "conj. inv. \n Gamma", "non-conj. inv. \n Gamma")
# colors1 <- sp::bpy.colors(6)[-c(1, 6)]
# 
# C <- matrix(sapply(paste("set1$C", 1:(set1$nclass*set1$D), sep=""), 
#                    function(s) {eval(parse(text=s))}), 
#             ncol=set1$nclass, nrow=set1$D)
# prior.means <- as.numeric((sapply(paste("set1$theta", 1:set1$D, sep=""),
#                                   function(s) {eval(parse(text=s))}) %*% C)/
#                             (set1$D/set1$nclass))
# temp1 <- data.frame(class=factor(rep(rep(c(1:set1$nclass), each=set1$D), 4)),
#                     method=rep(paste(1:4, labels1), each=set1$nclass*set1$D),
#                     prior.means=c((res1[, grepl("igauss.conj_theta", colnames(
#                       res1))] %*% C)/(set1$D/set1$nclass), 
#                       (res1[, grepl("igauss.non.conj_theta", colnames(res1))]
#                        %*% C)/(set1$D/set1$nclass), 
#                       ((res1[, grepl("igamma.conj_lambda", colnames(res1))]/
#                           (res1[, grepl("igamma.conj_eta", colnames(
#                             res1))] - 2)) %*% C)/(set1$D/set1$nclass), 
#                       ((res1[, grepl("igamma.non.conj_lambda", colnames(res1))]/
#                           (res1[, grepl("igamma.non.conj_eta", colnames(
#                             res1))] - 2)) %*% C)/(set1$D/set1$nclass)))
#   
# opar <- par(no.readonly=TRUE)
# at1 <- c(1:set1$nclass, 1:set1$nclass + set1$nclass + 1, 1:set1$nclass +
#            2*set1$nclass + 2, 1:set1$nclass + 3*set1$nclass + 3)
# at2 <- c(mean(1:set1$nclass), mean(1:set1$nclass + set1$nclass + 1),
#          mean(1:set1$nclass + 2*set1$nclass + 2),
#          mean(1:set1$nclass + 3*set1$nclass + 3))
# par(cex=1.2, mar=opar$mar*c(1, 1.3, 1/1.3, 1))
# boxplot(prior.means ~ class + method, data=temp1, outline=FALSE,
#         at=at1, col=colors1, names=NA, xaxt="n",
#         ylab=expression(hat(E)(gamma[d]^2)))
# abline(h=prior.means, col=colors1, lty=2, lwd=1.5)
# # text(at2, par("usr")[3] - 0.2, labels1, srt=45, pos=1, xpd=TRUE)
# axis(1, at=at2, labels=labels1, tick=FALSE)
# legend("topright", legend=c(paste("class", 1:set1$nclass), "true value"), 
#        lty=c(rep(NA, 4), 2), col=c(colors1, 1), seg.len=rep(1, 5), 
#        border=c(rep(1, 4), NA), fill=c(colors1, NA), merge=TRUE)
# par(opar)
# 
# # ---- boxplot_igaussian_res2_prior_mean ----  
# res2 <- as.matrix(read.table("results/simulations_igaussian_res2.csv"))
# set2 <- read.table("results/simulations_igaussian_set2.csv")
# labels2 <- c("conj. inv. \n Gaussian", "non-conj. inv. \n Gaussian",
#              "conj. inv. \n Gamma", "non-conj. inv. \n Gamma")
# colors2 <- sp::bpy.colors(6)[-c(1, 6)]
# 
# C <- matrix(sapply(paste("set2$C", 1:(set2$nclass*set2$D), sep=""), 
#                    function(s) {eval(parse(text=s))}), 
#             ncol=set2$nclass, nrow=set2$D)
# prior.means <- as.numeric((sapply(paste(
#   "set2$lambda", 1:set2$D, sep=""), function(s) {eval(parse(text=s))})/
#     (sapply(paste("set2$eta", 1:set2$D, sep=""), function(s) {
#       eval(parse(text=s))}) - 2)) %*% C)/(set1$D/set1$nclass)
# temp1 <- data.frame(class=factor(rep(rep(c(1:set2$nclass), each=set2$D), 4)),
#                     method=rep(paste(1:4, labels2), each=set2$nclass*set2$D),
#                     prior.means=c((res2[, grepl("igauss.conj_theta", colnames(
#                       res2))] %*% C)/(set2$D/set2$nclass), 
#                       (res2[, grepl("igauss.non.conj_theta", colnames(res2))]
#                        %*% C)/(set2$D/set2$nclass), 
#                       ((res2[, grepl("igamma.conj_lambda", colnames(res2))]/
#                           (res2[, grepl("igamma.conj_eta", colnames(
#                             res2))] - 2)) %*% C)/(set2$D/set2$nclass), 
#                       ((res2[, grepl("igamma.non.conj_lambda", colnames(res2))]/
#                           (res2[, grepl("igamma.non.conj_eta", colnames(
#                             res2))] - 2)) %*% C)/(set2$D/set2$nclass)))
# 
# opar <- par(no.readonly=TRUE)
# at1 <- c(1:set2$nclass, 1:set2$nclass + set2$nclass + 1, 1:set2$nclass +
#            2*set2$nclass + 2, 1:set2$nclass + 3*set2$nclass + 3)
# at2 <- c(mean(1:set2$nclass), mean(1:set2$nclass + set2$nclass + 1),
#          mean(1:set2$nclass + 2*set2$nclass + 2),
#          mean(1:set2$nclass + 3*set2$nclass + 3))
# par(cex=1.2, mar=opar$mar*c(1, 1.3, 1/1.3, 1))
# b2 <- boxplot(prior.means ~ class + method, data=temp1, outline=FALSE,
#               at=at1, col=colors2, names=NA, xaxt="n",
#               ylab=expression(hat(E)(gamma[d]^2)), 
#               subset=method!="2 nc. inv. \n Gaussian")
# r2 <- range(b2$stats, na.rm=TRUE)
# # arrows(mean(1:set2$nclass + set2$nclass + 1), r2[2] - (r2[2] - r2[1])/5, 
# #        mean(1:set2$nclass + set2$nclass + 1), r2[2])
# abline(h=prior.means, col=colors2, lty=2, lwd=1.5)
# # text(at2, par("usr")[3] - 0.2, labels1, srt=45, pos=1, xpd=TRUE)
# axis(1, at=at2, labels=labels2, tick=FALSE)
# legend("topright", legend=c(paste("class", 1:set2$nclass), "true value"), 
#        lty=c(rep(NA, 4), 2), col=c(colors2, 1), seg.len=rep(1, 5), 
#        border=c(rep(1, 4), NA), fill=c(colors2, NA), merge=TRUE)
# par(opar)

################################## supplement ##################################
# ---- hist_igaussian_res3_beta_post ---- 
methods <- c("MCMC", "ADVI", "VB")
labels <- c(methods, expression(hat(E)(beta[0]~"|"~bold(y))))
col <- 2:(length(methods) + 1)
which.beta <- c(1, sample(1:(p + 1), 5))

opar <- par(no.readonly=TRUE)
layout(matrix(c(rep(rep(c(1:3), each=2), 2), rep(rep(c(4:6), each=2), 2)), 
              byrow=TRUE, nrow=4, ncol=6), widths=1, heights=1, respect=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
for(j in 1:length(which.beta)) {
  h1 <- hist(samp.mcmc[, which.beta[j]], plot=FALSE, breaks=30)
  h2 <- hist(samp.advi[, which.beta[j]], plot=FALSE, breaks=30)
  ylim <- c(0, max(dnorm(best.vb[which.beta[j], ], 
                         mean=best.vb[which.beta[j], ], 
                         sd=sqrt(bvar.vb[which.beta[j], ])), 
                   h1$density, h2$density))
  plot(h1, freq=FALSE, 
       main=substitute(beta[i], list(i=which.beta[j] - 1)), 
       xlab=expression(beta),
       ylim=ylim, col=rgb(1, 0, 0, 0.5))
  plot(h2, freq=FALSE, add=TRUE, col=rgb(0, 1, 0, 0.5))
  curve(dnorm(x, mean=best.vb[which.beta[j], ], 
              sd=sqrt(bvar.vb[which.beta[j], ])), add=TRUE, col=col[3], lwd=2)
  abline(v=best.mcmc[which.beta[j], ], col=col[1], lty=2, lwd=2)
  abline(v=best.advi[which.beta[j], ], col=col[2], lty=2, lwd=2)
  abline(v=best.vb[which.beta[j], ], col=col[3], lty=2, lwd=2)
}
legend("topright", legend=labels, 
       lty=c(rep(NA, length(methods)), 2), col=c(col, 1), 
       seg.len=1, border=c(rep(1, length(methods)), NA), 
       fill=c(col, NA), merge=TRUE)
par(opar)

# ---- boxplot_igaussian_res1_mu_mse ----  
res1 <- as.matrix(read.table("results/simulations_igaussian_res1.csv"))
set1 <- read.table("results/simulations_igaussian_set1.csv")
labels1 <- c("conj. inv. \n Gaussian", "non-conj. inv. \n Gaussian",
             "conj. inv. \n Gamma", "non-conj. inv. \n Gamma")
colors1 <- sp::bpy.colors(6)[-c(1, 6)]

C <- matrix(sapply(paste("set1$C", 1:(set1$nclass*set1$D), sep=""), 
                   function(s) {eval(parse(text=s))}), 
            ncol=set1$nclass, nrow=set1$D)
temp1 <- data.frame(class=factor(rep(rep(c(1:set1$nclass), each=set1$D), 4)),
                    method=rep(paste(1:4, labels1), each=set1$nclass*set1$D),
                    prior.means=c((res1[, grepl("igauss.conj_mu", colnames(
                      res1))] %*% C)/(set1$D/set1$nclass), 
                      (res1[, grepl("igauss.non.conj_mu", colnames(res1))]
                       %*% C)/(set1$D/set1$nclass), 
                      (res1[, grepl("igamma.conj_mu", colnames(res1))]
                       %*% C)/(set1$D/set1$nclass),
                      (res1[, grepl("igamma.non.conj_mu", colnames(res1))]
                       %*% C)/(set1$D/set1$nclass)))

opar <- par(no.readonly=TRUE)
at1 <- c(1:set1$nclass, 1:set1$nclass + set1$nclass + 1, 1:set1$nclass +
           2*set1$nclass + 2, 1:set1$nclass + 3*set1$nclass + 3)
at2 <- c(mean(1:set1$nclass), mean(1:set1$nclass + set1$nclass + 1),
         mean(1:set1$nclass + 2*set1$nclass + 2),
         mean(1:set1$nclass + 3*set1$nclass + 3))
par(cex=1.2, mar=opar$mar*c(1, 1.3, 1/1.3, 1))
boxplot(prior.means ~ class + method, data=temp1, outline=FALSE,
        at=at1, col=colors1, names=NA, xaxt="n",
        ylab=expression("MSE" (hat(beta))))
axis(1, at=at2, labels=labels1, tick=FALSE)
legend("topleft", legend=paste("class", 1:set1$nclass), col=colors1, 
       seg.len=rep(1, 4), border=rep(1, 4), fill=colors1)
par(opar)

# ---- boxplot_igaussian_res2_mu_mse ----  
res2 <- as.matrix(read.table("results/simulations_igaussian_res2.csv"))
set2 <- read.table("results/simulations_igaussian_set2.csv")
labels2 <- c("conj. inv. \n Gaussian", "non-conj. inv. \n Gaussian",
             "conj. inv. \n Gamma", "non-conj. inv. \n Gamma")
colors2 <- sp::bpy.colors(6)[-c(1, 6)]

C <- matrix(sapply(paste("set2$C", 1:(set2$nclass*set2$D), sep=""), 
                   function(s) {eval(parse(text=s))}), 
            ncol=set2$nclass, nrow=set2$D)
temp1 <- data.frame(class=factor(rep(rep(c(1:set2$nclass), each=set2$D), 4)),
                    method=rep(paste(1:4, labels2), each=set2$nclass*set2$D),
                    prior.means=c((res2[, grepl("igauss.conj_mu", colnames(
                      res2))] %*% C)/(set2$D/set2$nclass), 
                      (res2[, grepl("igauss.non.conj_mu", colnames(res2))]
                       %*% C)/(set2$D/set2$nclass), 
                      (res2[, grepl("igamma.conj_mu", colnames(res2))]
                       %*% C)/(set2$D/set2$nclass),
                      (res2[, grepl("igamma.non.conj_mu", colnames(res2))]
                       %*% C)/(set2$D/set2$nclass)))

opar <- par(no.readonly=TRUE)
at1 <- c(1:set2$nclass, 1:set2$nclass + set2$nclass + 1, 1:set2$nclass +
           2*set2$nclass + 2, 1:set2$nclass + 3*set2$nclass + 3)
at2 <- c(mean(1:set2$nclass), mean(1:set2$nclass + set2$nclass + 1),
         mean(1:set2$nclass + 2*set2$nclass + 2),
         mean(1:set2$nclass + 3*set2$nclass + 3))
par(cex=1.2, mar=opar$mar*c(1, 1.3, 1/1.3, 1))
boxplot(prior.means ~ class + method, data=temp1, outline=FALSE,
        at=at1, col=colors2, names=NA, xaxt="n",
        ylab=expression("MSE" (hat(beta))))
axis(1, at=at2, labels=labels2, tick=FALSE)
legend("topleft", legend=paste("class", 1:set2$nclass), col=colors2, 
       seg.len=rep(1, 4), border=rep(1, 4), fill=colors2)
par(opar)


# ---- lines_igaussian_res2_igauss_conj_convergence ----  
fit2 <- read.table("results/simulations_igaussian_fit2.csv")
set2 <- read.table("results/simulations_igaussian_set2.csv")
colors2 <- sp::bpy.colors(6)[-c(1, 6)]

temp1 <- sapply(paste("fit2$igauss.conj.theta.", 1:set2$D, sep=""),
                function(s) {eval(parse(text=s))})
temp2 <- sapply(paste("fit2$igauss.conj.lambda.", 1:set2$D, sep=""),
                function(s) {eval(parse(text=s))})
temp3 <- temp1^3/temp2

opar <- par(no.readonly=TRUE)
layout(matrix(c(c(0, 0, 2, 2), c(1, 1, 2, 2), 
                c(1, 1, 3, 3), c(0, 0, 3, 3)), 4, 4), 
       widths=1, heights=1, respect=TRUE)
par(cex=1.3, mar=opar$mar*c(1/1.3, 1.3, 1/1.3, 1/1.3))
plot(temp2[, 1], type="l", ylim=range(temp2), ylab=expression(hat(lambda)),
     xlab="iteration", main="a)")
for(d in 2:set2$D) {
  lines(temp2[, d], col=d)
}
plot(temp1[, 1], type="l", ylim=range(temp1), 
     ylab=expression(hat(E)(gamma[d]^2)), xlab="iteration", main="b)")
for(d in 2:set2$D) {
  lines(temp1[, d], col=d)
}
plot(temp3[, 1], type="l", ylim=range(temp3), 
     ylab=expression(hat(V)(gamma[d]^2)), xlab="iteration", main="c)")
for(d in 2:set2$D) {
  lines(temp3[, d], col=d)
}
layout(matrix(1, 1, 1), widths=1, heights=1, respect=TRUE)
par(opar)

# ---- lines_igaussian_res2_conj_convergence ----
opar <- par(no.readonly=TRUE)
layout(matrix(c(rep(rep(c(1:3), each=2), 2), rep(c(0, 4, 4, 5, 5, 0), 2)), 
              byrow=TRUE, nrow=4, ncol=6), widths=1, heights=1, respect=TRUE)
par(cex=1.3, mar=opar$mar*c(1/2, 1.1, 1/2, 1/2))
plot(fit2.igauss.conj$seq.eb$alpha[, 1], type="l", xlab="iteration", main="a)",
     ylab=expression(hat(alpha)[c]), ylim=range(fit2.igauss.conj$seq.eb$alpha))
for(d in 2:ncol(fit2.igauss.conj$seq.eb$alpha)) {
  lines(fit2.igauss.conj$seq.eb$alpha[, d], col=d)
}
plot(fit2.igauss.conj$seq.eb$theta[, 1], type="l", xlab="iteration", main="b)",
     ylab=expression(hat(theta)[d]), ylim=range(fit2.igauss.conj$seq.eb$theta))
for(d in 2:ncol(fit2.igauss.conj$seq.eb$theta)) {
  lines(fit2.igauss.conj$seq.eb$theta[, d], col=d)
}
plot(fit2.igauss.conj$seq.eb$lambda[, 1], type="l", xlab="iteration", main="c)",
     ylab=expression(hat(lambda)[d]), ylim=range(fit2.igauss.conj$seq.eb$lambda))
for(d in 2:ncol(fit2.igauss.conj$seq.eb$lambda)) {
  lines(fit2.igauss.conj$seq.eb$lambda[, d], col=d)
}
plot(fit2.igamma.conj$seq.eb$eta[, 1], type="l", xlab="iteration", main="d)",
     ylab=expression(hat(eta)[d]), ylim=range(fit2.igamma.conj$seq.eb$eta))
for(d in 2:ncol(fit2.igamma.conj$seq.eb$eta)) {
  lines(fit2.igamma.conj$seq.eb$eta[, d], col=d)
}
plot(fit2.igamma.conj$seq.eb$lambda[, 1], type="l", xlab="iteration", main="e)",
     ylab=expression(hat(lambda)[d]), ylim=range(fit2.igamma.conj$seq.eb$lambda))
for(d in 2:ncol(fit2.igamma.conj$seq.eb$lambda)) {
  lines(fit2.igamma.conj$seq.eb$lambda[, d], col=d)
}
layout(matrix(1, 1, 1), widths=1, heights=1, respect=TRUE)
par(opar)

# ---- scatter_igaussian_res2_conj_prior ----
opar <- par(no.readonly=TRUE)
layout(matrix(c(rep(rep(c(1:3), each=2), 2), rep(c(4, 4, 5, 5, 6, 6), 2)), 
              byrow=TRUE, nrow=4, ncol=6), widths=1, heights=1, respect=TRUE)
par(cex=1.3, mar=opar$mar*c(1/2, 1.1, 1/2, 1/2))
plot(alpha, fit2.igauss.conj$seq.eb$alpha[fit2.igauss.conj$iter$eb, ], 
     ylab=expression(hat(alpha)), xlab=expression(alpha), main="a)")
plot(lambda, fit2.igauss.conj$seq.eb$lambda[fit2.igauss.conj$iter$eb, ],
     ylab=expression(hat(lambda)), xlab=expression(lambda), main="b)")
plot(lambda/(eta - 2), 
     fit2.igauss.conj$seq.eb$theta[fit2.igauss.conj$iter$eb, ],
     ylab=expression(hat(E)(gamma[d]^2)), xlab=expression(E(gamma[d]^2)), 
     main="c)")
plot(eta, fit2.igamma.conj$seq.eb$eta[fit2.igamma.conj$iter$eb, ], 
     ylab=expression(hat(eta)), xlab=expression(eta), main="d)")
plot(lambda, fit2.igamma.conj$seq.eb$lambda[fit2.igamma.conj$iter$eb, ],
     ylab=expression(hat(lambda)), xlab=expression(lambda), main="e)")
plot(lambda/(eta - 2), 
     fit2.igamma.conj$seq.eb$lambda[fit2.igamma.conj$iter$eb, ]/
       (fit2.igamma.conj$seq.eb$eta[fit2.igamma.conj$iter$eb, ] - 2),
     ylab=expression(hat(E)(gamma[d]^2)), xlab=expression(E(gamma[d]^2)), 
     main="f)")
layout(matrix(1, 1, 1), widths=1, heights=1, respect=TRUE)
par(opar)

# ---- scatter_igaussian_res2_conj_post ----
opar <- par(no.readonly=TRUE)
layout(matrix(c(rep(rep(c(1:3), each=2), 2), rep(c(4, 4, 5, 5, 6, 6), 2)), 
              byrow=TRUE, nrow=4, ncol=6), widths=1, heights=1, respect=TRUE)
par(cex=1.3, mar=opar$mar*c(1/2, 1.1, 1/2, 1/2))
plot(beta, fit2.igauss.conj$vb.post$mu, 
     ylab=expression(hat(E)(bold(beta)[d]~"|"~bold(y))),
     xlab=expression(bold(beta)[d]), main="a)")
plot(gamma^2, sqrt(
  fit2.igauss.conj$vb.post$delta/fit2.igauss.conj$seq.eb$lambda[
    fit2.igauss.conj$iter$eb, ])*fit2.igauss.conj$seq.eb$theta[
      fit2.igauss.conj$iter$eb, ]*sapply(sqrt(fit2.igauss.conj$seq.eb$lambda[
        fit2.igauss.conj$iter$eb, ]*fit2.igauss.conj$vb.post$delta)/
          fit2.igauss.conj$seq.eb$theta[fit2.igauss.conj$iter$eb, ], 
        ratio_besselK, nu=0.5*(p + 1)), 
  ylab=expression(hat(E)(gamma[d]^2~"|"~bold(y))),
  xlab=expression(gamma[d]^2), main="b)")
plot(sigma^2, 2*fit2.igauss.conj$vb.post$zeta/(p + n - 1), 
     ylab=expression(hat(E)(sigma[d]^2~"|"~bold(y))),
     xlab=expression(sigma[d]^2), main="c)")
plot(beta, fit2.igamma.conj$vb.post$mu, 
     ylab=expression(hat(E)(bold(beta)[d]~"|"~bold(y))),
     xlab=expression(bold(beta)[d]), main="d)")
plot(gamma^2, fit2.igamma.conj$vb.post$delta/
       (p + fit2.igamma.conj$seq.eb$eta[fit2.igamma.conj$iter$eb, ] - 2), 
     ylab=expression(hat(E)(gamma[d]^2~"|"~bold(y))),
     xlab=expression(gamma[d]^2), main="e)")
plot(sigma^2, 2*fit2.igamma.conj$vb.post$zeta/(p + n - 1), 
     ylab=expression(hat(E)(sigma[d]^2~"|"~bold(y))),
     xlab=expression(sigma[d]^2), main="f)")
layout(matrix(1, 1, 1), widths=1, heights=1, respect=TRUE)
par(opar)


# ---- scatter_igaussian_res3_beta_vs_best ----
opar <- par(no.readonly=TRUE)
layout(matrix(c(rep(rep(c(1:3), each=2), 2), rep(c(0, 4, 4, 5, 5, 0), 2)), 
              byrow=TRUE, nrow=4, ncol=6), widths=1, heights=1, respect=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
plot(rbind(0, beta), best.glmnet, ylab=expression(hat(beta)), 
     xlab=expression(beta),
     main=paste0("glmnet, MSE=", 
                 round(mean((best.glmnet - rbind(0, beta))^2), 4)))
abline(a=0, b=1, col=2, lty=2)
plot(rbind(0, beta), best.mcmc, ylab=expression(hat(beta)), 
     xlab=expression(beta),
     main=paste0("MCMC, MSE=", 
                 round(mean((best.mcmc - rbind(0, beta))^2), 4)))
abline(a=0, b=1, col=2, lty=2)
plot(rbind(0, beta), best.vb, ylab=expression(hat(beta)), 
     xlab=expression(beta),
     main=paste0("VB, MSE=", 
                 round(mean((best.vb - rbind(0, beta))^2), 4)))
abline(a=0, b=1, col=2, lty=2)
plot(rbind(0, beta), best.advi, ylab=expression(hat(beta)), 
     xlab=expression(beta),
     main=paste0("ADVI, MSE=", 
                 round(mean((best.advi - rbind(0, beta))^2), 4)))
abline(a=0, b=1, col=2, lty=2)
plot(rbind(0, beta), best.map, ylab=expression(hat(beta)), 
     xlab=expression(beta),
     main=paste0("MAP, MSE=", 
                 round(mean((best.map - rbind(0, beta))^2), 4)))
abline(a=0, b=1, col=2, lty=2)
par(opar)

# ---- scatter_igaussian_res3_mcmc_vs_best ----
opar <- par(no.readonly=TRUE)
layout(matrix(c(rep(rep(c(1:2), each=2), 2), rep(rep(c(3:4), each=2), 2)),
              byrow=TRUE, nrow=4, ncol=4), widths=1, heights=1, respect=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
plot(best.mcmc, best.vb, ylab=expression(hat(E)(beta["VB"]~"|"~bold(y))),
     xlab=expression(hat(E)(beta["MCMC"]~"|"~bold(y))))
abline(a=0, b=1, col=2, lty=2)
plot(best.mcmc, best.advi, ylab=expression(hat(beta)["ADVI"]),
     xlab=expression(hat(beta)["MCMC"]))
abline(a=0, b=1, col=2, lty=2)
plot(bvar.mcmc, bvar.vb, ylab=expression(hat(V)(beta["VB"]~"|"~bold(y))),
     xlab=expression(hat(V)(beta["MCMC"]~"|"~bold(y))))
abline(a=0, b=1, col=2, lty=2)
plot(bvar.mcmc, bvar.advi, ylab=expression(hat(V)(beta["ADVI"]~"|"~bold(y))),
     xlab=expression(hat(V)(beta["MCMC"]~"|"~bold(y))))
abline(a=0, b=1, col=2, lty=2)
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



        