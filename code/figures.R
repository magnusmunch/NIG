################################## supplement ##################################
# ---- figures ----
# ---- boxplot_igaussian_res1_prior_mean ----  
res1 <- as.matrix(read.table("results/simulations_igaussian_res1.csv"))
set1 <- read.table("results/simulations_igaussian_set1.csv")
labels1 <- c("c. inv. \n Gaussian", "nc. inv. \n Gaussian",
              "c. inv. \n Gamma", "nc. inv. \n Gamma")
colors1 <- sp::bpy.colors(6)[-c(1, 6)]

C <- matrix(sapply(paste("set1$C", 1:(set1$nclass*set1$D), sep=""), 
                   function(s) {eval(parse(text=s))}), 
            ncol=set1$nclass, nrow=set1$D)
prior.means <- as.numeric((sapply(paste("set1$theta", 1:set1$D, sep=""),
                                  function(s) {eval(parse(text=s))}) %*% C)/
                            (set1$D/set1$nclass))
temp1 <- data.frame(class=factor(rep(rep(c(1:set1$nclass), each=set1$D), 4)),
                    method=rep(paste(1:4, labels1), each=set1$nclass*set1$D),
                    prior.means=c((res1[, grepl("igauss.conj_theta", colnames(
                      res1))] %*% C)/(set1$D/set1$nclass), 
                      (res1[, grepl("igauss.non.conj_theta", colnames(res1))]
                       %*% C)/(set1$D/set1$nclass), 
                      ((res1[, grepl("igamma.conj_lambda", colnames(res1))]/
                          (res1[, grepl("igamma.conj_eta", colnames(
                            res1))] - 2)) %*% C)/(set1$D/set1$nclass), 
                      ((res1[, grepl("igamma.non.conj_lambda", colnames(res1))]/
                          (res1[, grepl("igamma.non.conj_eta", colnames(
                            res1))] - 2)) %*% C)/(set1$D/set1$nclass)))
  
opar <- par(no.readonly=TRUE)
at1 <- c(1:set1$nclass, 1:set1$nclass + set1$nclass + 1, 1:set1$nclass +
           2*set1$nclass + 2, 1:set1$nclass + 3*set1$nclass + 3)
at2 <- c(mean(1:set1$nclass), mean(1:set1$nclass + set1$nclass + 1),
         mean(1:set1$nclass + 2*set1$nclass + 2),
         mean(1:set1$nclass + 3*set1$nclass + 3))
par(cex=1.2, mar=opar$mar*c(1, 1.3, 1/1.3, 1))
boxplot(prior.means ~ class + method, data=temp1, outline=FALSE,
        at=at1, col=colors1, names=NA, xaxt="n",
        ylab=expression(hat(E)(gamma[d]^2)))
abline(h=prior.means, col=colors1, lty=2, lwd=1.5)
# text(at2, par("usr")[3] - 0.2, labels1, srt=45, pos=1, xpd=TRUE)
axis(1, at=at2, labels=labels1, tick=FALSE)
legend("topright", legend=c(paste("class", 1:set1$nclass), "true value"), 
       lty=c(rep(NA, 4), 2), col=c(colors1, 1), seg.len=rep(1, 5), 
       border=c(rep(1, 4), NA), fill=c(colors1, NA), merge=TRUE)
par(opar)

# ---- boxplot_igaussian_res1_mu_mse ----  
res1 <- as.matrix(read.table("results/simulations_igaussian_res1.csv"))
set1 <- read.table("results/simulations_igaussian_set1.csv")
labels1 <- c("c. inv. \n Gaussian", "nc. inv. \n Gaussian",
             "c. inv. \n Gamma", "nc. inv. \n Gamma")
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

# ---- boxplot_igaussian_res2_prior_mean ----  
res2 <- as.matrix(read.table("results/simulations_igaussian_res2.csv"))
set2 <- read.table("results/simulations_igaussian_set2.csv")
labels2 <- c("c. inv. \n Gaussian", "nc. inv. \n Gaussian",
             "c. inv. \n Gamma", "nc. inv. \n Gamma")
colors2 <- sp::bpy.colors(6)[-c(1, 6)]

C <- matrix(sapply(paste("set2$C", 1:(set2$nclass*set2$D), sep=""), 
                   function(s) {eval(parse(text=s))}), 
            ncol=set2$nclass, nrow=set2$D)
prior.means <- as.numeric((sapply(paste(
  "set2$lambda", 1:set2$D, sep=""), function(s) {eval(parse(text=s))})/
    (sapply(paste("set2$eta", 1:set2$D, sep=""), function(s) {
      eval(parse(text=s))}) - 2)) %*% C)/(set1$D/set1$nclass)
temp1 <- data.frame(class=factor(rep(rep(c(1:set2$nclass), each=set2$D), 4)),
                    method=rep(paste(1:4, labels2), each=set2$nclass*set2$D),
                    prior.means=c((res2[, grepl("igauss.conj_theta", colnames(
                      res2))] %*% C)/(set2$D/set2$nclass), 
                      (res2[, grepl("igauss.non.conj_theta", colnames(res2))]
                       %*% C)/(set2$D/set2$nclass), 
                      ((res2[, grepl("igamma.conj_lambda", colnames(res2))]/
                          (res2[, grepl("igamma.conj_eta", colnames(
                            res2))] - 2)) %*% C)/(set2$D/set2$nclass), 
                      ((res2[, grepl("igamma.non.conj_lambda", colnames(res2))]/
                          (res2[, grepl("igamma.non.conj_eta", colnames(
                            res2))] - 2)) %*% C)/(set2$D/set2$nclass)))

opar <- par(no.readonly=TRUE)
at1 <- c(1:set2$nclass, 1:set2$nclass + set2$nclass + 1, 1:set2$nclass +
           2*set2$nclass + 2, 1:set2$nclass + 3*set2$nclass + 3)
at2 <- c(mean(1:set2$nclass), mean(1:set2$nclass + set2$nclass + 1),
         mean(1:set2$nclass + 2*set2$nclass + 2),
         mean(1:set2$nclass + 3*set2$nclass + 3))
par(cex=1.2, mar=opar$mar*c(1, 1.3, 1/1.3, 1))
b2 <- boxplot(prior.means ~ class + method, data=temp1, outline=FALSE,
              at=at1, col=colors2, names=NA, xaxt="n",
              ylab=expression(hat(E)(gamma[d]^2)), 
              subset=method!="2 nc. inv. \n Gaussian")
r2 <- range(b2$stats, na.rm=TRUE)
arrows(mean(1:set2$nclass + set2$nclass + 1), r2[2] - (r2[2] - r2[1])/5, 
       mean(1:set2$nclass + set2$nclass + 1), r2[2])
abline(h=prior.means, col=colors2, lty=2, lwd=1.5)
# text(at2, par("usr")[3] - 0.2, labels1, srt=45, pos=1, xpd=TRUE)
axis(1, at=at2, labels=labels2, tick=FALSE)
legend("topright", legend=c(paste("class", 1:set2$nclass), "true value"), 
       lty=c(rep(NA, 4), 2), col=c(colors2, 1), seg.len=rep(1, 5), 
       border=c(rep(1, 4), NA), fill=c(colors2, NA), merge=TRUE)
par(opar)

# ---- boxplot_igaussian_res2_mu_mse ----  
res2 <- as.matrix(read.table("results/simulations_igaussian_res2.csv"))
set2 <- read.table("results/simulations_igaussian_set2.csv")
labels2 <- c("c. inv. \n Gaussian", "nc. inv. \n Gaussian",
             "c. inv. \n Gamma", "nc. inv. \n Gamma")
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
labels2 <- c("c. inv. \n Gaussian", "nc. inv. \n Gaussian",
             "c. inv. \n Gamma", "nc. inv. \n Gamma")
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











        