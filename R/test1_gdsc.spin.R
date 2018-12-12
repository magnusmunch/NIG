path.data <- "/Users/magnusmunch/Documents/OneDrive/PhD/project_cambridge/data/"


data <- read.table(paste(path.data, "v17.3_fitted_dose_response.csv", sep=""),
                   header=TRUE, sep=",", stringsAsFactors=FALSE)
data2 <- aggregate(LN_IC50 ~ CELL_LINE_NAME + DRUG_NAME, 
                   data=data[data$DRUG_ID > 1000, ], FUN=mean, na.rm=TRUE)
data3 <- reshape(data2, timevar="DRUG_NAME", idvar="CELL_LINE_NAME", 
                 direction="wide")


data4 <- na.omit(data3[, apply(data3, 2, function(x) {sum(is.na(x))}) <= 150])
str(data4)

plot(apply(data3[, -1], 2, function(x) {sum(is.na(x))}))
data4 <- na.omit(data3)

str(data4)


length(unique(data$DRUG_ID[data$DRUG_ID < 1000]))
length(unique(data$DRUG_ID[data$DRUG_ID > 1000]))
unique(data$DATASET_VERSION)




tau <- seq(-10.005, 10.005, 0.01)
T <- 20
plot(tau, sapply(tau, function(s) {sum(trigamma(s/2 + (1 - c(1:T))/2))}), type="l")



library(plot3Drgl)
rho <- seq(-0.99, 0.99, 0.01)
tau <- seq(0.01, 10, 0.01)
mll <- function(tau, rho, T, D, w1, w2, w3) {
  if(tau > max(rho - T*rho, rho, 0)) {
    out <- w1*tau - w2*tau^2  - w3*tau*rho - 
      D*sapply(tau, function(s) {sum(lgamma((s/2 + (1 - c(1:T))/2)))}) + 
      0.5*D*T*tau*log(tau) + 0.5*D*tau*log(tau - (1 - T)*rho) + 
      0.5*D*(T - 1)*tau*log(tau - rho)
    if(is.infinite(out)) {out <- NA}
  } else {
    out <- NA
  }
  return(out)
}

z <- matrix(NA, length(tau), length(rho))
for(i in 1:nrow(z)) {
  for(j in 1:ncol(z)) {
    z[i, j] <- mll(tau[i], rho[j], 21, 100, -10, 0.01, -10)
  }
}


persp3Drgl(tau, rho, z, xlab=expression(tau), ylab=expression(rho), zlab="mll")
par3d(windowRect=c(20, 30, 800, 800))


