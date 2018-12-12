########################################################################
# R script
# 
# title:
# empirical Bayes for general Bayesian SEM
#
# authors:
# Gwenael G.R. Leday
#
# date: 24/07/2018
#
########################################################################


# Work directory
setwd("somwhere in my computer")

# Clean environment
rm(list=ls());gc()

# Libraries
library(Rcpp)
library(RcppArmadillo)




#---------------------------#
#       Generate data
#---------------------------#

# sample size
n <- 100

# sizes of design matrices
p <- sample(100:200, 100, replace=TRUE)

# Generate outcomes
listy <- sapply(p, function(xx, nn){rnorm(nn)}, nn=n, simplify=FALSE)

# Generate design matrices
generateX <- function(myp, myn){
	theX <- matrix(rnorm(myn*myp), nrow = myn, ncol = myp)
	colnames(theX) <- paste("X", 1:ncol(theX), sep="")
	theX
}
listX <- sapply(p, generateX, myn = n)

# Generate matrix P (3 global shrinkage priors)
listP <- sapply(p, function(xx){sample(1:3, xx, replace=TRUE)})



#---------------------------#
#    run variational algo
#---------------------------#

# Compile cpp functions
sourceCpp("myfunctions.cpp")

# Input arguments
niter <- 20 # number of iterations
EBid <- 1:3 # indexes of gamma priors for which EB is desired. If no EB then non-informative. 

# Intialization
priors <- sort(unique(unlist(listP)))
nbpriors <- length(priors)
aprior <- bprior <- matrix(NA, niter+1, nbpriors)
colnames(aprior) <- paste("a", 1:nbpriors, sep="")
colnames(bprior) <- paste("b", 1:nbpriors, sep="")
aprior[1,] <- bprior[1,] <- 0.001
idxPriorList <- lapply(listP, function(x){sort(unique(x))})
allMLs <- matrix(NA, length(listX), niter)

for(j in 1:niter){
	cat("iteration ", j, "\n")

	# Prior as lists
	inputaList <- lapply(idxPriorList, function(x){aprior[j,x]})
	inputbList <- lapply(idxPriorList, function(x){bprior[j,x]})

	if(j==1){
		mydstarvec <- rep(0, length(inputaList))
		mybstarlist <- lapply(idxPriorList, function(x){rep(0, length(x))})
	}else{
		mydstarvec <- unlist(lapply(res$postSigList, function(x){x[2]}))
		mybstarlist <- lapply(res$postRandList, function(x){x[,2]})
	}

	# Fit BSEM
	res <- BSEMVarOneIter(ylist=listy, Xlist=listX, Plist=listP, alist=inputaList, blist=inputbList, bstarlist=mybstarlist, cSigma=0.001, dSigma=0.001, dstarvec=mydstarvec, lincomblist=list())

	# Marginal likelihood
	allMLs[,j] <- res$allmargs

	# Empirical Bayes
	for(ii in 1:nbpriors){
		if(ii%in%EBid){
			# Get posterior shape and rate parameters
			allaRandStar <- sapply(1:length(idxPriorList), function(x){res$postRandList[[x]][idxPriorList[[x]]==ii,1]}, simplify=TRUE)
			allbRandStar <- sapply(1:length(idxPriorList), function(x){res$postRandList[[x]][idxPriorList[[x]]==ii,2]}, simplify=TRUE)

			# Variational Empirical Bayes using fixed-point iteration as in Valpola and Honkela (2006)
			ab <- c(aprior[j, ii], bprior[j, ii])
			ab <- fixedPointIterEB(initab=ab, myallaRandStar=allaRandStar, myallbRandStar=allbRandStar, mymaxiter=20, myeps=1e-6)
ab
			aprior[j+1, ii] <- ab[1]
			bprior[j+1, ii] <- ab[2]
		}
	}
}


# Plot sum of marginal likelihood
plot(colSums(allMLs)[-1])

# Posterior expectation of gamma priors
aprior[niter,]/bprior[niter,]



