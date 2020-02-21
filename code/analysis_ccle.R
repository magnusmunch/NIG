#!/usr/bin/env Rscript

################################################################################
# regression GDSC log IC50 values on gene expression/mutations of cell lines
################################################################################

### libraries
packages <- c()
sapply(packages, library, character.only=TRUE)

### load data
load(file="data/data_ccle_dat1.Rdata")

# reassign data
resp <- resp.prep
expr <- expr.prep

################################## analysis 1 ##################################
### data preparation
# only retain the cell lines that are available in both response and expression
expr.prep <- expr[rownames(expr) %in% rownames(resp), ]
resp.prep <- resp[rownames(resp) %in% rownames(expr), ]

# select expressions that have non zero variance
expr.prep <- expr.prep[, apply(expr.prep, 2, sd)!=0]

# create data objects used in fitting
x <- scale(expr.prep, scale=FALSE)
y <- scale(log(resp.prep), scale=FALSE)
p <- ncol(x)
n <- nrow(y)

# create co-data
D <- ncol(resp.prep)
pvalues <- sapply(1:D, function(d) {
  apply(x, 2, function(s) {cor.test(s, y[, d])$p.value})})
feat.prep <- cbind(chifisher=-2*rowSums(log(pvalues)),
                   pfisher=-pchisq(-2*rowSums(log(pvalues)), 2*D, 
                                   lower.tail=FALSE, log.p=TRUE),
                   harmonic=log(rowSums(1/(pvalues*D))))
rownames(feat.prep) <- colnames(x)
write.table(feat.prep, file="results/analysis_ccle_res1.txt")

