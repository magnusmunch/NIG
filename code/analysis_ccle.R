#!/usr/bin/env Rscript

################################################################################
# CCLE log IC50 values and gene expression analyses
################################################################################

### libraries
packages <- c("glmnet", "cambridge")
sapply(packages, library, character.only=TRUE)

### load data
load(file="data/data_gdsc_dat1.Rdata")
expr.gdsc <- expr$expr
load(file="data/data_ccle_dat1.Rdata")

################################## analysis 1 ##################################
### data preparation
# reassign data
resp <- resp.prep
expr <- expr.prep

# select genes that are also in the GDSC data
expr.sel <- expr[, colnames(expr) %in% colnames(expr.gdsc)]

# select cell lines that are common to both response and expression
expr.sel <- sapply(resp, function(s) {
  log(expr.sel[rownames(expr.sel) %in% names(s), ] + 1)})
resp.sel <- sapply(1:length(resp), function(d) {
  resp[[d]][names(resp[[d]]) %in% rownames(expr.sel[[d]])]})

# select expressions that have non zero variance
expr.sel <- sapply(expr.sel, function(s) {
  s[, apply(s, 2, sd)!=0]})

# create data objects used in fitting
x <- sapply(expr.sel, scale)
y <- sapply(resp.sel, scale)

# calculate p-values for correlation
pvalues <- sapply(1:length(x), function(d) {
  apply(x[[d]], 2, function(s) {cor.test(s, y[[d]])$p.value})})
genenames <- unique(unlist(sapply(x, colnames)))
pvaluesvec <- unlist(pvalues)
harmonic <- sapply(genenames, function(s) {
  1/mean(1/pvaluesvec[names(pvaluesvec)==s])})
fisher <- sapply(genenames, function(s) {
  c(-2*sum(log(pvaluesvec[names(pvaluesvec)==s])), 
    2*sum(names(pvaluesvec)==s))})
pfisher <- apply(fisher, 2, function(s) {
  pchisq(s[1], df=s[2], lower.tail=FALSE)})

# create co-data
feat.prep <- cbind(chifisher=fisher[1, ], pfisher=pfisher, harmonic=harmonic)
rownames(feat.prep) <- genenames
write.table(feat.prep, file="results/analysis_ccle_res1.txt")

################################## analysis 2 ##################################
### data preparation
# reassign data
resp <- resp.prep
expr <- expr.prep

# select genes that are also in the GDSC data
expr.sel <- expr[, colnames(expr) %in% colnames(expr.gdsc)]

# select cell lines that are common to both response and expression
expr.sel <- sapply(resp, function(s) {
  log(expr.sel[rownames(expr.sel) %in% names(s), ] + 1)})
resp.sel <- sapply(1:length(resp), function(d) {
  resp[[d]][names(resp[[d]]) %in% rownames(expr.sel[[d]])]})

# create training and test data
set.seed(2020)
idtrain <- sapply(resp.sel, function(s) {
  sample(1:length(s), floor(length(s)/2))})
ytrain <- sapply(1:length(resp.sel), function(d) {
  resp.sel[[d]][idtrain[[d]]]})
idconstresp <- which(sapply(ytrain, sd)==0)
ytrain <- sapply(ytrain, function(s) {scale(s)[, 1]})[-idconstresp]
ytest <- sapply(1:length(resp.sel), function(d) {
  scale(resp.sel[[d]][-idtrain[[d]]])[, 1]})[-idconstresp]
xtrain <- sapply(1:length(expr.sel), function(d) {
  expr.sel[[d]][idtrain[[d]], ]})
idconst <- sapply(xtrain, function(s) {which(apply(s, 2, sd)==0)})
xtrain <- sapply(1:length(xtrain), function(d) {
  scale(xtrain[[d]][, -idconst[[d]]])})[-idconstresp]
xtest <- sapply(1:length(expr.sel), function(d) {
  dat <- expr.sel[[d]][-idtrain[[d]], -idconst[[d]]]
  apply(dat, 2, function(s) {
    if(sd(s)==0) {
      s
    } else {
      scale(s)
    }})})[-idconstresp]

### fit models
fit1.lasso <- lapply(1:length(ytrain), function(d) {
  cv.glmnet(xtrain[[d]], ytrain[[d]], alpha=1, intercept=FALSE)})
fit1.ridge <- lapply(1:length(ytrain), function(d) {
  cv.glmnet(xtrain[[d]], ytrain[[d]], alpha=0, intercept=FALSE)})

best <- list(lasso1=sapply(fit1.lasso, function(s) {
  coef(s, s="lambda.min")[-1, ]}),
  ridge1=sapply(fit1.ridge, function(s) {
    coef(s, s="lambda.min")[-1, ]}))
pred <- lapply(best, function(b) {sapply(1:length(ytrain), function(d) {
  as.numeric(xtest[[d]] %*% b[[d]])})})
pmse <- sapply(pred, function(s) {sapply(1:length(s), function(d) {
  mean((ytest[[d]] - s[[d]])^2)})})
colMeans(pmse)

sapply(1:length(xtrain), function(d) {
  sum(!(colnames(expr.gdsc) %in% colnames(xtrain[[d]])))})
mae <- sapply(best, function(s) {sapply(s, function(b) {mean(abs(b))})})
plot(mae[, 2], pmse[, 2])
str(best[[1]][[1]])

