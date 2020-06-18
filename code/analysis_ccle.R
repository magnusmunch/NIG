#!/usr/bin/env Rscript

### libraries
packages <- c()
sapply(packages, library, character.only=TRUE)

### load data
load(file="data/data_gdsc_dat1.Rdata")
expr.gdsc <- expr$expr
load(file="data/data_ccle_dat1.Rdata")

################################################################################
################################### analysis ###################################
################################################################################
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
