#!/usr/bin/env Rscript

### libraries
library(biomaRt)
# library(MASS)

### load data
load("data/P38MAPKpathwayKb100K.RData")
rm(G1, G2)

### feature data
distance <- sapply(1:length(Gn), function(g) {
  sapply(1:length(Gn[[g]][[4]]), function(s) {
    (1 - Gn[[g]][[4]][s])*min(abs(Gn[[g]][[5]][s] - Gn[[g]][[3]][c(1, 2)]))})})
inorout <- lapply(Gn, "[[", 4)
# View(listEnsembl())
geneid <- sapply(Gn, function(g) {
  strsplit(g[[7]]$gene_id, "\\.")[[1]][1]})
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",
                      GRCh=37)
# View(listFilters(ensembl))
gb <- getBM(attributes=c("ensembl_gene_id", "ensembl_exon_id", 
                         "exon_chrom_start", "exon_chrom_end",
                         "chromosome_name"), 
            filters=c("ensembl_gene_id"), 
            values=geneid, 
            mart=ensembl)
exon <- sapply(1:length(Gn), function(g) {
  vec <- sapply(Gn[[g]][[5]], function(s) {
    any(apply(gb[gb$ensembl_gene_id==geneid[[g]], c(3, 4)], 1, function(e) {
      (s >= e[1]) & (s <= e[2])}))}) + inorout[[g]] + 1; 
  ifelse(vec==1, "inter", ifelse(vec==2, "intron", "exon"))})
minoraf <- lapply(sapply(Gn, "[[", 2), function(g) {
  f <- apply(g, 2, function(s) {f <- mean(s); unname(min(f, 1 - f))})})
feat.prep <- list(distance=distance, exon=exon, minoraf=minoraf)
rm(distance, ensembl, exon, gb, inorout, minoraf)

### gene data
genelength <- sapply(1:length(Gn), function(g) {
  diff(as.numeric(Gn[[g]][[3]][c(1, 2)]))})
gene.prep <- data.frame(gene=geneid, length=genelength)
rm(genelength)

### gene expressions
expr <- sapply(Gn, "[[", 1)
colnames(expr) <- geneid
expr.prep <- expr
rm(expr)

### SNPs
snps <- sapply(Gn, "[[", 2)
names(snps) <- geneid
snps.prep <- snps
rm(snps, geneid, Gn)

### save data
save(gene.prep, snps.prep, expr.prep, feat.prep, 
     file="data/data_eqtl_dat1.Rdata")
