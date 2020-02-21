#!/usr/bin/env Rscript

################################ reading data ##################################
resp <- read.table("data/CCLE_NP24.2009_Drug_data_2015.02.24.csv", 
                   header=TRUE, sep=",", stringsAsFactors=FALSE)

expr <- read.table("data/CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz",
                   header=TRUE, stringsAsFactors=FALSE)

################################ data cleaning #################################
# expression data
expr1 <- t(expr[, -c(1, 2)])
colnames(expr1) <- sapply(strsplit(expr$gene_id, ".", fixed=TRUE), "[[", 1)
rownames(expr1)[substr(rownames(expr1), 1, 1)=="X"] <- 
  substr(rownames(expr1)[substr(rownames(expr1), 1, 1)=="X"], 2, 
         nchar(rownames(expr1)[substr(rownames(expr1), 1, 1)=="X"]))
expr2 <- expr1[, apply(expr1, 2, sd)!=0]
expr.prep <- expr2

# create a separate tissue variable for responses
resp$Tissue <- unlist(lapply(strsplit(resp$CCLE.Cell.Line.Name, "_"), 
                             function(l) {paste(l[-1], collapse=" ")}))

# select the sigmoid type fits (matches GDSC IC50)
resp2 <- resp[resp$FitType=="Sigmoid", c(1, 3, 11)]

# switch to wide format
resp3 <- reshape(resp2, timevar="Compound", idvar="CCLE.Cell.Line.Name",
                 direction="wide")
rownames(resp3) <- resp3$CCLE.Cell.Line.Name
resp4 <- resp3[, -1]
colnames(resp4) <- sapply(strsplit(colnames(resp4), "IC50..uM.."), "[[", 2)

# sequentially remove drug or cell line with highest proportion missings
nomiss <- FALSE
resp5 <- resp4
while(!nomiss) {
  pcell <- apply(resp5, 1, function(x) {sum(is.na(x))/length(x)})
  pdrug <- apply(resp5, 2, function(x) {sum(is.na(x))/length(x)})
  
  if(max(pcell) >= max(pdrug)) {
    resp5 <- resp5[-which.max(pcell), ]
  } else {
    resp5 <- resp5[, -which.max(pdrug)]
  }
  nomiss <- !any(is.na(resp5))
}
resp.prep <- resp5

rm(expr, expr1, expr2, nomiss, pcell, pdrug, resp, resp2, resp3, resp4, resp5)

save(expr.prep, resp.prep, file="data/data_ccle_dat1.Rdata")
