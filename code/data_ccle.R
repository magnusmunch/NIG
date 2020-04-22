#!/usr/bin/env Rscript

### libraries
library(gdata)
library(biomaRt)

### biomart
ensembl <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

################################# responses ####################################
resp1 <- read.table("data/CCLE_NP24.2009_Drug_data_2015.02.24.csv", 
                    header=TRUE, sep=",", stringsAsFactors=FALSE)
resp2 <- read.xls("data/CCLE_GNF_data_090613.xls", stringsAsFactors=FALSE,
                  strip.white=TRUE, skip=1, header=FALSE)
resp3 <- read.xls("data/CCLE_GNF_data_090613.xls", stringsAsFactors=FALSE,
                    strip.white=TRUE, nrows=1, header=FALSE)
resp4 <- read.table("data/CCLE_NP24.2009_profiling_2012.02.20.csv", header=TRUE,
                    stringsAsFactors=FALSE, sep=",")

# reassign headers
colnames(resp2) <- trimws(as.character(resp3))

# select the sigmoid type fits (matches GDSC IC50)
resp1.temp <- resp1[, c(1, 3, 11)]

# truncating IC50 to maximium dose value
resp1.temp$IC50..uM. <- apply(cbind(
  resp1.temp$IC50..uM., 
  sapply(strsplit(resp1$Doses..uM., ",", fixed=TRUE), function(s) {
    min(as.numeric(s))})), 1, max)
resp1.temp$IC50..uM. <- apply(cbind(
  resp1.temp$IC50..uM., 
  sapply(strsplit(resp1$Doses..uM., ",", fixed=TRUE), function(s) {
    max(as.numeric(s))})), 1, min)

# save censoring statusses
resp1.temp$cens <- ifelse(resp1.temp$IC50..uM. >= sapply(
  strsplit(resp1$Doses..uM., ",", fixed=TRUE), function(s) {
    max(as.numeric(s))}), 1, ifelse(resp1.temp$IC50..uM. <= sapply(
      strsplit(resp1$Doses..uM., ",", fixed=TRUE), function(s) {
        min(as.numeric(s))}), -1, 0))

# switch to wide format
resp1.cens <- reshape(resp1.temp[, -3], timevar="Compound", 
                      idvar="CCLE.Cell.Line.Name", direction="wide")
resp1.cens[is.na(resp1.cens)] <- 0
resp1.temp <- reshape(resp1.temp[, -4], timevar="Compound", 
                      idvar="CCLE.Cell.Line.Name", direction="wide")
rownames(resp1.temp) <- sapply(strsplit(resp1.temp$CCLE.Cell.Line.Name, "_"), 
                               "[", 1)
rownames(resp1.cens) <- sapply(strsplit(resp1.temp$CCLE.Cell.Line.Name, "_"), 
                               "[", 1)
resp1.temp <- resp1.temp[, -1]
resp1.cens <- resp1.cens[, -1]
colnames(resp1.temp) <- sapply(strsplit(colnames(resp1.temp), "IC50..uM.."), 
                               "[[", 2)
colnames(resp1.cens) <- sapply(strsplit(colnames(resp1.cens), "cens."), 
                               "[[", 2)
resp1.temp <- as.matrix(resp1.temp)
resp1.cens <- as.matrix(resp1.cens)

# more data
resp2.temp <- t(resp2[, -c(1:4)])
rownames(resp2.temp) <- toupper(gsub("/", "", gsub(".", "", gsub(
  "-", "", colnames(resp2)[-c(1:4)], fixed=TRUE), fixed=TRUE), fixed=TRUE))
colnames(resp2.temp) <- trimws(sapply(1:nrow(resp2), function(d) {
  s <- ifelse(resp2$"other name"[d]=="", resp2$"compound name"[d], 
              paste(resp2$"compound name"[d], resp2$"other name"[d], sep=","))
  gsub("\xa0", "", s)}))

# info on drugs
resp3$names <- sapply(1:nrow(resp3), function(d) {
  s <- ifelse(resp3$Compound..brand.name.[d]=="",
              resp3$Compound..code.or.generic.name.[d], 
              paste(resp3$Compound..code.or.generic.name.[d],
                    resp3$Compound..brand.name.[d], sep=","))
  gsub("\xa0", "", s)})
resp3$names <- gsub("PF2341066", "PF-2341066", resp3$names)
colnames(resp1.temp) <- gsub("PF2341066", "PF-2341066", colnames(resp1.temp))
colnames(resp1.cens) <- gsub("PF2341066", "PF-2341066", colnames(resp1.cens))
colnames(resp2.temp) <- gsub("PF2341066", "PF-2341066", colnames(resp2.temp))

# sunchronize all the names
names1 <- strsplit(colnames(resp1.temp), ",")
names2 <- strsplit(colnames(resp2.temp), ",")
names3 <- strsplit(resp3$names, ",")

fnames1 <- sapply(names1, function(s) {
  id2 <- which(sapply(names2, function(n2) {any(n2 %in% s)}))
  id3 <- which(sapply(names3, function(n3) {any(n3 %in% s)}))
  if(length(id2)==0) {id2 <- NA}
  if(length(id3)==0) {id3 <- NA}
  sort(unique(trimws(c(s, names2[[id2]], names3[[id3]]))))})
fnames2 <- sapply(names2, function(s) {
  id1 <- which(sapply(names1, function(n1) {any(n1 %in% s)}))
  id3 <- which(sapply(names3, function(n3) {any(n3 %in% s)}))
  if(length(id1)==0) {id1 <- NA}
  if(length(id3)==0) {id3 <- NA}
  sort(unique(trimws(c(s, names1[[id1]], names3[[id3]]))))})
fnames3 <- sapply(names3, function(s) {
  id2 <- which(sapply(names2, function(n2) {any(n2 %in% s)}))
  id1 <- which(sapply(names1, function(n1) {any(n1 %in% s)}))
  if(length(id2)==0) {id2 <- NA}
  if(length(id1)==0) {id1 <- NA}
  sort(unique(trimws(c(s, names2[[id2]], names1[[id1]]))))})
colnames(resp1.temp) <- sapply(fnames1, paste, collapse=",")
colnames(resp1.cens) <- sapply(fnames1, paste, collapse=",")
colnames(resp2.temp) <- sapply(fnames2, paste, collapse=",")

# cleaning data 2
resp2.temp[resp2.temp==""] <- NA
resp2.temp[resp2.temp==" "] <- NA
resp2.temp[resp2.temp=="NoFit"] <- NA
resp2.cens <- matrix(0, nrow=nrow(resp2.temp), ncol=ncol(resp2.temp))
dimnames(resp2.cens) <- dimnames(resp2.temp)
resp2.cens[substr(resp2.temp, 1, 1)==">" & !is.na(resp2.temp)] <- 1
resp2.temp[substr(resp2.temp, 1, 1)==">" & !is.na(resp2.temp)] <-
  substr(resp2.temp[substr(resp2.temp, 1, 1)==">" & !is.na(resp2.temp)],
         2, nchar(resp2.temp[substr(resp2.temp, 1, 1)==">" & 
                               !is.na(resp2.temp)]))
resp2.cens[substr(resp2.temp, 1, 1)=="<" & !is.na(resp2.temp)] <- -1
resp2.temp[substr(resp2.temp, 1, 1)=="<" & !is.na(resp2.temp)] <-
  substr(resp2.temp[substr(resp2.temp, 1, 1)=="<" & !is.na(resp2.temp)],
         2, nchar(resp2.temp[substr(resp2.temp, 1, 1)=="<" & 
                               !is.na(resp2.temp)]))
resp2.temp[substr(resp2.temp, 1, 1)=="~" & !is.na(resp2.temp)] <-
  substr(resp2.temp[substr(resp2.temp, 1, 1)=="~" & !is.na(resp2.temp)],
         2, nchar(resp2.temp[substr(resp2.temp, 1, 1)=="~" & 
                               !is.na(resp2.temp)]))
class(resp2.temp) <- "numeric"

# combing the data
alldrugs <- unique(c(colnames(resp1.temp), colnames(resp2.temp)))
resp1.temp <- cbind(resp1.temp, matrix(
  NA, ncol=sum(!(alldrugs %in% colnames(resp1.temp))), nrow=nrow(resp1.temp),
  dimnames=list(rownames(resp1.temp), 
                alldrugs[!(alldrugs %in% colnames(resp1.temp))])))
resp1.temp <- resp1.temp[, order(colnames(resp1.temp))]
resp1.cens <- cbind(resp1.cens, matrix(
  NA, ncol=sum(!(alldrugs %in% colnames(resp1.cens))), nrow=nrow(resp1.cens),
  dimnames=list(rownames(resp1.cens), 
                alldrugs[!(alldrugs %in% colnames(resp1.cens))])))
resp1.cens <- resp1.cens[, order(colnames(resp1.cens))]
resp2.temp <- cbind(resp2.temp, matrix(
  NA, ncol=sum(!(alldrugs %in% colnames(resp2.temp))), nrow=nrow(resp2.temp),
  dimnames=list(rownames(resp2.temp), 
                alldrugs[!(alldrugs %in% colnames(resp2.temp))])))
resp2.temp <- resp2.temp[, order(colnames(resp2.temp))]
resp2.cens <- cbind(resp2.cens, matrix(
  NA, ncol=sum(!(alldrugs %in% colnames(resp2.cens))), nrow=nrow(resp2.cens),
  dimnames=list(rownames(resp2.cens), 
                alldrugs[!(alldrugs %in% colnames(resp2.cens))])))
resp2.cens <- resp2.cens[, order(colnames(resp2.cens))]

allcelllines <- unique(c(rownames(resp1.temp), rownames(resp1.temp)))
resp.temp <- t(sapply(allcelllines, function(s) {
  if(!(s %in% rownames(resp1.temp))) {
    resp2.temp[s, ]
  } else if(!(s %in% rownames(resp2.temp))) {
    resp1.temp[s, ]
  } else {
    sapply(alldrugs, function(d) {
      if(is.na(resp1.cens[s, d]) & is.na(resp2.cens[s, d])) {
        NA
      } else if(is.na(resp1.cens[s, d]) | is.na(resp2.cens[s, d])) {
        mean(c(resp1.temp[s, d], resp2.temp[s, d]), na.rm=TRUE)
      } else if(resp1.cens[s, d]==0 & resp2.cens[s, d]==0) {
        mean(c(resp1.temp[s, d], resp2.temp[s, d]), na.rm=TRUE)
      } else if(resp1.cens[s, d]!=0 & resp2.cens[s, d]==0) {
        mean(c(resp1.temp[s, d], resp2.temp[s, d]), na.rm=TRUE)
      } else if(resp2.cens[s, d]!=0 & resp1.cens[s, d]==0) {
        mean(c(resp1.temp[s, d], resp2.temp[s, d]), na.rm=TRUE)
      } else if(resp1.cens[s, d]==-1 & resp2.cens[s, d]==1) {
        mean(c(resp1.temp[s, d], resp2.temp[s, d]), na.rm=TRUE)
      } else if(resp1.cens[s, d]==1 & resp2.cens[s, d]==-1) {
        mean(c(resp1.temp[s, d], resp2.temp[s, d]), na.rm=TRUE)
      } else if(resp1.cens[s, d]==1 & resp2.cens[s, d]==1) {
        max(c(resp1.temp[s, d], resp2.temp[s, d]), na.rm=TRUE)
      } else if(resp1.cens[s, d]==-1 & resp2.cens[s, d]==-1) {
        min(c(resp1.temp[s, d], resp2.temp[s, d]), na.rm=TRUE)
      }
    })
  }
}))

# create list of log IC50 values
resp.temp <- resp.temp[order(rownames(resp.temp)), ]
resp.prep <- apply(resp.temp, 2, function(s) {log(s[!is.na(s)])})
rm(allcelllines, alldrugs, fnames1, fnames2, fnames3, names1, names2, 
   names3, resp.temp, resp1.cens, resp1.temp, resp2, resp2.cens, 
   resp2.temp, resp1, resp3, resp4)

################################# expression ###################################
expr <- read.table("data/CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz",
                   header=TRUE, stringsAsFactors=FALSE)

# expression data
expr.temp <- t(expr[, -c(1, 2)])
colnames(expr.temp) <- sapply(strsplit(expr$gene_id, ".", fixed=TRUE), "[[", 1)
rownames(expr.temp)[substr(rownames(expr.temp), 1, 1)=="X"] <- 
  substr(rownames(expr.temp)[substr(rownames(expr.temp), 1, 1)=="X"], 2, 
         nchar(rownames(expr.temp)[substr(rownames(expr.temp), 1, 1)=="X"]))
rownames(expr.temp) <- sapply(strsplit(rownames(expr.temp), "_"), "[", 1)
did <- which(rownames(expr.temp)==rownames(expr.temp)[
  duplicated(rownames(expr.temp))])
expr.temp <- rbind(expr.temp[c(1:(min(did) - 1)), ], 
                   matrix(colMeans(expr.temp[did, ]), nrow=1, 
                          dimnames=list(rownames(expr.temp)[did][1], NULL)),
                   expr.temp[c((max(did) + 1):nrow(expr.temp)), ])
expr.temp <- expr.temp[, apply(expr.temp, 2, sd)!=0]
expr.temp <- expr.temp[order(rownames(expr.temp)), ]
expr.prep <- expr.temp
rm(did, expr, expr.temp)

################################# expression ###################################
rm(ensembl)
save(expr.prep, resp.prep, file="data/data_ccle_dat1.Rdata")
