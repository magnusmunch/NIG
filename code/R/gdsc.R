#!/usr/bin/env Rscript

# set working directory to source (this) file (only works with Rstudio)
if(any(grepl("rstudio", search()))) {
  setwd(gsub("(.*)/.*", "\\1", rstudioapi::getActiveDocumentContext()$path))  
}

### installation of package
# if(!("cambridge" %in% installed.packages())) {
#   if(!("devtools" %in% installed.packages())) {
#     install.packages("devtools")
#   }
#   library(devtools)
#   install_github("magnusmunch/cambridge/code", local=FALSE, 
#                  auth_token=Sys.getenv("GITHUB_PAT"))
# }
library(devtools)
install_github("magnusmunch/cambridge/code", local=FALSE, 
               auth_token=Sys.getenv("GITHUB_PAT"))

### libraries
library(cambridge)
library(gdata)
library(rcdk)

################################## analysis 1 ##################################
# loading the data either from a local file or from the GDSC website
if(file.exists("../../data/v17.3_fitted_dose_response.csv")) {
  resp1 <- read.table("../../data/v17.3_fitted_dose_response.csv", header=TRUE, 
                      sep=",", stringsAsFactors=FALSE)
} else {
  resp1 <- read.table("ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/v17.3_fitted_dose_response.xlsx",
                      header=TRUE, sep=",", stringsAsFactors=FALSE)
}
if(file.exists("../../data/sanger1018_brainarray_ensemblgene_rma.txt.gz")) {
  expr1 <- read.table("../../data/sanger1018_brainarray_ensemblgene_rma.txt.gz", 
                      header=TRUE, stringsAsFactors=FALSE)
} else {
  expr1 <- read.table(gcon(url("ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/sanger1018_brainarray_ensemblgene_rma.txt.gz"), 
                           text=TRUE), header=TRUE, stringsAsFactors=FALSE)
}
if(file.exists("../../data/Cell_Lines_Details.xlsx")) {
  cell1 <- read.xls("../../data/Cell_Lines_Details.xlsx", 
                    stringsAsFactors=FALSE)
} else {
  cell1 <- read.xls("ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/Cell_Lines_Details.xlsx", 
                    stringsAsFactors=FALSE)
}
if(file.exists("../../data/Screened_Compounds.xlsx")) {
  drug1 <- read.xls("../../data/Screened_Compounds.xlsx", 
                    stringsAsFactors=FALSE)
} else {
  drug1 <- read.xls("ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/Screened_Compounds.xlsx", 
                    stringsAsFactors=FALSE)
}

# load data from "A landscape of pharmacogenomic interaction in cancer" analysis
if(file.exists("../../data/TableS1F.xlsx")) {
  drug2 <- read.xls("../../data/TableS1F.xlsx", stringsAsFactors=FALSE, skip=1)
} else {
  drug2 <- read.xls("https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/suppData/TableS1F.xlsx", 
                    stringsAsFactors=FALSE, skip=1)
}
if(file.exists("../../data/TableS1G.xlsx")) {
  drug3 <- read.xls("../../data/TableS1G.xlsx", stringsAsFactors=FALSE, skip=1)
} else {
  drug3 <- read.xls("https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/suppData/TableS1G.xlsx", 
                    stringsAsFactors=FALSE, skip=1)
}

# fix mistakes with reading in excel file
drug4 <- drug3[, -c(6:25)]
drug4$Cluster.identifier <- sapply(1:nrow(drug4), function(s) {
  if(is.na(drug4$Cluster.identifier[s])) {
    drug4$Cluster.identifier[c(1:s)][tail(which(!is.na(drug4$Cluster.identifier[
      c(1:s)])), n=1)]
  } else {
    drug4$Cluster.identifier[s]
  }})

# take average IC50 if more than one measurement for drug cell line combination
resp2 <- aggregate(LN_IC50 ~ COSMIC_ID + CELL_LINE_NAME + DRUG_NAME, data=resp1, 
                   FUN=mean, na.rm=TRUE)

# put the data in the wide format (cell lines in rows and drugs in columns)
resp3 <- reshape(resp2, timevar="DRUG_NAME", 
                 idvar=c("COSMIC_ID", "CELL_LINE_NAME"), direction="wide")
rownames(resp3) <- resp3$COSMIC_ID
colnames(resp3) <- substr(colnames(resp3), 9, 10000L)
resp4 <- resp3[, -c(1, 2)]

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

# cell lines in rows and genes in columns
expr2 <- t(expr1[, -1])
colnames(expr2) <- expr1$ensembl_gene
rownames(expr2) <- substr(rownames(expr2), 2, 10000L)

# only retain the cell lines that are in both the response and expression set
resp <- resp5[rownames(resp5) %in% rownames(expr2), ]
expr <- expr2[rownames(expr2) %in% rownames(resp5), ]

# remove all in-between steps
rm(list=c("resp1", "resp2", "resp3", "resp4", "resp5", "expr1", "expr2"))

# scaling the response and expressions to mean 0 and variance 1
expr.scale <- scale(expr)
resp.scale <- scale(resp)

# creating some drug covariates
target <- drug1$TARGET[match(colnames(resp.scale), drug1$DRUG_NAME)]
pathway <- drug1$TARGET_PATHWAY[match(colnames(resp.scale), drug1$DRUG_NAME)]
stage <- drug2$Clinical.Stage
action <- drug2$Action
cluster <- drug4$Cluster.identifier

# retrieve CID and SMILES from PUG REST
smiles.list <- vector("list", nrow(drug1))
for(i in 1:nrow(drug1)) {
  cat("\r", "drug", i)
  
  # create a vector of all drug name synonyms
  syn <- c(drug1$DRUG_NAME[i], trimws(unlist(strsplit(drug1$SYNONYMS[i], ","))))
  
  # loop over the synonyms to find SMILES and CID
  sfound <- FALSE
  snum <- 0
  start <- -1
  while(!sfound & (snum < length(syn))) {
    snum <- snum + 1
    
    # see PUG REST website for url explanation
    url <- paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
                 paste(syn[snum], collapse=","), 
                 "/property/CanonicalSMILES/CSV", sep="")
    
    # PUG REST asks for not more than 5 requests per second
    elapsed <- proc.time()[3] - start
    Sys.sleep((0.2 - elapsed)*(elapsed < 0.2))
    start <- proc.time()[3]
    
###################### FIGURE OUT HOW TO READ ERROR CODES ######################
    # error handling if the name is not found in PUG REST (return NA's)
    smiles.list[[i]] <- suppressWarnings(tryCatch(read.table(
      url, sep=",", header=TRUE, stringsAsFactors=FALSE), error=function(e) {
        data.frame(CID=NA, CanonicalSMILES=NA)}))
    closeAllConnections()
    
    # indicator whether a CID and SMILES was found
    sfound <- !is.na(smiles.list[[i]]$CID[1])
    
  }
}
# just picking the first SMILES to use
smiles <- sapply(smiles.list, function(s) {s$CanonicalSMILES[1]})

# parse the smiles to IAtomContainer objects
mols <- parse.smiles(na.omit(smiles))

# fingerprint as binary string indicating occurence of structures
fps <- lapply(mols, get.fingerprint, type='circular')

fingerprint::bit.importance(list(fps[[1]]), fps[-1])
fingerprint::shannon(fps)

# creating clustering based on distance matrix between fingerprints
fp.sim <- fingerprint::fp.sim.matrix(fps, method='tanimoto')
fullclus <- hclust(as.dist(1 - fp.sim))

# choose clustering such that each cluster contains at least 10 drugs
clust <- replace(rep(9, nrow(drug1)), !is.na(smiles), cutree(fullclus, 8))

# trying to figure out how retrieve medicinal use information
library(XML)
url <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Erlotinib/description/XML"
test <- xmlParse(readLines(url), asText=TRUE)

url <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Erlotinib/assaysummary/CSV"
test <- read.table(url, sep=",", header=TRUE, stringsAsFactors=FALSE)
str(test)


# 1.2) Preselecting GEX features: we selected gene expression features 
# explaining 50% variations over all cell lines, since this results in not too 
# features and probably remains those important features. 
# 
# 1.3) Preselecting CNV: we used the same strategy as Garnett et al. (2012) where the GDSC data are from.
# 
# 1.4) Preselecting MUT: we deleted mutated genes with less than two cell lines.



