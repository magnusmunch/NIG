#!/usr/bin/env Rscript

### installation of package
# if(!("cambridge" %in% installed.packages())) {
#   if(!("devtools" %in% installed.packages())) {
#     install.packages("devtools")
#   }
#   library(devtools)
#   install_github("magnusmunch/cambridge/code", local=FALSE, 
#                  auth_token=Sys.getenv("GITHUB_PAT"))
# }

### libraries
detach("package:cambridge", unload=TRUE)
library(cambridge)
library(gdata)
# library(rcdk)
# library(fingerprint)
library(MASS)
library(biomaRt)
detach("package:pInc", unload=TRUE)
library(pInc)
library(glmnet)

################################## analysis 1 ##################################
# loading the data either from a local file or from the GDSC website
if(file.exists("data/v17.3_fitted_dose_response.csv")) {
  resp <- read.table("data/v17.3_fitted_dose_response.csv", header=TRUE,
                     sep=",", stringsAsFactors=FALSE)
} else {
  resp <- read.table("ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/v17.3_fitted_dose_response.xlsx",
                     header=TRUE, sep=",", stringsAsFactors=FALSE)
}
if(file.exists("data/sanger1018_brainarray_ensemblgene_rma.txt.gz")) {
  expr <- read.table("data/sanger1018_brainarray_ensemblgene_rma.txt.gz",
                     header=TRUE, stringsAsFactors=FALSE)
} else {
  expr <- read.table(gcon(url("ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/sanger1018_brainarray_ensemblgene_rma.txt.gz"),
                          text=TRUE), header=TRUE, stringsAsFactors=FALSE)
}
if(file.exists("data/Cell_Lines_Details.xlsx")) {
  cell <- read.xls("data/Cell_Lines_Details.xlsx",
                   stringsAsFactors=FALSE)
} else {
  cell <- read.xls("ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/Cell_Lines_Details.xlsx",
                   stringsAsFactors=FALSE)
}

# reading in the preformatted list of compounds
if(file.exists("data/Screened_Compounds.xlsx")) {
  drug1 <- read.xls("data/Screened_Compounds.xlsx",
                    stringsAsFactors=FALSE)
} else {
  drug1 <- read.xls("ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/Screened_Compounds.xlsx",
                    stringsAsFactors=FALSE)
}

# reading in the curent list of compounds
if(length(list.files(path="data", patt="Drug_list")) > 0) {
  drug2 <- read.table(paste0("data/", list.files(path="data",
                                                 patt="Drug_list")[1]),
                      header=TRUE, sep=",", quote='"', stringsAsFactors=FALSE)
} else {
  drug2 <- read.table("https://www.cancerrxgene.org/translation/drug_list?list=all&sEcho=6&iColumns=6&sColumns=&iDisplayStart=50&iDisplayLength=10&mDataProp_0=0&mDataProp_1=1&mDataProp_2=2&mDataProp_3=3&mDataProp_4=4&mDataProp_5=5&sSearch=&bRegex=false&sSearch_0=&bRegex_0=false&bSearchable_0=true&sSearch_1=&bRegex_1=false&bSearchable_1=true&sSearch_2=&bRegex_2=false&bSearchable_2=true&sSearch_3=&bRegex_3=false&bSearchable_3=true&sSearch_4=&bRegex_4=false&bSearchable_4=true&sSearch_5=&bRegex_5=false&bSearchable_5=true&iSortCol_0=0&sSortDir_0=asc&iSortingCols=1&bSortable_0=true&bSortable_1=true&bSortable_2=true&bSortable_3=true&bSortable_4=true&bSortable_5=true&export=csv",
                      header=TRUE, sep=",", quote='"', stringsAsFactors=FALSE)
}

# load data from "A landscape of pharmacogenomic interaction in cancer" analysis
if(file.exists("data/TableS1F.xlsx")) {
  drug3 <- read.xls("data/TableS1F.xlsx", stringsAsFactors=FALSE, skip=1)[, -9]
} else {
  drug3 <- read.xls("https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/suppData/TableS1F.xlsx",
                    stringsAsFactors=FALSE, skip=1)[, -9]
}
if(file.exists("data/TableS1G.xlsx")) {
  drug4 <- read.xls("data/TableS1G.xlsx", stringsAsFactors=FALSE, skip=1)[, 1:5]
} else {
  drug4 <- read.xls("https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/suppData/TableS1G.xlsx",
                    stringsAsFactors=FALSE, skip=1)[, 1:5]
}
drug4$Cluster.identifier <- sapply(1:nrow(drug4), function(s) {
  if(is.na(drug4$Cluster.identifier[s])) {
    drug4$Cluster.identifier[c(1:s)][tail(which(!is.na(drug4$Cluster.identifier[
      c(1:s)])), n=1)]
  } else {
    drug4$Cluster.identifier[s]
  }})

# combining and removing double drugs and rescreens
drug <- merge(merge(merge(drug1, drug2, by.x="DRUG_ID", by.y="drug_id",
                          all=TRUE), drug3, by.x="DRUG_ID", by.y="Identifier",
                    all=TRUE), drug4, by.x="Name.y", by.y="Drug", all=TRUE)
drug.temp <- drug[c(4:5, 12:13, 30:33, 34:37, 43:46, 49:50, 54:55,
                    56:59, 76:79, 108:109, 118:119, 128:129, 136:139, 178:179,
                    200:201, 211:212, 249:252, 16, 144), ]
rownames(drug.temp) <- 1:nrow(drug.temp)
drug.temp[2, c(20:23)] <- drug.temp[1, c(20:23)]
drug.temp[4, ] <- drug.temp[3, ]
drug.temp[8, ] <- drug.temp[5, ]
drug.temp[12, ] <- drug.temp[9, ]
drug.temp[12, 14] <- paste0(drug.temp[12, 14], ", ", drug.temp[11, 14])
drug.temp[16, ] <- drug.temp[13, ]
drug.temp[16, 14] <- paste0(drug.temp[16, 14], ", ", drug.temp[15, 14])
drug.temp[17, 11] <- drug.temp[18, 11]
drug.temp[18, ] <- drug.temp[17, ]
drug.temp[20, ] <- drug.temp[19, ]
drug.temp[24, ] <- drug.temp[21, ]
drug.temp[28, ] <- drug.temp[25, ]
drug.temp[29, 18] <- drug.temp[30, 18]
drug.temp[30, ] <- drug.temp[29, ]
drug.temp[31, 14] <- paste0(drug.temp[31, 14], ", ", drug.temp[32, 14])
drug.temp[32, ] <- drug.temp[31, ]
drug.temp[34, 4] <- paste0(drug.temp[34, 4], ", ", drug.temp[33, 1])
drug.temp[34, c(20:23)] <- drug.temp[33, c(20:23)]
drug.temp[38, ] <- drug.temp[35, ]
drug.temp[39, 17] <- drug.temp[40, 17]
drug.temp[40, ] <- drug.temp[39, ]
drug.temp[41, 15] <- drug.temp[42, 15]
drug.temp[42, ] <- drug.temp[41, ]
drug.temp[44, ] <- drug.temp[43, ]
drug.temp[50, c(2, 17, 18, 20, 23)] <- drug.temp[49, c(2, 17, 18, 20, 23)]
drug.temp <- drug.temp[-c(1, 3, 5:7, 9:11, 13:15, 17, 19, 21:23, 25:27, 29, 31,
                          33, 35:37, 39, 41, 43, 45:47, 49), ]
drug.temp <- rbind(drug[-c(4:5, 12:13, 30:33, 34:37, 43:46, 49:50, 54:55, 56:59,
                           76:79, 108:109, 118:119, 128:129, 136:139, 178:179,
                           200:201, 211:212, 249:252, 16, 144), ], drug.temp)
rownames(drug.temp) <- 1:nrow(drug.temp)

# create primary names
name <- drug.temp$DRUG_NAME
name[is.na(name)] <- drug.temp$Name.y[is.na(name)]

# combining synonyms
drug.temp$Brand.name[c(61, 71, 92)] <-
  gsub(")", "", gsub("(", ", ", drug.temp$Brand.name[c(61, 71, 92)],
                     fixed=TRUE), fixed=TRUE)
drug.temp$SYNONYMS[20] <- gsub(",", ", ", drug.temp$SYNONYMS[20])
drug.temp$SYNONYMS[92] <-
  gsub(")", "", gsub("(", ", ", drug.temp$SYNONYMS[92], fixed=TRUE), fixed=TRUE)
drug.temp$SYNONYMS[142] <- c("OSI 930, OSI930")
drug.temp$SYNONYMS[212] <- c("MK 0457, MK-0457, MK-045, VX-680, VX 680, VX-68")
drug.temp$Synonyms.x[92] <- "GW843682X, AN-13"
drug.temp$Synonyms.x[212] <- c("MK 0457, MK-0457, MK-045, VX-680, VX 680,
                               VX-68")
drug.temp$Synonyms.y[63] <- "KIN001-013, GNF-2, 3-(6-(4-(trifluoromethoxy)phenylamino)pyrimidin-4-yl)benzamide)"
drug.temp$Synonyms.y[101] <- "KIN001-205"
drug.temp$Synonyms.y[167] <- "2-(6,7-dimethoxyquinazolin-4-yl)-5-(pyridin-2-yl)-2H-1,2,4-triazol-3-amine"
drug.temp$Synonyms.y[-c(132, 167, 170, 176)] <-
  gsub(",", ", ", drug.temp$Synonyms.y[-c(132, 167, 170, 176)])
drug.temp$Brand.name <- gsub(",", ", ", drug.temp$Brand.name)
synonyms <- sapply(1:nrow(drug.temp), function(d) {
  vec <-c(strsplit(drug.temp$Name.y, ", ")[[d]],
          strsplit(drug.temp$DRUG_NAME, ", ")[[d]],
          strsplit(drug.temp$SYNONYMS, ", ")[[d]],
          strsplit(drug.temp$Name.x, ", ")[[d]],
          strsplit(drug.temp$Synonyms.x, ", ")[[d]],
          strsplit(drug.temp$Synonyms.y, ", ")[[d]],
          strsplit(drug.temp$Brand.name, ", ")[[d]])
  vec <- na.omit(vec)
  vec <- vec[vec!=""]
  vec <- trimws(vec)
  vec <- ifelse(substring(vec, 1, 1)==",",
                substring(vec, 2, nchar(vec)), vec)
  vec <- ifelse(substring(vec, nchar(vec), nchar(vec))==",",
                substring(vec, 1, nchar(vec) - 1), vec)
  vec <- sort(unique(vec))
  vec <- vec[vec!=name[d]]
  vec <- paste(vec, collapse=", ")
  vec[vec==""] <- NA
  return(vec)})

# combining targets
targets <- sapply(1:nrow(drug.temp), function(d) {
  vec <-c(strsplit(drug.temp$TARGET, ", ")[[d]],
          strsplit(drug.temp$Targets, ", ")[[d]],
          strsplit(drug.temp$Putative.Target, ", ")[[d]],
          strsplit(drug.temp$Target, ", ")[[d]])
  vec[vec=="unknown" | vec=="not defined" | vec=="others"] <- NA
  vec <- na.omit(vec)
  vec <- vec[vec!=""]
  vec <- trimws(vec)
  vec <- ifelse(substring(vec, 1, 1)==",",
                substring(vec, 2, nchar(vec)), vec)
  vec <- ifelse(substring(vec, nchar(vec), nchar(vec))==",",
                substring(vec, 1, nchar(vec) - 1), vec)
  vec <- sort(unique(vec))
  vec <- paste(vec, collapse=", ")
  vec[vec==""] <- NA
  return(vec)})

# combining pathways
pathways <- sapply(1:nrow(drug.temp), function(d) {
  vec <-c(strsplit(drug.temp$TARGET_PATHWAY, ", ")[[d]],
          strsplit(drug.temp$Target.pathway, ", ")[[d]],
          strsplit(drug.temp$Targeted.process.pathway, ", ")[[d]],
          strsplit(drug.temp$Targeted.process.or.pathway, ", ")[[d]])
  vec <- gsub("\\.", " ", vec)
  vec <- gsub("  ", " ", vec)
  vec[vec=="Other" | vec=="other"] <- NA
  vec <- na.omit(vec)
  vec <- vec[vec!=""]
  vec <- trimws(vec)
  vec <- ifelse(substring(vec, 1, 1)==",",
                substring(vec, 2, nchar(vec)), vec)
  vec <- ifelse(substring(vec, nchar(vec), nchar(vec))==",",
                substring(vec, 1, nchar(vec) - 1), vec)
  vec[vec=="chromain histone acetylation"] <- "Chromatin histone acetylation"
  vec[vec=="cytoskeleton"] <- "Cytoskeleton"
  vec[vec=="JNK and p38 signaling"] <- "JNK signaling, p38 signaling"
  vec[vec=="metabolism"] <- "Metabolism"
  vec[vec=="mitosis"] <- "Mitosis"
  vec[vec=="Unclassified"] <- NA
  vec <- vec[!duplicated(tolower(vec))]
  vec <- sort(vec)
  vec <- paste(vec, collapse=", ")
  vec[vec==""] <- NA
  return(vec)})

# making pubchem id into integer
pubchem <- drug.temp$PubCHEM
pubchem[pubchem=="none" | pubchem=="several"] <- NA
pubchem <- as.integer(pubchem)

# change action variable
action <- drug.temp$Action
action[action=="not defined"] <- NA

# combining the data
drug <- data.frame(name=name, synonyms=synonyms, id=drug.temp$DRUG_ID,
                   pubchem=pubchem, targets=targets, pathways=pathways,
                   nsamples=drug.temp$Sample.Size, count=drug.temp$Count,
                   stage=drug.temp$Clinical.Stage, action=action,
                   cluster=drug.temp$Cluster.identifier,
                   silhouette=drug.temp$Silhouette.Width,
                   stringsAsFactors=FALSE)
drug <- drug[order(drug$id), ]
rownames(drug) <- 1:nrow(drug)
rm(drug1, drug2, drug3, drug4, drug.temp, synonyms, name, targets, pathways,
   action, pubchem)

# take average IC50 if more than one measurement for drug cell line combination
resp.temp <- aggregate(LN_IC50 ~ COSMIC_ID + CELL_LINE_NAME + DRUG_NAME,
                       data=resp, FUN=mean, na.rm=TRUE)

# put the data in the wide format (cell lines in rows and drugs in columns)
resp.temp <- reshape(resp.temp, timevar="DRUG_NAME",
                     idvar=c("COSMIC_ID", "CELL_LINE_NAME"), direction="wide")
colnames(resp.temp) <- c("id", "cellline",
                         substr(colnames(resp.temp[-c(1, 2)]), 9, 10000L))
resp.temp <- resp.temp[order(resp.temp$id), ]
rownames(resp.temp) <- sort(resp.temp$id)
resp <- resp.temp[, -c(1, 2)]
rm(resp.temp)

# # sequentially remove drug or cell line with highest proportion missings
# nomiss <- FALSE
# while(!nomiss) {
#   pcell <- apply(resp, 1, function(x) {sum(is.na(x))/length(x)})
#   pdrug <- apply(resp, 2, function(x) {sum(is.na(x))/length(x)})
#   
#   if(max(pcell) >= max(pdrug)) {
#     resp <- resp[-which.max(pcell), ]
#   } else {
#     resp <- resp[, -which.max(pdrug)]
#   }
#   nomiss <- !any(is.na(resp))
# }
# rm(pcell, pdrug, nomiss)
# 
# # cell lines in rows and genes in columns
# expr.temp <- t(expr[, -1])
# colnames(expr.temp) <- expr$ensembl_gene
# rownames(expr.temp) <- substr(rownames(expr.temp), 2, 10000L)
# expr <- expr.temp[order(as.numeric(rownames(expr.temp))), ]
# rm(expr.temp)
# 
# # only retain the cell lines that are available in both response and expression
# expr.prep <- expr[rownames(expr) %in% rownames(resp), ]
# resp.prep <- resp[rownames(resp) %in% rownames(expr), ]
# 
# # keep drugs that are available in both response and drug data
# drug.prep <- drug[drug$name %in% colnames(resp.prep), ]
# drug.prep <- drug.prep[order(drug.prep$name), ]
# 
# # save(drug.prep, expr.prep, resp.prep, file="data/data_gdsc_dat1.Rdata")
load(file="data/data_gdsc_dat1.Rdata")

# selecting genes
# fit.enet <- apply(resp.scale, 2, function(y) {
#   glmnet(expr.scale, y, alpha=0.01)})
# idsel <- lapply(fit.enet, function(s) {
#   which(as.numeric(coef(s, s=s$lambda[
#     tail(which(s$df < 300), n=1)]))[-1]!=0)})
# save(fit.enet, idsel, file="results/data_gdsc_fit1.Rdata")
# load(file="results/data_gdsc_fit1.Rdata")
# expr.sel <- sapply(idsel, function(id) {expr.prep[, id]})

### creating extra external covariates
# retrieve CID and SMILES from PUG REST
# smiles.list <- vector("list", nrow(drug.prep))
# for(i in 1:nrow(drug.prep)) {
#   cat("\r", "drug", i)
#   
#   # create a vector of all drug name synonyms
#   syn <- na.omit(c(drug.prep$name[i], drug.prep$synonyms[i]))
#   
#   # loop over the synonyms to find SMILES and CID
#   sfound <- FALSE
#   snum <- 0
#   start <- -1
#   while(!sfound & (snum < length(syn))) {
#     snum <- snum + 1
#     
#     smiles.list[[i]] <- NA
#     
#     # see PUG REST website for url explanation
#     url <- paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
#                  paste(syn[snum], collapse=","), 
#                  "/property/CanonicalSMILES/CSV", sep="")
#     
#     
#     # PUG REST asks for not more than 5 requests per second
#     elapsed <- proc.time()[3] - start
#     Sys.sleep((0.2 - elapsed)*(elapsed < 0.2))
#     start <- proc.time()[3]
#     
#     # error handling if the name is not found in PUG REST (return NA's)
#     smiles.list[[i]] <- suppressWarnings(tryCatch(read.table(
#       url, sep=",", header=TRUE, stringsAsFactors=FALSE), error=function(e) {
#         data.frame(CID=NA, CanonicalSMILES=NA)}))
#     closeAllConnections()
#     
#     # indicator whether a CID and SMILES was found
#     sfound <- !is.na(smiles.list[[i]]$CID[1])
#     
#   }
#   if(all(is.na(smiles.list[[i]])) & !is.na(drug.prep$pubchem[i])) {
#     
#     url <- paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
#                  paste(drug.prep$pubchem[i], collapse=","), 
#                  "/property/CanonicalSMILES/CSV", sep="")
#     smiles.list[[i]] <- suppressWarnings(tryCatch(read.table(
#       url, sep=",", header=TRUE, stringsAsFactors=FALSE), error=function(e) {
#         data.frame(CID=NA, CanonicalSMILES=NA)}))
#   }
# }
# 
# # just picking the first SMILES to use
# smiles <- sapply(smiles.list, function(s) {s$CanonicalSMILES[1]})
# 
# # parse the smiles to IAtomContainer objects
# mols <- parse.smiles(smiles[!is.na(smiles)])
# 
# # fingerprint as binary string indicating occurence of structures
# fps <- lapply(mols, get.fingerprint, type='circular')
# 
# # sum of bit importances per drug
# bimp <- sapply(1:length(fps), function(s) {
#   sum(bit.importance(list(fps[[s]]), fps[-s]), na.rm=TRUE)})
# 
# # creating clustering based on distance matrix between fingerprints
# fp.sim <- fp.sim.matrix(fps, method='tanimoto')
# fullclus <- hclust(as.dist(1 - fp.sim))
# cluster <- replace(rep(5, nrow(drug.prep)), !is.na(smiles), 
#                    cutree(fullclus, k=4))

### find genes that are in the target pathway of the drug
# hand coded pathway names to KEGG and reactome ids
pathwayid <- data.frame(
  pathway=c('kinases', 'DNA replication', 'EGFR signaling', 'PI3K signaling', 
            'PI3K/MTOR signaling', 'RTK signaling', 
            'Chromatin histone acetylation', 
            'JNK signaling', 'p38 signaling', 'Cell cycle', 'Genome integrity', 
            'TOR signaling', 'Hormone-related', 
            'ERK MAPK signaling', 'IGFR signaling', 'ABL signaling', 
            'WNT signaling', 'Metabolism', 'Mitosis', 
            'Protein stability and degradation', 'Apoptosis regulation', 
            'Cytoskeleton', 'Chromatin other', 'p53 pathway', 
            'Chromatin histone methylation'),
  keggid=c(NA, "map03030", NA, "map04151", "map04150", NA, NA, "N00446", 
           "N00444", "map04110", NA, "map04150", NA, "map04010", "N00182", NA, 
           "map04310", NA, NA, NA, "map04210", "map04810", NA, "map04115", NA),
  reactomeid=c(NA, "R-HSA-69306", "R-HSA-177929", "R-HSA-109704", 
               "R-HSA-165159", "R-HSA-9006934", NA, NA, NA, "R-HSA-1640170", 
               NA, "R-HSA-165159", NA, "R-HSA-5683057", "R-HSA-2428924", NA, 
               "R-HSA-195721", "R-HSA-1430728", NA, NA, "R-HSA-169911", NA, 
               "R-HSA-4839726", "R-HSA-5633007", NA))
notes <- c("RTK is a family, many of which have identified pathways",
           "ERK is the same as MAPK",
           "IGFR is the same as IGF1R",
           "Cytoskeleton consists of actin and tubulin",
           "P53 is the same as TP53")

# retrieve gene lists of all reactome pathways
# mapid <- read.table("https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt",
#                     comment.char="", sep="\t", quote="",
#                     col.names=c("database", "reactomeid", "url", "name", 
#                                 "evidencecode", "species"))
mapid <- read.table("data/Ensembl2Reactome_All_Levels.txt",
                    comment.char="", sep="\t", quote="",
                    col.names=c("ensemblid", "reactomeid", "url", "name", 
                                "evidencecode", "species"))

# link drugs to ensembl symbols in their pathway
drug.reactomeid <- sapply(drug.prep$pathways, function(s) {
  pathwayid$reactomeid[match(strsplit(s, ", ")[[1]], pathwayid$pathway)]})
drug.ensemblid <- sapply(1:nrow(drug.prep), function(d) {
  vec <- mapid$ensemblid[mapid$reactomeid %in% drug.reactomeid[d]]
  if(length(vec)==0) {NA} else {vec}})

# determine for every feature whether it is in the target pathway or not
# inpathway <- sapply(1:nrow(drug.prep), function(d) {
#   as.numeric(colnames(expr.sel[[d]]) %in% drug.ensemblid[[d]])})
inpathway <- lapply(1:nrow(drug.prep), function(d) {
  as.numeric(colnames(expr.prep) %in% drug.ensemblid[[d]])})

psel <- 500
o <- order(-apply(expr.prep, 2, sd))
# idsel <- lapply(inpathway, function(s) {
#   m1 <- o[o %in% which(s==1)]; id1 <- head(m1, n=min(length(m1), psel/2))
#   m0 <- o[o %in% which(s==0)]; id0 <- head(m0, n=psel - length(id1))
#   c(id1, id0)})
# expr.sel <- lapply(idsel, function(s) {expr.prep[, s]})
idsel <- o[1:psel]
expr.sel <- lapply(1:D, function(d) {expr.prep[, idsel]})
inpathway <- lapply(inpathway, function(s) {s[idsel]})

# number of equations, features, and observations
D <- ncol(resp.prep)
p <- sapply(expr.sel, ncol)
n <- nrow(resp.prep)



### model fitting
# estimation settings
control <- list(conv.post=TRUE, trace=TRUE, epsilon.eb=1e-3, epsilon.vb=1e-3,
                maxit.eb=200, maxit.vb=1, maxit.post=100)

# estimation
x <- lapply(expr.sel, function(s) {scale(s)})
y <- scale(resp.prep)
Z <- matrix(1, nrow=D)
colnames(Z) <- c("intercept")
C <- lapply(inpathway, function(s) {
  s <- matrix(1, nrow=length(s)); colnames(s) <- c("intercept"); return(s)})
fit1.semnig <- semnig(x=x, y=y, C=C, Z=Z, unpenalized=NULL,
                      standardize=FALSE, intercept=FALSE, fixed.eb="none",
                      full.post=TRUE, init=NULL, control=control)
# create co-data matrices
Z <- model.matrix(~ replace(drug.prep$stage, is.na(drug.prep$stage),
                            "experimental") + replace(drug.prep$action, is.na(
                              drug.prep$action), "unknown"))
colnames(Z) <- c("intercept", "experimental", "in clinical development",
                 "targeted", "unknown")
C <- lapply(inpathway, function(s) {
  s <- cbind(1, s); colnames(s) <- c("intercept", "in pathway"); return(s)})
fit2.semnig <- semnig(x=x, y=y, C=C, Z=Z, unpenalized=NULL,
                      standardize=FALSE, intercept=FALSE, fixed.eb=FALSE,
                      full.post=TRUE, init=NULL, control=control)
fit3.semnig <- semnig(x=x, y=y, C=C, Z=NULL, unpenalized=NULL,
                      standardize=FALSE, intercept=FALSE, fixed.eb=FALSE,
                      full.post=TRUE, init=NULL, control=control)
# fit1.bSEM <- bSEM(rep(list(x), D), split(y, 1:ncol(y)),
#                   lapply(inpathway, "+", 1),
#                   control=list(maxit=200, trace=TRUE, epsilon=1e-3))

# save(C, Z, fit1.semnig, fit2.semnig, fit1.bSEM,
#      file="results/data_gdsc_fit2.Rdata")
load("results/data_gdsc_fit2.Rdata")

### cross-validation
# estimation settings
control <- list(conv.post=TRUE, trace=FALSE, epsilon.eb=1e-3, epsilon.vb=1e-3, 
                maxit.eb=200, maxit.vb=1, maxit.post=100)

# multiple splits
nreps <- 50
ntrain <- floor(nrow(resp.prep)/2)
# apmse <- matrix(NA, nrow=nreps, ncol=5)
# colnames(apmse) <- c("NIG1", "NIG2", "bSEM", "lasso", "ridge")
# brank <- matrix(NA, nrow=nreps, ncol=sum(p)*5)
# colnames(brank) <- paste0(rep(c("NIG1", "NIG2", "bSEM", "lasso", "ridge"), 
#                               each=sum(p)), ".", rep(paste0(
#                                 rep(1:D, times=p), ".",
#                                 unlist(sapply(p, function(s) {
#                                   return(1:s)}))), 5))
# elbo <- matrix(NA, ncol=3, nrow=nreps)
# colnames(elbo) <- c("NIG1", "NIG2", "bSEM")
apmse <- matrix(NA, nrow=nreps, ncol=4)
colnames(apmse) <- c("NIG1", "NIG2", "NIG3", "NIG4")
brank <- matrix(NA, nrow=nreps, ncol=sum(p)*4)
colnames(brank) <- paste0(rep(c("NIG1", "NIG2", "NIG3", "NIG4"), 
                              each=sum(p)), ".", rep(paste0(
                                rep(1:D, times=p), ".",
                                unlist(sapply(p, function(s) {
                                  return(1:s)}))), 4))
elbo <- matrix(NA, ncol=4, nrow=nreps)
colnames(elbo) <- c("NIG1", "NIG2", "NIG3", "NIG4")

for(r in 1:nreps) {
  cat("\r", "replication", r)
  set.seed(2019 + r)
  
  idtrain <- sample(1:nrow(resp.prep), ntrain)
  xtrain <- lapply(expr.sel, function(s) {scale(s[idtrain, ])})
  ytrain <- scale(resp.prep[idtrain, ])
  xtest <- lapply(expr.sel, function(s) {scale(s[-idtrain, ])})
  ytest <- scale(resp.prep[-idtrain, ])
  
  Z <- matrix(1, nrow=D)
  colnames(Z) <- c("intercept")
  C <- lapply(inpathway, function(s) {
    s <- matrix(1, nrow=length(s)); colnames(s) <- c("intercept"); return(s)})
  cv1.semnig <- semnig(x=xtrain, y=ytrain, C=C, Z=Z, unpenalized=NULL, 
                        standardize=FALSE, intercept=FALSE, fixed.eb="none", 
                        full.post=TRUE, init=NULL, control=control)
  # create co-data matrices
  Z <- model.matrix(~ replace(drug.prep$stage, is.na(drug.prep$stage),
                              "experimental") + replace(drug.prep$action, is.na(
                                drug.prep$action), "unknown"))
  colnames(Z) <- c("intercept", "experimental", "in clinical development", 
                   "targeted", "unknown")
  C <- lapply(inpathway, function(s) {
    s <- cbind(1, s); colnames(s) <- c("intercept", "in pathway"); return(s)})
  cv2.semnig <- semnig(x=xtrain, y=ytrain, C=C, Z=Z, unpenalized=NULL, 
                        standardize=FALSE, intercept=FALSE, fixed.eb="none", 
                        full.post=TRUE, init=NULL, control=control)
  cv3.semnig <- semnig(x=xtrain, y=ytrain, C=C, Z=NULL, unpenalized=NULL, 
                       standardize=FALSE, intercept=FALSE, fixed.eb="none", 
                       full.post=TRUE, init=NULL, control=control)
  cv4.semnig <- semnig(x=xtrain, y=ytrain, C=NULL, Z=Z, unpenalized=NULL, 
                       standardize=FALSE, intercept=FALSE, fixed.eb="none", 
                       full.post=TRUE, init=NULL, control=control)
  # cv1.bSEM <- bSEM(xtrain, split(ytrain, 1:ncol(ytrain)), 
  #                  lapply(inpathway, "+", 1), 
  #                   control=list(maxit=200, trace=FALSE, epsilon=1e-3))
  
  # cv1.lasso <- lapply(1:D, function(d) {
  #   cv.glmnet(xtrain[[d]], ytrain[, d], intercept=FALSE, standardize=FALSE)})
  # cv1.ridge <- lapply(1:D, function(d) {
  #   cv.glmnet(xtrain[[d]], ytrain[, d], alpha=0, intercept=FALSE, 
  #             standardize=FALSE)})
  
  # apmse[r, ] <- c(mean(sapply(1:D, function(d) {
  #                   mean((ytest[, d] - xtest[[d]] %*%
  #                           cv1.semnig$vb$mpost$beta[[d]])^2)})),
  #                 mean(sapply(1:D, function(d) {
  #                   mean((ytest[, d] - xtest[[d]] %*%
  #                         cv2.semnig$vb$mpost$beta[[d]])^2)})),
  #                 mean(sapply(1:D, function(d) {
  #                   mean((ytest[, d] - xtest[[d]] %*%
  #                         cv1.bSEM$vb$beta[[d]][, 1])^2)})),
  #                 mean((ytest - sapply(1:D, function(d) {
  #                   predict(cv1.lasso[[d]], xtest[[d]], s="lambda.min")}))^2),
  #                 mean((ytest - sapply(1:D, function(d) {
  #                   predict(cv1.ridge[[d]], xtest[[d]], s="lambda.min")}))^2))
  # 
  # brank[r, 1:sum(p)] <- unlist(lapply(cv1.semnig$vb$mpost$beta, function(s) {
  #   rank(abs(s), ties.method="average")}))
  # brank[r, (sum(p) + 1):(2*sum(p))] <-
  #   unlist(lapply(cv2.semnig$vb$mpost$beta, function(s) {
  #     rank(abs(s), ties.method="average")}))
  # brank[r, (2*sum(p) + 1):(3*sum(p))] <-
  #   unlist(lapply(cv1.bSEM$vb$beta, function(s) {
  #     rank(abs(s[, 1]), ties.method="average")}))
  # brank[r, (3*sum(p) + 1):(4*sum(p))] <- unlist(sapply(cv1.lasso, function(s) {
  #   rank(abs(as.numeric(coef(s, s="lambda.min"))[-1]), ties.method="average")}))
  # brank[r, (4*sum(p) + 1):(5*sum(p))] <- unlist(sapply(cv1.ridge, function(s) {
  #   rank(abs(as.numeric(coef(s, s="lambda.min"))[-1]), ties.method="average")}))
  # 
  # elbo[r, ] <- c(mean(cv1.semnig$seq.elbo[nrow(cv1.semnig$seq.elbo), ]),
  #                mean(cv2.semnig$seq.elbo[nrow(cv2.semnig$seq.elbo), ]),
  #                mean(cv1.bSEM$mll[, ncol(cv1.bSEM$mll)]))
  apmse[r, ] <- c(mean(sapply(1:D, function(d) {
    mean((ytest[, d] - xtest[[d]] %*%
            cv1.semnig$vb$mpost$beta[[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((ytest[, d] - xtest[[d]] %*%
              cv2.semnig$vb$mpost$beta[[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((ytest[, d] - xtest[[d]] %*%
              cv3.semnig$vb$mpost$beta[[d]])^2)})),
    mean(sapply(1:D, function(d) {
      mean((ytest[, d] - xtest[[d]] %*%
              cv4.semnig$vb$mpost$beta[[d]])^2)})))
  
  brank[r, 1:sum(p)] <- unlist(lapply(cv1.semnig$vb$mpost$beta, function(s) {
    rank(abs(s), ties.method="average")}))
  brank[r, (sum(p) + 1):(2*sum(p))] <-
    unlist(lapply(cv2.semnig$vb$mpost$beta, function(s) {
      rank(abs(s), ties.method="average")}))
  brank[r, (2*sum(p) + 1):(3*sum(p))] <-
    unlist(lapply(cv3.semnig$vb$mpost$beta, function(s) {
      rank(abs(s), ties.method="average")}))
  brank[r, (3*sum(p) + 1):(4*sum(p))] <- 
    unlist(lapply(cv4.semnig$vb$mpost$beta, function(s) {
      rank(abs(s), ties.method="average")}))
  
  elbo[r, ] <- c(mean(cv1.semnig$seq.elbo[nrow(cv1.semnig$seq.elbo), ]),
                 mean(cv2.semnig$seq.elbo[nrow(cv2.semnig$seq.elbo), ]),
                 mean(cv3.semnig$seq.elbo[nrow(cv3.semnig$seq.elbo), ]),
                 mean(cv4.semnig$seq.elbo[nrow(cv4.semnig$seq.elbo), ]))
                 
}

# save(apmse, brank, elbo, file="results/data_gdsc_cv1.Rdata")
save(apmse, brank, elbo, file="results/data_gdsc_cv2.Rdata")
# load("results/data_gdsc_cv1.Rdata")

# EB estimates
tab <- rbind(c(fit1.semnig$eb$alphad, rep(NA, 4), fit1.semnig$eb$alphaf, NA),
             c(fit2.semnig$eb$alphad, fit2.semnig$eb$alphaf),
             c(rep(NA, 5), fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 1]/
                 fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 1],
               fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 2]/
                 fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 2] -
                 fit1.bSEM$eb$a[nrow(fit1.bSEM$eb$a), 1]/
                 fit1.bSEM$eb$b[nrow(fit1.bSEM$eb$b), 1]))
colnames(tab) <- c(colnames(Z), colnames(C[[1]]))
rownames(tab) <- c("NIG1", "NIG2", "bSEM")

# Spearman correlations between ranks
brankcor <- cbind(NIG1=cor(t(brank[, substr(colnames(brank), 1, 4)=="NIG1"]), 
                           method="spearman")[
                             upper.tri(matrix(0, nrow(brank), nrow(brank)))],
                  NIG2=cor(t(brank[, substr(colnames(brank), 1, 4)=="NIG2"]), 
                           method="spearman")[
                             upper.tri(matrix(0, nrow(brank), nrow(brank)))],
                  bSEM=cor(t(brank[, substr(colnames(brank), 1, 4)=="bSEM"]), 
                           method="spearman")[
                             upper.tri(matrix(0, nrow(brank), nrow(brank)))])

save(apmse, elbo, tab, brankcor, file="results/data_gdsc_res1.Rdata")




x=xtrain
y=split(ytrain, 1:ncol(ytrain))
P=lapply(inpathway, "+", 1)
EBid=NULL
control=list(maxit=200, trace=FALSE, epsilon=1e-3)
bSEM <- function(x, y, P, EBid=NULL,
                 control=list(maxit=100, trace=TRUE, epsilon=1e-3)) {
  
  if(is.null(EBid)) {
    EBid <- sort(unique(unlist(P)))
  }
  
  niter <- control$maxit
  listX <- x
  listy <- y
  listP <- P
  
  # Intialization
  priors <- sort(unique(unlist(listP)))
  nbpriors <- length(priors)
  aprior <- bprior <- matrix(NA, 1, nbpriors)
  colnames(aprior) <- paste("a", 1:nbpriors, sep="")
  colnames(bprior) <- paste("b", 1:nbpriors, sep="")
  aprior[1,] <- bprior[1,] <- 0.001
  idxPriorList <- lapply(listP, function(x){sort(unique(x))})
  allMLs <- matrix(NA, length(listX), 1)
  allMLs[, 1] <- -Inf
  
  conv <- FALSE
  iter <- 0
  while(!conv & (iter < control$maxit)) {
    iter <- j <- iter + 1
    if(control$trace) {cat("iteration ", j, "\n")}
    
    # Prior as lists
    inputaList <- lapply(idxPriorList, function(x){aprior[j,x]})
    inputbList <- lapply(idxPriorList, function(x){bprior[j,x]})
    
    if(j==1){
      mydstarvec <- rep(0, length(inputaList))
      mybstarlist <- lapply(idxPriorList, function(x){rep(0, length(x))})
    } else {
      mydstarvec <- unlist(lapply(res$postSigList, function(x){x[2]}))
      mybstarlist <- lapply(res$postRandList, function(x){x[,2]})
    }
    
    # Fit BSEM
    res <- BSEMVarOneIter(ylist=listy, Xlist=listX, Plist=listP, alist=inputaList, blist=inputbList, bstarlist=mybstarlist, cSigma=0.001, dSigma=0.001, dstarvec=mydstarvec, lincomblist=list())
    
    # Marginal likelihood
    allMLs <- cbind(allMLs, res$allmargs)
    
    # Empirical Bayes
    aprior <- rbind(aprior, NA)
    bprior <- rbind(bprior, NA)
    for(ii in 1:nbpriors){
      if(ii%in%EBid){
        # Get posterior shape and rate parameters
        allaRandStar <- sapply(1:length(idxPriorList), function(x){res$postRandList[[x]][idxPriorList[[x]]==ii,1]}, simplify=TRUE)
        if(is.list(allaRandStar)){
          allaRandStar <- unlist(allaRandStar)
        }
        allbRandStar <- sapply(1:length(idxPriorList), function(x){res$postRandList[[x]][idxPriorList[[x]]==ii,2]}, simplify=TRUE)
        if(is.list(allbRandStar)){
          allbRandStar <- unlist(allbRandStar)
        }
        
        # Variational Empirical Bayes using fixed-point iteration as in Valpola and Honkela (2006)
        ab <- c(aprior[j, ii], bprior[j, ii])
        ab <- fixedPointIterEB(initab=ab, myallaRandStar=allaRandStar, myallbRandStar=allbRandStar, mymaxiter=20, myeps=1e-6)
        if(control$trace) {ab}
        aprior[j+1, ii] <- ab[1]
        bprior[j+1, ii] <- ab[2]
      }
    }
    conv <- all(abs(allMLs[, j] - allMLs[, j + 1]) <= control$epsilon)
  }
  
  out <- list(eb=list(a=aprior, b=bprior),
              vb=list(gamma.sq=lapply(res$postRandList, function(s) {
                c(a=s[1, 1], b=s[1, 2])}),
                sigma.sq=lapply(res$postSigList, function(s) {
                  c(c=s[1, 1], d=s[2, 1])}),
                beta=lapply(res$postBetaList, function(s) {
                  cbind(mu=s[, 1], dSigma=s[, 2]^2)})),
              mll=allMLs, conv=conv, iter=iter)
  return(out)
}