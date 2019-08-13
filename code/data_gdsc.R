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
library(cambridge)
library(gdata)
library(rcdk)
library(fingerprint)
library(glmnet)
library(MASS)
library(biomaRt)
library(rstan)

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
# if(file.exists("data/Cell_Lines_Details.xlsx")) {
#   cell1 <- read.xls("data/Cell_Lines_Details.xlsx", 
#                     stringsAsFactors=FALSE)
# } else {
#   cell1 <- read.xls("ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/Cell_Lines_Details.xlsx", 
#                     stringsAsFactors=FALSE)
# }

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

# sequentially remove drug or cell line with highest proportion missings
nomiss <- FALSE
while(!nomiss) {
  pcell <- apply(resp, 1, function(x) {sum(is.na(x))/length(x)})
  pdrug <- apply(resp, 2, function(x) {sum(is.na(x))/length(x)})
  
  if(max(pcell) >= max(pdrug)) {
    resp <- resp[-which.max(pcell), ]
  } else {
    resp <- resp[, -which.max(pdrug)]
  }
  nomiss <- !any(is.na(resp))
}
rm(pcell, pdrug, nomiss)

# cell lines in rows and genes in columns
expr.temp <- t(expr[, -1])
colnames(expr.temp) <- expr$ensembl_gene
rownames(expr.temp) <- substr(rownames(expr.temp), 2, 10000L)
expr <- expr.temp[order(as.numeric(rownames(expr.temp))), ]
rm(expr.temp)

# only retain the cell lines that are available in both response and expression
expr.prep <- expr[rownames(expr) %in% rownames(resp), ]
resp.prep <- resp[rownames(resp) %in% rownames(expr), ]

# keep drugs that are available in both response and drug data
drug.prep <- drug[drug$name %in% colnames(resp.prep), ]
drug.prep <- drug.prep[order(drug.prep$name), ]

# save(drug.prep, expr.prep, resp.prep, file="data/data_gdsc_dat1.Rdata")
load(file="data/data_gdsc_dat1.Rdata")

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
drug.ensemblid <- sapply(1:D, function(d) {
  vec <- mapid$ensemblid[mapid$reactomeid %in% drug.reactomeid[d]]
  if(length(vec)==0) {NA} else {vec}})

# determine for every feature whether it is in the target pathway or not
inpathway <- sapply(1:D, function(d) {
  as.numeric(colnames(expr.sel[[d]]) %in% drug.ensemblid[[d]])})

# selecting genes
# fit.enet <- apply(resp.scale, 2, function(y) {
#   glmnet(expr.scale, y, alpha=0.01)})
# s <- fit.enet[[1]]
# idsel <- lapply(fit.enet, function(s) {
#   which(as.numeric(coef(s, s=s$lambda[
#     tail(which(s$df < 300), n=1)]))[-1]!=0)})
# save(fit.enet, idsel, file="results/data_gdsc_res1.Rdata")
load(file="results/data_gdsc_res1.Rdata")
expr.sel <- sapply(idsel, function(id) {expr.prep[, id]})

# number of equations, features, and observations
d <- ncol(resp.prep)
p <- sapply(idsel, length)
n <- nrow(resp.prep)

# create co-data matrix
C <- lapply(1:ncol(resp.prep), function(s) {
  stage <- replace(drug.prep$stage, is.na(drug.prep$stage), "experimental")
  stage <- matrix(model.matrix(~ stage)[s, -1], nrow=p[s], ncol=2, byrow=TRUE)
  action <- matrix(model.matrix(~ replace(drug.prep$action, is.na(
    drug.prep$action), "unknown"))[s, -1], nrow=p[s], ncol=2, byrow=TRUE)
  mat <- cbind(stage, action, inpathway[[s]])
  colnames(mat) <- c("stage.experimental", "stage.in clinical development",
                     "action.targeted", "action.unknown", "inpathway")
  return(mat)})

# analysis
D <- ncol(resp.prep)
p <- sapply(expr.sel, ncol)
n <- nrow(resp.prep)
set.seed(2019)
ntrain <- floor(nrow(resp.prep)/2)
idtrain <- sample(1:nrow(resp.prep), ntrain)
xtrain <- lapply(expr.sel, function(s) {scale(s[idtrain, ])})
ytrain <- scale(resp.prep[idtrain, ])
xtest <- lapply(expr.sel, function(s) {scale(s[-idtrain, ])})
ytest <- scale(resp.prep[-idtrain, ])

control <- list(conv.post=TRUE, trace=TRUE,
                epsilon.eb=1e-3, epsilon.vb=1e-3,
                epsilon.opt=sqrt(.Machine$double.eps),
                maxit.eb=300, maxit.vb=2, maxit.opt=100,
                maxit.post=100)
fit1.enig <- enig(xtrain, ytrain, NULL, mult.lambda=TRUE, intercept.eb=TRUE,
                  fixed.eb="none", full.post=FALSE, init=NULL, control=control)
fit2.enig <- enig(xtrain, ytrain, C, mult.lambda=TRUE, intercept.eb=TRUE,
                  fixed.eb="none", full.post=FALSE, init=NULL, control=control)
save(fit1.enig, fit2.enig, file="results/data_gdsc_fit1.Rdata")
load(file="results/data_gdsc_fit1.Rdata")

# MCMC samples
# stan.enig <- stan_model("code/enig.stan")
# nsamples <- 1000
# nwarmup <- 1000
# mcmc1.enig <- sapply(1:D, function(d) {
#   sampling(stan.enig, chains=1, warmup=nwarmup, iter=nsamples + nwarmup, 
#            cores=1, refresh=0, control=list(adapt_delta=0.8, max_treedepth=12),
#            data=list(n=ntrain, p=p[d], y=ytrain[, d], x=xtrain[[d]], 
#                      ctalphainv=fit1.enig$eb$mprior[[d]], 
#                      lambda=fit1.enig$eb$lambda[d]))})
# 
# mcmc2.enig <- sapply(1:D, function(d) {
#   sampling(stan.enig, chains=1, warmup=nwarmup, iter=nsamples + nwarmup, 
#            cores=1, refresh=0, control=list(adapt_delta=0.8, max_treedepth=12),
#            data=list(n=ntrain, p=p[d], y=ytrain[, d], x=xtrain[[d]], 
#                      ctalphainv=fit2.enig$eb$mprior[[d]], 
#                      lambda=fit2.enig$eb$lambda[d]))})
# bmpost1.enig <- sapply(1:D, function(d) {
#   get_posterior_mean(mcmc1.enig[[d]])[1:(p[d] + 1), ]})
# bmpost2.enig <- sapply(1:D, function(d) {
#   get_posterior_mean(mcmc2.enig[[d]])[1:(p[d] + 1), ]})

mse1.enig <- sapply(1:D, function(d) {
  mean((ytest[, d] - xtest[[d]] %*% fit1.enig$vb$mpost$beta[[d]])^2)})
mse2.enig <- sapply(1:D, function(d) {
  mean((ytest[, d] - xtest[[d]] %*% fit2.enig$vb$mpost$beta[[d]])^2)})

mean(mse1.enig)
mean(mse2.enig)
plot(mse1.enig, mse2.enig)
abline(a=0, b=1, col=2, lty=2)

tab <- rbind(c(fit1.enig$eb$alpha, rep(NA, 5)),
             fit2.enig$eb$alpha)
colnames(tab) <- c("intercept", colnames(C[[1]]))
rownames(tab) <- c("intercept only", "model 1")
tab

plot(fit2.enig$seq.eb$alpha[, 1], ylim=range(fit2.enig$seq.eb$alpha), type="l", 
     col=1)
lines(fit2.enig$seq.eb$alpha[, 2], ylim=range(fit2.enig$seq.eb$alpha), col=2)
lines(fit2.enig$seq.eb$alpha[, 3], ylim=range(fit2.enig$seq.eb$alpha), col=3)
lines(fit2.enig$seq.eb$alpha[, 4], ylim=range(fit2.enig$seq.eb$alpha), col=4)
lines(fit2.enig$seq.eb$alpha[, 5], ylim=range(fit2.enig$seq.eb$alpha), col=5)

plot(fit1.enig$seq.eb$alpha, type="l")

barplot(unique(unlist(fit2.enig$eb$mprior)))
