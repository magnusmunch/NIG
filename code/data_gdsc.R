#!/usr/bin/env Rscript

### libraries
library(gdata)
library(biomaRt)

################################ reading data ##################################
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
  # same data, different file
  # expr <- read.csv("data/Cell_line_RMA_proc_basalExp.txt", 
  #                  header=TRUE, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  
} else {
  expr <- read.table("ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/sanger1018_brainarray_ensemblgene_rma.txt.gz", 
                     header=TRUE, stringsAsFactors=FALSE)
  expr <- read.table("https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/preprocessed/Cell_line_RMA_proc_basalExp.txt.zip",
                     header=TRUE, stringsAsFactors=FALSE, sep="\t", dec=",", 
                     colClasses="character", nrows=1)
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

# retrieve gene lists of all reactome pathways
if(file.exists("data/Ensembl2Reactome_All_Levels.txt")) {
  mapid <- read.table("data/Ensembl2Reactome_All_Levels.txt",
                      comment.char="", sep="\t", quote="",
                      col.names=c("ensemblid", "reactomeid", "url", "name", 
                                  "evidencecode", "species"),
                      header=TRUE, stringsAsFactors=FALSE)
} else {
  mapid <- read.table("https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt",
                      comment.char="", sep="\t", quote="",
                      col.names=c("ensemblid", "reactomeid", "url", "name",
                                  "evidencecode", "species"),
                      header=TRUE, stringsAsFactors=FALSE)
}

# mutation data
mut1 <- read.table(unz("data/CellLines_CG_BEMs.zip", 
                       "CellLines_CG_BEMs/PANCAN_SEQ_BEM.txt"), header=TRUE)

# methylation data
meth1 <- read.table(unz("data/METH_CELL_DATA.txt.zip", "F2_METH_CELL_DATA.txt"), 
                    header=TRUE)
meth2 <- read.xls("data/methSampleId_2_cosmicIds.xlsx", stringsAsFactors=FALSE)

################################ data cleaning #################################
# prepping mutation data
mut <- t(mut1[, -1])
rownames(mut) <- substr(rownames(mut), 2, nchar(rownames(mut)))
mut <- mut[order(as.numeric(rownames(mut))), ]
colnames(mut) <- mut1$CG

# prepping methylation data
meth <- t(meth1)
rownames(meth) <- meth2$cosmic_id[match(
  paste0("X", meth2$Sentrix_ID, "_", meth2$Sentrix_Position), colnames(meth1))]
meth <- meth[order(as.numeric(rownames(meth))), ]

# library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# annid <- sapply(1:nrow(meth1), function(s) {
#   which(ann450k$Islands_Name==rownames(meth1)[s])})
# geneid <- sapply(annid, function(s) {
#   unique(unlist(strsplit(ann450k$UCSC_RefGene_Name[s], ";")))})
# island <- sapply(annid, function(s) {
#   unique(unlist(strsplit(ann450k$Relation_to_Island[s], ";")))})
# 
# # View(listEnsembl())
# # View(listFilters(ensembl))
# # View(listAttributes(ensembl))
# ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# gb <- getBM(attributes=c("external_gene_name", "ensembl_gene_id", 
#                          "hgnc_symbol"), 
#             filters=c("external_gene_name"), 
#             values=unique(unlist(geneid)), 
#             mart=ensembl)

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

# prepped data
mut.prep <- mut
expr.prep <- expr
resp.prep <- as.matrix(resp)
meth.prep <- meth

# keep drugs that are available in both response and drug data
drug.prep <- drug[drug$name %in% colnames(resp.prep), ]
drug.prep <- drug.prep[order(drug.prep$name), ]

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

# link drugs to ensembl symbols in their pathway
drug.reactomeid <- sapply(drug.prep$pathways, function(s) {
  pathwayid$reactomeid[match(strsplit(s, ", ")[[1]], pathwayid$pathway)]})
drug.ensemblid <- sapply(1:nrow(drug.prep), function(d) {
  vec <- mapid$ensemblid[mapid$reactomeid %in% drug.reactomeid[d]]
  if(length(vec)==0) {NA} else {vec}})

# determine for every feature whether it is in the target pathway or not
inpathway <- lapply(1:nrow(drug.prep), function(d) {
  as.numeric(colnames(expr.prep) %in% drug.ensemblid[[d]])})

# create object that contains feature information
feat.prep <- list(inpathway=inpathway)
rm(drug.ensemblid, drug.reactomeid, notes, pathwayid, mapid, inpathway,
   drug, expr, resp, cell, meth, meth1, meth2, mut, mut1)

save(drug.prep, expr.prep, resp.prep, feat.prep, meth.prep, mut.prep,
     file="data/data_gdsc_dat1.Rdata")
