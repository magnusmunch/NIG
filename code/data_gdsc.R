#!/usr/bin/env Rscript

### libraries
library(biomaRt)
library(gdata)

### biomart
ensembl <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

################################ response data #################################
# loading the data
resp <- read.table("data/v17.3_fitted_dose_response.csv", header=TRUE,
                   sep=",", stringsAsFactors=FALSE)

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
resp.temp <- resp.temp[, -c(1, 2)]
resp.temp <- lapply(resp.temp, function(s) {
  names(s) <- rownames(resp.temp); s[!is.na(s)]})

resp.prep <- resp.temp
rm(resp, resp.temp)

#################################### drugs #####################################
# reading in the preformatted list of compounds
drug1 <- read.xls("data/Screened_Compounds.xlsx", stringsAsFactors=FALSE)

# reading in the curent list of compounds
drug2 <- read.table("data/Drug_listTue Jun 25 15_19_34 2019.csv",
                    header=TRUE, sep=",", quote='"', stringsAsFactors=FALSE)

# load data from "A landscape of pharmacogenomic interaction in cancer" analysis
drug3 <- read.xls("data/TableS1F.xlsx", stringsAsFactors=FALSE, skip=1)[, -9]
drug4 <- read.xls("data/TableS1G.xlsx", stringsAsFactors=FALSE, skip=1)[, 1:5]
drug4$Cluster.identifier <- sapply(1:nrow(drug4), function(s) {
  if(is.na(drug4$Cluster.identifier[s])) {
    drug4$Cluster.identifier[c(1:s)][tail(which(!is.na(drug4$Cluster.identifier[
      c(1:s)])), n=1)]
  } else {
    drug4$Cluster.identifier[s]
  }})

# retrieve gene lists of all reactome pathways
mapid <- read.table("data/Ensembl2Reactome_All_Levels.txt",
                    comment.char="", sep="\t", quote="",
                    col.names=c("ensemblid", "reactomeid", "url", "name", 
                                "evidencecode", "species"),
                    header=TRUE, stringsAsFactors=FALSE)

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
                   pubchem=pubchem, target=targets, pathway=pathways,
                   nsamples=drug.temp$Sample.Size, count=drug.temp$Count,
                   stage=drug.temp$Clinical.Stage, action=action,
                   cluster=drug.temp$Cluster.identifier,
                   silhouette=drug.temp$Silhouette.Width,
                   stringsAsFactors=FALSE)
drug <- drug[order(drug$id), ]
rownames(drug) <- 1:nrow(drug)

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
reactomeid <- sapply(drug$pathway, function(s) {
  pathwayid$reactomeid[match(strsplit(s, ", ")[[1]], pathwayid$pathway)]}, 
  USE.NAMES=FALSE)
ensemblid <- sapply(1:nrow(drug), function(d) {
  vec <- mapid$ensemblid[mapid$reactomeid %in% reactomeid[d]]
  if(length(vec)==0) {NA} else {vec}}, USE.NAMES=FALSE)

# cleaning targets
targets <- sapply(drug$target, function(s) {
  s <- gsub(",", ", ", s, fixed=TRUE)
  s <- gsub("  ", " ", s, fixed=TRUE)
  s <- gsub("ABL", "ABL1", s)
  s <- gsub("ABL11", "ABL1", s)
  s <- gsub("BCR-ABL1", "BCR, ABL1", s, fixed=TRUE)
  s <- gsub(" [", ", ", s, fixed=TRUE)
  s <- gsub("]", "", s, fixed=TRUE)
  s <- gsub("HDAC inhibitor Class I, IIa, IIb, IV", 
            paste(paste0("HDAC inhibitor Class ", c("I", "IIa", "IIb", "IV")), 
                  collapse=", "), s, fixed=TRUE)
  s <- gsub("Antimetabolite (DNA & RNA)", "RNA antimetabolite", s, fixed=TRUE)
  s <- gsub("Retinioic X receptor (RXR) agonist", 
            "Retinioic X receptor agonist", s, fixed=TRUE)
  s <- gsub("AAPK1 (AMPK) agonist", "AAPK1 agonist, AMPK agonist", s, 
            fixed=TRUE)
  s <- gsub("PI3K (Class 1)", "PI3K class 1, PI3K", s, fixed=TRUE)
  s <- gsub("PI3K (class 1)", "PI3K class 1", s, fixed=TRUE)
  s <- gsub(" and ", ", ", s, fixed=TRUE)
  s <- gsub("VEGFR3/FLT4", "VEGFR3, FLT4", s, fixed=TRUE)
  s <- gsub("PPARdelta, PPARgamma, unknown", "PPARdelta, PPARgamma", s,
            fixed=TRUE)
  s <- gsub(" (", ", ", s, fixed=TRUE)
  s <- gsub("(", ", ", s, fixed=TRUE)
  s <- gsub(")", "", s, fixed=TRUE)
  s <- paste(unique(strsplit(s, ", ")[[1]]), collapse=", ")
  return(s)}, USE.NAMES=FALSE)

# matching targets to ensembl ids
# View(listFilters(ensembl))
# View(listAttributes((ensembl)))
query1 <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"), 
                filters="external_gene_name", 
                values=unique(na.omit(unlist(strsplit(targets,", ")))), 
                mart=ensembl)  
query2 <- getBM(attributes=c("ensembl_gene_id", "external_synonym"), 
                filters="external_synonym", 
                values=unique(na.omit(unlist(strsplit(targets,", ")))), 
                mart=ensembl)  
query3 <- getBM(attributes=c("ensembl_gene_id", "entrezgene_accession"), 
                filters="entrezgene_accession", 
                values=unique(na.omit(unlist(strsplit(targets,", ")))), 
                mart=ensembl)  
query4 <- getBM(attributes=c("ensembl_gene_id", "wikigene_name"), 
                filters="wikigene_name", 
                values=unique(na.omit(unlist(strsplit(targets,", ")))), 
                mart=ensembl)  
colnames(query1) <- colnames(query2) <- colnames(query3) <- colnames(query4) <-
  c("ensembl_gene_id", "gene_name")
query <- rbind(query1, query2, query3, query4)
query <- query[!duplicated(query), ]
query <- query[!is.na(query$gene_name), ]
targets_ensembl_gene_id <- sapply(targets, function(s) {
  vec <- unique(unlist(sapply(strsplit(s, ", ")[[1]], function(m) {
    query$ensembl_gene_id[query$gene_name==m]}, USE.NAMES=FALSE)))
  vec <- paste(vec[length(vec)!=0], collapse=", ")
  ifelse(vec=="", NA, vec)}, USE.NAMES=FALSE)

drug$target <- targets
drug$pathway_reactome_id <- sapply(reactomeid, function(s) {
  ifelse(all(is.na(s)), NA, paste(s, collapse=", "))})
drug$pathway_ensembl_gene_id <- sapply(ensemblid, function(s) {
  ifelse(all(is.na(s)), NA, paste(s, collapse=", "))})
drug$target_ensembl_gene_id <- targets_ensembl_gene_id

drug.prep <- drug

rm(action, drug, drug.temp, drug1, drug2, drug3, drug4, ensemblid, mapid, name,
   notes, pathwayid, pathways, pubchem, query, query1, query2, query3, query4, 
   reactomeid, synonyms, targets, targets_ensembl_gene_id)

############################### gene expression ################################
expr1 <- read.table("data/sanger1018_brainarray_ensemblgene_rma.txt.gz",
                    header=FALSE, stringsAsFactors=FALSE, skip=1)
expr2 <- read.table("data/sanger1018_brainarray_ensemblgene_rma.txt.gz",
                    header=FALSE, stringsAsFactors=FALSE, nrows=1)

# Cancer Gene Census (https://cancer.sanger.ac.uk/census)
census <- read.table("data/Census_allFri Feb 21 14_49_10 2020.csv", header=TRUE,
                     stringsAsFactors=FALSE, sep=",", quote='"')

# cell lines in rows and genes in columns
expr.temp <- t(expr1[, -1])
colnames(expr.temp) <- expr1$V1
rownames(expr.temp) <- as.character(expr2)[-1]
expr.temp <- expr.temp[order(as.numeric(rownames(expr.temp))), ]
expr.temp <- expr.temp[, order(colnames(expr.temp))]

# average duplicate rows
duprows <- rownames(expr.temp)[which(duplicated(rownames(expr.temp)))]
mat <- t(sapply(duprows, function(s) {
  apply(expr.temp[rownames(expr.temp)==s, ], 2, mean)}))
expr.temp <- rbind(expr.temp[!(rownames(expr.temp) %in% duprows), ], mat)
expr.temp <- expr.temp[order(rownames(expr.temp)), ]
expr <- expr.temp

# determine for every feature whether it is in the target pathway or not
in_pathway <- lapply(1:nrow(drug.prep), function(d) {
  s <- as.numeric(colnames(expr) %in% strsplit(
    drug.prep$pathway_ensembl_gene_id[[d]], ", ")[[1]])
  names(s) <- colnames(expr)
  s})
names(in_pathway) <- drug.prep$name
is_target <- lapply(1:nrow(drug.prep), function(d) {
  s <- as.numeric(colnames(expr) %in% strsplit(
    drug.prep$target_ensembl_gene_id[[d]], ", ")[[1]])
  names(s) <- colnames(expr)
  s})
names(is_target) <- drug.prep$name

# cancer gene census cleaning
cancergenes <- strsplit(census$Synonyms, ",")
cancergenes <- sapply(cancergenes, function(s) {
  s <- unlist(strsplit(s[which(substring(s, 1, 4)=="ENSG")], ".", 
                       fixed=TRUE))[1];
  ifelse(is.null(s), NA, s)})

# View(listFilters(ensembl))
# View(listAttributes((ensembl)))
query <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"), 
               filters="external_gene_name", 
               values=na.omit(census$Gene.Symbol[is.na(cancergenes)]), 
               mart=ensembl)  
cancergenes[match(query$external_gene_name, census$Gene.Symbol)] <- 
  query$ensembl_gene_id
census$ensembl_gene_id <- cancergenes
is_cancergene <- as.numeric(colnames(expr) %in% cancergenes)
names(is_cancergene) <- colnames(expr)

# create object that contains feature information
expr.prep <- list(expr=expr, in_pathway=in_pathway, is_target=is_target, 
                  is_cancergene=is_cancergene)
rm(cancergenes, census, expr1, expr2, expr, expr.temp, in_pathway, 
   is_cancergene, is_target, query, duprows, mat)

################################# cell lines ###################################
cell <- read.xls("data/Cell_Lines_Details.xlsx",
                 stringsAsFactors=FALSE)
rm(cell)

################################## mutations ###################################
mut1 <- read.table(unz("data/CellLines_CG_BEMs.zip", 
                       "CellLines_CG_BEMs/PANCAN_SEQ_BEM.txt"), header=TRUE)

# prepping mutation data
mut <- t(mut1[, -1])
rownames(mut) <- substr(rownames(mut), 2, nchar(rownames(mut)))
mut <- mut[order(as.numeric(rownames(mut))), ]

# translating gene names to ensembl ids
query1 <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"), 
                filters="external_gene_name", 
                values=mut1$CG, 
                mart=ensembl) 
query2 <- getBM(attributes=c("ensembl_gene_id", "external_synonym"), 
                filters="external_synonym", 
                values=mut1$CG, 
                mart=ensembl)  
colnames(query1) <- colnames(query2) <- c("ensembl_gene_id", "gene_name")
query <- rbind(query1, query2)
query <- query[!duplicated(query), ]
query <- query[!is.na(query$gene_name), ]
ensemblids <- sapply(mut1$CG, function(s) {
  query$ensembl_gene_id[query$gene_name==s]})
ensemblids[sapply(ensemblids, length)!=1] <- 
  lapply(ensemblids[sapply(ensemblids, length)!=1], function(s) {
    if(any(s %in% colnames(expr.prep$expr))) {
      paste(s[s %in% colnames(expr.prep$expr)], collapse=",")
      } else {
        s[1]
      }})
colnames(mut) <- unlist(ensemblids)
mut.prep <- mut
rm(ensemblids, mut, mut1, query, query1, query2)

################################## methylation #################################
meth1 <- read.table(unz("data/METH_CELL_DATA.txt.zip", "F2_METH_CELL_DATA.txt"), 
                    header=TRUE)
meth2 <- read.xls("data/methSampleId_2_cosmicIds.xlsx", stringsAsFactors=FALSE)

# prepping methylation data
meth <- t(meth1)
rownames(meth) <- meth2$cosmic_id[match(
  paste0("X", meth2$Sentrix_ID, "_", meth2$Sentrix_Position), colnames(meth1))]
meth <- meth[order(as.numeric(rownames(meth))), ]
meth.prep <- meth
rm(meth, meth1, meth2)

################################### save data ##################################
drug <- drug.prep
expr <- expr.prep
meth <- meth.prep
mut <- mut.prep
resp <- resp.prep
rm(ensembl, drug.prep, expr.prep, meth.prep, mut.prep, resp.prep)
save(drug, expr, resp, meth, mut, file="data/data_gdsc_dat1.Rdata")
