## Script to extract P38MAPK pathway and snps from Geuvadis data
rm(list=ls());gc()
library(R.utils)
library(refGenome)
library(snpStats)
library(rtracklayer)

## Define the working directory
path <- "/Users/gino/Documents/project3/data_plink"
setwd(path)

## Import annotation
gtf_dir <- "gencode.v12.annotation.gtf"
gtf <- import(gtf_dir)

## Keep annotation for protein coding genes only
index <- mcols(gtf)$gene_type == "protein_coding" & mcols(gtf)$type == "gene"
gtf2  <- gtf[index]
P38MAPK <- c("ATF1", "ATF2", "ATF6", "CEBPB", "CREB1", "CSNK2A1", "CSNK2A2", "CSNK2B", "DDIT3", "EIF4E",
             "EIF4EBP1", "ELK4", "ESR1", "GDI1", "HBP1", "HSPB1", "JUN", "GADD45G", "MAP3K1", "MAP3K10", "MAP3K4",
             "MAP3K5", "MAP3K6", "MAP3K7", "TAB2", "TAOK1", "TAOK2", "TAOK3", "MAPKAPK5", "MEF2A", "MEF2C", "MITF",
             "MKNK1", "NOS2","PLA2G4A", "PPARGC1A", "PTGS2", "RAB5A", "RPS6KA5", 
             "SLC9A1", "TP53", "USF1", "BLK", "CCM2", "CDC42", "DUSP1", "DUSP10", "DUSP16", "DUSP8", "FGR", "FYN", 
             "RAC1", "RALA", "RALB", "RIPK1","HCK", "LCK", "LYN", "MAP2K3", "MAP2K4", "MAP2K6", "MAP3K12",
             "MAP3K3", "PAK1", "PAK2", "TRAF2", "SRC", "TAB1", "YES1", "PRKG1", "RPS6KA4", "TRAF6", "EEF2K",
             "CCND1", "MAPK12", "MAPK13", "MAPT", "PKN1", "SNTA1", "STMN1", "AC013461.1", "CDC25B",
             "LSP1", "RAF1", "SFN", "SRF", "TCF3", "TSC2", "YWHAB", "YWHAE", "YWHAG", "YWHAH", "YWHAQ", "YWHAZ",
             "ATM", "CAMK2B", "GADD45A", "GADD45B", "KRT8", "MAPK11", "MAPK14", "MAPKAPK2", "MAPKAPK3", "TXN")

## "TXN" could not be found on any chromosome 1-22 thus I removed it (oops chr X)
## These genes "PRKG1", "RPS6KA4", "TRAF6", "EEF2K" have missing data, namely sample size of length one instead of 373

## Restrict gtf annotation to P38MAPK pathway genes
gtf_P38MAPK <- gtf2[mcols(gtf2)$gene_name %in% P38MAPK]
P38MAPK_id  <- substr(mcols(gtf_P38MAPK)$gene_id, start=1,stop=15)

## Load and order expression data
Y <- list()
Y$data <- utils::read.table(file=file.path(path,"GD462.GeneQuantRPKM.50FN.samplename.resk10.txt"),header=TRUE,nrows=-1)
Y$data <- Y$data[order(Y$data$Coord),]
Y$data <- Y$data[order(Y$data$Chr),]
Y$data <- Y$data[(Y$data$Chr!="X"),]


## Restrict genes to P38MAPK pathway genes and order gene positions with corresponding annotation
Y$data         <- Y$data[substr(as.character(Y$data$TargetID), start=1, stop=15) %in% P38MAPK_id,]
index2         <- match(substr(as.character(Y$data$TargetID), start=1, stop=15), P38MAPK_id)
P38MAPK_id     <- P38MAPK_id[index2]
P38MAPK        <- P38MAPK[index2]
df_gtf_P38MAPK <- mcols(gtf_P38MAPK)[index2,]
Y$map_gtf      <- df_gtf_P38MAPK
Y$map          <- Y$data[,c(1:4)]
Y$eData        <- t(as.matrix(Y$data[,-c(1:4)]))

## Obtain pathway genes positions and convert into dataframe
df_iranges  <- data.frame(iranges = gtf_P38MAPK@ranges)
df_iranges  <- df_iranges[index2,]

G   <- vector("list", length(index2))
Kb  <- 10^5
#################################################
for(i in 1:22){
  print(i)
  cat("\n")
  indices <- which(Y$map$Chr==i)
  if(length(indices)>=1){
    cat("Yes\n")
    chr <- i
    Z         <- list()
    Z$data    <- Y$data[indices,]
    Z$map     <- Z$data[,c(1:4)]
    Z$eData   <- t(as.matrix(Z$data[,-c(1:4)]))
    Z$map_gtf <- Y$map_gtf[indices,]
    iranges   <- df_iranges[indices,]
    ## Load genotype data
    bed <- file.path(path,paste0("snps",chr,".bed"))
    bim <- file.path(path,paste0("snps",chr,".bim"))
    fam <- file.path(path,paste0("snps",chr,".fam"))
    X   <- snpStats::read.plink(bed=bed,bim=bim,fam=fam)
    X$fam <- NULL
    names(X)[1] <- "data"
    all(diff(X$map$position) > 0)
    
    ##Exclude samples without expression and genotype data
    intersect  <- intersect(rownames(Z$eData),rownames(X$data))
    Z$eData    <- Z$eData[intersect,]
    X$data     <- X$data[intersect,]
    
    ## Exclude samples from population YRI
    samples  <- utils::read.delim(file=file.path(path,"E-GEUV-1.sdrf.txt"))
    if(length(indices)>1){Match <- match(rownames(Z$eData),samples$Source.Name)
    cond     <- samples$Factor.Value.population[Match]!="YRI"
    Z$eData  <- Z$eData[cond,]
    X$data   <- X$data[cond,]
    }else{Match <- match(names(Z$eData),samples$Source.Name)
    cond     <- samples$Factor.Value.population[Match]!="YRI"
    Z$eData  <- Z$eData[cond]
    X$data   <- X$data[cond,]
    }
    
    ## Exclude SNPs with $\text{maf}<5\%$
    maf    <- snpStats::col.summary(X$data)$MAF
    cond   <- maf >= 0.05
    X$data <- X$data[,cond]
    X$map  <- X$map[cond,]
    
    ## Exclude genes with zero-variance
    if(length(indices)>1){
    cond      <- apply(Z$eData,2,function(x) sd(x)>0)
    Z$eData   <- Z$eData[,cond]
    Z$map     <- Z$map[cond,]
    Z$map_gtf <- Z$map_gtf[cond,]}
    
###################
    ## Identify SNPs within two thousand base (2Kb) pairs around gene start
    match <- data.frame(from=integer(length(indices)),med=integer(length(indices)),to=integer(length(indices)))
    for(j in 1:length(indices)){
      ## conditions: same chromosome, within window
      cond1 <- X$map$position >= Z$map$Coord[j] - Kb # change!
      cond2 <- X$map$position <= Z$map$Coord[j] + Kb # change!
      boundInd <- as.integer(which(cond1 & cond2))
      if(!is.finite(max(which(X$map$position <= Z$map$Coord[j])))){next}
      media    <- max(which(X$map$position <= Z$map$Coord[j]))
        
      ## compression: first and last entry
      if(length(boundInd)==0){next}
      match$from[j] <- min(boundInd)
      match$med[j]  <- media
      match$to[j]   <- max(boundInd)
      if(length(boundInd)==1){next}
      if(!all(diff(boundInd)==1)){ stop("Wrong order!")}
      
      Xj            <- X$data[,seq(from=match$from[j],to=match$to[j],by=1)]
      Xj            <- abs(as(Xj,"numeric")-2) # major/minor allele
      Xj[Xj==2]     <- 1 ## to obtain dominant effect
      XjMap         <- X$map[seq(from=match$from[j],to=match$to[j],by=1),]
      snpsPos       <- XjMap$position
      if(length(indices)>1){Zj <- (Z$eData[,j] - mean(Z$eData[,j]))/sd(Z$eData[,j])}else{Zj <- (Z$eData[j] - mean(Z$eData[j]))/sd(Z$eData[j])}
      ZjMap         <- Z$map[j,]
      ZjMap_gtf     <- Z$map_gtf[j,]
      iranges_j     <- iranges[j,]
      Pj            <- vector(,length(snpsPos))
      for(s in 1:length(snpsPos)){
        if(snpsPos[s]>=iranges_j$iranges.start & snpsPos[s]<=iranges_j$iranges.end){Pj[s]=1}else{Pj[s]=0}
      }
      matchRanges   <- c(match$from[j],match$med[j],match$to[j])
      names(matchRanges) <- c("from","media","to")
      Gj            <- list(Zj,Xj,iranges_j,Pj,snpsPos,ZjMap,ZjMap_gtf,XjMap,matchRanges)
      G[indices[j]] <- list(Gj)
    }
  }
} 

cond <- logical(length(G))
for(r in 1:length(G)){
  cond[r] <- length(G[[r]][[1]])==373
}
Gn <- G[cond]


save(Gn, file = "/Users/gino/Documents/project3/data_plink/P38MAPKpathwayKb100K.RData")

