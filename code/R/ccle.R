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
library(devtools)
install_github("magnusmunch/cambridge/code", local=FALSE, 
               auth_token=Sys.getenv("GITHUB_PAT"))

### libraries
library(cambridge)

################################## analysis 1 ##################################
# loading the data either from a local file or from the GDSC website
# resp <- read.table("https://data.broadinstitute.org/ccle_legacy_data/pharmacological_profiling/CCLE_NP24.2009_Drug_data_2015.02.24.csv", 
#                    header=TRUE, sep=",", stringsAsFactors=FALSE)
resp <- read.table("../../data/CCLE_NP24.2009_Drug_data_2015.02.24.csv", 
                   header=TRUE, sep=",", stringsAsFactors=FALSE)

# create a separate tissue variable
resp$Tissue <- unlist(lapply(strsplit(resp$CCLE.Cell.Line.Name, "_"), 
                             function(l) {paste(l[-1], collapse=" ")}))

# select the sigmoid type fits
resp2 <- resp[resp$FitType=="Sigmoid", c(2, 3, 11)]

# switch to wide format
resp3 <- reshape(resp2, timevar="Compound", idvar="Primary.Cell.Line.Name",
                 direction="wide")
rownames(resp3) <- resp3$Primary.Cell.Line.Name
resp4 <- resp3[, -1]

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






