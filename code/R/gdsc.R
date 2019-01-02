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
# data <- read.table("ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/v17.3_fitted_dose_response.xlsx",
#                    header=TRUE, sep=",", stringsAsFactors=FALSE)
resp <- read.table("../../data/v17.3_fitted_dose_response.csv", header=TRUE, 
                   sep=",", stringsAsFactors=FALSE)


# take average if more than one measurement for drug cell line combination
resp2 <- aggregate(LN_IC50 ~ CELL_LINE_NAME + DRUG_NAME, data=resp, FUN=mean, 
                   na.rm=TRUE)

# put the data in the wide format (cell lines in rows and drugs in columns)
resp3 <- reshape(resp2, timevar="DRUG_NAME", idvar="CELL_LINE_NAME",
                 direction="wide")
rownames(resp3) <- resp3$CELL_LINE_NAME
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

str(resp5)

# 1.2) Preselecting GEX features: we selected gene expression features 
# explaining 50% variations over all cell lines, since this results in not too 
# features and probably remains those important features. 
# 
# 1.3) Preselecting CNV: we used the same strategy as Garnett et al. (2012) where the GDSC data are from.
# 
# 1.4) Preselecting MUT: we deleted mutated genes with less than two cell lines.



