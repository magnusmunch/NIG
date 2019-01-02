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
install_github("cancerrxgene/gdscIC50", build_vignettes=TRUE)

### libraries
library(cambridge)
library(gdscIC50)
vignette("gdscIC50")

test <- gdsc_example


path.data <- "/Users/magnusmunch/Documents/OneDrive/PhD/project_cambridge/data/"


data <- read.table(paste(path.data, "v17.3_fitted_dose_response.csv", sep=""),
                   header=TRUE, sep=",", stringsAsFactors=FALSE)
data2 <- aggregate(LN_IC50 ~ CELL_LINE_NAME + DRUG_NAME, 
                   data=data[data$DRUG_ID > 1000, ], FUN=mean, na.rm=TRUE)
data3 <- reshape(data2, timevar="DRUG_NAME", idvar="CELL_LINE_NAME", 
                 direction="wide")

hist(colMeans(is.na(data)), xlab="Proportion missing", 
     main="Proportion missing in features")
hist(rowMeans(is.na(data)), xlab="Proportion missing", 
     main="Proportion missing in observations")

data4 <- na.omit(data3[, apply(data3, 2, function(x) {sum(is.na(x))}) <= 150])

# 1.1) Deleting missing drug response data: we sequentially deleted missing drugs and missing cell lines from bigger proportions to small proportions, to be more possible to get as big nonmissing matrix as possible. But we still cannot guarantee that is the biggest nonmissing matrix.
# 
# 1.2) Preselecting GEX features: we selected gene expression features explaining 50% variations over all cell lines, since this results in not too features and probably remains those important features. 
# 
# 1.3) Preselecting CNV: we used the same strategy as Garnett et al. (2012) where the GDSC data are from.
# 
# 1.4) Preselecting MUT: we deleted mutated genes with less than two cell lines.