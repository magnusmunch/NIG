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