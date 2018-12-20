# Cambridge project
## Manuscript
* compile docs/manuscript.rnw for manuscript (e.g., in R studio using knitr)
* citation style is in docs/author_short3.bst file
* bib style references are docs/refs.bib

## Code
* R functions are in code/R/functions.R file
* cpp functions are in code/src/functions.cpp file
* Gwen's functions and simple estimation loop are in code/src/myfunction.cpp and code/R/script.R, respectively

## R package
* R CMD build code for R package compilation
* R CMD check cambridge_0.0.0.9000.tar.gz to check compiled package
* installation: 
  1. install.packages("cambridge_0.0.0.9000.tar.gz", repos=NULL) 
  2. library(devtools); install_github("magnusmunch/cambridge/code", local=FALSE, auth_token="yourtoken")