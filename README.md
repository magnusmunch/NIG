# Cambridge project
## Manuscript
* compile manuscript.rnw for manuscript (e.g., in R studio using knitr)
* citation style is in author_short3.bst file
* bib style references are refs.bib

## Code
* R functions are in functions.R file
* cpp functions are in functions.cpp file
* Gwen's functions and simple estimation loop are in myfunction.cpp and script.R, respectively. Changes to script.R are in script_MM.R.

## R package
* R CMD build code for R package compilation
* R CMD check cambridge_0.0.0.9000.tar.gz to check compiled package
* installation: 
  1. install.packages("path/cambridge_0.0.0.9000.tar.gz", repos=NULL) 
  OR
  2. library(devtools)
  install_github("magnusmunch/cambridge/code", local=FALSE,
                 auth_token="da10f2b37c513e3383c5b2e0aa1300288329c636")