# Normal inverse Gaussian regression
All commands below assume that the working directory is the NIG project folder.

## packages
required R packages, available from CRAN: 
* knitr 
* kableExtra
* gdata
* sp
* GeneralizedHyperbolic
available from bioconductor:
* biomaRt

## Code
* R code files are in code folder
* GDSC data based simulations are in code/simulations_gdsc.R
* data analysis scripts are in code/data_ccle.R and code/data_gdsc.R
* figure creation script is code/figures.R

## Manuscript
* compile docs/manuscript.rnw for manuscript (e.g., in R studio using knitr)
* compile docs/supplement.rnw for supplementary material (e.g., in R studio using knitr)
* citation style is in docs/author_short3.bst file
* bib references are in docs/refs.bib

## R package
* `R CMD build rpackage` for R package compilation
* `R CMD check NIG_0.0.0.9000.tar.gz` to check compiled package
* installation: 
  1. `install.packages("NIG_0.0.0.9000.tar.gz", repos=NULL)` 
  2. `library(devtools); install_github("magnusmunch/NIG/rpackage", local=FALSE, auth_token="yourtoken")`, where `"yourtoken"` is exchanged with a personal Github token, created with Settings > Developer settings > Personal access tokens. Alternatively, one (1) creates a personal token on Github as explained before, (2) copies the token and (3) adds the following line to the environment settings in .Rprofile: `Sys.setenv(GITHUB_PAT="yourtoken")`, where `"yourtoken"` is exchanged with the copied token. (4) Installation is then `library(devtools); install_github("magnusmunch/NIG/code", local=FALSE, auth_token=Sys.getenv("GITHUB_PAT"))`.
  
## Session info
```r
sessionInfo()
```
