# Cambridge project
## Manuscript
* compile docs/manuscript.rnw for manuscript (e.g., in R studio using knitr)
* compile docs/supplement.rnw for supplementary material (e.g., in R studio using knitr)
* citation style is in docs/author_short3.bst file
* bib references are in docs/refs.bib

## Code
* R code files are in code
* Gwen's original simple estimation loop is in code/script.R
* simple model simulations (without tissue effects and omics feature groups) are in code/simulations_igaussian.R. Re-run them with `Rscript --no-save --no-restore --verbose simulations_igaussian.R > simulations_igaussian.Rout 2>&1`
* data analysis scripts are in code/data_ccle.R and code/data_gdsc.R
* figure creation script is code/figures.R

## R package
* `R CMD build code` for R package compilation
* `R CMD check cambridge_0.0.0.9000.tar.gz` to check compiled package
* installation: 
  1. `install.packages("path/cambridge_0.0.0.9000.tar.gz", repos=NULL)` 
  2. `library(devtools); install_github("magnusmunch/cambridge/code", local=FALSE, auth_token="yourtoken")`. Best practice would be to (1) create a personal token on Github with Settings > Developer settings > Personal access tokens, (2) copy the token and (3) add the following line to your environment settings in .Rprofile: `Sys.setenv(GITHUB_PAT="yourtoken")`, where you exchange `"yourtoken"` with the copied token.
  
## Session info
```r
sessionInfo()
```