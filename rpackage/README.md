[![Build Status](https://travis-ci.com/magnusmunch/cambridge.svg?token=5MrUgcz2TnCF4zpaKnRM&branch=master)](https://travis-ci.com/magnusmunch/cambridge)

## R package
* `R CMD build code` for R package compilation
* `R CMD check cambridge_0.0.0.9000.tar.gz` to check compiled package
* installation: 
  1. `install.packages("path/cambridge_0.0.0.9000.tar.gz", repos=NULL)` 
  2. `library(devtools); install_github("magnusmunch/cambridge/code", local=FALSE, auth_token="yourtoken")`. Best practice would be to (1) create a personal token on Github with Settings > Developer settings > Personal access tokens, (2) copy the token and (3) add the following line to your environment settings in .Rprofile: `Sys.setenv(GITHUB_PAT="yourtoken")`, where you exchange `"yourtoken"` with the copied token.

## Simulations
* simple model simulations (without tissue effects and omics feature groups) are in R/simulations_igaussian.R. Re-run them with `Rscript --no-save --no-restore --verbose simulations_igaussian.R > simulations_igaussian.Rout 2>&1`

## Data