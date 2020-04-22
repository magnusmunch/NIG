[![Build Status](https://travis-ci.com/magnusmunch/cambridge.svg?token=5MrUgcz2TnCF4zpaKnRM&branch=master)](https://travis-ci.com/magnusmunch/cambridge)

## R package
* `R CMD build rpackage` for R package compilation
* `R CMD check NIG_0.0.0.9000.tar.gz` to check compiled package
* installation: 
  1. `install.packages("NIG_0.0.0.9000.tar.gz", repos=NULL)` 
  2. `library(devtools); install_github("magnusmunch/NIG/rpackage", local=FALSE, auth_token="yourtoken")`, where `"yourtoken"` is exchanged with a personal Github token, created with Settings > Developer settings > Personal access tokens. Alternatively, one (1) creates a personal token on Github as explained before, (2) copies the token and (3) adds the following line to the environment settings in .Rprofile: `Sys.setenv(GITHUB_PAT="yourtoken")`, where `"yourtoken"` is exchanged with the copied token. (4) Installation is then `library(devtools); install_github("magnusmunch/NIG/code", local=FALSE, auth_token=Sys.getenv("GITHUB_PAT"))`.