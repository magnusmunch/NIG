## Project info
* All commands below assume that the working directory is the NIG project folder.
* Title: Drug sensitivity prediction with normal inverse Gaussian shrinkage 
informed by external data
* Authors: Magnus M. M&uuml;nch (responsible for writing code, email: m.munch@amsterdamumc.nl), Mark A. van de Wiel, Sylvia Richardson, and Gwena&euml;l G. R. Leday

* configuration:
```
─ Session info ──────────────────────────────────────────────────────────────────────────
 setting  value                       
 version  R version 3.6.3 (2020-02-29)
 os       macOS Catalina 10.15.3      
 system   x86_64, darwin15.6.0        
 ui       RStudio                     
 language (EN)                        
 collate  en_US.UTF-8                 
 ctype    en_US.UTF-8                 
 tz       Europe/Brussels             
 date     2020-04-28                  

─ Packages ──────────────────────────────────────────────────────────────────────────────
 package       * version    date       lib source        
 AnnotationDbi   1.48.0     2019-10-29 [1] Bioconductor  
 askpass         1.1        2019-01-13 [1] CRAN (R 3.6.0)
 assertthat      0.2.1      2019-03-21 [1] CRAN (R 3.6.0)
 backports       1.1.6      2020-04-05 [1] CRAN (R 3.6.2)
 Biobase         2.46.0     2019-10-29 [1] Bioconductor  
 BiocFileCache   1.10.2     2019-11-08 [1] Bioconductor  
 BiocGenerics    0.32.0     2019-10-29 [1] Bioconductor  
 biomaRt       * 2.42.1     2020-03-26 [1] Bioconductor  
 bit             1.1-15.2   2020-02-10 [1] CRAN (R 3.6.0)
 bit64           0.9-7      2017-05-08 [1] CRAN (R 3.6.0)
 blob            1.2.1      2020-01-20 [1] CRAN (R 3.6.0)
 callr           3.4.3      2020-03-28 [1] CRAN (R 3.6.2)
 cli             2.0.2      2020-02-28 [1] CRAN (R 3.6.0)
 codetools       0.2-16     2018-12-24 [1] CRAN (R 3.6.3)
 colorspace      1.4-1      2019-03-18 [1] CRAN (R 3.6.0)
 crayon          1.3.4      2017-09-16 [1] CRAN (R 3.6.0)
 curl            4.3        2019-12-02 [1] CRAN (R 3.6.0)
 DBI             1.1.0      2019-12-15 [1] CRAN (R 3.6.0)
 dbplyr          1.4.2      2019-06-17 [1] CRAN (R 3.6.0)
 desc            1.2.0      2018-05-01 [1] CRAN (R 3.6.0)
 devtools        2.3.0      2020-04-10 [1] CRAN (R 3.6.3)
 digest          0.6.25     2020-02-23 [1] CRAN (R 3.6.0)
 doParallel    * 1.0.15     2019-08-02 [1] CRAN (R 3.6.0)
 dplyr           0.8.5      2020-03-07 [1] CRAN (R 3.6.0)
 ellipsis        0.3.0      2019-09-20 [1] CRAN (R 3.6.0)
 fansi           0.4.1      2020-01-08 [1] CRAN (R 3.6.0)
 foreach       * 1.5.0      2020-03-30 [1] CRAN (R 3.6.2)
 fs              1.4.1      2020-04-04 [1] CRAN (R 3.6.2)
 gdata         * 2.18.0     2017-06-06 [1] CRAN (R 3.6.0)
 ggplot2       * 3.3.0      2020-03-05 [1] CRAN (R 3.6.0)
 glmnet        * 3.0-2      2019-12-11 [1] CRAN (R 3.6.0)
 glue            1.4.0      2020-04-03 [1] CRAN (R 3.6.2)
 gridExtra       2.3        2017-09-09 [1] CRAN (R 3.6.0)
 gsl             2.1-6      2019-03-25 [1] CRAN (R 3.6.3)
 gtable          0.3.0      2019-03-25 [1] CRAN (R 3.6.0)
 gtools          3.8.2      2020-03-31 [1] CRAN (R 3.6.2)
 hms             0.5.3      2020-01-08 [1] CRAN (R 3.6.0)
 httr            1.4.1      2019-08-05 [1] CRAN (R 3.6.0)
 inline          0.3.15     2018-05-18 [1] CRAN (R 3.6.0)
 IRanges         2.20.2     2020-01-13 [1] Bioconductor  
 iterators     * 1.0.12     2019-07-26 [1] CRAN (R 3.6.0)
 lattice         0.20-38    2018-11-04 [1] CRAN (R 3.6.3)
 lifecycle       0.2.0      2020-03-06 [1] CRAN (R 3.6.0)
 loo             2.2.0      2019-12-19 [1] CRAN (R 3.6.0)
 magrittr        1.5        2014-11-22 [1] CRAN (R 3.6.0)
 Matrix        * 1.2-18     2019-11-27 [1] CRAN (R 3.6.3)
 matrixStats     0.56.0     2020-03-13 [1] CRAN (R 3.6.0)
 memoise         1.1.0      2017-04-21 [1] CRAN (R 3.6.0)
 munsell         0.5.0      2018-06-12 [1] CRAN (R 3.6.0)
 NIG           * 0.0.0.9000 2020-04-21 [1] local         
 openssl         1.4.1      2019-07-18 [1] CRAN (R 3.6.0)
 pillar          1.4.3      2019-12-20 [1] CRAN (R 3.6.0)
 pkgbuild        1.0.6      2019-10-09 [1] CRAN (R 3.6.0)
 pkgconfig       2.0.3      2019-09-22 [1] CRAN (R 3.6.0)
 pkgload         1.0.2      2018-10-29 [1] CRAN (R 3.6.0)
 prettyunits     1.1.1      2020-01-24 [1] CRAN (R 3.6.0)
 processx        3.4.2      2020-02-09 [1] CRAN (R 3.6.0)
 progress        1.2.2      2019-05-16 [1] CRAN (R 3.6.0)
 ps              1.3.2      2020-02-13 [1] CRAN (R 3.6.0)
 purrr           0.3.3      2019-10-18 [1] CRAN (R 3.6.0)
 R6              2.4.1      2019-11-12 [1] CRAN (R 3.6.0)
 rappdirs        0.3.1      2016-03-28 [1] CRAN (R 3.6.0)
 Rcpp            1.0.4.6    2020-04-09 [1] CRAN (R 3.6.3)
 remotes         2.1.1      2020-02-15 [1] CRAN (R 3.6.0)
 rlang           0.4.5      2020-03-01 [1] CRAN (R 3.6.0)
 rprojroot       1.3-2      2018-01-03 [1] CRAN (R 3.6.0)
 RSQLite         2.2.0      2020-01-07 [1] CRAN (R 3.6.0)
 rstan         * 2.19.3     2020-02-11 [1] CRAN (R 3.6.0)
 rstudioapi      0.11       2020-02-07 [1] CRAN (R 3.6.0)
 S4Vectors       0.24.4     2020-04-09 [1] Bioconductor  
 scales          1.1.0      2019-11-18 [1] CRAN (R 3.6.0)
 sessioninfo     1.1.1      2018-11-05 [1] CRAN (R 3.6.0)
 shape           1.4.4      2018-02-07 [1] CRAN (R 3.6.0)
 StanHeaders   * 2.21.0-1   2020-01-19 [1] CRAN (R 3.6.0)
 statmod       * 1.4.34     2020-02-17 [1] CRAN (R 3.6.0)
 stringi         1.4.6      2020-02-17 [1] CRAN (R 3.6.0)
 stringr         1.4.0      2019-02-10 [1] CRAN (R 3.6.0)
 testthat        2.3.2      2020-03-02 [1] CRAN (R 3.6.0)
 tibble          3.0.0      2020-03-30 [1] CRAN (R 3.6.2)
 tidyselect      1.0.0      2020-01-27 [1] CRAN (R 3.6.0)
 usethis         1.6.0      2020-04-09 [1] CRAN (R 3.6.3)
 vctrs           0.2.4      2020-03-10 [1] CRAN (R 3.6.0)
 withr           2.1.2      2018-03-15 [1] CRAN (R 3.6.0)
 XML             3.99-0.3   2020-01-20 [1] CRAN (R 3.6.0)
 yaml            2.2.1      2020-02-01 [1] CRAN (R 3.6.0)
```

## packages
required R packages, available from CRAN: 
* knitr 
* kableExtra
* gdata
* sp
* GeneralizedHyperbolic

available from bioconductor:
* biomaRt

## install R package
* `R CMD build rpackage` for R package compilation
* `R CMD check NIG_0.0.0.9000.tar.gz` to check compiled package
* installation via one of options: 
  1. `install.packages("NIG_0.0.0.9000.tar.gz", repos=NULL)` 
  2. `library(devtools); install_github("magnusmunch/NIG/rpackage", local=FALSE, auth_token="yourtoken")`, where `"yourtoken"` is exchanged with a personal Github token, created with Settings > Developer settings > Personal access tokens. Alternatively, one (1) creates a personal token on Github as explained before, (2) copies the token and (3) adds the following line to the environment settings in .Rprofile: `Sys.setenv(GITHUB_PAT="yourtoken")`, where `"yourtoken"` is exchanged with the copied token. (4) Installation is then `library(devtools); install_github("magnusmunch/NIG/code", local=FALSE, auth_token=Sys.getenv("GITHUB_PAT"))`.

## Results
* R code files are in code folder
* data preparation scripts are in code/data_ccle.R and code/data_gdsc.R
* GDSC data based simulations are in code/simulations_gdsc.R
* data analysis scripts are in code/analysis_ccle.R and code/analysis_gdsc.R
* figure creation script is code/figures.R
* to generate the results, run the code in code/all_results.R. This may take
a few days. To speed up, consider lowering the number of cross validation 
folds `nfolds`, increasing the number of cores to compute with `ncores`, or
decreasing the number of simulation reps `nreps` at the top of the file.

## Manuscript
* compile docs/manuscript.rnw for manuscript (e.g., in R studio using knitr)
* compile docs/supplement.rnw for supplementary material (e.g., in R studio using knitr)
* citation style is in docs/author_short3.bst file
* bib references are in docs/refs.bib
  
## Data
* CCLE data retrieved from: https://portals.broadinstitute.org/ccle
* GDSC data retrieved from: https://www.cancerrxgene.org
