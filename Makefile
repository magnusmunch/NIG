PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" code/DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" code/DESCRIPTION)
RDIR := code/R
R_OPTS := --no-save --no-restore --verbose

# all: check clean

build: 
	R CMD build --no-manual code

check: build
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz

code/R/simulations_igaussian.Rout: $(PKGNAME)_$(PKGVERS).tar.gz code/R/simulations_igaussian.R
	Rscript $(R_OPTS) code/R/simulations_igaussian.R > code/R/simulations_igaussian.Rout 2>&1

docs/supplement.pdf: supplement.Rnw code/R/figures.R refs.bib author_short3.bst results/simulations_igaussian_fit1.csv results/simulations_igaussian_fit2.csv results/simulations_igaussian_res1.csv results/simulations_igaussian_res2.csv results/simulations_igaussian_set1.csv results/simulations_igaussian_set2.csv
	Rscript -e "library(knitr); knit('docs/supplement.Rnw')"
	pdflatex docs/supplement.tex

docs/manuscript.pdf: manuscript.Rnw code/R/figures.R refs.bib author_short3.bst results/simulations_igaussian_fit1.csv results/simulations_igaussian_fit2.csv results/simulations_igaussian_res1.csv results/simulations_igaussian_res2.csv results/simulations_igaussian_set1.csv results/simulations_igaussian_set2.csv
	Rscript -e "library(knitr); knit('docs/manuscript.Rnw')"
	pdflatex docs/manuscript.tex


# results/%.csv: code/R/simulations_igaussian.R 
# 	Rscript $(R_OPTS) simulations_igaussian.R > simulations_igaussian.Rout 2>&1
#
# docs/supplement.pdf: supplement.Rnw code/R/figures.R refs.bib author_short3.bst
# 	Rscript -e "library(knitr); knit('docs/supplement.Rnw')"
# 	pdflatex supplement.tex
#
# list R files
# RFILES := $(wildcard $(RDIR)/simulations_*.R)

# run simulations and data analysis scripts
# $(RDIR)/%.Rout: $(RDIR)/%.R $(RDIR)/functions.R
# 	Rscript $(R_OPTS) simulations_igaussian.R > simulations_igaussian.Rout 2>&1

