################################### PREAMBLE ###################################

# set directories
CDIR := code
DDIR := docs
PKGDIR := rpackage

# set R package compile options
ROPTS := --no-save --no-restore --verbose

# retrieve package name and version
PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" $(PKGDIR)/DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" $(PKGDIR)/DESCRIPTION)

#################################### RULES #####################################

all: $(DDIR)/supplement.pdf $(DDIR)/manuscript.pdf build

# build the package if contents changed
#$(PKGNAME)_$(PKGVERS).tar.gz: $(shell find $(PKGDIR))
#	R CMD build $(PKGDIR)
build:
	R CMD build $(PKGDIR)

# check the package if the compiled package changed
#$(PKGNAME).Rcheck: $(PKGNAME)_$(PKGVERS).tar.gz
#	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz
check: build
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz

install: build
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

# run the simulations if the script or package changed
$(CDIR)/simulations_igaussian.Rout: $(CDIR)/simulations_igaussian.R
	cd $(CDIR);\
	Rscript $(ROPTS) simulations_igaussian.R > simulations_igaussian.Rout 2>&1

# run knitr if simulations output, or the knitr file changed
$(DDIR)/supplement.tex: $(CDIR)/simulations_igaussian.Rout $(DDIR)/supplement.Rnw
	Rscript -e "library(knitr); knit('$(DDIR)/supplement.Rnw', '$(DDIR)/supplement.tex')"

# compile the supplement if the tex file or refs file changed
$(DDIR)/supplement.pdf: $(DDIR)/supplement.tex $(DDIR)/refs.bib
	pdflatex -output-directory=$(DDIR) $(DDIR)/supplement.tex

# run knitr if the knitr file changed
# $(DDIR)/manuscript.tex: $(DDIR)/manuscript.Rnw
# 	Rscript -e "library(knitr); knit('$(DDIR)/manuscript.Rnw', '$(DDIR)/supplement.tex')"

# compile the manuscript if the tex file or refs file changed
# $(DDIR)/manuscript.pdf: $(DDIR)/manuscript.tex $(DDIR)/refs.bib
# 	pdflatex -output-directory=$(DDIR) $(DDIR)/manuscript.tex

# have to check these two:
$(DDIR)/%.tex: $(DDIR)/%.Rnw
	Rscript -e "library(knitr); knit('$<', '$(DDIR)/%.tex')"

$(DDIR)/%.pdf: $(DDIR)/%.tex $(DDIR)/refs.bib
	pdflatex -output-directory=$(DDIR) $<

