# Get the version info for later
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)

all: docs check clean

docs:
	Rscript -e 'library("roxygen2"); roxygenise(".")'

build: docs
	cd ..;\
	R CMD build coenocliner

check: build
	cd ..;\
	R CMD check coenocliner_$(PKGVERS).tar.gz

check-cran: build
	cd ..;\
	R CMD check --as-cran coencoliner_$(PKGVERS).tar.gz

install: build
	cd ..;\
	R CMD INSTALL coencoliner_$(PKGVERS).tar.gz

move:
	cp coenocliner.Rcheck/coenocliner-Ex.Rout ./tests/Examples/coenocliner-Ex.Rout.save

clean:
	cd ..;\
	rm -r coenocliner.Rcheck/
