# Intention is to be able to build from the source distribution into binary
# distributions in one step.  Ideally including testing.

#################################################################
#################################################################
# VARIABLES
#################################################################
#################################################################

# Current version:
ver = 4


#################################################################
# Variables to let 'make' function.
#################################################################

# The directory names to use in the binary package
packagedir = kSNP$(ver)_Linux_package
binarydir = $(packagedir)/kSNP$(ver)

all_products = kSNP$(ver)_Source.zip kSNP$(ver).zip Examples.zip


# Define the pythonbin target as the compiled binary for ever .py file,
# stored in the binaries directory
python = $(wildcard *.py)
pythonbin = $(patsubst %.py,binaries/%,$(wildcard *.py))
pythonlib = ksnpConfig.py ksnpCache.py


# Likewise perlbin
perl = $(wildcard *.pl)
perlbin = $(patsubst %.pl,binaries/%,$(wildcard *.pl))


shellscripts = binaries/kSNP4 binaries/buildtree binaries/extractNthLocus4 binaries/selectNodeAnnotations4 binaries/installkSNP


dependencies = binaries/FastTreeMP binaries/parsimonator binaries/mummer binaries/consense binaries/jellyfish


installer = binaries/installkSNP



#################################################################
# Which targets are "phony" and should be built even if the
# target exists?
#################################################################
.PHONY: clean distclean all testrun testclean


#################################################################
# Default build target
#################################################################

all: $(all_products)



#################################################################
# Other build target
#################################################################


testrun: RunExamples.sh Examples all
	PATH=`pwd`/binaries:"${PATH}" ./RunExamples.sh

testclean:
	rm -r tmp.*



#################################################################
#################################################################
# BUILD PACKAGES FOR DISTRIBUTION
#################################################################
#################################################################



# The documents that are to be included in the binary and source distributions.
docs: binaries/THE\ BSD\ OPENSOURCE\ LICENSE.pdf binaries/kSNP3.021\ User\ Guide\ .pdf


# Note that packaging the source depends on the build and packaging of the binaries
# being successful.
# I think this defaults to the revision currently checked out, minus local
# revisions, and excludes the .hg directory to avoid sharing extra detail
kSNP$(ver)_Source.zip: kSNP$(ver)_Source
	zip -r $@ $<


kSNP$(ver)_Source: kSNP$(ver).zip
	hg archive --exclude ".hg*" $@  || echo "No mercurial available, switching to git"
	git clone . $@ && rm -fr $@/.git || echo "No git available.  Did mercurial work?"


$(packagedir): $(docs) kSNP$(ver) $(perlbin) $(pythonbin) $(shellscripts) $(dependencies)
	mkdir -p $(packagedir)
	mkdir -p $(binarydir)
	for doc in $(docs) ; do cp $$doc $(packagedir) ; done
	for bin in $(perlbin) $(pythonbin) kSNP$(ver) $(shellscripts) ; do cp $$bin $(binarydir) ; done
	for dep in $(dependencies) ; do cp $$dep $(binarydir) ; done
	cp $(installer) $(packagedir)

kSNP$(ver).zip: $(packagedir)
	cd $(packagedir) && zip --symlinks -r ../$@ kSNP$(ver) installkSNP

#################################################################
#################################################################
# BUILD BINARIES / DISTRIBUTION-READY DOCS
#################################################################
#################################################################

# This section is for the building of the binaries for distribution.
# May need separate sections for Mac and Linux?  PC builds later?


#################################################################
# Perl
#################################################################


# Build perl scripts into binaries using pp, a tool from CPAN in the PAR::Packer module
binaries/%: %.pl
	pp -o $@ $<

#################################################################
# Python
#################################################################

# Exception for ksnpConfig.py, which we don't need to build yet.
binaries/ksnpConfig: ksnpConfig.py
	echo "Skipping for now."

# Build python binaries using pyinstaller.  Wrapped with a temporary directory
# creation and removal to avoid cluttering the working directory.  May be more
# efficient to create it and re-use it?
binaries/%: %.py $(pythonlib)
	pybuilddir=`mktemp -d` ; \
		   echo "Build dir: $$pybuilddir" ; \
		   pyinstaller --workpath "$$pybuilddir" --specpath "$$pybuilddir" --clean --onefile --distpath binaries $< ; \
		   echo "Removing build directory: $$pybuilddir" ;\
		   rm -r "$$pybuilddir"


#################################################################
# Shell
#################################################################

$(shellscripts): binaries/%: %
	cp $< $@



#################################################################
# Docs
#################################################################

# How to prepare the PDF files to be packaged.  Note that we can generate PDF
# from other source if desirable.
binaries/%.pdf : %.pdf
	cp "$<" "$@"



#################################################################
# DEPENDENCIES
#################################################################

#
# NOTENOTENOTENOTENOTENOTENOTENOTENOTENOTENOTENOTE
# Files from this section need to be generally cleaned up
# to avoid cluttering the source directory.
# NOTENOTENOTENOTENOTENOTENOTENOTENOTENOTENOTENOTE
#


deps: $(dependencies)

# FastTreeMP
# See: http://www.microbesonline.org/fasttree/
binaries/FastTreeMP: FastTree.c
	gcc -DOPENMP -fopenmp -O3 -finline-functions -funroll-loops -Wall -o $@ $< -lm

FastTree.c:
	curl -L "http://www.microbesonline.org/fasttree/FastTree.c" > $@

# Jellyfish
# This is available as a package from Debian with the version available on 
# GitHub here: https://github.com/gmarcais/Jellyfish/releases/tag/v2.3.0
# It is also available in binary form from GitHub, as well as source.  The
# source doesn't build for me, so punting on rebuilding something that we can
# get more easily two other ways.

# Consense
# More info available here: https://evolution.genetics.washington.edu/phylip.html
# See testing here: https://evolution.gs.washington.edu/phylip/doc/consense.html
consense_build = phylip-3.697
binaries/consense: phylip-3.697.tar.gz
	tar -xzf $<
	# Patch consense to permit build and to permit longer (200 char in this case) names.
	patch -p 0 < consense.patch
	cd $(consense_build)/src && make -f Makefile.unx consense
	cp $(consense_build)/src/consense binaries/

phylip-3.697.tar.gz:
	curl -L "http://evolution.gs.washington.edu/phylip/download/phylip-3.697.tar.gz" > $@


# Mummer
mummer_build = mummer-4.0.0rc1
binaries/mummer: mummer-src.tgz
	tar -xzf $<
	cd $(mummer_build) && autoreconf -fi && ./configure && make mummer
	cp $(mummer_build)/.libs/mummer binaries/mummer-bin
	cp $(mummer_build)/.libs/libumdmummer.so.0.0.0 binaries/libumdmummer.so.0
	echo '#!/bin/bash\nLD_LIBRARY_PATH=`dirname "$${0}"`:$${LD_LIBRARY_PATH} `dirname "$${0}"`/mummer-bin "$${@}"' > "$@"
	chmod a+x "$@"
	# This is a release candidate, guessing that's why the build doesn't make the standalone binary?

mummer-src.tgz:
	curl -L "https://github.com/mummer4/mummer/archive/refs/tags/v4.0.0rc1.tar.gz" > $@

# Parsimonator
parsimonator_build=Parsimonator-1.0.2-master
binaries/parsimonator: parsimonator-src.zip
	unzip $<
	cd $(parsimonator_build)/ && make -f Makefile.gcc
	cp $(parsimonator_build)/parsimonator binaries/
	rm -r $(parsimonator_build)

parsimonator-src.zip:
	curl -L "https://github.com/stamatak/Parsimonator-1.0.2/archive/refs/heads/master.zip" >$@


jellyfish_build=jellyfish-2.2.6
binaries/jellyfish: jellyfish-src.tgz
	tar -xf $<
	cd $(jellyfish_build) && autoreconf -i && ./configure && make -j 4
	cp $(jellyfish_build)/bin/.libs/jellyfish $@


jellyfish-src.tgz:
	curl -L "https://github.com/gmarcais/Jellyfish/releases/download/v2.2.6/jellyfish-2.2.6.tar.gz" >$@


Examples: Examples.zip
	unzip "$<"
	rm -r __MACOSX/

Examples.zip:
	curl -L "https://downloads.sourceforge.net/project/ksnp/Examples.zip" >$@


clean:
	rm binaries/* || true
	rm -r $(packagedir) || true
	rm $(all_products) || true
	rm $(dependencies) || true

distclean: clean
	rm parsimonator-src.zip mummer-src.tgz mummer-bin libumdmummer.so.0 mummer phylip-3.697.tar.gz jellyfish-src.tgz FastTree.c Examples.zip  || true
	rm -r  $(mummer_build) $(parsimonator_build) $(consense_build) $(jellyfish_build) __MACOSX Examples __pycache__ || true


