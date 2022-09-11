# Intention is to be able to build from the source distribution into binary
# distributions in one step.  Ideally including testing.

#################################################################
#################################################################
# VARIABLES
#################################################################
#################################################################

# Current version:
ver = 4

# Whether to build dependencies from source or not
deps_from = download
# deps_from = source

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

download_cache = ~/kSNP/Build



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

# This is a process to acquire all of the below dependencies into the current
# directory through another means.  Note that this breaks dependency checking,
# which means that if the 'download' value is set for $(deps_from), the
# binaries will be rebuilt _every_time_ because deps_download does not exist.
deps_download: $(download_cache)/kSNP3.1_Linux_package.zip $(download_cache)/jellyfish $(download_cache)/FastTreeMP
	cp $(download_cache)/jellyfish .
	cp $(download_cache)/FastTreeMP .
	unzip -j -o $(download_cache)/kSNP3.1_Linux_package.zip kSNP3.1_Linux_package/kSNP3/mummer kSNP3.1_Linux_package/kSNP3/consense \
		kSNP3.1_Linux_package/kSNP3/ kSNP3.1_Linux_package/kSNP3/parsimonator

# This triggers rules below.
deps_source: consense mummer jellyfish FastTreeMP parsimonator
	echo "Built dependent binaries from source"

deps: $(dependencies)
	echo "Dependencies are in place."

$(dependencies): deps_$(deps_from)
	cp `basename $@` binaries/
#
# Some (all?) of these dependencies can come from the distro the user is using.
#
# See:  apt-get install phylip jellyfish fasttree mummer
#
# Note that consense is not available.
#


# FastTreeMP
# See: http://www.microbesonline.org/fasttree/
FastTreeMP: $(download_cache)/FastTree.c
	gcc -DOPENMP -fopenmp -O3 -finline-functions -funroll-loops -Wall -o $@ $< -lm

$(download_cache)/FastTree.c:
	mkdir -p $(download_cache)
	curl -L "http://www.microbesonline.org/fasttree/FastTree.c" > $@

$(download_cache)/FastTreeMP:
	mkdir -p $(download_cache)
	curl -L "http://www.microbesonline.org/fasttree/FastTreeMP" > $@

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
consense: $(download_cache)/phylip-3.697.tar.gz
	tar -xzf $<
	# Patch consense to permit build and to permit longer (200 char in this case) names.
	patch -p 0 < consense.patch
	cd $(consense_build)/src && make -f Makefile.unx consense
	cp $(consense_build)/src/consense binaries/

$(download_cache)/phylip-3.697.tar.gz:
	mkdir -p $(download_cache)
	curl -L "http://evolution.gs.washington.edu/phylip/download/phylip-3.697.tar.gz" > $@


# Mummer
mummer_build = mummer-4.0.0rc1
mummer: $(download_cache)/mummer-src.tgz
	tar -xzf $<
	cd $(mummer_build) && autoreconf -fi && ./configure && make mummer
	cp $(mummer_build)/.libs/mummer binaries/mummer-bin
	cp $(mummer_build)/.libs/libumdmummer.so.0.0.0 binaries/libumdmummer.so.0
	echo '#!/bin/bash\nLD_LIBRARY_PATH=`dirname "$${0}"`:$${LD_LIBRARY_PATH} `dirname "$${0}"`/mummer-bin "$${@}"' > "$@"
	chmod a+x "$@"
	# This is a release candidate, guessing that's why the build doesn't make the standalone binary?

$(download_cache)/mummer-src.tgz:
	mkdir -p $(download_cache)
	curl -L "https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz" > $@

# Parsimonator
parsimonator_build=Parsimonator-1.0.2-master
parsimonator: $(download_cache)/parsimonator-src.zip
	unzip $<
	cd $(parsimonator_build)/ && make -f Makefile.gcc
	cp $(parsimonator_build)/parsimonator binaries/
	rm -r $(parsimonator_build)

$(download_cache)/parsimonator-src.zip:
	mkdir -p $(download_cache)
	curl -L "https://github.com/stamatak/Parsimonator-1.0.2/archive/refs/heads/master.zip" >$@


jellyfish_build=Jellyfish-2.3.0
jellyfish: $(download_cache)/jellyfish-src.tgz
	tar -xf $<
	cd $(jellyfish_build) && autoreconf -i && ./configure && make -j 4
	cp $(jellyfish_build)/bin/.libs/jellyfish $@


$(download_cache)/jellyfish-src.tgz:
	mkdir -p $(download_cache)
	curl -L "https://github.com/gmarcais/Jellyfish/releases/download/v2.3.0/jellyfish-2.3.0.tar.gz" >$@

$(download_cache)/jellyfish:
	mkdir -p $(download_cache)
	curl -L "https://github.com/gmarcais/Jellyfish/releases/download/v2.3.0/jellyfish-linux" >$@


Examples: $(download_cache)/Examples.zip
	unzip "$<"
	rm -r __MACOSX/

Examples.zip: $(download_cache)/Examples.zip
	cp $< $@

$(download_cache)/Examples.zip:
	mkdir -p $(download_cache)
	curl -L "https://downloads.sourceforge.net/project/ksnp/Examples.zip" >$@

$(download_cache)/kSNP3.1_Linux_package.zip:
	mkdir -p $(download_cache)
	curl -L "https://downloads.sourceforge.net/project/ksnp/kSNP3.1_Linux_package.zip" >$@






#################################################################
# CLEANLINESS
#################################################################




clean:
	rm binaries/* || true
	rm -r $(packagedir) || true
	rm $(all_products) || true
	rm $(dependencies) || true

distclean: clean
	rm parsimonator-src.zip mummer-src.tgz mummer-bin parsimonator consense jellyfish libumdmummer.so.0 mummer phylip-3.697.tar.gz jellyfish-src.tgz FastTree.c Examples.zip || true
	rm -r  $(mummer_build) $(parsimonator_build) $(consense_build) $(jellyfish_build) __MACOSX Examples __pycache__ || true


