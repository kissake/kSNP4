# Intention is to be able to build from the source distribution into binary
# distributions in one step.  Ideally including testing.

#################################################################
#################################################################
# VARIABLES
#################################################################
#################################################################

# Current version:
ver = 4


# Set this to 'true' to generate a smaller binary package, but note that
# the individual scripts are no longer standalone; they reference another
# built file named 'perlscripts' that is required for them to function.
perlmonolith = false
# perlmonolith = true


#################################################################
# Variables to let 'make' function.
#################################################################

# The directory names to use in the binary package
packagedir = kSNP$(ver)_Linux_package
binarydir = $(packagedir)/kSNP$(ver)

all_products = kSNP$(ver)_Source.zip kSNP$(ver).zip Examples.zip

# All of the perl scripts.  Used to generate binaries in the perlmonolith case.
perl = add_paths3 annotate_SNPs_from_genbankFiles3 \
       CheckFileNames core_SNPs3 delete_allele_conflicts3 distance_tree3 \
       find_allele3 \
       find_unresolved_clusters3 force_binary_tree genome_names3 \
       get_genbank_file3 get_quantile3 labelTree_AlleleCount-new3 \
       label_tree_nodes3 LE2Unix \
       merge_fasta_reads3 NodeChiSquare2tree3 \
       parallel_commands3 \
       parse_mummer4kSNP3 parse_protein_annotation_counts3 parse_SNPs2VCF3 \
       pick_snps_from_kmer_genome_counts3 rc_kmer_freqs3 \
       rename_from_table3 renumber_probes3 \
       SNP_matrix2dist_matrix3 SNPs2fastaQuery3 SNPs2nodes-new3 \
       subset_mer_list3 subset_mers3 \
       subset_SNPs_all3 tree_nodes3

# All of the perl binaries that need to be created.
perlbin = binaries/add_paths3 binaries/CheckFileNames binaries/core_SNPs3 \
	binaries/distance_tree3 binaries/force_binary_tree binaries/genome_names3 \
	binaries/labelTree_AlleleCount-new3 binaries/label_tree_nodes3 binaries/LE2Unix \
	binaries/MakeKSNP4infile binaries/merge_fasta_reads3 \
	binaries/NodeChiSquare2tree3 binaries/NodeChiSquare2tree4 \
	binaries/parse_mummer4kSNP3 binaries/parse_SNPs2VCF3 binaries/rc_kmer_freqs3 \
	binaries/renumber_probes3 binaries/rmNodeNamesFromTree4 \
	binaries/SNPs2fastaQuery3 binaries/SNPs2nodes-new3 \
	binaries/SNPs_all_2_fasta_matrix3 binaries/subset_SNPs_all3 \
	binaries/tree_nodes3 binaries/find_unresolved_clusters3 \
	binaries/rename_from_table3

# Ideally we can rename this, or rename all of the other perl scripts to match
# the .pl suffix for this one.  Consistency is the name of the game and
# permits us to treat them all the same.
perlpl = SNPs_all_2_fasta_matrix3.pl

# Names of the python binaries that need to be created.  The source files are
# determined from these filenames.
pythonbin = binaries/FTPgenomes binaries/number_SNPs_all3 binaries/ParAnn binaries/parse_assembly_summary \
	binaries/find_snps binaries/get_filtered_kmers binaries/guessPartition \
	binaries/inline_frequency_check binaries/partitionKmers binaries/Kchooser4

# Define the pythonbin target as the compiled binary for ever .py file,
# stored in the binaries directory
python = $(wildcard *.py)
pythonbin = $(patsubst %.py,binaries/%,$(wildcard *.py))

# Likewise perlbin
perl = $(wildcard *.pl)
perlbin = $(patsubst %.pl,binaries/%,$(wildcard *.pl))


shellscripts = binaries/kSNP4 binaries/buildtree binaries/extractNthLocus4 binaries/selectNodeAnnotations4


dependencies = binaries/FastTreeMP binaries/parsimonator binaries/mummer binaries/consense binaries/jellyfish



#################################################################
# Which targets are "phony" and should be built even if the
# target exists?
#################################################################
.PHONY: clean distclean all testrun


#################################################################
# Default build target
#################################################################

all: $(all_products)



#################################################################
# Other build target
#################################################################


testrun: RunExamples.sh Examples all
	PATH=`pwd`/binaries:"${PATH}" ./RunExamples.sh



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
kSNP$(ver)_Source.zip: kSNP$(ver).zip
	hg archive --exclude ".hg*" --prefix kSNP$(ver)_Source $@


# Build a different zip file depending on whether we are using a monolithic
# perl binary.
ifeq ($(perlmonolith),true)
$(packagedir): $(docs) kSNP$(ver) $(perlbin) binaries/perlscripts $(pythonbin) $(shellscripts) $(dependencies)
	mkdir -p $(packagedir)
	mkdir -p $(binarydir)
	for doc in $(docs) ; do cp $$doc $(packagedir) ; done
	for bin in $(perlbin) $(pythonbin) binaries/perlscripts kSNP$(ver) ; do cp $$bin $(binarydir) ; done
	for dep in $(dependencies) ; do cp $$dep $(binarydir) ; done
else
$(packagedir): $(docs) kSNP$(ver) $(perlbin) $(pythonbin) $(shellscripts) $(dependencies)
	mkdir -p $(packagedir)
	mkdir -p $(binarydir)
	for doc in $(docs) ; do cp $$doc $(packagedir) ; done
	for bin in $(perlbin) $(pythonbin) kSNP$(ver) ; do cp $$bin $(binarydir) ; done
	for dep in $(dependencies) ; do cp $$dep $(binarydir) ; done
endif

kSNP$(ver).zip: $(packagedir)
	zip --symlinks -r $@ $(packagedir)

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

# Defined above, determines whether we build a single monolithic perl
# executable that contains all scripts, or multiple binaries.  Significant
# space savings to build one, at the cost of risking confusing the user.
ifeq ($(perlmonolith),true)	

# Builds all perl sripts into 'perlscripts' executable that determines
# which one to execute based on the name used when calling it (thus the
# symbolic links)
binaries/perlscripts $(perlbin) &: $(perl)
	echo "NERF"
	pp -o binaries/perlscripts $(perl)
	for script in $(perl); do ln -s perlscripts binaries/$$script; done

else

# Build perl scripts into binaries using pp, a tool from CPAN in the PAR::Packer module
binaries/%: %.pl
	pp -o $@ $<

endif

#################################################################
# Python
#################################################################

# Build python binaries using pyinstaller.  Wrapped with a temporary directory
# creation and removal to avoid cluttering the working directory.  May be more
# efficient to create it and re-use it?
binaries/%: %.py
	pybuilddir=`mktemp -d` ; \
		   echo "Build dir: $$pybuilddir" ; \
		   pyinstaller --workpath "$$pybuilddir" --specpath "$$pybuilddir" --clean --onefile --distpath binaries $< ; \
		   echo "Removing build directory: $$pybuilddir" ;\
		   rm -r "$$pybuilddir"


#################################################################
# Shell
#################################################################

$(shellscripts): $(@F)
	cp $(@F) $@



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


jellyfish_build=Jellyfish-2.3.0
binaries/jellyfish: jellyfish-src.tgz
	tar -xf $<
	cd $(jellyfish_build) && autoreconf -i && ./configure && make -j 4
	cp $(jellyfish_build)/bin/.libs/jellyfish $@


jellyfish-src.tgz:
	curl -L "https://github.com/gmarcais/Jellyfish/archive/refs/tags/v2.3.0.tar.gz" >$@


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


