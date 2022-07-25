# Intention is to be able to build from the source distribution into binary
# distributions in one step.  Ideally including testing.

#################################################################
#################################################################
# VARIABLES
#################################################################
#################################################################

# Current version:
ver = 3.1


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
binarydir = $(packagedir)/kSNP3

all_products = kSNP3.1_Source.zip kSNP3.zip Examples.zip

# All of the perl scripts.  Used to generate binaries in the perlmonolith case.
perl = add_paths3 annotate_SNPs_from_genbankFiles3 \
       CheckFileNames core_SNPs3 delete_allele_conflicts3 distance_tree3 \
       divide_input3 fasta_remove_new_lines3 FG2IF find_allele3 \
       find_unresolved_clusters3 force_binary_tree genome_names3 \
       get_genbank_file3 get_quantile3 Kchooser labelTree_AlleleCount-new3 \
       labelTree_HGT3 label_tree_nodes3 LE2Unix MakeFasta \
       MakeKSNP3infile merge_fasta_reads3 NodeChiSquare2tree3 \
       NodeChiSquare2tree31 parallel_commands3 \
       parse_mummer4kSNP3 parse_protein_annotation_counts3 parse_SNPs2VCF3 \
       pick_snps_from_kmer_genome_counts3 rc_kmer_freqs3 \
       rename_from_table3 renumber_probes3 rm_node_names_from_tree3 \
       SNP_matrix2dist_matrix3 SNPs2fastaQuery3 SNPs2nodes-new3 \
       split_by_fasta_entry3 subset_mer_list3 subset_mers3 \
       subset_SNPs_all3 tree_bioperl3 tree_nodes3 tree_plotter3

# All of the perl binaries that need to be created.
perlbin = binaries/add_paths3 binaries/annotate_SNPs_from_genbankFiles3 \
	  binaries/CheckFileNames binaries/core_SNPs3 \
	  binaries/delete_allele_conflicts3 binaries/distance_tree3 \
	  binaries/divide_input3 binaries/fasta_remove_new_lines3 \
	  binaries/FG2IF binaries/find_allele3 \
	  binaries/find_unresolved_clusters3 binaries/force_binary_tree \
	  binaries/genome_names3 binaries/get_genbank_file3 \
	  binaries/get_quantile3 binaries/Kchooser \
	  binaries/labelTree_AlleleCount-new3 binaries/labelTree_HGT3 \
	  binaries/label_tree_nodes3 binaries/LE2Unix binaries/MakeFasta \
	  binaries/MakeKSNP3infile binaries/merge_fasta_reads3 \
	  binaries/NodeChiSquare2tree3 binaries/NodeChiSquare2tree31 \
	  binaries/parallel_commands3 binaries/parse_mummer4kSNP3 \
	  binaries/parse_protein_annotation_counts3 binaries/parse_SNPs2VCF3 \
	  binaries/pick_snps_from_kmer_genome_counts3 binaries/rc_kmer_freqs3 \
	  binaries/rename_from_table3 binaries/renumber_probes3 \
	  binaries/rm_node_names_from_tree3 binaries/SNP_matrix2dist_matrix3 \
	  binaries/SNPs2fastaQuery3 binaries/SNPs2nodes-new3 \
	  binaries/split_by_fasta_entry3 binaries/subset_mer_list3 \
	  binaries/subset_mers3 binaries/subset_SNPs_all3 \
	  binaries/tree_bioperl3 binaries/tree_nodes3 binaries/tree_plotter3

# Ideally we can rename this, or rename all of the other perl scripts to match
# the .pl suffix for this one.  Consistency is the name of the game and
# permits us to treat them all the same.
perlpl = SNPs_all_2_fasta_matrix3.pl

# Names of the python binaries that need to be created.  The source files are
# determined from these filenames.
pythonbin = binaries/FTPgenomes binaries/number_SNPs_all3 binaries/ParAnn binaries/parse_assembly_summary

dependencies = FastTreeMP parsimonator mummer consense #jellyfish


#################################################################
# Default build target
#################################################################

all: $(all_products)

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
kSNP3.1_Source.zip: kSNP3.zip
	hg archive --exclude ".hg*" --prefix kSNP3.1_Source $@


# Build a different zip file depending on whether we are using a monolithic
# perl binary.
ifeq ($(perlmonolith),true)
$(packagedir): $(docs) kSNP3 $(perlbin) binaries/perlscripts $(pythonbin) $(dependencies)
	mkdir -p $(packagedir)
	mkdir -p $(binarydir)
	for doc in $(docs) ; do cp $$doc $(packagedir) ; done
	for bin in $(perlbin) $(pythonbin) binaries/perlscripts kSNP3 ; do cp $$bin $(binarydir) ; done
	for dep in $(dependencies) ; do cp $$dep $(binarydir) ; done
else
$(packagedir): $(docs) kSNP3 $(perlbin) $(pythonbin) $(dependencies)
	mkdir -p $(packagedir)
	mkdir -p $(binarydir)
	for doc in $(docs) ; do cp $$doc $(packagedir) ; done
	for bin in $(perlbin) $(pythonbin) kSNP3 ; do cp $$bin $(binarydir) ; done
	for dep in $(dependencies) ; do cp $$dep $(binarydir) ; done
endif

kSNP3.zip: $(packagedir)
	zip -y --symlinks -r $@ $(packagedir)

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
	pp -o binaries/perlscripts $(perl)
	for script in $(perl); do ln -s perlscripts binaries/$$script; done

else

# Build perl scripts into binaries using pp, a tool from CPAN in the PAR::Packer module
$(perlbin): binaries/%: %
	pp -o $@ $<

endif

#################################################################
# Python
#################################################################

# Build python binaries using pyinstaller.  Wrapped with a temporary directory
# creation and removal to avoid cluttering the working directory.  May be more
# efficient to create it and re-use it?
$(pythonbin):
	pybuilddir=`mktemp -d` ; \
		   echo "Build dir: $$pybuilddir" ; \
		   pyinstaller --workpath "$$pybuilddir" --specpath "$$pybuilddir" --clean --onefile --distpath binaries $(@F).py ; \
		   echo "Removing build directory: $$pybuilddir" ;\
		   rm -r "$$pybuilddir"



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
FastTreeMP: FastTree.c
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
consense: phylip-3.697.tar.gz
	tar -xzf $<
	# Patch consense to permit build and to permit longer (200 char in this case) names.
	patch -p 0 < consense.patch
	cd $(consense_build)/src && make -f Makefile.unx consense
	cp $(consense_build)/src/consense .

phylip-3.697.tar.gz:
	curl -L "http://evolution.gs.washington.edu/phylip/download/phylip-3.697.tar.gz" > $@


# Mummer
mummer_build = mummer-4.0.0rc1
mummer: mummer-src.tgz
	tar -xzf $<
	cd $(mummer_build) && autoreconf -fi && ./configure && make mummer
	cp $(mummer_build)/mummer .
	rm -r $(mummer_build)

mummer-src.tgz:
	curl -L "https://github.com/mummer4/mummer/archive/refs/tags/v4.0.0rc1.tar.gz" > $@

# Parsimonator
parsimonator_build=Parsimonator-1.0.2-master
parsimonator: parsimonator-src.zip
	unzip $<
	cd $(parsimonator_build)/ && make -f Makefile.gcc
	cp $(parsimonator_build)/parsimonator .
	rm -r $(parsimonator_build)

parsimonator-src.zip:
	curl -L "https://github.com/stamatak/Parsimonator-1.0.2/archive/refs/heads/master.zip" >$@


# Automated testing for the build

Examples.zip: Example1 Example2

test: example1 example2

example1: binaries kSNP3

example2: binaries kSNP3

clean:
	rm binaries/* || echo "Clean"
	rm -r $(packagedir) || echo "Clean"
	rm $(all_products) || echo "Clean"
	rm $(dependencies) || echo "Clean"

distclean: clean
	rm parsimonator-src.zip || echo "Clean"
	rm mummer-src.tgz || echo "Clean"
	rm phylip-3.697.tar.gz || echo "Clean"
	rm FastTree.c || echo "Clean"
	rm -r $(mummer_build) || echo "Clean"
	rm -r $(parsimonator_build) || echo "Clean"
	rm -r $(consense_build) || echo "Clean"
