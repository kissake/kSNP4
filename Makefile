# Intention is to be able to build from the source distribution into binary
# distributions in one step.  Ideally including testing.

python = FTPgenomes number_SNPs_all3 ParAnn parse_assembly_summary

perl = SNPs_all_2_fasta_matrix3

all: kSNP3-source.zip kSNP3.zip Examples.zip 

# Note that packaging the source depends on the build and packaging of the binaries
# being successful.
kSNP3-source.zip: kSNP3.zip

kSNP3.zip: docs binaries kSNP3 test

Examples.zip: Example1 Example2

docs: THE\ BSD\ OPENSOURCE\ LICENSE.pdf kSNP3.1.2\ User\ Guide.pdf

# This section is for the building of the binaries for distribution.
# May need separate sections for Mac and Linux?  PC builds later?
binaries: $(python) $(perl) 

$(python): %: %.py

$(perl): %: %.pl


# Automated testing for the build
test: example1 example2

example1: binaries kSNP3

example2: binaries kSNP3
