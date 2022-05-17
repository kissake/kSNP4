#!/usr/bin/python3

# Cannot find this command: number_SNPs_all3 - JN
# It exists in the binary distribution. - JN
# Additional observations about this executable:  Originally written in perl, compiled to
# binary using ActiveState's ActivePerl.  Takes a single commandline argument
# (will bail, possibly after creating a spurious output file if missing), and
# will create a file named <argument>_labelLoci as well as a COUNT_SNPs file.
# Output in testingwas the same as the content of the COUNT_SNPs file.
# all_SNPs sorted to all_SNPs_sorted_labelLoci doesn't add any lines, does
# seem to maintain sorting by the first column in the original file (kmers?), but
# does NOT maintain sorting by subsequent columns.  Appears to:
# separate distinct kmers with newlines, and apply a single sequential number
# to each distinct kmer, prepending that number on each line that kmer appears

# This version functions differently.  Specifically:
# It does not take any commandline arguments, but instead transforms its input to its output.
# It will function similarly to the original if the input is sorted, but not identically, as
#  the original will reorder lines within a given SNP section.
# It does not require the input to be sorted.
# It does require enough memory to contain all of the SNPs.  Another version that does not
#   require that is possible, assuming the input is sorted.

# Example original usage:
#   renumber_SNPs_all3 all_SNPs_sorted
# Revised usage:
#   renumber_SNPs_all3 < all_SNPs_sorted > all_SNPs_sorted_labelLoci

import sys as sys

snpFile = sys.stdin
seenSNPs = {}
lastSnp = ""
maxSnpNumber = 0
for line in snpFile.readlines():
    args = line.strip().split("\t")

    snp = args[0]

    snpNumber = seenSNPs.setdefault(snp,maxSnpNumber)
    
    # Ensure that we can use maxSnpNumber to initialize the next new SNP.
    if snpNumber == maxSnpNumber:
        maxSnpNumber = maxSnpNumber+1
        print("") # And output a divider.

    print("\t".join([str(snpNumber)]+args))
