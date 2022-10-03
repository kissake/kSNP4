#!/usr/bin/python3

# We will be run with a file that should fit entirely in memory, and which should not
# break across a SNP boundary.  We don't need to know more than that.
# We will take filtered individual kmers, and output SNPs
#
# Performance note:  We use about 350 bytes/kmer discovered, and from a dataset with
# about 3.9million kmers as a median of five genomes, there were 8.7million total
# kmers found (in this dataset).  Need to assume that the overlap is worse than that
# (however, the worst case, completely divergent genomes, doesn't give us any data,
# so how important is it to plan for that case?)
#
# Revised estimate is closer to 115bytes or even 33bytes per kmer.  Keep testing.
# May be worth trying to crunch memory usage.  Should be doable.

# Import standard Python libraries for argument parsing, interacting with the
# filesystem, and environment and time / date processing.
import argparse as argparse

# Permit use of stdout
import sys

# Help with debug data:
import logging as logging


# Import tools for parsing kSNP configuration data.
# TODO - JN
# import ksnpConfig


# The override option to this function can be used to cause the parsing of a provided list
# instead of the programs passed options, for testing.
def parseCommandline(override=None):

    description = '''Identify SNPs from filtered kmers, and prepare files for input into mummer to 
determine SNP location.  Will work on partitioned data as long as the partition 
kmer is the same and isn't on a loci boundary.'''
    
    example = '''The term "partitioned" refers to creating files that have only a subset of the 
kmers (e.g. those beginning with AA).  All files should contain the same subset 
for proper operation.

To output the SNPs for the files fsplit0.kmer.part0, fsplit1.kmer.part0, and 
fsplit2.kmer.part0:

%(prog)s SNPs.part0 fsplit0.kmer.part0 fsplit1.kmer.part0 fsplit2.kmer.part0

The outputs will be:
  fsplit0.kmer.part0.SNPs.fasta 
  fsplit1.kmer.part0.SNPs.fasta
  fsplit2.kmer.part0.SNPs.fasta
and also (with the overall list of SNPs):
  SNPs.part0
'''
    
    parser = argparse.ArgumentParser(description= description, epilog=example,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    # parser.add_argument('SNPsFile', metavar='SNPsFile', nargs=1,
    #                     help='The filename to use for SNPs discovered in this input set')
    
    parser.add_argument('filteredGenomes', metavar='filteredGenome', nargs='+',
                        help='Filtered genomes (possibly partitioned) from which to extract SNPs.')

    parser.add_argument('--fasta-suffix', action='store', help='The suffix to use for fasta files to output to mummer for position data', default='SNPs.fasta')
    parser.add_argument('--snp-suffix', action='store', help='The suffix to use for SNP files', default='SNPs')
    parser.add_argument('--snps-all', action='store', help='A set of SNPs output by a previous run, for comparing to this run.')
    parser.add_argument('--debug', action='store_true', help='Output diagnostic data to the screen using STDERR')
    parser.add_argument('--info', action='store_true', help='Output status information on STDERR', default=True)
    parser.add_argument('--quiet', action='store_true', help='Silence debug and other messages. (warnings only)')

    if override is None:
        return parser.parse_args()
    else:
        return parser.parse_args(override)





    

def loadSNPsFromSNPsAll(bucket, snpsAllFile, startGenomeBit=1):
    '''
    This function loads a set of loci from an existing SNPs_all file.  The input
    file represents the sum total of the loci that can involve our previous data.

    If we only have one new genome, this will be all of the loci ever.  Otherwise,
    this will ensure that all of these loci are represented, including where they
    appear in the new genomes, in addition to all SNPs between new genomes.
    '''
    # read buffer (trying this out)
    bufferSize = 32 * 1024

    tab = '\t'

    # Column identifiers
    center = 2
    loci = 1
    genome = 4

    addedSNPsCount = 0
    lociAddedCount = 0
    lastSnpLocus = ''

    maxGenomeBit = startGenomeBit
    # Mapping of genomes to the associated bit in the bit field.
    genomeBits = {}
    
    # Default / initial value for a given locus.
    noSNPsNoGenomes = {
        'A': 0,
        'C': 0,
        'G': 0,
        'T': 0,
        }
    
    snps = open(snpsAllFile, 'r')

    nextLines = snps.readlines(bufferSize)
    while nextLines != []:
        for kmer in nextLines:
            snpData = kmer.split(tab)
            snpCenter = snpData[center]
            snpLocus = snpData[loci]
            snpGenome = snpData[genome]

            thisBit = genomeBits.setdefault(snpGenome, maxGenomeBit)
            if thisBit == maxGenomeBit:
                maxGenomeBit = maxGenomeBit * 2
            
            # Set the bit corresponding to this central nucleotide for this genome
            bucket[snpLocus][snpCenter] = bucket.setdefault(snpLocus, noSNPsNoGenomes.copy())[snpCenter] | thisBit

            # Overhead:
            if snpLocus != lastSnpLocus:
                lociAddedCount = lociAddedCount + 1
                # logging.debug("Added SNP: %s: %s", lastSnpLocus, bucket[lastSnpLocus])
            lastSnpLocus = snpLocus
            addedSNPsCount = addedSNPsCount + 1
        nextLines = snps.readlines(bufferSize)

    # logging.debug("Added SNP: %s: %s", lastSnpLocus, bucket[lastSnpLocus])

    snps.close()
    
    logging.info("Added %s SNP kmers to bucket (%s loci), new bucket size: %s loci", addedSNPsCount, lociAddedCount, len(bucket))
    
    return (bucket, genomeBits)
    

    
def fillBucket(genomeList, k ):
    # read buffer (trying this out)
    bufferSize = 32 * 1024
    # To go faster, might need to switch to asynchronous file I/O.

    centerpos = int((k-1)/2)
    maxGenomeBit = 1
    # Mapping of genomes to the associated bit in the bit field.
    genomeBits = {}
    
    thisBucket = {}
    # Default / initial value for a given locus.
    noSNPsNoGenomes = {
        'A': 0,
        'C': 0,
        'G': 0,
        'T': 0,
        }

    for genome in genomeList:
        logging.info("Loading genome: %s", genome)
        genomeBits[genome] = maxGenomeBit
        thisGenome = open(genome,'r')
        nextLines = thisGenome.readlines(bufferSize)
        while nextLines != []:
            for kmer in nextLines:
                center = kmer[centerpos]
                locus=kmer[:centerpos] + '.' + kmer[centerpos+1:k]
                # Set the bit corresponding to this central nucleotide for this genome
                thisBucket[locus][center] = thisBucket.setdefault(locus, noSNPsNoGenomes.copy())[center] | maxGenomeBit
                # logging.debug("%s: %s", locus, thisBucket[locus])
            nextLines = thisGenome.readlines(bufferSize)
        maxGenomeBit = maxGenomeBit * 2 # Move to a more significant bit.
        thisGenome.close()
        logging.debug("bucket size: %s loci", len(thisBucket))

    return (thisBucket, genomeBits)


def findConflicts(genomeSets):
    conflictingGenomes = 0
    remainingGenomeSets = genomeSets.copy()
    
    while remainingGenomeSets:
        (nucleotideA, setA) = remainingGenomeSets.pop()
        for (nucleotideB, setB) in remainingGenomeSets:

            # This should run for every pair of nucleotide integers.
            # It will set any bits in conflictingGenomes that correspond to a genome
            # that has more than one nucleotide associated with this locus.
            conflictingGenomes = conflictingGenomes | ( setA & setB )

    return conflictingGenomes


def deconflict(genomeSets):
    # Figure out which genomes have conflicts in this locus.
    conflicts = findConflicts(genomeSets)
    # Remove all references to conflicted genomes from the genome sets, trim
    # unused center nucleotides.
    return [ (nucleotide, genomeSet & ~ conflicts ) for (nucleotide, genomeSet) in genomeSets if genomeSet & ~ conflicts != 0]


# def writeSNP(locus, nucleotide, genome, SNPsFile, k):
def writeSNP(locus, nucleotide, genome, centerPos):
    newline = '\n'
    tab = '\t'
    
    genome['fasta'].write('>' + locus + '_' + nucleotide + newline + locus[:centerPos] + nucleotide + locus[centerPos+1:] + newline)
    genome['snp'].write(tab.join([locus, nucleotide, locus[:centerPos] + nucleotide + locus[centerPos+1:], 'FILENAME']) + newline)
    
    # If we also want to create a central list of SNPs as well NOW, we can do it here.
    # SNPsFile.write(locus[:centerpos] + nucleotide + locus[centerpos+1:] + newline)



# def dumpBucket(bucket, genomeBits, genomeList, SNPsFile, k):
def dumpBucket(bucket, genomeBits, genomeList, k):
    # We then browse through values looking for:
    #   Only one center char in a given mask (not a SNP)
    #   The same genome in the same mask twice, with different center chars, and if not: SNP??

    # Identify the index of the center nucleotide (the polymorphic one)
    centerPos = int((k-1)/2)

    # Collect statistics for debugging
    totalSNPs = 0
    writtenKmers = 0
    
    # create reverse map from genomeBit to genome:
    bitGenomes = [ (bitValue, genome) for (genome, bitValue) in genomeBits.items() ]
    
    # Now the thisBucket hash table is fully populated...
    for (locus, nucleotideGenomes) in bucket.items():
        # First, deconflict; remove from consideration genomes that show up for more than
        # one nucleotide for this locus:
        kmerGenomeMap = list(nucleotideGenomes.items())

        # Remove all references to conflicted genomes from the genome sets, trim
        # unused center nucleotides.
        deconflicted = deconflict(kmerGenomeMap)

        if len(deconflicted) > 1:
            # logging.debug("deconflicted left something to play with: %s", deconflicted)
            # Deconfliction removes any nucleotides with no corresponding genomes after
            # the deconfliction is done.  The possible outcomes are:
            #  0 nucleotides: All of the hits were from genomes with conflicts.
            #  1 nucleotide: After all conflicts removed, there is no polymorphism (only
            #      one nucleotide represented means no mutation)
            #  2-4 nucleotides: There is an SNP here.
            
            totalSNPs = totalSNPs + 1 # Count the number of SNPs found.
            # logging.debug("Found SNP: %s", deconflicted)

            # These two loops should be v. fast; small lists, bitwise operations, at most
            # 4*len(genomeList) tests and len(genomeList)*2 writes.
            for (nucleotide, genomeSet) in deconflicted:
                for (bitValue, genome) in bitGenomes:
                    if bitValue & genomeSet != 0:
                        # This should be called at most len(genomeList) times.
                        # writeSNP(locus, nucleotide, genome, SNPsFile)
                        writtenKmers = writtenKmers + 1
                        writeSNP(locus, nucleotide, genomeList[genome], centerPos)
    return (totalSNPs, writtenKmers)

                        

def getK(genomeList):
    # Don't do this a _lot_, but it is low cost to do once in a while.
    sampleFile = open(genomeList[0],'r')
    # Grab the first space delimited "word", which should be the kmer
    k = len(sampleFile.readline().strip().split(" ")[0])
    sampleFile.close()
    return k
    

def findSNPs(genomeList, fastaSuffix, snpSuffix, snpsAll=None):
    k = getK(genomeList)

    logging.debug('Looking for SNPs in %s genomes', len(genomeList))
    (thisBucket, genomeBits) = fillBucket(genomeList, k)
    logging.debug('Found %s loci in the %s genomes', len(thisBucket), len(genomeList))

    if snpsAll is not None:
        logging.info("snpsAll: %s -> Therefore we are going to load SNPs from that file", snpsAll)
        # Start genome bit for new genomes:
        nextGenomeBit = max(list(genomeBits.values())) * 2
        (nextBucket, newGenomeBits) = loadSNPsFromSNPsAll(thisBucket, snpsAll, nextGenomeBit)

        # Swap in our new bucket. (may not be necessary if we modified the one we passed in-place?)
        thisBucket = nextBucket

        # Nothing else to do;  dumpBucket will only dump the new (non-SNPs_all related) genomes, and
        # those only if they either intersect with the SNPs_all data after deconflicting, or they
        # internally (within new genomes) create new SNPs.
        
        
    outputFiles = {}
    for genome in genomeList:
        outputFiles[genome] = {
            'fasta': open(genome + fastaSuffix, 'w'),
            'snp':   open(genome + snpSuffix, 'w'),
            }
        # Avoid the cost of flushing at the end of every line.
        outputFiles[genome]['fasta'].reconfigure(line_buffering=False)
        outputFiles[genome]['snp'].reconfigure(line_buffering=False)
    
    # dumpBucket(thisBucket, genomeBits, genomeList, SNPsFile, k)
    (totalSNPs, writtenKmers) = dumpBucket(thisBucket, genomeBits, outputFiles, k)
    logging.info('Found %s SNPs in this partition, wrote %s SNP kmers', totalSNPs, writtenKmers)

    # Clean up, flush data to disk.
    for outputFiles in outputFiles.values():
        for fileType in outputFiles.values():
            fileType.close()


     
            
####################################
####################################
###
###
### Main
###
###
####################################
####################################


if __name__ == "__main__":

    options = parseCommandline()

    if options.debug:
        logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
    elif options.info:
        logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
    else:
        logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.WARN)
            

    logging.debug('Commandline options, as parsed: %s', str(options))

    # Primary function of this script:
    #  - Load kmers into memory,
    #  - Traverse kmers, filtering out conflicts, 


    genomeFiles = options.filteredGenomes # Looks like fsplitX.SNPs.partY.mers AND fsplitX.SNPs.partY.mers.fasta (last for mummer)

    # This output will be piped into mummer to get position information to link us with annotation data.
    findSNPs(genomeFiles, options.fasta_suffix, options.snp_suffix, options.snps_all)

