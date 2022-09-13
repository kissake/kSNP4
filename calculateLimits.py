#!/usr/bin/python3

import psutil as psutil
import resource as resource
import sys as sys
import os as os
import logging as logging
import math as math
import argparse as argparse


def parseCommandline(override=None):

    description = '''Make a best estimate as to number of partitions to use, number of open files that will
be needed, and number of physical CPUs available'''
    
    example = '''Typical usage is:

NUMPARTS=`$(prog)s in_file-Eco12 19 0 | cut -f 1`
'''
    
    parser = argparse.ArgumentParser(description= description, epilog=example,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('fileToGenome', metavar='input_file', help='Input file as you would pass it to kSNP4 -in <>')
    parser.add_argument('k', help='Value of k (kmer length) as you would pass it to kSNP4.')
    parser.add_argument('parallelism', help='Number of threads you expect to run at once.  Value of 0 means one per physical CPU.')
    parser.add_argument('--kmers', help='For a better estimate as to number of needed partitions, you can provide a more accurate number of kmers (otherwise number of kmers will be estimated from size of the fasta file, possibly guessing 1k times too high in some cases.)', default=None)

    parser.add_argument('--debug', action='store_true', help='Output diagnostic data to the screen using STDERR')
    parser.add_argument('--info', action='store_true', help='Output status information on STDERR', default=True)
    parser.add_argument('--quiet', action='store_true', help='Silence debug and other messages. (warnings only)')

    if override is None:
        return parser.parse_args()
    else:
        return parser.parse_args(override)





def parseFilesToGenome(filename):
    input = open(filename,'r')
    columns = [ 'filename', 'genome' ]
    
    genomeData = [ {
        key : value for (key, value) in zip(columns, line.split())
    } for line in input.readlines() ]

    input.close()

    return genomeData


def estimateResourceRequirements(kmerEstimate, numberGenomes, k, parallelism=0):
    # Estimated bytes per kmer.  Get this from find_snps.py if possible.
    bytesPerKmer = 115
    
    # Get the external reality that we are working within.
    numberCPUs = psutil.cpu_count() # This could be None.
    physicalRAM = psutil.virtual_memory().total
    currentOpenFiles, maxOpenFiles = resource.getrlimit(resource.RLIMIT_NOFILE)

    # Derate; assume a significant fraction is going into useful other things.
    availableRAM = physicalRAM * 0.75

    # Default to number of CPUs if parallelism isn't set (0 is an indicator value)
    if parallelism == 0:
        parallelism = numberCPUs
    
    memoryRequired = kmerEstimate * bytesPerKmer
    
    multiplier = 1
    if memoryRequired > int(availableRAM):
        multiplier = memoryRequired / availableRAM

    # The number of partitions we need is the multiplier above times the amount of
    # parallelism.  This ensures that when we have parallelism number of processes running,
    # we still aren't using more than the availableRAM (derated).
    numberParts = math.ceil( multiplier * parallelism )

    # Our processes use a lot of open files.  One holds one open file for every partition, and
    # another holds two for every genome in the list.  Then we add 3 to account for stdin, stdout
    # and stderr.
    openFilesNeeded = max( numberGenomes * 2, numberParts ) + 3

    openFilesRequested = openFilesNeeded + 20

    if openFilesRequested > maxOpenFiles:
        logging.warn("WARNING: we are likely to run out of open file handles; try increasing the resource limits per: https://medium.com/mindful-technology/too-many-open-files-limit-ulimit-on-mac-os-x-add0f1bfddde or man limits.conf")

    logging.debug('Based on a total of %s kmers of length %s from %s genomes, split over %s threads at once,', kmerEstimate, k, numberGenomes, parallelism)
    logging.debug('Recommend using %s partitions (%s max files) over %s genomes (%s max files), which means maxing out at %s files.', numberParts, numberParts + 3, numberGenomes, (numberGenomes * 2)+3, openFilesRequested)

    return numberParts, openFilesRequested, numberCPUs
    

def guessFromInFile(inFile):


    # Get the provided genome data.
    genomes = parseFilesToGenome( inFile )

    # Calculate the total genome size (roughly; this should be an over-estimate)
    
    totalSize = 0

    for file in [ genome['filename'] for genome in genomes ]:
        totalSize = totalSize + os.stat(file).st_size

    # Roughly one kmer per byte.  Only becomes trivially non-true when k becomes crazy-big.
    kmerEstimate = totalSize * k

    return len(genomes), kmerEstimate



if __name__ == '__main__':
    options = parseCommandline()

    if options.debug:
        logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
    elif options.info:
        logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
    else:
        logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.WARN)
            

    logging.debug('Commandline options, as parsed: %s', str(options))

    
    # Estimated bytes per kmer.
    fileToGenome = options.fileToGenome
    k = int(options.k)
    parallelism = int(options.parallelism)
    kmers = options.kmers

    genomeCount, kmerEstimate = guessFromInFile(fileToGenome)

    if kmers is None:
        logging.info("No estimate of # of kmers provided.  Using guess from input file: %s", kmerEstimate)
        kmers = kmerEstimate
    else:
        kmers = int(kmers)
    
    parts, openFiles, numberCPUs = estimateResourceRequirements(kmers, genomeCount, k, parallelism)

    
    print("\t".join([str(x) for x in [ parts, openFiles, numberCPUs ]]))
    
