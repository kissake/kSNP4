#!/usr/bin/python

'''
usage: guessPartition.py [-h] [--all] [--samples SAMPLES] [--input INPUT] [--debug] partitions

Guesses good values for partitioning lists of kmers based on a subset of example data.

positional arguments:
  partitions         The number of buckets of data you want to generate. The output will be a series of kmer prefixes one fewer than the number of partitions.

optional arguments:
  -h, --help         show this help message and exit
  --all              Parse _all_ values in the source file. Useful if the source file is sorted.
  --samples SAMPLES  The number of samples to use (overrides the defailt)
  --input INPUT      The input file to read from (overrides the default: STDIN)
  --debug            Output diagnostic data to STDERR

NOTE: This program requires UNSORTED input for it to function with sampling. Otherwise you should use --all to ensure that you don't get worst-case behavior.


This program's function is to guess a good place to partition the kmer space
based on:
 - How many partitions are needed, and
 - The variation in the first hundred thousand records in the largest (any?) genome.

NOTE NOTE NOTE:
If the input is sorted, this partitioning will be quite wrong.  If that is
the case for you, you should use --all.

The input is a number of partitions, specified on the commandline, and the 
filename of a file that lists kmers as the first element on the line 
(starting at the first column) with minimal repetition of kmers.

While we do sort(), we limit the cost by limiting the number of lines.  On
the author's system, this takes under a tenth of a second, and lets us maximize
our use of memory for speed in dependent scripts.

Samples values as low as 1000 are likely fine, but don't offer significant time
savings.  The more samples, the more evenly distributed the kmers should be, 
with significantly diminishing returns.
'''

import sys as sys
import argparse as argparse
import logging as logging


def getGuessPartitionArgs(override=None):
    parser = argparse.ArgumentParser(description='Guesses good values for partitioning lists of kmers based on a subset of example data.', epilog="NOTE: This program requires UNSORTED input for it to function with sampling.  Otherwise you should use --all to ensure that you don't get worst-case behavior.")
    parser.add_argument('partitions', action='store', help='The number of buckets of data you want to generate.  The output will be a series of kmer prefixes one fewer than the number of partitions.')
    parser.add_argument('--all', action='store_true', help='Parse _all_ values in the source file.  Useful if the source file is sorted.')
    parser.add_argument('--samples', action='store', help='The number of samples to use (overrides the defailt)', default=100000)
    parser.add_argument('--input', action='store', help='The input file to read from (overrides the default: STDIN)', default=sys.stdin)
    parser.add_argument('--debug', action='store_true', help='Output diagnostic data to STDERR')

    if override is None:
        return parser.parse_args()
    else:
        return parser.parse_args(override)
    
def moreCleverSortKey(kmer):
    # UNUSED
    # Will sort SNPs to be adjacent.
    output = ''
    centerPos = int(len(kmer) / 2) - 1
    for count in range( centerPos ):
        output = output + kmer[count] + kmer[-count]
    output = output + kmer[centerPos]
    return output

def getKeygenFunction(N):
    '''
    Produces a function to quickly generate a key for a string
    where the Nth character is the center character.  Maybe faster than
    the above?  Worth testing.  Yay lambda and recursion!
    '''
    if N == 1:
        return lambda a: a
    else:
        return lambda a: a[0] + a[-1] + getKeygen(N-1)(a[1:-1])



if __name__ == '__main__':

    options = getGuessPartitionArgs()

    if options.debug:
        logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
    else:
        logging.basicConfig(format='%(levelname)s:%(message)s')

    logging.debug('Commandline options, as parsed: %s', str(options))

    if not options.all:
        
        maxLines = int(options.samples)

    
    

    # Performance drops off significantly below 4k and above ~128k
    bufferSize=32*1024

    # Ensure input file is open.
    inputFile = options.input
    if inputFile == sys.stdin:
        input = inputFile
    else:
        input = open(inputFile,'r')

    array = []

    batch = input.readlines(bufferSize)
    while batch != []:
        array.extend(batch)
        if not options.all and len(array) > maxLines:
            break

        batch = input.readlines(bufferSize)

    input.close()

    # This or the I/O is going to take the bulk of the time.  If we have
    # to use swap to hold the data, this will probably overwhelm.
    array.sort()

    # Assume that the kmer is the first item on the line, and that
    # the line is space-separated.
    logging.debug("Length of the sorted array: %s", len(array))
    k = len(array[0].split(' ')[0])

    # Identify the center position of the kmer (this is the location
    # of the SNP)
    centerPos = int((k-1)/2)

    buckets = int(options.partitions)

    partitionStep = int(len(array)/buckets)

    for X in range(1,buckets):
        # We need to make sure we don't establish a partition in the middle of a locus.
        # To ensure this, we truncate the output to avoid including the center nucleotide
        # If more granular sorting is required (more granular than the first half of the
        # kmer permits), you'll need a more clever sort key (but it is possible)
        logging.debug("X: %s, partitionStep: %s, index: %s, array[index]: %s ",X, partitionStep, X*partitionStep, array[X*partitionStep])
        print( array[X * partitionStep].strip()[:centerPos] )
