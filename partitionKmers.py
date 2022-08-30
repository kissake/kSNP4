#!/usr/bin/python3

import sys
import bisect as bisect
import argparse as argparse
import logging as logging


def getPartitionKmersArgs( override = None ):

    description = '''Take input files containing lists of kmers and distribute their contents line-by-line to different sub-files 
(buckets) based on the kmer value, such that like kmers are in the same bucket.'''

    example = '''To partition the file fsplit0.kmer into fsplit0.kmer.part0 fsplit0.kmer.part1, given the file partition.txt that consists of a single line with the value 'CACCGTACAT':

%(prog)s --partition partition.txt --suffix .part fsplit0.kmer'''
    
    parser = argparse.ArgumentParser(description=description, epilog = example,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('--partition', help='The file containing (one-per-line) the kmer prefix to use as the divider between files', default=sys.stdin)
    parser.add_argument('--suffix', help='The suffix to use before the number indicating which part is being output.', default='.part')
    parser.add_argument('kmers', metavar='kmerFile', nargs='+', help='A list of kmers to be filtered into multiple partitins (probably for multiprocessing)')

    parser.add_argument('--debug', action='store_true', help='Output diagnostic data to the screen using STDERR')

    if override is None:
        return parser.parse_args()
    else:
        return parser.parse_args(override)


if __name__ == '__main__':

    
    options = getPartitionKmersArgs()

    if options.debug:
        logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
    else:
        logging.basicConfig(format='%(levelname)s:%(message)s')

    logging.debug('Commandline options, as parsed: %s', str(options))

    # Buffer to improve processing speed for large files.
    # Note that hasn't been tested on gigabyte-sized files.  Could be
    # some room for performance improvement at that level.
    bufferSize = 32 * 1024

    fileSuffix = options.suffix


    if options.partition == sys.stdin:
        partitionsFile = options.partition
    else:
        partitionsFile = open(options.partition, 'r')
        
    # Note that 'sorted' is likely redundant; the output of guessPartitions.py
    # will be sorted.
    partitions = sorted([ line.strip() for line in partitionsFile.readlines() ])

    # Our input is kmers; all kmers are less than ZZZZ.  This value is
    # intended to provide a terminal bucket for kmers larger than the last
    # partition value.
    partitions.append('ZZZZ') 
    sys.stdin.close()

    for file in options.kmers:
        outFiles = {}
        input = open(file, 'r')

        count = 0

        for partition in partitions:
            outFiles[partition] = open(file + fileSuffix + str(count),'w', bufferSize)
            logging.debug("Opened partition file: %s", file + fileSuffix + str(count))
            # Set the output file to NOT use line buffering because it results
            # in WAY too many writes.
            outFiles[partition].reconfigure(line_buffering=False)
            count = count + 1

        block = input.readlines(bufferSize)
        while len(block) > 0:
            for line in block:
                # This line is... obscure looking.  It will write the line to
                # the correct output file based on its content and the partitions
                # that were input.
                outFiles[partitions[bisect.bisect_right(partitions,line)]].write(line)
            block = input.readlines(bufferSize)


        for outputFile in outFiles.values():
            outputFile.close()
        input.close()
